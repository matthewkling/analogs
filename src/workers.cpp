// src/workers.cpp
#include "workers.h"
#include "ridge_solve.h"
#include "profiling.h"

namespace analogs {

// Define thread-local storage
thread_local PairWorker::ThreadLocalStorage PairWorker::tls;
thread_local AggWorker::ThreadLocalStorage AggWorker::tls;

// ---------------------------------------------------------------------------
// Helper: Check if focal point has NA in coordinates or any environment variable
// Returns true if any value is NA/NaN, meaning this focal should be skipped.
// ---------------------------------------------------------------------------
inline bool focal_has_na(const double* focal_ptr, std::size_t i,
                         int stride_f, int n_env) {
      // Check coordinates
      if (ISNAN(focal_ptr[i]) || ISNAN(focal_ptr[i + stride_f])) {
            return true;
      }
      // Check environment values (strided access)
      for (int d = 0; d < n_env; ++d) {
            if (ISNAN(focal_ptr[i + (2 + d) * stride_f])) {
                  return true;
            }
      }
      return false;
}

// ---------------------------------------------------------------------------
// Helpers for per-point pool weights (cell-area and user-supplied).
// Each returns 1.0 when the corresponding mechanism is inactive, otherwise
// the per-point value at 0-based pool index j.
// ---------------------------------------------------------------------------
static inline double get_area_weight(bool has, const double* ptr, std::size_t j) {
      return has ? ptr[j] : 1.0;
}
static inline double get_user_weight(bool has, const double* ptr, std::size_t j) {
      return has ? ptr[j] : 1.0;
}

void PairWorker::operator()(std::size_t begin, std::size_t end) {
      const bool knn_mode      = (scode == SelectCode::KNN_ENV || scode == SelectCode::KNN_GEOG);
      const bool rank_by_env  = (scode == SelectCode::KNN_ENV);
      const bool rank_by_geog  = (scode == SelectCode::KNN_GEOG);

      // Get thread-local storage
      auto& fgeo_vec = tls.fgeo_vec;
      auto& fenv_vec = tls.fenv_vec;
      auto& q_geo = tls.q_geo;
      auto& q_env = tls.q_env;
      auto& cand = tls.cand;
      auto& cand_weights = tls.cand_weights;

      // Pre-compute inverse covariance matrices if using Mahalanobis
      std::vector< std::vector<double> > inv_cov_matrices;
      if (use_mahalanobis) {
            inv_cov_matrices.resize(end - begin);

            for (std::size_t i = begin; i < end; ++i) {
                  const size_t local_idx = i - begin;

                  // Skip NA focals - leave inv_cov empty (won't be used)
                  if (focal_has_na(focal_ptr, i, stride_f, n_env)) {
                        continue;
                  }

                  // Extract covariance components for focal i from column-major x_cov matrix
                  // x_cov is n_focal × n_cov_components stored column-major
                  // Row i, col j is at: x_cov_ptr[i + j * x_cov_stride]
                  std::vector<double> cov_vec(n_cov_components);
                  for (int comp = 0; comp < n_cov_components; ++comp) {
                        cov_vec[comp] = x_cov_ptr[i + comp * x_cov_stride];
                  }

                  // Reconstruct covariance matrix
                  std::vector<double> cov_matrix;
                  reconstruct_cov_matrix(cov_vec.data(), n_env, cov_matrix);

                  // Invert covariance matrix
                  std::vector<double>& inv_cov = inv_cov_matrices[local_idx];
                  if (!invert_cov_matrix(cov_matrix, n_env, inv_cov)) {
                        // Matrix not positive definite - use identity (Euclidean)
                        inv_cov.resize(n_env * n_env, 0.0);
                        for (int k = 0; k < n_env; ++k) {
                              inv_cov[k * n_env + k] = 1.0;
                        }
                  }
            }
      }

      for (std::size_t i = begin; i < end; ++i) {
            const double fx = focal_ptr[i];
            const double fy = focal_ptr[i + stride_f];
            const double* f_env_col = focal_ptr + i + 2 * stride_f;

            // -----------------------------------------------------------------
            // Check for NA in focal coordinates or environment values.
            // If NA, leave out_indices[i] and out_weights[i] empty.
            // Post-processing in core.cpp handles k=1 case by inserting 0 sentinel.
            // -----------------------------------------------------------------
            if (focal_has_na(focal_ptr, i, stride_f, n_env)) {
                  // out_indices[i] and out_weights[i] are already empty
                  continue;
            }

            // -----------------------------------------------------------------
            // Fast path: KNN modes + lattice + (planar OR ECEF)
            // Note: Mahalanobis doesn't use the fast path because lattice
            // needs bounding box adjustments. Falls through to regular path.
            // -----------------------------------------------------------------
            if (knn_mode &&
                use_lattice &&
                lattice_ptr != nullptr &&
                !use_haversine &&
                !use_mahalanobis)  // Mahalanobis uses regular path
            {
                  // Use lattice KNN expansion (Euclidean only)
                  const size_tu n_geo = lattice_ptr->n_geo_dims;
                  fgeo_vec.resize(n_geo);
                  fenv_vec.resize(n_env);

                  // Geographic coords
                  if (use_ecef && n_geo == 3) {
                        double X, Y, Z;
                        lonlat_to_ecef(fx, fy, R_earth, X, Y, Z);
                        fgeo_vec[0] = X;
                        fgeo_vec[1] = Y;
                        fgeo_vec[2] = Z;
                  } else {
                        fgeo_vec[0] = fx;
                        fgeo_vec[1] = fy;
                  }

                  // Environment coords
                  for (int d = 0; d < n_env; ++d) {
                        fenv_vec[d] = f_env_col[d * stride_f];
                  }

                  // Expanded KNN search (Euclidean)
                  // Reuse cand vector for results
                  cand.clear();
                  cand_weights.clear();

                  const double geog_thresh = use_ecef ?
                  max_geog_chord : max_geog;

                  // Matching inner-radius threshold for the annulus. Uses the
                  // same metric as geog_thresh: chord distance on the ECEF path,
                  // planar distance otherwise. 0 when inactive (min_geog check
                  // is skipped inside knn_query).
                  const double geog_min_thresh = use_geog_min
                  ? (use_ecef ? min_geog_chord : min_geog)
                        : 0.0;

                        {
                              ANALOGS_PROFILE_SCOPE(GATHER);
                              lattice_ptr->knn_query(
                                          fgeo_vec.data(),
                                          fenv_vec.data(),
                                          ref_ptr,
                                          static_cast<size_tu>(stride_r),
                                          static_cast<size_tu>(n_env),
                                          rank_by_geog,
                                          geog_thresh,
                                          geog_min_thresh,
                                          use_scalar_env,
                                          max_env_pervar,
                                          max_env_scalar,
                                          k,
                                          cand,
                                          cand_weights
                              );
                        }

                        // Convert 0-based indices to 1-based and emit three weight streams.
                        const int m = static_cast<int>(cand.size());
                        std::vector<int>    keep(m);
                        std::vector<double> sw(m), aw(m), uw(m);
                        for (int t = 0; t < m; ++t) {
                              const std::size_t j0 = static_cast<std::size_t>(cand[t]);
                              keep[t] = static_cast<int>(j0) + 1;
                              sw[t]   = cand_weights[t];
                              aw[t]   = get_area_weight(has_area_weight, area_weight_ptr, j0);
                              uw[t]   = get_user_weight(has_user_weight, user_weight_ptr, j0);
                        }
                        out_indices[i]        = std::move(keep);
                        out_sample_weights[i] = std::move(sw);
                        out_area_weights[i]   = std::move(aw);
                        out_user_weights[i]   = std::move(uw);
                        continue;
            }

            // -----------------------------------------------------------------
            // Fallback path: gather candidates
            // (Used for: Haversine, Mahalanobis, or ALL mode)
            // -----------------------------------------------------------------
            cand.clear();
            cand_weights.clear();

            if (use_lattice && lattice_ptr != nullptr) {
                  const size_tu n_geo = lattice_ptr->n_geo_dims;
                  q_geo.resize(n_geo);

                  if (use_ecef && n_geo == 3) {
                        double X, Y, Z;
                        lonlat_to_ecef(fx, fy, R_earth, X, Y, Z);
                        q_geo[0] = X;
                        q_geo[1] = Y;
                        q_geo[2] = Z;
                  } else {
                        q_geo[0] = fx;
                        q_geo[1] = fy;
                  }

                  q_env.resize(n_env);
                  for (int d = 0; d < n_env; ++d) {
                        q_env[d] = f_env_col[d * stride_f];
                  }

                  // Adjust environment bounds if using Mahalanobis
                  double effective_max_env = max_env_scalar;
                  std::vector<double> mahal_bounds;

                  if (use_mahalanobis && std::isfinite(max_env_scalar)) {
                        const std::vector<double>& inv_cov = inv_cov_matrices[i - begin];
                        mahalanobis_bounding_box(q_env.data(), inv_cov, n_env,
                                                 max_env_scalar, mahal_bounds);
                        // Note: We could use mahal_bounds to adjust lattice query,
                        // but for simplicity we just use a conservative max_env.
                        // The actual Mahalanobis filtering happens below.
                  }

                  {
                        ANALOGS_PROFILE_SCOPE(GATHER);
                        lattice_ptr->query(
                                    q_geo.data(),
                                    q_env.data(),
                                    max_geog,
                                    use_scalar_env,
                                    max_env_pervar,
                                    effective_max_env,
                                    cand,
                                    cand_weights
                        );
                  }
            } else {
                  cand.resize(n_ref);
                  cand_weights.resize(n_ref, 1.0);
                  for (size_t j = 0; j < static_cast<size_t>(n_ref); ++j) {
                        cand[j] = static_cast<index_t>(j);
                  }
            }

            // Filter and collect results from candidates
            if (scode == SelectCode::ALL) {
                  ANALOGS_PROFILE_SCOPE(EXACT);
                  std::vector<int>    keep;
                  std::vector<double> sw, aw, uw;
                  keep.reserve(cand.size());
                  sw.reserve(cand.size());
                  aw.reserve(cand.size());
                  uw.reserve(cand.size());

                  double fx_ecef = 0, fy_ecef = 0, fz_ecef = 0;
                  if (use_ecef) {
                        lonlat_to_ecef(fx, fy, R_earth, fx_ecef, fy_ecef, fz_ecef);
                  }

                  // Get inverse covariance matrix for this focal (if using Mahalanobis)
                  const std::vector<double>* inv_cov_ptr = nullptr;
                  if (use_mahalanobis) {
                        inv_cov_ptr = &inv_cov_matrices[i - begin];
                  }

                  for (size_t t = 0; t < cand.size(); ++t) {
                        const index_t j = cand[t];
                        const double weight = cand_weights[t];
                        const double rx = ref_ptr[j];
                        const double ry = ref_ptr[j + stride_r];

                        // Geog distance & filter (annulus: min <= gdist <= max).
                        // gdist is needed when either bound is active.
                        if (use_geog_filter || use_geog_min) {
                              double gdist;
                              if (use_ecef) {
                                    const double rx_ecef = ref_latt_ptr[j];
                                    const double ry_ecef = ref_latt_ptr[j + stride_latt_r];
                                    const double rz_ecef = ref_latt_ptr[j + 2 * stride_latt_r];
                                    const double dx = fx_ecef - rx_ecef;
                                    const double dy = fy_ecef - ry_ecef;
                                    const double dz = fz_ecef - rz_ecef;
                                    gdist = std::sqrt(dx*dx + dy*dy + dz*dz);
                                    if (use_geog_filter && gdist > max_geog_chord) continue;
                                    if (use_geog_min && gdist < min_geog_chord) continue;
                              } else {
                                    gdist = geo_distance_km(fx, fy, rx, ry, use_haversine);
                                    if (use_geog_filter && gdist > max_geog) continue;
                                    if (use_geog_min && gdist < min_geog) continue;
                              }
                        }

                        // Environment checks
                        const double* r_env_col = ref_ptr + j + 2 * stride_r;
                        bool ok;

                        if (use_mahalanobis) {
                              // Use Mahalanobis distance
                              auto okd = mahalanobis_ok_and_dist(
                                    f_env_col, r_env_col,
                                    *inv_cov_ptr,
                                    n_env, stride_f, stride_r,
                                    max_env_scalar, false
                              );
                              ok = okd.first;
                        } else {
                              // Use Euclidean distance
                              auto okd = env_ok_and_dist(
                                    f_env_col, r_env_col,
                                    n_env, stride_f, stride_r,
                                    use_pervar_env, max_env_pervar,
                                    use_scalar_env, max_env_scalar
                              );
                              ok = okd.first;
                        }

                        if (!ok) continue;

                        const std::size_t j0 = static_cast<std::size_t>(j);
                        keep.push_back(static_cast<int>(j0) + 1); // 1-based
                        sw.push_back(weight);
                        aw.push_back(get_area_weight(has_area_weight, area_weight_ptr, j0));
                        uw.push_back(get_user_weight(has_user_weight, user_weight_ptr, j0));
                  }

                  out_indices[i]        = std::move(keep);
                  out_sample_weights[i] = std::move(sw);
                  out_area_weights[i]   = std::move(aw);
                  out_user_weights[i]   = std::move(uw);
                  continue;
            }

            // kNN modes (Environment or Geog) without expanding-search lattice path
            ANALOGS_PROFILE_SCOPE(EXACT);
            auto cmp = [](const Neighbor& a, const Neighbor& b) {
                  return a.first < b.first; // max-heap: top has largest distance
            };
            std::priority_queue<Neighbor, std::vector<Neighbor>, decltype(cmp)> pq(cmp);

            // Build map from candidate index to weight for later lookup
            std::unordered_map<index_t, double> idx_to_weight;
            for (size_t t = 0; t < cand.size(); ++t) {
                  idx_to_weight[cand[t]] = cand_weights[t];
            }

            // ECEF focal coords (if needed)
            double fx_ecef = 0, fy_ecef = 0, fz_ecef = 0;
            if (use_ecef) {
                  lonlat_to_ecef(fx, fy, R_earth, fx_ecef, fy_ecef, fz_ecef);
            }

            // Get inverse covariance matrix for this focal (if using Mahalanobis)
            const std::vector<double>* inv_cov_ptr = nullptr;
            if (use_mahalanobis) {
                  inv_cov_ptr = &inv_cov_matrices[i - begin];
            }

            for (size_t t = 0; t < cand.size(); ++t) {
                  const index_t j = cand[t];
                  const double rx = ref_ptr[j];
                  const double ry = ref_ptr[j + stride_r];

                  // Geog distance & filter (annulus: min <= gdist <= max)
                  double gdist;
                  if (use_ecef) {
                        // Use ECEF chord distance
                        const double rx_ecef = ref_latt_ptr[j];
                        const double ry_ecef = ref_latt_ptr[j + stride_latt_r];
                        const double rz_ecef = ref_latt_ptr[j + 2 * stride_latt_r];
                        const double dx = fx_ecef - rx_ecef;
                        const double dy = fy_ecef - ry_ecef;
                        const double dz = fz_ecef - rz_ecef;
                        gdist = std::sqrt(dx*dx + dy*dy + dz*dz);
                        if (use_geog_filter && gdist > max_geog_chord) continue;
                        if (use_geog_min && gdist < min_geog_chord) continue;
                  } else {
                        gdist = geo_distance_km(fx, fy, rx, ry, use_haversine);
                        if (use_geog_filter && gdist > max_geog) continue;
                        if (use_geog_min && gdist < min_geog) continue;
                  }

                  // Environment checks & distance
                  const double* r_env_col = ref_ptr + j + 2 * stride_r;
                  double env_dist;
                  bool ok;

                  if (use_mahalanobis) {
                        // Use Mahalanobis distance
                        auto okd = mahalanobis_ok_and_dist(
                              f_env_col, r_env_col,
                              *inv_cov_ptr,
                              n_env, stride_f, stride_r,
                              max_env_scalar, true
                        );
                        ok = okd.first;
                        env_dist = okd.second;
                  } else {
                        // Use Euclidean distance
                        auto okd = env_ok_and_dist(
                              f_env_col, r_env_col,
                              n_env, stride_f, stride_r,
                              use_pervar_env, max_env_pervar,
                              use_scalar_env, max_env_scalar
                        );
                        ok = okd.first;
                        env_dist = okd.second;
                  }

                  if (!ok) continue;

                  const double key = rank_by_env ?
                  env_dist : gdist;
                  const index_t ref_index_1based = j + 1;

                  if (static_cast<int>(pq.size()) < k) {
                        pq.emplace(key, ref_index_1based);
                  } else if (!pq.empty() && key < pq.top().first) {
                        pq.pop();
                        pq.emplace(key, ref_index_1based);
                  }
            }

            const int m = static_cast<int>(pq.size());
            std::vector<int>    idx_vec(m);
            std::vector<double> sw_vec(m), aw_vec(m), uw_vec(m);
            for (int pos = m - 1; pos >= 0; --pos) {
                  const Neighbor& nb = pq.top();
                  const int idx_1based = nb.second;
                  const std::size_t j0 = static_cast<std::size_t>(idx_1based - 1);

                  idx_vec[pos] = idx_1based;

                  // sample_weight from the candidate map (bin-level downsampling)
                  auto it = idx_to_weight.find(static_cast<index_t>(j0));
                  sw_vec[pos] = (it != idx_to_weight.end()) ? it->second : 1.0;
                  aw_vec[pos] = get_area_weight(has_area_weight, area_weight_ptr, j0);
                  uw_vec[pos] = get_user_weight(has_user_weight, user_weight_ptr, j0);

                  pq.pop();
            }
            out_indices[i]        = std::move(idx_vec);
            out_sample_weights[i] = std::move(sw_vec);
            out_area_weights[i]   = std::move(aw_vec);
            out_user_weights[i]   = std::move(uw_vec);
      }
}


void AggWorker::operator()(std::size_t begin, std::size_t end) {
      // Get thread-local storage
      auto& q_geo = tls.q_geo;
      auto& q_env = tls.q_env;
      auto& cand = tls.cand;
      auto& cand_weights = tls.cand_weights;

      // Determine which stats are regular vs value-based vs regression vs tabulate.
      // Also record the position of the WEIGHTED_MEAN stat within the value
      // stats (if any), since that's the only value stat with SE support.
      std::vector<int> regular_stat_positions;
      std::vector<int> value_stat_positions;
      int wm_value_idx = -1;            // index into value_stat_positions of WEIGHTED_MEAN, or -1
      bool has_regression = false;
      bool has_tabulate = false;

      for (int s = 0; s < static_cast<int>(acodes.size()); ++s) {
            if (acodes[s] == AggregateCode::SUM ||
                acodes[s] == AggregateCode::MEAN ||
                acodes[s] == AggregateCode::WEIGHTED_SUM ||
                acodes[s] == AggregateCode::WEIGHTED_MEAN) {
                  if (acodes[s] == AggregateCode::WEIGHTED_MEAN) {
                        wm_value_idx = static_cast<int>(value_stat_positions.size());
                  }
                  value_stat_positions.push_back(s);
            } else if (acodes[s] == AggregateCode::REGRESSION) {
                  // REGRESSION: regression is handled separately
                  has_regression = true;
            } else if (acodes[s] == AggregateCode::TABULATE) {
                  // TABULATE: handled separately, K_v columns per y variable
                  has_tabulate = true;
            } else {
                  regular_stat_positions.push_back(s);  // COUNT, SUM_WEIGHTS, MEAN_WEIGHTS, ESS
            }
      }

      const int n_regular = static_cast<int>(regular_stat_positions.size());
      const int n_value = static_cast<int>(value_stat_positions.size());

      // SE bookkeeping: do we need weighted_mean SE accumulators this call?
      const bool want_se = (se_code != SeCode::NONE);
      const bool wm_se = want_se && (wm_value_idx >= 0) && has_values;
      const bool wm_se_design = wm_se && (se_code == SeCode::DESIGN);

      // Check if any stats need weights
      bool need_weights = false;
      for (size_t s = 0; s < acodes.size(); ++s) {
            if (acodes[s] == AggregateCode::SUM_WEIGHTS ||
                acodes[s] == AggregateCode::MEAN_WEIGHTS ||
                acodes[s] == AggregateCode::WEIGHTED_SUM ||
                acodes[s] == AggregateCode::WEIGHTED_MEAN ||
                acodes[s] == AggregateCode::ESS ||
                acodes[s] == AggregateCode::REGRESSION ||
                acodes[s] == AggregateCode::TABULATE) {
                  need_weights = true;
                  break;
            }
      }

      // TABULATE: pre-compute per-variable bin offsets so we can index into a
      // single flat tabulate accumulator. tab_offsets[v] is the starting index
      // in tabulate_accum for variable v's K_v bins; tab_total is the total
      // number of tabulate bins across all y variables.
      std::vector<int> tab_offsets;
      int tab_total = 0;
      if (has_tabulate) {
            tab_offsets.resize(n_vars, 0);
            for (int v = 0; v < n_vars; ++v) {
                  tab_offsets[v] = tab_total;
                  if (v < static_cast<int>(n_classes_per_var.size())) {
                        tab_total += n_classes_per_var[v];
                  }
            }
      }

      // Pre-compute inverse covariance matrices if using Mahalanobis
      std::vector< std::vector<double> > inv_cov_matrices;
      if (use_mahalanobis) {
            inv_cov_matrices.resize(end - begin);

            for (std::size_t i = begin; i < end; ++i) {
                  const size_t local_idx = i - begin;

                  // Skip NA focals - leave inv_cov empty (won't be used)
                  if (focal_has_na(focal_ptr, i, stride_f, n_env)) {
                        continue;
                  }

                  // Extract covariance components for focal i from column-major x_cov matrix
                  std::vector<double> cov_vec(n_cov_components);
                  for (int comp = 0; comp < n_cov_components; ++comp) {
                        cov_vec[comp] = x_cov_ptr[i + comp * x_cov_stride];
                  }

                  // Reconstruct covariance matrix
                  std::vector<double> cov_matrix;
                  reconstruct_cov_matrix(cov_vec.data(), n_env, cov_matrix);

                  // Invert covariance matrix
                  std::vector<double>& inv_cov = inv_cov_matrices[local_idx];
                  if (!invert_cov_matrix(cov_matrix, n_env, inv_cov)) {
                        // Matrix not positive definite - use identity (Euclidean)
                        inv_cov.resize(n_env * n_env, 0.0);
                        for (int k = 0; k < n_env; ++k) {
                              inv_cov[k * n_env + k] = 1.0;
                        }
                  }
            }
      }

      // REGRESSION: pre-compute regression output dimension.
      // Each variable contributes reg_dim columns for coefficients, and
      // (when se_code != NONE) an additional reg_dim columns for SEs.
      const int reg_dim = has_covariates ? (n_covs + 1) : 0;
      const bool reg_se = want_se && has_regression;

      for (std::size_t i = begin; i < end; ++i) {
            // Skip focal points with NA
            if (focal_has_na(focal_ptr, i, stride_f, n_env)) {
                  // Leave all output as NA_REAL (pre-initialized)
                  continue;
            }

            const double fx = focal_ptr[i];
            const double fy = focal_ptr[i + stride_f];
            const double* f_env_col = focal_ptr + i + 2 * stride_f;

            // Get candidates from lattice
            if (use_lattice) {
                  q_geo.resize(use_ecef ? 3 : 2);
                  q_env.resize(n_env);

                  if (use_ecef) {
                        lonlat_to_ecef(fx, fy, R_earth, q_geo[0], q_geo[1], q_geo[2]);
                  } else {
                        q_geo[0] = fx;
                        q_geo[1] = fy;
                  }

                  for (int d = 0; d < n_env; ++d) {
                        q_env[d] = f_env_col[d * stride_f];
                  }

                  {
                        ANALOGS_PROFILE_SCOPE(GATHER);
                        lattice_ptr->query(q_geo.data(), q_env.data(),
                                           use_geog_filter ? max_geog : -1.0,
                                           use_scalar_env,
                                           max_env_pervar,
                                           max_env_scalar,
                                           cand, cand_weights);
                  }
            } else {
                  cand.resize(n_ref);
                  cand_weights.resize(n_ref, 1.0);
                  for (size_t j = 0; j < static_cast<size_t>(n_ref); ++j) {
                        cand[j] = static_cast<index_t>(j);
                  }
            }

            // Initialize accumulators for regular and value stats
            std::vector<double> regular_accum(n_regular, 0.0);
            std::vector<double> value_accum;
            if (has_values && n_value > 0) {
                  value_accum.resize(n_vars * n_value, 0.0);
            }

            // TABULATE: per-focal accumulator, length = sum of K_v across vars.
            // Layout is contiguous per-variable: var 0's bins, then var 1's, etc.
            std::vector<double> tabulate_accum;
            if (has_tabulate && tab_total > 0) {
                  tabulate_accum.assign(tab_total, 0.0);
            }

            // weighted_mean SE accumulators (per y-var).
            // Allocated only when actually needed.
            //   Common (ESS + DESIGN): sum_wy2
            //   DESIGN extra: sum_w2y, sum_w2y2, sum_w2 (scalar, not per-var)
            std::vector<double> wm_sum_wy2;
            std::vector<double> wm_sum_w2y;
            std::vector<double> wm_sum_w2y2;
            double wm_sum_w2 = 0.0;
            if (wm_se) {
                  wm_sum_wy2.assign(n_vars, 0.0);
                  if (wm_se_design) {
                        wm_sum_w2y.assign(n_vars, 0.0);
                        wm_sum_w2y2.assign(n_vars, 0.0);
                  }
            }

            int count = 0;
            double sum_weights = 0.0;
            double sum_weights_squared = 0.0;

            double fx_ecef = 0, fy_ecef = 0, fz_ecef = 0;
            if (use_ecef) {
                  lonlat_to_ecef(fx, fy, R_earth, fx_ecef, fy_ecef, fz_ecef);
            }

            // Get inverse covariance matrix for this focal (if using Mahalanobis)
            const std::vector<double>* inv_cov_ptr = nullptr;
            if (use_mahalanobis) {
                  inv_cov_ptr = &inv_cov_matrices[i - begin];
            }

            // REGRESSION: collect analog indices and combined weights for regression solve
            std::vector<std::size_t> reg_analog_indices;
            std::vector<double> reg_combined_weights;
            if (has_regression) {
                  reg_analog_indices.reserve(cand.size());
                  reg_combined_weights.reserve(cand.size());
            }

            // Iterate over candidates once
            {
                  ANALOGS_PROFILE_SCOPE(EXACT);
                  for (size_t t = 0; t < cand.size(); ++t) {
                        const index_t j = cand[t];

                        // Self-exclusion: skip pool row that corresponds to this focal.
                        // Only active when exclude_self is set; the R layer enforces
                        // identical(x, pool) so j == i is the correct identity check.
                        if (exclude_self && static_cast<std::size_t>(j) == i) {
                              continue;
                        }

                        const double sample_weight = cand_weights[t];
                        const double rx = ref_ptr[j];
                        const double ry = ref_ptr[j + stride_r];

                        // Geog distance & filter (annulus: min <= gdist <= max)
                        double gdist = 0.0;
                        if (use_geog_filter || use_geog_min || need_weights) {
                              if (use_ecef) {
                                    const double rx_ecef = ref_latt_ptr[j];
                                    const double ry_ecef = ref_latt_ptr[j + stride_latt_r];
                                    const double rz_ecef = ref_latt_ptr[j + 2 * stride_latt_r];
                                    const double dx = fx_ecef - rx_ecef;
                                    const double dy = fy_ecef - ry_ecef;
                                    const double dz = fz_ecef - rz_ecef;
                                    gdist = std::sqrt(dx*dx + dy*dy + dz*dz);
                                    if (use_geog_filter && gdist > max_geog_chord) continue;
                                    if (use_geog_min && gdist < min_geog_chord) continue;
                              } else {
                                    gdist = geo_distance_km(fx, fy, rx, ry, use_haversine);
                                    if (use_geog_filter && gdist > max_geog) continue;
                                    if (use_geog_min && gdist < min_geog) continue;
                              }
                        }

                        // Environment checks & distance
                        const double* r_env_col = ref_ptr + j + 2 * stride_r;
                        double env_dist = 0.0;
                        bool ok;

                        if (use_mahalanobis) {
                              // Use Mahalanobis distance
                              const bool need_dist = need_weights;
                              auto okd = mahalanobis_ok_and_dist(
                                    f_env_col, r_env_col,
                                    *inv_cov_ptr,
                                    n_env, stride_f, stride_r,
                                    max_env_scalar, need_dist
                              );
                              ok = okd.first;
                              env_dist = okd.second;
                        } else {
                              // Use Euclidean distance
                              auto okd = env_ok_and_dist(
                                    f_env_col, r_env_col,
                                    n_env, stride_f, stride_r,
                                    use_pervar_env, max_env_pervar,
                                    use_scalar_env, max_env_scalar
                              );
                              ok = okd.first;
                              env_dist = okd.second;
                        }

                        if (!ok) continue;

                        // Analog passed all filters - compute distance weight if needed.
                        // Combined weight is the product of the two per-family
                        // kernels (each UNIFORM short-circuits to 1.0).
                        const double dist_weight = need_weights
                        ? weight_from_families(env_kernel, env_dist, env_wparam,
                                               geog_kernel, gdist, geog_wparam)
                              : 1.0;

                        // Per-point pool weights (cell-area and user-supplied).
                        // Each defaults to 1.0 when its mechanism is inactive.
                        const std::size_t j0 = static_cast<std::size_t>(j);
                        const double area_w = get_area_weight(has_area_weight, area_weight_ptr, j0);
                        const double user_w = get_user_weight(has_user_weight, user_weight_ptr, j0);

                        // COUNT preserves existing semantics: it adjusts for downsampling
                        // (sample_weight) but is intentionally NOT affected by area_w or
                        // user_w. Users wanting area- or user-weighted counts should
                        // request stat = "sum_weights" with kernel = "uniform".
                        count += static_cast<int>(sample_weight);

                        // Combined weight folds in all four streams: kernel,
                        // downsampling, cell area, and user-provided.
                        const double combined_weight = dist_weight * sample_weight * area_w * user_w;
                        sum_weights += combined_weight;
                        sum_weights_squared += (combined_weight * combined_weight);

                        // Update regular stat accumulators
                        for (int idx = 0; idx < n_regular; ++idx) {
                              int s = regular_stat_positions[idx];

                              if (acodes[s] == AggregateCode::COUNT) {
                                    // Count tracked separately
                              } else if (acodes[s] == AggregateCode::SUM_WEIGHTS) {
                                    regular_accum[idx] += combined_weight;
                              } else if (acodes[s] == AggregateCode::MEAN_WEIGHTS) {
                                    regular_accum[idx] += combined_weight;
                              } else if (acodes[s] == AggregateCode::ESS) {
                                    // ESS accumulates sum_weights and sum_weights_squared
                                    // Will compute final ESS in finalization step
                              }
                        }

                        // Update value stat accumulators.
                        // SUM and MEAN treat area/user weights as case weights (no
                        // kernel involvement); WEIGHTED_SUM and WEIGHTED_MEAN use
                        // the full combined_weight.
                        const double point_weight = sample_weight * area_w * user_w;
                        if (has_values && n_value > 0) {
                              for (int v = 0; v < n_vars; ++v) {
                                    const double val = values_ptr[j + v * values_stride];

                                    for (int idx = 0; idx < n_value; ++idx) {
                                          int s = value_stat_positions[idx];
                                          int accum_idx = v * n_value + idx;

                                          if (acodes[s] == AggregateCode::SUM) {
                                                value_accum[accum_idx] += (val * point_weight);
                                          } else if (acodes[s] == AggregateCode::MEAN) {
                                                value_accum[accum_idx] += (val * point_weight);
                                          } else if (acodes[s] == AggregateCode::WEIGHTED_SUM) {
                                                value_accum[accum_idx] += (val * combined_weight);
                                          } else if (acodes[s] == AggregateCode::WEIGHTED_MEAN) {
                                                value_accum[accum_idx] += (val * combined_weight);
                                          }
                                    }

                                    // weighted_mean SE accumulators (extra work only
                                    // when weighted_mean SE is requested)
                                    if (wm_se) {
                                          // NA y is skipped for SE too; ISNAN test
                                          if (!ISNAN(val)) {
                                                const double w = combined_weight;
                                                wm_sum_wy2[v] += w * val * val;
                                                if (wm_se_design) {
                                                      const double w2 = w * w;
                                                      wm_sum_w2y[v]  += w2 * val;
                                                      wm_sum_w2y2[v] += w2 * val * val;
                                                }
                                          }
                                    }
                              }
                              if (wm_se_design) {
                                    wm_sum_w2 += combined_weight * combined_weight;
                              }
                        }

                        // TABULATE: add this analog's combined weight to the bin
                        // indexed by its (1-based) class code, for each y variable.
                        // NA codes (ISNAN) are skipped; out-of-range codes are
                        // ignored defensively (R-side validation should prevent them).
                        // Note: y values are passed as doubles holding 1-based integer
                        // class codes (or NA_REAL for missing); int round-trip is exact.
                        if (has_tabulate && has_values && tab_total > 0) {
                              for (int v = 0; v < n_vars; ++v) {
                                    const int K_v = (v < static_cast<int>(n_classes_per_var.size()))
                                    ? n_classes_per_var[v] : 0;
                                    if (K_v <= 0) continue;

                                    const double val = values_ptr[j + v * values_stride];
                                    if (ISNAN(val)) continue;

                                    const int code = static_cast<int>(val) - 1; // 1-based -> 0-based
                                    if (code < 0 || code >= K_v) continue;

                                    tabulate_accum[tab_offsets[v] + code] += combined_weight;
                              }
                        }

                        // REGRESSION: collect this analog for post-loop regression solve
                        if (has_regression) {
                              reg_analog_indices.push_back(static_cast<std::size_t>(j));
                              reg_combined_weights.push_back(combined_weight);
                        }
                  }
            } // end EXACT scope

            // -----------------------------------------------------------------
            // Finalize and store results
            // -----------------------------------------------------------------
            ANALOGS_PROFILE_SCOPE(AGG);
            int col_idx = 0;

            // Write regular stats
            for (int idx = 0; idx < n_regular; ++idx) {
                  int s = regular_stat_positions[idx];
                  double result = NA_REAL;

                  if (acodes[s] == AggregateCode::COUNT) {
                        result = static_cast<double>(count);
                  } else if (acodes[s] == AggregateCode::SUM_WEIGHTS) {
                        result = regular_accum[idx];
                  } else if (acodes[s] == AggregateCode::MEAN_WEIGHTS) {
                        result = (count > 0) ?
                        (regular_accum[idx] / count) : NA_REAL;
                  } else if (acodes[s] == AggregateCode::ESS) {
                        if (sum_weights_squared > 0.0 && sum_weights > 0.0) {
                              result = (sum_weights * sum_weights) / sum_weights_squared;
                        } else {
                              result = 0.0;
                        }
                  }

                  agg[i * n_stats + col_idx] = result;
                  col_idx++;
            }

            // Write value stats.
            // Layout: for each value var v, for each requested value stat s (in
            // order they appear in acodes), write one column. When weighted_mean
            // SE is requested, the SE column is emitted immediately after its
            // weighted_mean column (still within the v-block).
            if (has_values && n_value > 0) {
                  for (int v = 0; v < n_vars; ++v) {
                        for (int idx = 0; idx < n_value; ++idx) {
                              int s = value_stat_positions[idx];
                              int accum_idx = v * n_value + idx;
                              double result = NA_REAL;

                              if (acodes[s] == AggregateCode::SUM) {
                                    result = value_accum[accum_idx];
                              } else if (acodes[s] == AggregateCode::MEAN) {
                                    result = (count > 0) ? (value_accum[accum_idx] / count) : NA_REAL;
                              } else if (acodes[s] == AggregateCode::WEIGHTED_SUM) {
                                    result = value_accum[accum_idx];
                              } else if (acodes[s] == AggregateCode::WEIGHTED_MEAN) {
                                    result = (sum_weights > 0) ? (value_accum[accum_idx] / sum_weights) : NA_REAL;
                              }

                              agg[i * n_stats + col_idx] = result;
                              col_idx++;

                              // Immediately after weighted_mean, emit SE if requested.
                              if (wm_se && idx == wm_value_idx) {
                                    double se_val = NA_REAL;
                                    if (sum_weights > 0.0 && !ISNAN(result)) {
                                          const double ybar = result;  // weighted mean for var v
                                          const double sum_wy2 = wm_sum_wy2[v];

                                          if (se_code == SeCode::ESS) {
                                                // n_eff = (Σw)^2 / Σw^2
                                                // var_w(y) = Σwy^2 / Σw - ybar^2
                                                // SE = sqrt(var_w / n_eff)
                                                if (sum_weights_squared > 0.0) {
                                                      const double n_eff =
                                                            (sum_weights * sum_weights) / sum_weights_squared;
                                                      double var_w =
                                                            (sum_wy2 / sum_weights) - (ybar * ybar);
                                                      if (var_w < 0.0) var_w = 0.0;  // fp guard
                                                      if (n_eff > 0.0 && std::isfinite(var_w)) {
                                                            se_val = std::sqrt(var_w / n_eff);
                                                      }
                                                }
                                          } else {
                                                // DESIGN: SE = sqrt(Σ w^2 (y - ybar)^2) / Σw
                                                //   numerator = Σw^2 y^2 - 2*ybar*Σw^2 y + ybar^2 * Σw^2
                                                const double sum_w2y  = wm_sum_w2y[v];
                                                const double sum_w2y2 = wm_sum_w2y2[v];
                                                double num = sum_w2y2
                                                - 2.0 * ybar * sum_w2y
                                                + ybar * ybar * wm_sum_w2;
                                                if (num < 0.0) num = 0.0;  // fp guard
                                                if (std::isfinite(num)) {
                                                      se_val = std::sqrt(num) / sum_weights;
                                                }
                                          }
                                    }
                                    agg[i * n_stats + col_idx] = se_val;
                                    col_idx++;
                              }
                        }
                  }
            }

            // TABULATE: write per-class weighted vote totals.
            // Layout: for each y variable v (in column order), write K_v
            // contiguous columns, one per class.
            if (has_tabulate && tab_total > 0) {
                  for (int b = 0; b < tab_total; ++b) {
                        agg[i * n_stats + col_idx] = tabulate_accum[b];
                        col_idx++;
                  }
            }

            // REGRESSION: solve and write regression coefficients (and SEs if requested).
            if (has_regression && has_covariates && has_values) {
                  if (reg_analog_indices.empty()) {
                        // Zero analogs: write NAs for all coefficients (and SEs if requested)
                        for (int v = 0; v < n_vars; ++v) {
                              for (int d = 0; d < reg_dim; ++d) {
                                    agg[i * n_stats + col_idx] = NA_REAL;
                                    col_idx++;
                              }
                        }
                        if (reg_se) {
                              for (int v = 0; v < n_vars; ++v) {
                                    for (int d = 0; d < reg_dim; ++d) {
                                          agg[i * n_stats + col_idx] = NA_REAL;
                                          col_idx++;
                                    }
                              }
                        }
                  } else {
                        // Solve ridge regression
                        std::vector<double> coeffs(n_vars * reg_dim);
                        // Allocate SE buffer only if requested, but pass a
                        // harmless pointer either way so solve_ridge can be
                        // uniform. When se_code == NONE, solve_ridge does not
                        // touch out_se.
                        std::vector<double> se_buf;
                        double* se_ptr = nullptr;
                        if (reg_se) {
                              se_buf.assign(n_vars * reg_dim, 0.0);
                              se_ptr = se_buf.data();
                        }

                        solve_ridge(reg_analog_indices,
                                    reg_combined_weights,
                                    values_ptr,
                                    values_stride,
                                    n_vars,
                                    covariates_ptr,
                                    covariates_stride,
                                    n_covs,
                                    lambda,
                                    se_code,
                                    coeffs.data(),
                                    se_ptr);

                        // Write coefficients: for each variable, intercept then covariate slopes
                        for (int v = 0; v < n_vars; ++v) {
                              for (int d = 0; d < reg_dim; ++d) {
                                    agg[i * n_stats + col_idx] = coeffs[v * reg_dim + d];
                                    col_idx++;
                              }
                        }
                        // Write standard errors (same layout) if requested.
                        if (reg_se) {
                              for (int v = 0; v < n_vars; ++v) {
                                    for (int d = 0; d < reg_dim; ++d) {
                                          agg[i * n_stats + col_idx] = se_buf[v * reg_dim + d];
                                          col_idx++;
                                    }
                              }
                        }
                  }
            }
      }
}

} // namespace analogs
