// src/workers.cpp
#include "workers.hpp"

namespace analogs {

// Define thread-local storage
thread_local PairWorker::ThreadLocalStorage PairWorker::tls;
thread_local AggWorker::ThreadLocalStorage AggWorker::tls;

void PairWorker::operator()(std::size_t begin, std::size_t end) {
      const bool knn_mode      = (scode == SelectCode::KNN_CLIM || scode == SelectCode::KNN_GEOG);
      const bool rank_by_clim  = (scode == SelectCode::KNN_CLIM);
      const bool rank_by_geog  = (scode == SelectCode::KNN_GEOG);

      // Get thread-local storage
      auto& fgeo_vec = tls.fgeo_vec;
      auto& fclim_vec = tls.fclim_vec;
      auto& q_geo = tls.q_geo;
      auto& q_clim = tls.q_clim;
      auto& cand = tls.cand;

      // Pre-compute inverse covariance matrices if using Mahalanobis
      std::vector< std::vector<double> > inv_cov_matrices;
      if (use_mahalanobis) {
            inv_cov_matrices.resize(end - begin);

            for (std::size_t i = begin; i < end; ++i) {
                  const size_t local_idx = i - begin;

                  // Extract covariance components for focal i from column-major x_cov matrix
                  // x_cov is n_focal × n_cov_components stored column-major
                  // Row i, col j is at: x_cov_ptr[i + j * x_cov_stride]
                  std::vector<double> cov_vec(n_cov_components);
                  for (int comp = 0; comp < n_cov_components; ++comp) {
                        cov_vec[comp] = x_cov_ptr[i + comp * x_cov_stride];
                  }

                  // Reconstruct covariance matrix
                  std::vector<double> cov_matrix;
                  reconstruct_cov_matrix(cov_vec.data(), n_clim, cov_matrix);

                  // Invert covariance matrix
                  std::vector<double>& inv_cov = inv_cov_matrices[local_idx];
                  if (!invert_cov_matrix(cov_matrix, n_clim, inv_cov)) {
                        // Matrix not positive definite - use identity (Euclidean)
                        inv_cov.resize(n_clim * n_clim, 0.0);
                        for (int k = 0; k < n_clim; ++k) {
                              inv_cov[k * n_clim + k] = 1.0;
                        }
                  }
            }
      }

      for (std::size_t i = begin; i < end; ++i) {
            const double fx = focal_ptr[i];
            const double fy = focal_ptr[i + stride_f];
            const double* f_clim_col = focal_ptr + i + 2 * stride_f;

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
                  fclim_vec.resize(n_clim);

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

                  // Climate coords
                  for (int d = 0; d < n_clim; ++d) {
                        fclim_vec[d] = f_clim_col[d * stride_f];
                  }

                  // Expanded KNN search (Euclidean)
                  // Reuse cand vector for results
                  cand.clear();

                  const double geog_thresh = use_ecef ? max_geog_chord : max_geog;

                  lattice_ptr->knn_query(
                              fgeo_vec.data(),
                              fclim_vec.data(),
                              ref_ptr,
                              static_cast<size_tu>(stride_r),
                              static_cast<size_tu>(n_clim),
                              rank_by_geog,
                              geog_thresh,
                              use_scalar_clim,
                              max_clim_pervar,
                              max_clim_scalar,
                              k,
                              cand
                  );

                  // Convert 0-based indices to 1-based for R
                  std::vector<int> keep(cand.size());
                  for (std::size_t t = 0; t < cand.size(); ++t) {
                        keep[t] = static_cast<int>(cand[t]) + 1;
                  }
                  out_indices[i] = std::move(keep);
                  continue;
            }

            // -----------------------------------------------------------------
            // Fallback path: gather candidates
            // (Used for: Haversine, Mahalanobis, or ALL mode)
            // -----------------------------------------------------------------
            cand.clear();

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

                  q_clim.resize(n_clim);
                  for (int d = 0; d < n_clim; ++d) {
                        q_clim[d] = f_clim_col[d * stride_f];
                  }

                  // Adjust climate bounds if using Mahalanobis
                  double effective_max_clim = max_clim_scalar;
                  std::vector<double> mahal_bounds;

                  if (use_mahalanobis && std::isfinite(max_clim_scalar)) {
                        const std::vector<double>& inv_cov = inv_cov_matrices[i - begin];
                        mahalanobis_bounding_box(q_clim.data(), inv_cov, n_clim,
                                                 max_clim_scalar, mahal_bounds);
                        // Note: We could use mahal_bounds to adjust lattice query,
                        // but for simplicity we just use a conservative max_clim.
                        // The actual Mahalanobis filtering happens below.
                  }

                  lattice_ptr->query(
                              q_geo.data(),
                              q_clim.data(),
                              max_geog,
                              use_scalar_clim,
                              max_clim_pervar,
                              effective_max_clim,
                              cand
                  );
            } else {
                  cand.resize(n_ref);
                  for (size_t j = 0; j < static_cast<size_t>(n_ref); ++j) {
                        cand[j] = static_cast<index_t>(j);
                  }
            }

            // Filter and collect results from candidates
            if (scode == SelectCode::ALL) {
                  std::vector<int> keep;
                  keep.reserve(cand.size());

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

                        // Geog distance & filter
                        if (use_geog_filter) {
                              double gdist;
                              if (use_ecef) {
                                    const double rx_ecef = ref_latt_ptr[j];
                                    const double ry_ecef = ref_latt_ptr[j + stride_latt_r];
                                    const double rz_ecef = ref_latt_ptr[j + 2 * stride_latt_r];
                                    const double dx = fx_ecef - rx_ecef;
                                    const double dy = fy_ecef - ry_ecef;
                                    const double dz = fz_ecef - rz_ecef;
                                    gdist = std::sqrt(dx*dx + dy*dy + dz*dz);
                                    if (gdist > max_geog_chord) continue;
                              } else {
                                    gdist = geo_distance_km(fx, fy, rx, ry, use_haversine);
                                    if (gdist > max_geog) continue;
                              }
                        }

                        // Climate checks
                        const double* r_clim_col = ref_ptr + j + 2 * stride_r;
                        bool ok;

                        if (use_mahalanobis) {
                              // Use Mahalanobis distance
                              auto okd = mahalanobis_ok_and_dist(
                                    f_clim_col, r_clim_col,
                                    *inv_cov_ptr,
                                    n_clim, stride_f, stride_r,
                                    max_clim_scalar, false
                              );
                              ok = okd.first;
                        } else {
                              // Use Euclidean distance
                              auto okd = clim_ok_and_dist(
                                    f_clim_col, r_clim_col,
                                    n_clim, stride_f, stride_r,
                                    use_pervar_clim, max_clim_pervar,
                                    use_scalar_clim, max_clim_scalar
                              );
                              ok = okd.first;
                        }

                        if (!ok) continue;

                        keep.push_back(j + 1); // 1-based
                  }

                  out_indices[i] = std::move(keep);
                  continue;
            }

            // kNN modes (Climate or Geog) without expanding-search lattice path
            auto cmp = [](const Neighbor& a, const Neighbor& b) {
                  return a.first < b.first; // max-heap: top has largest distance
            };
            std::priority_queue<Neighbor, std::vector<Neighbor>, decltype(cmp)> pq(cmp);

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

                  // Geog distance & filter
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
                  } else {
                        gdist = geo_distance_km(fx, fy, rx, ry, use_haversine);
                        if (use_geog_filter && gdist > max_geog) continue;
                  }

                  // Climate checks & distance
                  const double* r_clim_col = ref_ptr + j + 2 * stride_r;
                  double clim_dist;
                  bool ok;

                  if (use_mahalanobis) {
                        // Use Mahalanobis distance
                        auto okd = mahalanobis_ok_and_dist(
                              f_clim_col, r_clim_col,
                              *inv_cov_ptr,
                              n_clim, stride_f, stride_r,
                              max_clim_scalar, true
                        );
                        ok = okd.first;
                        clim_dist = okd.second;
                  } else {
                        // Use Euclidean distance
                        auto okd = clim_ok_and_dist(
                              f_clim_col, r_clim_col,
                              n_clim, stride_f, stride_r,
                              use_pervar_clim, max_clim_pervar,
                              use_scalar_clim, max_clim_scalar
                        );
                        ok = okd.first;
                        clim_dist = okd.second;
                  }

                  if (!ok) continue;

                  const double key = rank_by_clim ? clim_dist : gdist;
                  const index_t ref_index_1based = j + 1;

                  if (static_cast<int>(pq.size()) < k) {
                        pq.emplace(key, ref_index_1based);
                  } else if (!pq.empty() && key < pq.top().first) {
                        pq.pop();
                        pq.emplace(key, ref_index_1based);
                  }
            }

            const int m = static_cast<int>(pq.size());
            std::vector<int> idx_vec(m);
            for (int pos = m - 1; pos >= 0; --pos) {
                  const Neighbor& nb = pq.top();
                  idx_vec[pos] = nb.second;
                  pq.pop();
            }
            out_indices[i] = std::move(idx_vec);
      }
}


void AggWorker::operator()(std::size_t begin, std::size_t end) {
      // Get thread-local storage
      auto& q_geo = tls.q_geo;
      auto& q_clim = tls.q_clim;
      auto& cand = tls.cand;

      // Determine which stats are regular vs value-based
      std::vector<int> regular_stat_positions;
      std::vector<int> value_stat_positions;

      for (int s = 0; s < static_cast<int>(acodes.size()); ++s) {
            if (acodes[s] == AggregateCode::SUM ||
                acodes[s] == AggregateCode::MEAN ||
                acodes[s] == AggregateCode::WEIGHTED_SUM ||
                acodes[s] == AggregateCode::WEIGHTED_MEAN) {
                  value_stat_positions.push_back(s);
            } else {
                  regular_stat_positions.push_back(s);
            }
      }

      const int n_regular = static_cast<int>(regular_stat_positions.size());
      const int n_value = static_cast<int>(value_stat_positions.size());

      // Check if any stats need weights
      bool need_weights = false;
      for (size_t s = 0; s < acodes.size(); ++s) {
            if (acodes[s] == AggregateCode::SUM_WEIGHTS ||
                acodes[s] == AggregateCode::MEAN_WEIGHTS ||
                acodes[s] == AggregateCode::WEIGHTED_SUM ||
                acodes[s] == AggregateCode::WEIGHTED_MEAN) {
                  need_weights = true;
                  break;
            }
      }

      // Pre-compute inverse covariance matrices if using Mahalanobis
      std::vector< std::vector<double> > inv_cov_matrices;
      if (use_mahalanobis) {
            inv_cov_matrices.resize(end - begin);

            for (std::size_t i = begin; i < end; ++i) {
                  const size_t local_idx = i - begin;

                  // Extract covariance components for focal i from column-major x_cov matrix
                  // x_cov is n_focal × n_cov_components stored column-major
                  // Row i, col j is at: x_cov_ptr[i + j * x_cov_stride]
                  std::vector<double> cov_vec(n_cov_components);
                  for (int comp = 0; comp < n_cov_components; ++comp) {
                        cov_vec[comp] = x_cov_ptr[i + comp * x_cov_stride];
                  }

                  // Reconstruct covariance matrix
                  std::vector<double> cov_matrix;
                  reconstruct_cov_matrix(cov_vec.data(), n_clim, cov_matrix);

                  // Invert covariance matrix
                  std::vector<double>& inv_cov = inv_cov_matrices[local_idx];
                  if (!invert_cov_matrix(cov_matrix, n_clim, inv_cov)) {
                        // Matrix not positive definite - use identity (Euclidean)
                        inv_cov.resize(n_clim * n_clim, 0.0);
                        for (int k = 0; k < n_clim; ++k) {
                              inv_cov[k * n_clim + k] = 1.0;
                        }
                  }
            }
      }

      for (std::size_t i = begin; i < end; ++i) {
            const double fx = focal_ptr[i];
            const double fy = focal_ptr[i + stride_f];
            const double* f_clim_col = focal_ptr + i + 2 * stride_f;

            cand.clear();

            if (use_lattice && lattice_ptr != nullptr) {
                  const size_tu n_geo = lattice_ptr->n_geo_dims;
                  q_geo.resize(n_geo);

                  if (use_ecef && n_geo == 3) {
                        // Lattice was built in ECEF; query in ECEF.
                        double X, Y, Z;
                        lonlat_to_ecef(fx, fy, R_earth, X, Y, Z);
                        q_geo[0] = X;
                        q_geo[1] = Y;
                        q_geo[2] = Z;
                  } else {
                        // Planar / non-ECEF lattice: use raw x,y.
                        q_geo[0] = fx;
                        q_geo[1] = fy;
                  }

                  q_clim.resize(n_clim);
                  for (int d = 0; d < n_clim; ++d) {
                        q_clim[d] = f_clim_col[d * stride_f];
                  }

                  // Adjust climate bounds if using Mahalanobis
                  double effective_max_clim = max_clim_scalar;
                  std::vector<double> mahal_bounds;

                  if (use_mahalanobis && std::isfinite(max_clim_scalar)) {
                        const std::vector<double>& inv_cov = inv_cov_matrices[i - begin];
                        mahalanobis_bounding_box(q_clim.data(), inv_cov, n_clim,
                                                 max_clim_scalar, mahal_bounds);
                  }

                  lattice_ptr->query(
                              q_geo.data(),
                              q_clim.data(),
                              max_geog,
                              use_scalar_clim,
                              max_clim_pervar,
                              effective_max_clim,
                              cand
                  );
            } else {
                  cand.resize(n_ref);
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

            int count = 0;
            double sum_weights = 0.0;

            double fx_ecef = 0, fy_ecef = 0, fz_ecef = 0;
            if (use_ecef) {
                  lonlat_to_ecef(fx, fy, R_earth, fx_ecef, fy_ecef, fz_ecef);
            }

            // Get inverse covariance matrix for this focal (if using Mahalanobis)
            const std::vector<double>* inv_cov_ptr = nullptr;
            if (use_mahalanobis) {
                  inv_cov_ptr = &inv_cov_matrices[i - begin];
            }

            // Iterate over candidates once
            for (size_t t = 0; t < cand.size(); ++t) {
                  const index_t j = cand[t];
                  const double rx = ref_ptr[j];
                  const double ry = ref_ptr[j + stride_r];

                  // Geog distance & filter
                  double gdist = 0.0;
                  if (use_geog_filter || need_weights) {
                        if (use_ecef) {
                              const double rx_ecef = ref_latt_ptr[j];
                              const double ry_ecef = ref_latt_ptr[j + stride_latt_r];
                              const double rz_ecef = ref_latt_ptr[j + 2 * stride_latt_r];
                              const double dx = fx_ecef - rx_ecef;
                              const double dy = fy_ecef - ry_ecef;
                              const double dz = fz_ecef - rz_ecef;
                              gdist = std::sqrt(dx*dx + dy*dy + dz*dz);
                              if (use_geog_filter && gdist > max_geog_chord) continue;
                        } else {
                              gdist = geo_distance_km(fx, fy, rx, ry, use_haversine);
                              if (use_geog_filter && gdist > max_geog) continue;
                        }
                  }

                  // Climate checks & distance
                  const double* r_clim_col = ref_ptr + j + 2 * stride_r;
                  double clim_dist = 0.0;
                  bool ok;

                  if (use_mahalanobis) {
                        // Use Mahalanobis distance
                        const bool need_dist = need_weights;
                        auto okd = mahalanobis_ok_and_dist(
                              f_clim_col, r_clim_col,
                              *inv_cov_ptr,
                              n_clim, stride_f, stride_r,
                              max_clim_scalar, need_dist
                        );
                        ok = okd.first;
                        clim_dist = okd.second;
                  } else {
                        // Use Euclidean distance
                        auto okd = clim_ok_and_dist(
                              f_clim_col, r_clim_col,
                              n_clim, stride_f, stride_r,
                              use_pervar_clim, max_clim_pervar,
                              use_scalar_clim, max_clim_scalar
                        );
                        ok = okd.first;
                        clim_dist = okd.second;
                  }

                  if (!ok) continue;

                  // Analog passed all filters - compute weight if needed
                  const double w = need_weights
                  ? weight_from_codes(wcode, clim_dist, gdist, weight_param1, weight_param2)
                        : 1.0;

                  count++;
                  sum_weights += w;

                  // Update regular stat accumulators
                  for (int idx = 0; idx < n_regular; ++idx) {
                        int s = regular_stat_positions[idx];

                        if (acodes[s] == AggregateCode::COUNT) {
                              // Count tracked separately
                        } else if (acodes[s] == AggregateCode::SUM_WEIGHTS) {
                              regular_accum[idx] += w;
                        } else if (acodes[s] == AggregateCode::MEAN_WEIGHTS) {
                              regular_accum[idx] += w;
                        }
                  }

                  // Update value stat accumulators
                  if (has_values && n_value > 0) {
                        for (int v = 0; v < n_vars; ++v) {
                              const double val = values_ptr[j + v * values_stride];

                              for (int idx = 0; idx < n_value; ++idx) {
                                    int s = value_stat_positions[idx];
                                    int accum_idx = v * n_value + idx;

                                    if (acodes[s] == AggregateCode::SUM) {
                                          value_accum[accum_idx] += val;
                                    } else if (acodes[s] == AggregateCode::MEAN) {
                                          value_accum[accum_idx] += val;
                                    } else if (acodes[s] == AggregateCode::WEIGHTED_SUM) {
                                          value_accum[accum_idx] += val * w;
                                    } else if (acodes[s] == AggregateCode::WEIGHTED_MEAN) {
                                          value_accum[accum_idx] += val * w;
                                    }
                              }
                        }
                  }
            }

            // Finalize and store results
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
                        result = (count > 0) ? (regular_accum[idx] / count) : NA_REAL;
                  }

                  agg[i * n_stats + col_idx] = result;
                  col_idx++;
            }

            // Write value stats (grouped by variable)
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
                        }
                  }
            }
      }
}

} // namespace analogs
