#include "workers.hpp"

namespace analogs {

// Define thread-local storage
thread_local PairWorker::ThreadLocalStorage PairWorker::tls;
thread_local AggWorker::ThreadLocalStorage AggWorker::tls;

void PairWorker::operator()(std::size_t begin, std::size_t end) {
      const bool knn_mode      = (mcode == ModeCode::KNN_CLIM || mcode == ModeCode::KNN_GEOG);
      const bool rank_by_clim  = (mcode == ModeCode::KNN_CLIM);
      const bool rank_by_geog  = (mcode == ModeCode::KNN_GEOG);

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
                  const double* cov_vec = x_cov_ptr + i;  // Column-major access

                  // Reconstruct covariance matrix
                  std::vector<double> cov_matrix;
                  reconstruct_cov_matrix(cov_vec, n_clim, cov_matrix);

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
            // -----------------------------------------------------------------
            if (knn_mode &&
                use_lattice &&
                lattice_ptr != nullptr &&
                (!use_haversine || use_ecef)) {

                  const size_tu n_geo = lattice_ptr->n_geo_dims;

                  // Reuse fgeo_vec
                  fgeo_vec.resize(n_geo);

                  if (use_ecef) {
                        // Convert lon/lat (degrees) to ECEF (km)
                        double X, Y, Z;
                        lonlat_to_ecef(fx, fy, R_earth, X, Y, Z);
                        fgeo_vec[0] = X;
                        if (n_geo > 1) fgeo_vec[1] = Y;
                        if (n_geo > 2) fgeo_vec[2] = Z;
                  } else {
                        // Planar projected
                        fgeo_vec[0] = fx;
                        if (n_geo > 1) fgeo_vec[1] = fy;
                  }

                  // Reuse fclim_vec
                  fclim_vec.resize(n_clim);
                  for (int kdim = 0; kdim < n_clim; ++kdim) {
                        fclim_vec[kdim] = f_clim_col[kdim * stride_f];
                  }

                  // Use chord distance threshold for ECEF, regular for planar
                  const double geog_thresh = use_ecef ? max_geog_chord : max_geog;

                  // Reuse cand vector
                  cand.clear();

                  lattice_ptr->knn_query(
                              fgeo_vec.data(),
                              fclim_vec.data(),
                              ref_latt_ptr,
                              static_cast<size_tu>(stride_latt_r),
                              static_cast<size_tu>(n_clim),
                              /*rank_by_geog*/ rank_by_geog,
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
            // Otherwise: use existing candidate generation (lattice query or brute)
            // -----------------------------------------------------------------
            cand.clear();

            if (use_lattice && lattice_ptr != nullptr) {
                  // Query lattice - use ECEF coords if applicable
                  q_geo.resize(3); // max 3 for ECEF
                  if (use_ecef) {
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
                  for (int kdim = 0; kdim < n_clim; ++kdim) {
                        q_clim[kdim] = f_clim_col[kdim * stride_f];
                  }

                  const double geog_thresh = use_ecef ? max_geog_chord : max_geog;

                  lattice_ptr->query(q_geo.data(),
                                     q_clim.data(),
                                     geog_thresh,
                                     use_scalar_clim,
                                     max_clim_pervar,
                                     use_scalar_clim ? max_clim_scalar
                                           : std::numeric_limits<double>::infinity(),
                                             cand);
            } else {
                  cand.reserve(n_ref);
                  for (int j = 0; j < n_ref; ++j) {
                        cand.push_back(static_cast<index_t>(j));
                  }
            }

            // ALL mode: keep all matches passing filters
            if (mcode == ModeCode::ALL) {
                  std::vector<int> keep;
                  keep.reserve(cand.size());

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

                        // Geog filter
                        if (use_geog_filter) {
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

      // Pre-compute inverse covariance matrices if using Mahalanobis
      std::vector< std::vector<double> > inv_cov_matrices;
      if (use_mahalanobis) {
            inv_cov_matrices.resize(end - begin);

            for (std::size_t i = begin; i < end; ++i) {
                  const size_t local_idx = i - begin;
                  const double* cov_vec = x_cov_ptr + i;  // Column-major access

                  // Reconstruct covariance matrix
                  std::vector<double> cov_matrix;
                  reconstruct_cov_matrix(cov_vec, n_clim, cov_matrix);

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
                        if (n_geo > 1) q_geo[1] = fy;
                  }

                  q_clim.resize(n_clim);
                  for (int kdim = 0; kdim < n_clim; ++kdim) {
                        q_clim[kdim] = f_clim_col[kdim * stride_f];
                  }

                  const double geog_thresh = use_ecef ? max_geog_chord : max_geog;

                  lattice_ptr->query(
                              q_geo.data(),
                              q_clim.data(),
                              geog_thresh,
                              use_scalar_clim,
                              max_clim_pervar,
                              use_scalar_clim
                              ? max_clim_scalar
                  : std::numeric_limits<double>::infinity(),
                    cand
                  );
            } else {
                  cand.reserve(n_ref);
                  for (int j = 0; j < n_ref; ++j) {
                        cand.push_back(static_cast<index_t>(j));
                  }
            }

            // Aggregate over candidates
            double sum = 0.0;
            double sum_w = 0.0;
            int    count = 0;

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

                  // Geog filter (in original space: planar, Haversine, or ECEF)
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

                  ++count;

                  if (mcode == ModeCode::SUM || mcode == ModeCode::MEAN) {
                        double w = weight_from_codes(wcode, clim_dist, gdist,
                                                     weight_param1, weight_param2);
                        sum   += w * 1.0;
                        sum_w += w;
                  }

            }

            double val = NA_REAL;
            if (mcode == ModeCode::COUNT) {
                  val = static_cast<double>(count);
            } else if ((mcode == ModeCode::SUM || mcode == ModeCode::MEAN) &&
                  sum_w > 0.0) {
                  double mu = sum / sum_w;
                  if (mcode == ModeCode::SUM) {
                        val = sum;
                  } else { // MEAN
                        val = mu;
                  }
            }

            agg[i] = val;
      }
}

} // namespace analogs
