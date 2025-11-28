// src/emit_pairs.cpp
#include <Rcpp.h>
#include <cmath>
#include "mahalanobis.hpp"

using namespace Rcpp;

inline double haversine_km(double lon1, double lat1,
                           double lon2, double lat2) {
      const double R = 6371.0088; // Earth's mean radius in km
      const double to_rad = M_PI / 180.0;

      lon1 *= to_rad;
      lat1 *= to_rad;
      lon2 *= to_rad;
      lat2 *= to_rad;

      const double dlon = lon2 - lon1;
      const double dlat = lat2 - lat1;

      const double sdlat = std::sin(0.5 * dlat);
      const double sdlon = std::sin(0.5 * dlon);

      const double a = sdlat * sdlat +
            std::cos(lat1) * std::cos(lat2) * sdlon * sdlon;

      const double c = 2.0 * std::asin(std::min(1.0, std::sqrt(a)));
      return R * c;
}

// [[Rcpp::export(.emit_pairs_cpp)]]
SEXP emit_pairs_cpp(List res,
                    NumericMatrix focal_mm,
                    NumericMatrix ref_mm,
                    bool report_dist,
                    std::string geo_mode,
                    Nullable<NumericMatrix> x_cov = R_NilValue) {

      const int n_f = focal_mm.nrow();
      const int n_ref = ref_mm.nrow();

      if (res.size() != n_f) {
            stop("Internal error: length(res) != nrow(focal_mm).");
      }

      const int ncol_focal = focal_mm.ncol();
      const int ncol_ref   = ref_mm.ncol();

      if (ncol_focal < 3 || ncol_ref < 3) {
            stop("Internal error: focal_mm and ref_mm must have >= 3 columns (x, y, climate...).");
      }
      if (ncol_focal != ncol_ref) {
            stop("Internal error: focal_mm and ref_mm must have the same number of columns.");
      }

      // First pass: count total number of pairs
      std::size_t total_pairs = 0;
      for (int i = 0; i < n_f; ++i) {
            IntegerVector idx = res[i];
            total_pairs += static_cast<std::size_t>(idx.size());
      }

      // Handle "no matches anywhere" case
      if (total_pairs == 0) {
            if (report_dist) {
                  return DataFrame::create(
                        _["focal_index"] = IntegerVector(0),
                        _["focal_x"]     = NumericVector(0),
                        _["focal_y"]     = NumericVector(0),
                        _["analog_index"]= IntegerVector(0),
                        _["analog_x"]    = NumericVector(0),
                        _["analog_y"]    = NumericVector(0),
                        _["clim_dist"]   = NumericVector(0),
                        _["geog_dist"]   = NumericVector(0),
                        _["stringsAsFactors"] = false
                  );
            } else {
                  return DataFrame::create(
                        _["focal_index"] = IntegerVector(0),
                        _["focal_x"]     = NumericVector(0),
                        _["focal_y"]     = NumericVector(0),
                        _["analog_index"]= IntegerVector(0),
                        _["analog_x"]    = NumericVector(0),
                        _["analog_y"]    = NumericVector(0),
                        _["stringsAsFactors"] = false
                  );
            }
      }

      // Allocate output vectors
      IntegerVector focal_index(total_pairs);
      NumericVector focal_x(total_pairs), focal_y(total_pairs);
      IntegerVector analog_index(total_pairs);
      NumericVector analog_x(total_pairs), analog_y(total_pairs);
      NumericVector clim_dist, geog_dist;

      if (report_dist) {
            clim_dist = NumericVector(total_pairs);
            geog_dist = NumericVector(total_pairs);
      }

      const bool use_lonlat = (geo_mode == "lonlat");
      const int clim_start_col = 2; // 0-based index: col 2 == third column
      const int n_clim = ncol_focal - clim_start_col;

      // Parse x_cov if provided for Mahalanobis distance calculation
      bool use_mahalanobis = false;
      const double* x_cov_ptr = nullptr;
      int x_cov_stride = 0;
      int n_cov_components = 0;

      if (x_cov.isNotNull()) {
            NumericMatrix x_cov_mat = x_cov.get();

            if (x_cov_mat.nrow() != n_f) {
                  stop("Internal error: x_cov must have same number of rows as focal data");
            }

            n_cov_components = n_clim * (n_clim + 1) / 2;
            if (x_cov_mat.ncol() != n_cov_components) {
                  stop("Internal error: x_cov must have n_clim * (n_clim + 1) / 2 columns");
            }

            use_mahalanobis = true;
            x_cov_ptr = REAL(x_cov_mat);
            x_cov_stride = n_f;  // Column-major stride
      }

      std::size_t pos = 0;

      for (int i = 0; i < n_f; ++i) {
            IntegerVector idx = res[i];
            const int m = idx.size();
            if (m == 0)
                  continue;

            // Focal coords (constant for this i)
            const double fx = focal_mm(i, 0);
            const double fy = focal_mm(i, 1);

            // Cache focal climate row
            std::vector<double> f_clim(n_clim);
            for (int k = 0; k < n_clim; ++k) {
                  f_clim[k] = focal_mm(i, clim_start_col + k);
            }

            // Pre-compute inverse covariance matrix for this focal if using Mahalanobis
            std::vector<double> inv_cov;
            if (use_mahalanobis && report_dist) {
                  // Extract covariance components for focal i from column-major x_cov matrix
                  // Row i, col j is at: x_cov_ptr[i + j * x_cov_stride]
                  std::vector<double> cov_vec(n_cov_components);
                  for (int comp = 0; comp < n_cov_components; ++comp) {
                        cov_vec[comp] = x_cov_ptr[i + comp * x_cov_stride];
                  }

                  // Reconstruct covariance matrix
                  std::vector<double> cov_matrix;
                  analogs::reconstruct_cov_matrix(cov_vec.data(), n_clim, cov_matrix);

                  // Invert covariance matrix
                  if (!analogs::invert_cov_matrix(cov_matrix, n_clim, inv_cov)) {
                        // Matrix not positive definite - use identity (Euclidean)
                        inv_cov.resize(n_clim * n_clim, 0.0);
                        for (int k = 0; k < n_clim; ++k) {
                              inv_cov[k * n_clim + k] = 1.0;
                        }
                  }
            }

            for (int j = 0; j < m; ++j, ++pos) {
                  int ref_idx = idx[j];

                  if (ref_idx < 1 || ref_idx > n_ref) {
                        stop("Internal error: analog index out of bounds.");
                  }
                  const int ref_row = ref_idx - 1; // 0-based

                  const double ax = ref_mm(ref_row, 0);
                  const double ay = ref_mm(ref_row, 1);

                  focal_index[pos]  = i + 1;   // 1-based
                  focal_x[pos]      = fx;
                  focal_y[pos]      = fy;
                  analog_index[pos] = ref_idx; // keep 1-based index
                  analog_x[pos]     = ax;
                  analog_y[pos]     = ay;

                  if (report_dist) {
                        // Climate distance
                        if (use_mahalanobis) {
                              // Compute Mahalanobis distance
                              std::vector<double> diff(n_clim);
                              for (int k = 0; k < n_clim; ++k) {
                                    diff[k] = ref_mm(ref_row, clim_start_col + k) - f_clim[k];
                              }

                              // d^T * inv_cov * d
                              double result = 0.0;
                              for (int ii = 0; ii < n_clim; ++ii) {
                                    double sum = 0.0;
                                    for (int jj = 0; jj < n_clim; ++jj) {
                                          sum += inv_cov[ii * n_clim + jj] * diff[jj];
                                    }
                                    result += diff[ii] * sum;
                              }
                              clim_dist[pos] = std::sqrt(std::max(0.0, result));
                        } else {
                              // Compute Euclidean distance in climate space
                              double sum_sq = 0.0;
                              for (int k = 0; k < n_clim; ++k) {
                                    const double diff = ref_mm(ref_row, clim_start_col + k) - f_clim[k];
                                    sum_sq += diff * diff;
                              }
                              clim_dist[pos] = std::sqrt(sum_sq);
                        }

                        // Geographic distance
                        if (use_lonlat) {
                              geog_dist[pos] = haversine_km(fx, fy, ax, ay);
                        } else {
                              const double dx = ax - fx;
                              const double dy = ay - fy;
                              geog_dist[pos]  = std::sqrt(dx * dx + dy * dy);
                        }
                  }
            }
      }

      if (pos != total_pairs) {
            stop("Internal error: mismatch between allocated and filled pair counts.");
      }

      if (report_dist) {
            return DataFrame::create(
                  _["focal_index"] = focal_index,
                  _["focal_x"]     = focal_x,
                  _["focal_y"]     = focal_y,
                  _["analog_index"]= analog_index,
                  _["analog_x"]    = analog_x,
                  _["analog_y"]    = analog_y,
                  _["clim_dist"]   = clim_dist,
                  _["geog_dist"]   = geog_dist,
                  _["stringsAsFactors"] = false
            );
      } else {
            return DataFrame::create(
                  _["focal_index"] = focal_index,
                  _["focal_x"]     = focal_x,
                  _["focal_y"]     = focal_y,
                  _["analog_index"]= analog_index,
                  _["analog_x"]    = analog_x,
                  _["analog_y"]    = analog_y,
                  _["stringsAsFactors"] = false
            );
      }
}
