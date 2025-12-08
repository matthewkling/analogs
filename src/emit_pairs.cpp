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
                    Nullable<NumericMatrix> x_cov = R_NilValue,
                    Nullable<NumericMatrix> values = R_NilValue,
                    Nullable<CharacterVector> values_names = R_NilValue) {

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

      // Parse values if provided
      bool has_values = values.isNotNull();
      const double* values_ptr = nullptr;
      int n_vars = 0;
      int values_stride = 0;
      std::vector<std::string> var_names;

      if (has_values) {
            NumericMatrix values_mat = values.get();

            if (values_mat.nrow() != n_ref) {
                  stop("Internal error: values must have same number of rows as reference data");
            }

            n_vars = values_mat.ncol();
            values_ptr = REAL(values_mat);
            values_stride = n_ref;

            // Get variable names if provided
            if (values_names.isNotNull()) {
                  CharacterVector names_vec = values_names.get();
                  for (int i = 0; i < n_vars; ++i) {
                        var_names.push_back(as<std::string>(names_vec[i]));
                  }
            } else {
                  // Generate default names
                  for (int i = 0; i < n_vars; ++i) {
                        var_names.push_back("value_" + std::to_string(i + 1));
                  }
            }
      }

      // First pass: count total number of pairs
      std::size_t total_pairs = 0;
      for (int i = 0; i < n_f; ++i) {
            IntegerVector idx = res[i];
            total_pairs += static_cast<std::size_t>(idx.size());
      }

      // Handle "no matches anywhere" case
      if (total_pairs == 0) {
            // Build empty data.frame with appropriate columns
            List df_cols;
            df_cols["focal_index"] = IntegerVector(0);
            df_cols["focal_x"] = NumericVector(0);
            df_cols["focal_y"] = NumericVector(0);
            df_cols["analog_index"] = IntegerVector(0);
            df_cols["analog_x"] = NumericVector(0);
            df_cols["analog_y"] = NumericVector(0);

            if (report_dist) {
                  df_cols["clim_dist"] = NumericVector(0);
                  df_cols["geog_dist"] = NumericVector(0);
            }

            if (has_values) {
                  for (int v = 0; v < n_vars; ++v) {
                        df_cols[var_names[v]] = NumericVector(0);
                  }
            }

            DataFrame df(df_cols);
            df.attr("stringsAsFactors") = false;
            return df;
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

      // Allocate value vectors if needed
      std::vector<NumericVector> value_cols;
      if (has_values) {
            value_cols.resize(n_vars);
            for (int v = 0; v < n_vars; ++v) {
                  value_cols[v] = NumericVector(total_pairs);
            }
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

                  // Copy values if provided
                  if (has_values) {
                        for (int v = 0; v < n_vars; ++v) {
                              value_cols[v][pos] = values_ptr[ref_row + v * values_stride];
                        }
                  }
            }
      }

      if (pos != total_pairs) {
            stop("Internal error: mismatch between allocated and filled pair counts.");
      }

      // Build output DataFrame
      List df_cols;
      df_cols["focal_index"] = focal_index;
      df_cols["focal_x"] = focal_x;
      df_cols["focal_y"] = focal_y;
      df_cols["analog_index"] = analog_index;
      df_cols["analog_x"] = analog_x;
      df_cols["analog_y"] = analog_y;

      if (report_dist) {
            df_cols["clim_dist"] = clim_dist;
            df_cols["geog_dist"] = geog_dist;
      }

      if (has_values) {
            for (int v = 0; v < n_vars; ++v) {
                  df_cols[var_names[v]] = value_cols[v];
            }
      }

      DataFrame df(df_cols);
      df.attr("stringsAsFactors") = false;
      return df;
}
