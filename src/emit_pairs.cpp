// src/emit_pairs.cpp
#include <Rcpp.h>
#include <cmath>
#include "mahalanobis.hpp"

using namespace Rcpp;
using namespace analogs;

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

      const int n_clim = ncol_focal - 2;

      // Check if we have x_cov for Mahalanobis distance
      bool use_mahalanobis = false;
      NumericMatrix x_cov_mat;
      int n_cov_components = 0;
      if (x_cov.isNotNull()) {
            use_mahalanobis = true;
            x_cov_mat = x_cov.get();
            n_cov_components = x_cov_mat.ncol();
            if (x_cov_mat.nrow() != n_f) {
                  stop("x_cov must have same number of rows as focal_mm");
            }
      }

      // Check if we have values
      bool has_values = false;
      NumericMatrix values_mat;
      int n_vars = 0;
      CharacterVector var_names;
      if (values.isNotNull()) {
            has_values = true;
            values_mat = values.get();
            n_vars = values_mat.ncol();
            if (values_mat.nrow() != n_ref) {
                  stop("values must have same number of rows as ref_mm");
            }
            if (values_names.isNotNull()) {
                  var_names = values_names.get();
                  if (var_names.size() != n_vars) {
                        stop("values_names length must match number of value columns");
                  }
            }
      }

      // Count total pairs
      int total_pairs = 0;
      for (int i = 0; i < n_f; ++i) {
            IntegerVector analog_idx = res[i];
            total_pairs += analog_idx.size();
      }

      // Allocate output vectors
      IntegerVector focal_index(total_pairs);
      NumericVector focal_x(total_pairs);
      NumericVector focal_y(total_pairs);
      IntegerVector analog_index(total_pairs);
      NumericVector analog_x(total_pairs);
      NumericVector analog_y(total_pairs);
      NumericVector clim_dist(total_pairs);
      NumericVector geog_dist(total_pairs);

      // Value columns if needed
      std::vector<NumericVector> value_cols;
      if (has_values) {
            value_cols.resize(n_vars);
            for (int v = 0; v < n_vars; ++v) {
                  value_cols[v] = NumericVector(total_pairs);
            }
      }

      // Fill output vectors
      int out_idx = 0;
      for (int i = 0; i < n_f; ++i) {
            IntegerVector analog_idx_vec = res[i];
            const int n_analogs = analog_idx_vec.size();

            const double fx = focal_mm(i, 0);
            const double fy = focal_mm(i, 1);

            for (int j = 0; j < n_analogs; ++j) {
                  int ref_idx = analog_idx_vec[j];

                  focal_index[out_idx] = i + 1;  // 1-based for R
                  focal_x[out_idx] = fx;
                  focal_y[out_idx] = fy;

                  // Handle NA case (ref_idx == 0, since we use 1-based indexing from C++)
                  if (ref_idx < 1) {
                        analog_index[out_idx] = NA_INTEGER;
                        analog_x[out_idx] = NA_REAL;
                        analog_y[out_idx] = NA_REAL;
                        clim_dist[out_idx] = NA_REAL;
                        geog_dist[out_idx] = NA_REAL;

                        // NA for values too
                        if (has_values) {
                              for (int v = 0; v < n_vars; ++v) {
                                    value_cols[v][out_idx] = NA_REAL;
                              }
                        }
                  } else if (ref_idx > n_ref) {
                        stop("Internal error: analog index out of bounds");
                  } else {
                        // Valid analog - compute distances
                        // Convert from 1-based (R) to 0-based (C++) indexing
                        int ref_row = ref_idx - 1;
                        analog_index[out_idx] = ref_idx;  // Keep 1-based for R
                        analog_x[out_idx] = ref_mm(ref_row, 0);
                        analog_y[out_idx] = ref_mm(ref_row, 1);

                        // Compute climate distance
                        if (use_mahalanobis) {
                              // Mahalanobis distance calculation
                              // Extract covariance components for this focal
                              std::vector<double> cov_vec(n_cov_components);
                              for (int c = 0; c < n_cov_components; ++c) {
                                    cov_vec[c] = x_cov_mat(i, c);
                              }

                              // Reconstruct covariance matrix using helper function
                              std::vector<double> cov_mat;
                              reconstruct_cov_matrix(cov_vec.data(), n_clim, cov_mat);

                              // Invert covariance matrix
                              std::vector<double> inv_cov;
                              bool success = invert_cov_matrix(cov_mat, n_clim, inv_cov);

                              if (!success) {
                                    // Matrix not positive definite - use Euclidean distance
                                    double sum_sq = 0.0;
                                    for (int d = 0; d < n_clim; ++d) {
                                          double diff = focal_mm(i, 2 + d) - ref_mm(ref_row, 2 + d);
                                          sum_sq += diff * diff;
                                    }
                                    clim_dist[out_idx] = std::sqrt(sum_sq);
                              } else {
                                    // Compute Mahalanobis distance
                                    // Note: mahalanobis_distance expects pointers with stride access
                                    const double* f_clim = &focal_mm(i, 2);
                                    const double* r_clim = &ref_mm(ref_row, 2);

                                    // stride is number of rows (for column-major access)
                                    int stride_f = focal_mm.nrow();
                                    int stride_r = ref_mm.nrow();

                                    clim_dist[out_idx] = mahalanobis_distance(
                                          f_clim, r_clim, inv_cov, n_clim, stride_f, stride_r
                                    );
                              }
                        } else {
                              // Euclidean distance
                              double sum_sq = 0.0;
                              for (int d = 0; d < n_clim; ++d) {
                                    double diff = focal_mm(i, 2 + d) - ref_mm(ref_row, 2 + d);
                                    sum_sq += diff * diff;
                              }
                              clim_dist[out_idx] = std::sqrt(sum_sq);
                        }

                        // Compute geographic distance
                        if (report_dist) {
                              if (geo_mode == "lonlat") {
                                    geog_dist[out_idx] = haversine_km(
                                          fx, fy,
                                          analog_x[out_idx], analog_y[out_idx]
                                    );
                              } else {
                                    double dx = analog_x[out_idx] - fx;
                                    double dy = analog_y[out_idx] - fy;
                                    geog_dist[out_idx] = std::sqrt(dx * dx + dy * dy);
                              }
                        } else {
                              geog_dist[out_idx] = NA_REAL;
                        }

                        // Copy values if present
                        if (has_values) {
                              for (int v = 0; v < n_vars; ++v) {
                                    value_cols[v][out_idx] = values_mat(ref_row, v);
                              }
                        }
                  }

                  out_idx++;
            }
      }

      // Build output data.frame
      List out_list;
      out_list["focal_index"] = focal_index;
      out_list["focal_x"] = focal_x;
      out_list["focal_y"] = focal_y;
      out_list["analog_index"] = analog_index;
      out_list["analog_x"] = analog_x;
      out_list["analog_y"] = analog_y;
      out_list["clim_dist"] = clim_dist;
      out_list["geog_dist"] = geog_dist;

      // Add value columns
      if (has_values) {
            for (int v = 0; v < n_vars; ++v) {
                  std::string col_name = var_names.size() > 0 ?
                  Rcpp::as<std::string>(var_names[v]) :
                  ("value_" + std::to_string(v + 1));
                  out_list[col_name] = value_cols[v];
            }
      }

      DataFrame df(out_list);
      df.attr("class") = CharacterVector::create("data.frame");

      return df;
}
