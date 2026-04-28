// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>

#include "lattice.h"
#include "metrics.h"
#include "types.h"
// #include "profiling.h"
#include "geometry.h"
#include "climate.h"
#include "weights.h"
#include "workers.h"
#include "mahalanobis.h"
#include "se_code.h"

#include <vector>
#include <limits>

using namespace Rcpp;
using namespace analogs;
using namespace RcppParallel;

// namespace analogs { thread_local ProfileTimer g_profiler; }


// =========================================================================
// Build analog index (lattice) and return as Xptr
// =========================================================================

// [[Rcpp::export]]
SEXP build_analog_index_cpp(const NumericMatrix& ref_mm,
                            const std::string& coord_type,
                            int index_res,
                            double downsample,
                            unsigned int seed)
{
      const int n_ref = ref_mm.nrow();
      const int ncol_ref = ref_mm.ncol();

      if (ncol_ref < 3) {
            stop("Reference data must have at least 2 coordinate columns and 1 climate variable");
      }

      const int n_clim = ncol_ref - 2;
      const bool use_haversine = (coord_type == "lonlat");
      const bool use_ecef = use_haversine; // Use ECEF for lonlat coords
      const double R_earth = 6371.0088;

      // Create lattice object
      Lattice* lattice = new Lattice();

      // Storage for coordinate ranges (for metadata)
      std::vector<double> coord_mins(2), coord_maxs(2);
      std::vector<double> clim_mins(n_clim), clim_maxs(n_clim);

      // Compute ranges
      const double* ref_ptr = REAL(ref_mm);
      const int stride = n_ref;

      for (int i = 0; i < n_ref; ++i) {
            // Coordinate ranges
            for (int d = 0; d < 2; ++d) {
                  double v = ref_ptr[i + d * stride];
                  if (i == 0) {
                        coord_mins[d] = coord_maxs[d] = v;
                  } else {
                        if (v < coord_mins[d]) coord_mins[d] = v;
                        if (v > coord_maxs[d]) coord_maxs[d] = v;
                  }
            }

            // Climate ranges
            for (int k = 0; k < n_clim; ++k) {
                  double v = ref_ptr[i + (2 + k) * stride];
                  if (i == 0) {
                        clim_mins[k] = clim_maxs[k] = v;
                  } else {
                        if (v < clim_mins[k]) clim_mins[k] = v;
                        if (v > clim_maxs[k]) clim_maxs[k] = v;
                  }
            }
      }

      // Prepare data for lattice building
      const double* lattice_data_ptr;
      std::vector<double> ecef_storage;
      int lattice_stride;
      size_t n_geo_dims;

      if (use_ecef) {
            // Build ECEF (X,Y,Z) + climate lattice
            n_geo_dims = 3;
            const size_t stride_e = static_cast<size_t>(n_ref);
            const size_t n_cols = n_geo_dims + static_cast<size_t>(n_clim);

            ecef_storage.assign(n_cols * n_ref, 0.0);

            for (size_t j = 0; j < static_cast<size_t>(n_ref); ++j) {
                  const double lon = ref_ptr[j];
                  const double lat = ref_ptr[j + stride_e];

                  double X, Y, Z;
                  lonlat_to_ecef(lon, lat, R_earth, X, Y, Z);

                  ecef_storage[j] = X;
                  ecef_storage[j + stride_e] = Y;
                  ecef_storage[j + 2 * stride_e] = Z;
            }

            // Copy climate variables
            for (size_t k = 0; k < static_cast<size_t>(n_clim); ++k) {
                  const double* src_col = ref_ptr + (2 + k) * stride_e;
                  double* dst_col = ecef_storage.data() + (3 + k) * stride_e;
                  std::copy(src_col, src_col + n_ref, dst_col);
            }

            lattice_data_ptr = ecef_storage.data();
            lattice_stride = n_ref;
      } else {
            // Use original projected coordinates
            n_geo_dims = 2;
            lattice_data_ptr = ref_ptr;
            lattice_stride = n_ref;
      }

      // Build lattice (always index both geo and climate dimensions)
      MetricType metric = use_ecef ? MetricType::Chord3D : MetricType::Planar;

      lattice->build(
                  lattice_data_ptr,
                  static_cast<size_tu>(n_ref),
                  static_cast<size_tu>(n_geo_dims),
                  static_cast<size_tu>(n_clim),
                  static_cast<size_tu>(lattice_stride),
                  metric,
                  std::numeric_limits<double>::infinity(),
                  false,
                  std::vector<double>(),
                  std::numeric_limits<double>::infinity(),
                  index_res,
                  true,       // use_geo_lattice
                  true,       // use_clim_lattice
                  downsample, // NEW
                  seed        // NEW
      );

      // Create Xptr
      XPtr<Lattice> lattice_xptr(lattice, true);

      // Store ECEF data if needed (for queries)
      SEXP ecef_xptr = R_NilValue;
      if (use_ecef) {
            std::vector<double>* ecef_copy = new std::vector<double>(ecef_storage);
            ecef_xptr = XPtr< std::vector<double> >(ecef_copy, true);
      }

      // Return list with lattice Xptr and metadata
      List result = List::create(
            Named("lattice_xptr") = lattice_xptr,
            Named("ecef_xptr") = ecef_xptr,
            Named("coord_type") = coord_type,
            Named("n_ref") = n_ref,
            Named("n_clim") = n_clim,
            Named("index_res") = index_res,
            Named("coord_mins") = coord_mins,
            Named("coord_maxs") = coord_maxs,
            Named("clim_mins") = clim_mins,
            Named("clim_maxs") = clim_maxs,
            Named("use_ecef") = use_ecef,
            Named("total_bins") = static_cast<double>(lattice->total_bins),
            Named("n_cells_nonempty") = static_cast<double>(lattice->n_cells_nonempty),
            Named("min_cell_occ") = static_cast<double>(lattice->min_cell_occ),
            Named("max_cell_occ") = static_cast<double>(lattice->max_cell_occ),
            Named("downsample_actual") = lattice->downsample_actual  // NEW
      );

      return result;
}


// =========================================================================
// Query analog index with focal points
// =========================================================================

// [[Rcpp::export]]
SEXP query_analog_index_cpp(SEXP index_list,
                            const NumericMatrix& focal_mm,
                            const NumericMatrix& ref_mm,
                            int k,
                            const NumericVector& max_clim,
                            double max_geog,
                            int select_code,
                            const IntegerVector& aggregate_codes,
                            int weight_code,
                            const NumericVector& theta,
                            SEXP x_cov_sexp,
                            SEXP values_sexp,
                            SEXP covariates_sexp,
                            double lambda,
                            int se_code,
                            const IntegerVector& n_classes_per_var,
                            bool exclude_self = false)
{
      // Extract lattice and metadata from index
      List idx = as<List>(index_list);

      XPtr<Lattice> lattice_xptr = idx["lattice_xptr"];
      Lattice* lattice_ptr = lattice_xptr.get();

      if (lattice_ptr == nullptr) {
            stop("Invalid lattice pointer in index");
      }

      std::string coord_type = as<std::string>(idx["coord_type"]);
      int n_ref = as<int>(idx["n_ref"]);
      int n_clim = as<int>(idx["n_clim"]);
      bool use_ecef = as<bool>(idx["use_ecef"]);

      // Validate dimensions
      if (ref_mm.nrow() != n_ref) {
            stop("Reference data dimensions don't match index");
      }
      if (ref_mm.ncol() - 2 != n_clim) {
            stop("Number of climate variables doesn't match index");
      }
      if (focal_mm.ncol() != ref_mm.ncol()) {
            stop("Focal and reference data must have same number of columns");
      }

      const int n_focal = focal_mm.nrow();
      const bool use_haversine = (coord_type == "lonlat");
      const double R_earth = 6371.0088;

      // Parse climate constraints
      bool use_scalar_clim = false;
      bool use_pervar_clim = false;
      double max_clim_scalar = std::numeric_limits<double>::infinity();
      std::vector<double> max_clim_pervar_std(n_clim, std::numeric_limits<double>::infinity());

      if (max_clim.size() == 1) {
            double v = max_clim[0];
            if (std::isfinite(v)) {
                  use_scalar_clim = true;
                  max_clim_scalar = v;
            }
      } else if (max_clim.size() == n_clim) {
            for (int i = 0; i < n_clim; ++i) {
                  max_clim_pervar_std[i] = max_clim[i];
            }
            use_pervar_clim = true;
      } else if (max_clim.size() > 1) {
            stop("max_clim must be length 1 or equal to number of climate variables");
      }

      const bool use_geog_filter = std::isfinite(max_geog);
      const double max_geog_chord = (use_ecef && std::isfinite(max_geog))
            ? km_to_chord(max_geog, R_earth)
                  : max_geog;

      const SelectCode scode = static_cast<SelectCode>(select_code);

      // Parse aggregate_codes vector
      const int n_stats = aggregate_codes.size();
      if (n_stats == 0) {
            stop("aggregate_codes must have at least one element");
      }

      std::vector<AggregateCode> acodes(n_stats);
      for (int i = 0; i < n_stats; ++i) {
            acodes[i] = static_cast<AggregateCode>(aggregate_codes[i]);
      }

      // Check for "none" (pairs mode)
      const bool return_pairs = (n_stats == 1 && acodes[0] == AggregateCode::NONE);

      const WeightCode wcode = static_cast<WeightCode>(weight_code);

      // Validate and parse SE code
      if (se_code < 0 || se_code > 2) {
            stop("se_code must be 0 (none), 1 (ess), or 2 (design)");
      }
      const SeCode scode_se = static_cast<SeCode>(se_code);

      // Pre-compute weight parameters for efficiency
      auto weight_params = precompute_weight_params(wcode, theta);
      double weight_param1 = weight_params.first;
      double weight_param2 = weight_params.second;

      // Parse x_cov parameter
      bool use_mahalanobis = false;
      const double* x_cov_ptr = nullptr;
      int x_cov_stride = 0;
      int n_cov_components = 0;

      if (!Rf_isNull(x_cov_sexp) && x_cov_sexp != R_NilValue) {
            NumericMatrix x_cov_mat = as<NumericMatrix>(x_cov_sexp);

            // Validate dimensions
            if (x_cov_mat.nrow() != n_focal) {
                  stop("x_cov must have same number of rows as focal data");
            }

            n_cov_components = n_clim * (n_clim + 1) / 2;
            if (x_cov_mat.ncol() != n_cov_components) {
                  stop("x_cov must have n_clim * (n_clim + 1) / 2 columns");
            }

            use_mahalanobis = true;
            x_cov_ptr = REAL(x_cov_mat);
            x_cov_stride = n_focal;  // Column-major stride
      }

      // Parse values parameter
      bool has_values = false;
      const double* values_ptr = nullptr;
      int values_stride = 0;
      int n_vars = 0;

      if (!Rf_isNull(values_sexp) && values_sexp != R_NilValue) {
            NumericMatrix values_mat = as<NumericMatrix>(values_sexp);

            // Validate dimensions
            if (values_mat.nrow() != n_ref) {
                  stop("values must have same number of rows as reference data");
            }

            has_values = true;
            values_ptr = REAL(values_mat);
            values_stride = n_ref;  // Column-major stride
            n_vars = values_mat.ncol();
      }

      // Parse covariates parameter (for regression stat)
      bool has_covariates = false;
      const double* covariates_ptr = nullptr;
      int covariates_stride = 0;
      int n_covs = 0;

      if (!Rf_isNull(covariates_sexp) && covariates_sexp != R_NilValue) {
            NumericMatrix covariates_mat = as<NumericMatrix>(covariates_sexp);

            if (covariates_mat.nrow() != n_ref) {
                  stop("covariates must have same number of rows as reference data");
            }

            has_covariates = true;
            covariates_ptr = REAL(covariates_mat);
            covariates_stride = n_ref;
            n_covs = covariates_mat.ncol();
      }

      // Convert n_classes_per_var (R IntegerVector) to std::vector<int>.
      // Length 0 means tabulate is not requested. When tabulate IS requested,
      // R must pass length == n_vars.
      std::vector<int> n_classes_per_var_std;
      if (n_classes_per_var.size() > 0) {
            if (n_classes_per_var.size() != n_vars) {
                  stop("n_classes_per_var length must equal number of y columns");
            }
            n_classes_per_var_std.assign(n_classes_per_var.begin(),
                                         n_classes_per_var.end());
      }

      // Get ECEF data pointer if applicable
      const double* ref_latt_ptr;
      int stride_latt_r;

      if (use_ecef) {
            SEXP ecef_sexp = idx["ecef_xptr"];
            if (ecef_sexp == R_NilValue) {
                  stop("ECEF data missing from index");
            }
            XPtr< std::vector<double> > ecef_xptr(ecef_sexp);
            if (ecef_xptr.get() == nullptr) {
                  stop("ECEF data pointer is null");
            }
            ref_latt_ptr = ecef_xptr->data();
            stride_latt_r = n_ref;
      } else {
            ref_latt_ptr = REAL(ref_mm);
            stride_latt_r = n_ref;
      }

      // Execute query using workers
      if (return_pairs) {
            const int k_knn = (scode == SelectCode::ALL ? 0 : k);
            // Pair mode
            std::vector< std::vector<int> > out_indices(n_focal);
            std::vector< std::vector<double> > out_weights(n_focal);

            PairWorker worker(focal_mm,
                              ref_mm,
                              ref_latt_ptr,
                              stride_latt_r,
                              true,  // use_lattice
                              use_geog_filter,
                              use_haversine,
                              use_scalar_clim,
                              use_pervar_clim,
                              max_clim_scalar,
                              max_geog,
                              max_geog_chord,
                              max_clim_pervar_std,
                              scode,
                              k_knn,
                              lattice_ptr,
                              use_ecef,
                              R_earth,
                              use_mahalanobis,
                              x_cov_ptr,
                              x_cov_stride,
                              n_cov_components,
                              out_indices,
                              out_weights);

            parallelFor(0, static_cast<std::size_t>(n_focal), worker);

            // Handle k=1 NA cases - need to add weight too
            if (k == 1 && (scode == SelectCode::KNN_CLIM || scode == SelectCode::KNN_GEOG)) {
                  for (int i = 0; i < n_focal; ++i) {
                        if (out_indices[i].empty()) {
                              out_indices[i].push_back(0);
                              out_weights[i].push_back(1.0);
                        }
                  }
            }

            // Convert to List with both indices and weights
            // Structure: List with two named elements: "indices" and "weights"
            List out;

            // Create indices list
            List indices_list(n_focal);
            for (int i = 0; i < n_focal; ++i) {
                  const std::vector<int>& v = out_indices[i];
                  IntegerVector idx_vec(v.size());
                  for (std::size_t j = 0; j < v.size(); ++j) {
                        idx_vec[j] = v[j];
                  }
                  indices_list[i] = idx_vec;
            }

            // Create weights list
            List weights_list(n_focal);
            for (int i = 0; i < n_focal; ++i) {
                  const std::vector<double>& w = out_weights[i];
                  NumericVector wgt_vec(w.size());
                  for (std::size_t j = 0; j < w.size(); ++j) {
                        wgt_vec[j] = w[j];
                  }
                  weights_list[i] = wgt_vec;
            }

            // Package both into the output list
            out = List::create(
                  Named("indices") = indices_list,
                  Named("weights") = weights_list
            );

            // Attach diagnostics
            out.attr("n_focal") = n_focal;
            out.attr("n_ref") = n_ref;
            out.attr("n_clim") = n_clim;
            out.attr("max_dist") = max_geog;
            out.attr("max_clim") = max_clim;
            out.attr("geo_mode") = coord_type;
            out.attr("binning_method") = use_ecef ? "lattice_ecef" : "lattice";
            out.attr("total_bins") = static_cast<double>(lattice_ptr->total_bins);
            out.attr("n_bins_nonempty") = static_cast<double>(lattice_ptr->n_cells_nonempty);
            out.attr("min_bin_occupancy") = static_cast<double>(lattice_ptr->min_cell_occ);
            out.attr("max_bin_occupancy") = static_cast<double>(lattice_ptr->max_cell_occ);

            double avg_bin_occupancy = 0.0;
            if (lattice_ptr->total_bins > 0) {
                  avg_bin_occupancy = static_cast<double>(lattice_ptr->n_points) /
                        static_cast<double>(lattice_ptr->total_bins);
            }
            out.attr("avg_bin_occupancy") = avg_bin_occupancy;

            double avg_nonempty_bin_occupancy = 0.0;
            if (lattice_ptr->n_cells_nonempty > 0) {
                  avg_nonempty_bin_occupancy = static_cast<double>(lattice_ptr->n_points) /
                        static_cast<double>(lattice_ptr->n_cells_nonempty);
            }
            out.attr("avg_nonempty_bin_occupancy") = avg_nonempty_bin_occupancy;

            return out;
      }

      // Aggregate modes
      // Count regular vs value stats vs tabulate to calculate total columns
      int n_regular_stats = 0;
      int n_value_stats = 0;
      bool has_weighted_mean_stat = false;
      bool has_regression_stat = false;
      bool has_tabulate_stat = false;

      for (int i = 0; i < n_stats; ++i) {
            if (acodes[i] == AggregateCode::SUM ||
                acodes[i] == AggregateCode::MEAN ||
                acodes[i] == AggregateCode::WEIGHTED_SUM ||
                acodes[i] == AggregateCode::WEIGHTED_MEAN) {
                  n_value_stats++;
                  if (acodes[i] == AggregateCode::WEIGHTED_MEAN) {
                        has_weighted_mean_stat = true;
                  }
            } else if (acodes[i] == AggregateCode::REGRESSION) {
                  has_regression_stat = true;
            } else if (acodes[i] == AggregateCode::TABULATE) {
                  has_tabulate_stat = true;
            } else {
                  n_regular_stats++;
            }
      }

      // Column counts
      // - Value stats: n_value_stats × n_vars columns
      // - Regression: reg_dim columns per var for coefficients,
      //               plus reg_dim columns per var for SEs if se != NONE
      // - weighted_mean SE: one additional column per var if se != NONE AND
      //                    weighted_mean is among requested stats
      // - Tabulate: sum of K_v across y variables (one column per class)
      const bool want_se = (scode_se != SeCode::NONE);

      const int reg_dim = has_regression_stat ? (n_covs + 1) : 0;
      const int reg_cols_per_var = reg_dim + (want_se ? reg_dim : 0);
      const int n_regression_cols = has_regression_stat
      ? (n_vars * reg_cols_per_var)
            : 0;

      const int wm_se_cols = (want_se && has_weighted_mean_stat) ? n_vars : 0;

      int n_tabulate_cols = 0;
      if (has_tabulate_stat) {
            for (int v = 0; v < static_cast<int>(n_classes_per_var_std.size()); ++v) {
                  n_tabulate_cols += n_classes_per_var_std[v];
            }
      }

      // Total columns = regular stats + (value stats × n_vars) + wm SE cols
      //                 + tabulate cols + regression cols
      const int n_total_cols = n_regular_stats
      + (n_value_stats * n_vars)
            + wm_se_cols
            + n_tabulate_cols
            + n_regression_cols;

            // Allocate flat vector for all stats
            // Layout: [focal0_stat0, focal0_stat1, ..., focal1_stat0, focal1_stat1, ...]
            std::vector<double> agg_vals(n_focal * n_total_cols, NA_REAL);

            AggWorker aworker(focal_mm,
                              ref_mm,
                              ref_latt_ptr,
                              stride_latt_r,
                              true,  // use_lattice
                              use_geog_filter,
                              use_haversine,
                              use_scalar_clim,
                              use_pervar_clim,
                              max_clim_scalar,
                              max_geog,
                              max_geog_chord,
                              max_clim_pervar_std,
                              scode,
                              acodes,
                              wcode,
                              weight_param1,
                              weight_param2,
                              lattice_ptr,
                              use_ecef,
                              R_earth,
                              use_mahalanobis,
                              x_cov_ptr,
                              x_cov_stride,
                              n_cov_components,
                              has_values,
                              values_ptr,
                              values_stride,
                              n_vars,
                              has_covariates,
                              covariates_ptr,
                              covariates_stride,
                              n_covs,
                              lambda,
                              scode_se,
                              n_classes_per_var_std,
                              exclude_self,
                              agg_vals,
                              n_total_cols);

            parallelFor(0, static_cast<std::size_t>(n_focal), aworker);

            // Convert flat vector to matrix: n_focal rows x n_total_cols columns
            NumericMatrix agg(n_focal, n_total_cols);
            for (int i = 0; i < n_focal; ++i) {
                  for (int s = 0; s < n_total_cols; ++s) {
                        agg(i, s) = agg_vals[i * n_total_cols + s];
                  }
            }

            // Add diagnostics as attributes
            agg.attr("n_focal") = n_focal;
            agg.attr("n_ref") = n_ref;
            agg.attr("n_clim") = n_clim;
            agg.attr("max_dist") = max_geog;
            agg.attr("max_clim") = max_clim;
            agg.attr("geo_mode") = coord_type;
            agg.attr("binning_method") = use_ecef ? "lattice_ecef" : "lattice";
            agg.attr("total_bins") = static_cast<double>(lattice_ptr->total_bins);
            agg.attr("n_bins_nonempty") = static_cast<double>(lattice_ptr->n_cells_nonempty);
            agg.attr("min_bin_occupancy") = static_cast<double>(lattice_ptr->min_cell_occ);
            agg.attr("max_bin_occupancy") = static_cast<double>(lattice_ptr->max_cell_occ);

            double avg_bin_occupancy = 0.0;
            if (lattice_ptr->total_bins > 0) {
                  avg_bin_occupancy = static_cast<double>(lattice_ptr->n_points) /
                        static_cast<double>(lattice_ptr->total_bins);
            }
            agg.attr("avg_bin_occupancy") = avg_bin_occupancy;

            double avg_nonempty_bin_occupancy = 0.0;
            if (lattice_ptr->n_cells_nonempty > 0) {
                  avg_nonempty_bin_occupancy = static_cast<double>(lattice_ptr->n_points) /
                        static_cast<double>(lattice_ptr->n_cells_nonempty);
            }
            agg.attr("avg_nonempty_bin_occupancy") = avg_nonempty_bin_occupancy;

            return agg;
}
