#include <RcppArmadillo.h>
#include "analog_utils.cpp"
#include "lattice_index.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// Aggregation function types
enum AggType { NMAX, SUM, MEAN, COUNT, ALL };

inline AggType parse_agg_type(std::string fun_str) {
      if (fun_str == "nmax") return NMAX;
      if (fun_str == "sum") return SUM;
      if (fun_str == "mean") return MEAN;
      if (fun_str == "count") return COUNT;
      if (fun_str == "all") return ALL;
      return NMAX;  // default
}


// 1. Matrix focal, Raster ref
// [[Rcpp::export]]
DataFrame find_analogs_matrix_raster(
            IntegerVector focal_cells,
            arma::mat focal_coords,
            arma::mat focal_climate,
            arma::cube ref_array,
            arma::mat ref_coords,
            double xres,
            double yres,
            int nrow,
            int ncol,
            Nullable<double> max_dist,
            Nullable<double> max_clim,
            std::string weight,
            std::string fun,
            Nullable<int> n
) {

      int n_focal = focal_cells.size();
      int n_vars = ref_array.n_slices;

      WeightType weight_type = parse_weight_type(weight);
      AggType agg_type = parse_agg_type(fun);

      int n_analogs = 0;
      if (agg_type == NMAX) {
            if (n.isNull()) {
                  stop("n must be specified when fun = 'nmax'");
            }
            n_analogs = as<int>(n);
      }

      double radius_km = max_dist.isNotNull() ? as<double>(max_dist) : R_PosInf;
      double clim_threshold = max_clim.isNotNull() ? as<double>(max_clim) : R_PosInf;

      // Calculate bounding box radius in cells
      int radius_cells_x, radius_cells_y;
      if (max_dist.isNull()) {
            // No distance constraint - search entire raster
            radius_cells_x = ncol;
            radius_cells_y = nrow;
      } else {
            double km_per_cell_x = xres * 111.0;
            double km_per_cell_y = yres * 111.0;
            radius_cells_x = std::ceil(radius_km / km_per_cell_x);
            radius_cells_y = std::ceil(radius_km / km_per_cell_y);
      }

      // Results storage
      std::vector<int> result_focal_idx;
      std::vector<double> result_focal_x;
      std::vector<double> result_focal_y;
      std::vector<int> result_analog_idx;
      std::vector<double> result_analog_x;
      std::vector<double> result_analog_y;
      std::vector<double> result_clim_dist;
      std::vector<double> result_geo_dist;
      std::vector<double> result_weight;
      std::vector<double> result_value;

      // Loop over focal sites
      for (int f = 0; f < n_focal; f++) {

            int focal_cell = focal_cells[f] - 1;
            int focal_row = focal_cell / ncol;
            int focal_col = focal_cell % ncol;

            double focal_x = focal_coords(f, 0);
            double focal_y = focal_coords(f, 1);

            arma::rowvec focal_clim = focal_climate.row(f);

            if (focal_clim.has_nan()) {
                  Rf_warning("Focal site %d has NA climate values, skipping", f + 1);
                  continue;
            }

            // Bounding box
            int row_min = std::max(0, focal_row - radius_cells_y);
            int row_max = std::min(nrow - 1, focal_row + radius_cells_y);
            int col_min = std::max(0, focal_col - radius_cells_x);
            int col_max = std::min(ncol - 1, focal_col + radius_cells_x);

            // Collect candidates
            std::vector<int> candidate_cells;
            std::vector<double> candidate_x;
            std::vector<double> candidate_y;
            std::vector<double> geo_distances;
            std::vector<double> clim_distances;
            std::vector<double> weights;

            for (int r = row_min; r <= row_max; r++) {
                  for (int c = col_min; c <= col_max; c++) {

                        // Extract candidate climate
                        arma::rowvec cand_clim(n_vars);
                        bool has_na = false;
                        for (int v = 0; v < n_vars; v++) {
                              cand_clim(v) = ref_array(r, c, v);
                              if (!arma::is_finite(cand_clim(v))) {
                                    has_na = true;
                                    break;
                              }
                        }

                        if (has_na) continue;

                        // Get candidate coordinates
                        int cand_cell = r * ncol + c;
                        double cand_x = ref_coords(cand_cell, 0);
                        double cand_y = ref_coords(cand_cell, 1);

                        // Calculate distances
                        double geo_dist = calc_geo_distance(focal_x, focal_y, cand_x, cand_y);
                        double clim_dist = calc_climate_distance(focal_clim, cand_clim);

                        if (geo_dist > radius_km) continue;
                        if (clim_dist > clim_threshold) continue;

                        double w = calc_weight(clim_dist, geo_dist, weight_type);

                        candidate_cells.push_back(cand_cell);
                        candidate_x.push_back(cand_x);
                        candidate_y.push_back(cand_y);
                        geo_distances.push_back(geo_dist);
                        clim_distances.push_back(clim_dist);
                        weights.push_back(w);
                  }
            }

            int n_candidates = candidate_cells.size();

            // For aggregation functions, handle no candidates case
            if (n_candidates == 0) {
                  if (agg_type == COUNT) {
                        result_focal_idx.push_back(f + 1);
                        result_focal_x.push_back(focal_x);
                        result_focal_y.push_back(focal_y);
                        result_value.push_back(0.0);
                  } else if (agg_type == SUM || agg_type == MEAN) {
                        result_focal_idx.push_back(f + 1);
                        result_focal_x.push_back(focal_x);
                        result_focal_y.push_back(focal_y);
                        result_value.push_back(NA_REAL);
                  }
                  continue;
            }

            // Apply aggregation
            if (agg_type == NMAX) {
                  arma::vec weight_vec(weights);
                  arma::uvec top_idx = find_top_n_weights(weight_vec, n_analogs);

                  for (size_t i = 0; i < top_idx.n_elem; i++) {
                        int idx = top_idx[i];
                        result_focal_idx.push_back(f + 1);
                        result_focal_x.push_back(focal_x);
                        result_focal_y.push_back(focal_y);
                        result_analog_idx.push_back(candidate_cells[idx] + 1);
                        result_analog_x.push_back(candidate_x[idx]);
                        result_analog_y.push_back(candidate_y[idx]);
                        result_clim_dist.push_back(clim_distances[idx]);
                        result_geo_dist.push_back(geo_distances[idx]);
                        result_weight.push_back(weights[idx]);
                  }

            } else if (agg_type == SUM) {
                  double sum_weights = 0.0;
                  for (double w : weights) sum_weights += w;
                  result_focal_idx.push_back(f + 1);
                  result_focal_x.push_back(focal_x);
                  result_focal_y.push_back(focal_y);
                  result_value.push_back(sum_weights);

            } else if (agg_type == MEAN) {
                  double sum_weights = 0.0;
                  for (double w : weights) sum_weights += w;
                  result_focal_idx.push_back(f + 1);
                  result_focal_x.push_back(focal_x);
                  result_focal_y.push_back(focal_y);
                  result_value.push_back(sum_weights / n_candidates);

            } else if (agg_type == COUNT) {
                  result_focal_idx.push_back(f + 1);
                  result_focal_x.push_back(focal_x);
                  result_focal_y.push_back(focal_y);
                  result_value.push_back(n_candidates);

            } else if (agg_type == ALL) {
                  for (int i = 0; i < n_candidates; i++) {
                        result_focal_idx.push_back(f + 1);
                        result_focal_x.push_back(focal_x);
                        result_focal_y.push_back(focal_y);
                        result_analog_idx.push_back(candidate_cells[i] + 1);
                        result_analog_x.push_back(candidate_x[i]);
                        result_analog_y.push_back(candidate_y[i]);
                        result_clim_dist.push_back(clim_distances[i]);
                        result_geo_dist.push_back(geo_distances[i]);
                        result_weight.push_back(weights[i]);
                  }
            }
      }

      if (agg_type == SUM || agg_type == MEAN || agg_type == COUNT) {
            return DataFrame::create(
                  Named("focal_index") = result_focal_idx,
                  Named("focal_x") = result_focal_x,
                  Named("focal_y") = result_focal_y,
                  Named("value") = result_value
            );
      } else {
            return DataFrame::create(
                  Named("focal_index") = result_focal_idx,
                  Named("focal_x") = result_focal_x,
                  Named("focal_y") = result_focal_y,
                  Named("analog_index") = result_analog_idx,
                  Named("analog_x") = result_analog_x,
                  Named("analog_y") = result_analog_y,
                  Named("clim_dist") = result_clim_dist,
                  Named("geog_dist") = result_geo_dist,
                  Named("weight") = result_weight
            );
      }
}


// 2. Matrix focal, Matrix ref
// [[Rcpp::export]]
DataFrame find_analogs_matrix_matrix(
            arma::mat focal_coords,
            arma::mat focal_climate,
            arma::mat ref_coords,
            arma::mat ref_climate,
            Nullable<double> max_dist,
            Nullable<double> max_clim,
            std::string weight,
            std::string fun,
            Nullable<int> n
) {

      int n_focal = focal_coords.n_rows;
      int n_ref = ref_coords.n_rows;
      int n_vars = focal_climate.n_cols;

      WeightType weight_type = parse_weight_type(weight);
      AggType agg_type = parse_agg_type(fun);

      int n_analogs = 0;
      if (agg_type == NMAX) {
            if (n.isNull()) {
                  stop("n must be specified when fun = 'nmax'");
            }
            n_analogs = as<int>(n);
      }

      double radius_km = max_dist.isNotNull() ? as<double>(max_dist) : R_PosInf;
      double clim_threshold = max_clim.isNotNull() ? as<double>(max_clim) : R_PosInf;

      // Build indices if beneficial (n_ref large enough)
      LatticeIndex<2> geo_idx;
      bool use_geo_idx = false;
      if (max_dist.isNotNull() && n_ref > 5000) {
            geo_idx.build(ref_coords, 20);
            use_geo_idx = true;
      }

      LatticeIndex<10> clim_idx;  // Support up to 10 climate variables
      bool use_clim_idx = false;
      if (max_clim.isNotNull() && n_ref > 5000 && n_vars <= 10) {
            clim_idx.build(ref_climate, 15);
            use_clim_idx = true;
      }

      std::vector<int> result_focal_idx;
      std::vector<double> result_focal_x;
      std::vector<double> result_focal_y;
      std::vector<int> result_analog_idx;
      std::vector<double> result_analog_x;
      std::vector<double> result_analog_y;
      std::vector<double> result_clim_dist;
      std::vector<double> result_geo_dist;
      std::vector<double> result_weight;
      std::vector<double> result_value;

      for (int f = 0; f < n_focal; f++) {

            double focal_x = focal_coords(f, 0);
            double focal_y = focal_coords(f, 1);
            arma::rowvec focal_clim = focal_climate.row(f);

            if (focal_clim.has_nan()) {
                  Rf_warning("Focal site %d has NA climate values, skipping", f + 1);
                  continue;
            }

            // Get candidate set using indices
            std::vector<int> candidate_pool;

            if (use_geo_idx && use_clim_idx) {
                  // Use both indices - intersect results
                  auto geo_candidates = geo_idx.query_box(focal_coords.row(f), radius_km);
                  auto clim_candidates = clim_idx.query_box(focal_clim, clim_threshold);

                  // Intersect (candidates in both sets)
                  std::set<int> geo_set(geo_candidates.begin(), geo_candidates.end());
                  for (int c : clim_candidates) {
                        if (geo_set.count(c)) {
                              candidate_pool.push_back(c);
                        }
                  }

            } else if (use_geo_idx) {
                  candidate_pool = geo_idx.query_box(focal_coords.row(f), radius_km);

            } else if (use_clim_idx) {
                  candidate_pool = clim_idx.query_box(focal_clim, clim_threshold);

            } else {
                  // No indices - use all ref points
                  candidate_pool.resize(n_ref);
                  std::iota(candidate_pool.begin(), candidate_pool.end(), 0);
            }

            // Now process candidates (same as before, but on filtered set)
            std::vector<int> candidate_idx;
            std::vector<double> geo_distances;
            std::vector<double> clim_distances;
            std::vector<double> weights;

            for (int r : candidate_pool) {

                  arma::rowvec ref_clim = ref_climate.row(r);

                  if (ref_clim.has_nan()) continue;

                  double geo_dist = calc_geo_distance(focal_x, focal_y,
                                                      ref_coords(r, 0), ref_coords(r, 1));
                  double clim_dist = calc_climate_distance(focal_clim, ref_clim);

                  if (geo_dist > radius_km) continue;
                  if (clim_dist > clim_threshold) continue;

                  double w = calc_weight(clim_dist, geo_dist, weight_type);

                  candidate_idx.push_back(r);
                  geo_distances.push_back(geo_dist);
                  clim_distances.push_back(clim_dist);
                  weights.push_back(w);
            }

            int n_candidates = candidate_idx.size();

            // For aggregation functions, handle no candidates case
            if (n_candidates == 0) {
                  if (agg_type == COUNT) {
                        result_focal_idx.push_back(f + 1);
                        result_focal_x.push_back(focal_x);
                        result_focal_y.push_back(focal_y);
                        result_value.push_back(0.0);
                  } else if (agg_type == SUM || agg_type == MEAN) {
                        result_focal_idx.push_back(f + 1);
                        result_focal_x.push_back(focal_x);
                        result_focal_y.push_back(focal_y);
                        result_value.push_back(NA_REAL);
                  }
                  continue;
            }

            if (agg_type == NMAX) {
                  arma::vec weight_vec(weights);
                  arma::uvec top_idx = find_top_n_weights(weight_vec, n_analogs);

                  for (size_t i = 0; i < top_idx.n_elem; i++) {
                        int idx = top_idx[i];
                        int ref_idx = candidate_idx[idx];
                        result_focal_idx.push_back(f + 1);
                        result_focal_x.push_back(focal_x);
                        result_focal_y.push_back(focal_y);
                        result_analog_idx.push_back(ref_idx + 1);
                        result_analog_x.push_back(ref_coords(ref_idx, 0));
                        result_analog_y.push_back(ref_coords(ref_idx, 1));
                        result_clim_dist.push_back(clim_distances[idx]);
                        result_geo_dist.push_back(geo_distances[idx]);
                        result_weight.push_back(weights[idx]);
                  }

            } else if (agg_type == SUM) {
                  double sum_weights = 0.0;
                  for (double w : weights) sum_weights += w;
                  result_focal_idx.push_back(f + 1);
                  result_focal_x.push_back(focal_x);
                  result_focal_y.push_back(focal_y);
                  result_value.push_back(sum_weights);

            } else if (agg_type == MEAN) {
                  double sum_weights = 0.0;
                  for (double w : weights) sum_weights += w;
                  result_focal_idx.push_back(f + 1);
                  result_focal_x.push_back(focal_x);
                  result_focal_y.push_back(focal_y);
                  result_value.push_back(sum_weights / n_candidates);

            } else if (agg_type == COUNT) {
                  result_focal_idx.push_back(f + 1);
                  result_focal_x.push_back(focal_x);
                  result_focal_y.push_back(focal_y);
                  result_value.push_back(n_candidates);

            } else if (agg_type == ALL) {
                  for (int i = 0; i < n_candidates; i++) {
                        int ref_idx = candidate_idx[i];
                        result_focal_idx.push_back(f + 1);
                        result_focal_x.push_back(focal_x);
                        result_focal_y.push_back(focal_y);
                        result_analog_idx.push_back(ref_idx + 1);
                        result_analog_x.push_back(ref_coords(ref_idx, 0));
                        result_analog_y.push_back(ref_coords(ref_idx, 1));
                        result_clim_dist.push_back(clim_distances[i]);
                        result_geo_dist.push_back(geo_distances[i]);
                        result_weight.push_back(weights[i]);
                  }
            }
      }

      if (agg_type == SUM || agg_type == MEAN || agg_type == COUNT) {
            return DataFrame::create(
                  Named("focal_index") = result_focal_idx,
                  Named("focal_x") = result_focal_x,
                  Named("focal_y") = result_focal_y,
                  Named("value") = result_value
            );
      } else {
            return DataFrame::create(
                  Named("focal_index") = result_focal_idx,
                  Named("focal_x") = result_focal_x,
                  Named("focal_y") = result_focal_y,
                  Named("analog_index") = result_analog_idx,
                  Named("analog_x") = result_analog_x,
                  Named("analog_y") = result_analog_y,
                  Named("clim_dist") = result_clim_dist,
                  Named("geog_dist") = result_geo_dist,
                  Named("weight") = result_weight
            );
      }
}


// 3. Raster focal, Matrix ref
// [[Rcpp::export]]
DataFrame find_analogs_raster_matrix(
            arma::cube focal_array,
            arma::mat focal_coords,
            arma::mat ref_coords,
            arma::mat ref_climate,
            Nullable<double> max_dist,
            Nullable<double> max_clim,
            std::string weight,
            std::string fun,
            Nullable<int> n
) {

      int n_focal = focal_array.n_rows * focal_array.n_cols;
      int n_ref = ref_coords.n_rows;
      int n_vars = focal_array.n_slices;

      WeightType weight_type = parse_weight_type(weight);
      AggType agg_type = parse_agg_type(fun);

      int n_analogs = 0;
      if (agg_type == NMAX) {
            if (n.isNull()) {
                  stop("n must be specified when fun = 'nmax'");
            }
            n_analogs = as<int>(n);
      }

      double radius_km = max_dist.isNotNull() ? as<double>(max_dist) : R_PosInf;
      double clim_threshold = max_clim.isNotNull() ? as<double>(max_clim) : R_PosInf;

      std::vector<int> result_focal_idx;
      std::vector<double> result_focal_x;
      std::vector<double> result_focal_y;
      std::vector<int> result_analog_idx;
      std::vector<double> result_analog_x;
      std::vector<double> result_analog_y;
      std::vector<double> result_clim_dist;
      std::vector<double> result_geo_dist;
      std::vector<double> result_weight;
      std::vector<double> result_value;

      for (int cell = 0; cell < n_focal; cell++) {

            double focal_x = focal_coords(cell, 0);
            double focal_y = focal_coords(cell, 1);

            int focal_row = cell / focal_array.n_cols;
            int focal_col = cell % focal_array.n_cols;
            arma::rowvec focal_clim(n_vars);
            bool has_na = false;
            for (int v = 0; v < n_vars; v++) {
                  focal_clim(v) = focal_array(focal_row, focal_col, v);
                  if (!arma::is_finite(focal_clim(v))) {
                        has_na = true;
                        break;
                  }
            }

            if (has_na) continue;

            std::vector<int> candidate_idx;
            std::vector<double> geo_distances;
            std::vector<double> clim_distances;
            std::vector<double> weights;

            for (int r = 0; r < n_ref; r++) {

                  arma::rowvec ref_clim = ref_climate.row(r);

                  if (ref_clim.has_nan()) continue;

                  double geo_dist = calc_geo_distance(focal_x, focal_y,
                                                      ref_coords(r, 0), ref_coords(r, 1));
                  double clim_dist = calc_climate_distance(focal_clim, ref_clim);

                  if (geo_dist > radius_km) continue;
                  if (clim_dist > clim_threshold) continue;

                  double w = calc_weight(clim_dist, geo_dist, weight_type);

                  candidate_idx.push_back(r);
                  geo_distances.push_back(geo_dist);
                  clim_distances.push_back(clim_dist);
                  weights.push_back(w);
            }

            int n_candidates = candidate_idx.size();

            if (n_candidates == 0) continue;

            if (agg_type == NMAX) {
                  arma::vec weight_vec(weights);
                  arma::uvec top_idx = find_top_n_weights(weight_vec, n_analogs);

                  for (size_t i = 0; i < top_idx.n_elem; i++) {
                        int idx = top_idx[i];
                        int ref_idx = candidate_idx[idx];
                        result_focal_idx.push_back(cell + 1);
                        result_focal_x.push_back(focal_x);
                        result_focal_y.push_back(focal_y);
                        result_analog_idx.push_back(ref_idx + 1);
                        result_analog_x.push_back(ref_coords(ref_idx, 0));
                        result_analog_y.push_back(ref_coords(ref_idx, 1));
                        result_clim_dist.push_back(clim_distances[idx]);
                        result_geo_dist.push_back(geo_distances[idx]);
                        result_weight.push_back(weights[idx]);
                  }

            } else if (agg_type == SUM) {
                  double sum_weights = 0.0;
                  for (double w : weights) sum_weights += w;
                  result_focal_idx.push_back(cell + 1);
                  result_focal_x.push_back(focal_x);
                  result_focal_y.push_back(focal_y);
                  result_value.push_back(sum_weights);

            } else if (agg_type == MEAN) {
                  double sum_weights = 0.0;
                  for (double w : weights) sum_weights += w;
                  result_focal_idx.push_back(cell + 1);
                  result_focal_x.push_back(focal_x);
                  result_focal_y.push_back(focal_y);
                  result_value.push_back(sum_weights / n_candidates);

            } else if (agg_type == COUNT) {
                  result_focal_idx.push_back(cell + 1);
                  result_focal_x.push_back(focal_x);
                  result_focal_y.push_back(focal_y);
                  result_value.push_back(n_candidates);

            } else if (agg_type == ALL) {
                  for (int i = 0; i < n_candidates; i++) {
                        int ref_idx = candidate_idx[i];
                        result_focal_idx.push_back(cell + 1);
                        result_focal_x.push_back(focal_x);
                        result_focal_y.push_back(focal_y);
                        result_analog_idx.push_back(ref_idx + 1);
                        result_analog_x.push_back(ref_coords(ref_idx, 0));
                        result_analog_y.push_back(ref_coords(ref_idx, 1));
                        result_clim_dist.push_back(clim_distances[i]);
                        result_geo_dist.push_back(geo_distances[i]);
                        result_weight.push_back(weights[i]);
                  }
            }
      }

      if (agg_type == SUM || agg_type == MEAN || agg_type == COUNT) {
            return DataFrame::create(
                  Named("focal_index") = result_focal_idx,
                  Named("focal_x") = result_focal_x,
                  Named("focal_y") = result_focal_y,
                  Named("value") = result_value
            );
      } else {
            return DataFrame::create(
                  Named("focal_index") = result_focal_idx,
                  Named("focal_x") = result_focal_x,
                  Named("focal_y") = result_focal_y,
                  Named("analog_index") = result_analog_idx,
                  Named("analog_x") = result_analog_x,
                  Named("analog_y") = result_analog_y,
                  Named("clim_dist") = result_clim_dist,
                  Named("geog_dist") = result_geo_dist,
                  Named("weight") = result_weight
            );
      }
}
