#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Calculate Euclidean climate distance between two climate vectors
inline double calc_climate_distance(const arma::rowvec& focal_clim,
                                    const arma::rowvec& ref_clim) {
      arma::rowvec diff = focal_clim - ref_clim;
      return std::sqrt(arma::dot(diff, diff));
}

// Calculate geographic distance between two points (simple Euclidean * 111 km/deg)
inline double calc_geo_distance(double x1, double y1, double x2, double y2) {
      double dx = (x2 - x1) * 111.0;  // approximate km per degree longitude
      double dy = (y2 - y1) * 111.0;  // approximate km per degree latitude
      return std::sqrt(dx * dx + dy * dy);
}

// Calculate geographic distance on a regular grid (in cells)
inline double calc_grid_distance(int row1, int col1, int row2, int col2,
                                 double yres, double xres) {
      double dy = (row2 - row1) * yres * 111.0;
      double dx = (col2 - col1) * xres * 111.0;
      return std::sqrt(dx * dx + dy * dy);
}

// Weight calculation
enum WeightType { UNIFORM, INVERSE_CLIM, INVERSE_DIST };

inline double calc_weight(double clim_dist, double geo_dist, WeightType type) {
      switch(type) {
      case UNIFORM:
            return 1.0;
      case INVERSE_CLIM:
            return 1.0 / (1.0 + clim_dist);
      case INVERSE_DIST:
            return 1.0 / (1.0 + geo_dist);
      default:
            return 1.0;
      }
}

// Convert weight string to enum
inline WeightType parse_weight_type(std::string weight_str) {
      if (weight_str == "uniform") return UNIFORM;
      if (weight_str == "inverse_clim") return INVERSE_CLIM;
      if (weight_str == "inverse_dist") return INVERSE_DIST;
      return UNIFORM;  // default
}

// Find indices of top N largest values (for weights, larger is better)
inline arma::uvec find_top_n_weights(const arma::vec& weights, int n) {
      arma::uvec sorted_indices = arma::sort_index(weights, "descend");
      int n_return = std::min(n, (int)weights.n_elem);
      return sorted_indices.head(n_return);
}

// Apply geographic radius filter, return indices of values <= radius
inline arma::uvec filter_by_radius(const arma::vec& geo_distances, double radius_km) {
      return arma::find(geo_distances <= radius_km);
}

// Apply climate threshold filter, return indices of values <= threshold
inline arma::uvec filter_by_clim(const arma::vec& clim_distances, double clim_threshold) {
      return arma::find(clim_distances <= clim_threshold);
}
