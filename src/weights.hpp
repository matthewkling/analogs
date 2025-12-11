#pragma once

#include <Rcpp.h>
#include <cmath>
#include <limits>

namespace analogs {

// Weight code enums (consistent with R wrapper mapping)
enum class WeightCode : int {
      NONE            = 0,
            UNIFORM         = 1,
            INVERSE_CLIM    = 2,
            INVERSE_GEOG    = 3,
            GAUSSIAN_CLIM   = 4,
            GAUSSIAN_GEOG   = 5,
            GAUSSIAN_JOINT  = 6,
            INVERSE_JOINT   = 7
};

// Pre-compute weight parameters for efficiency
// Returns (param1, param2) for use in weight_from_codes
inline std::pair<double, double> precompute_weight_params(
            WeightCode wc,
            const Rcpp::NumericVector& theta_vec)
{
      double param1 = 0.0;
      double param2 = 0.0;

      switch (wc) {
      case WeightCode::UNIFORM:
      case WeightCode::NONE:
            break;

      case WeightCode::INVERSE_CLIM:
      case WeightCode::INVERSE_GEOG: {
            // param1 = epsilon
            if (theta_vec.size() > 0 && std::isfinite(theta_vec[0]) && theta_vec[0] > 0.0) {
            param1 = theta_vec[0];
      } else {
            // Default epsilon
            param1 = (wc == WeightCode::INVERSE_CLIM) ? 1e-12 : 1e-6;
      }
      break;
      }

      case WeightCode::GAUSSIAN_CLIM:
      case WeightCode::GAUSSIAN_GEOG: {
            // param1 = -1/(2*sigma^2)
            double sigma = 1.0;
            if (theta_vec.size() > 0 && std::isfinite(theta_vec[0]) && theta_vec[0] > 0.0) {
                  sigma = theta_vec[0];
            }
            param1 = -1.0 / (2.0 * sigma * sigma);
            break;
      }

      case WeightCode::GAUSSIAN_JOINT: {
            // param1 = -1/(2*sigma_clim^2), param2 = -1/(2*sigma_geog^2)
            double sigma_clim = 1.0;
            double sigma_geog = 1.0;
            if (theta_vec.size() >= 2) {
                  if (std::isfinite(theta_vec[0]) && theta_vec[0] > 0.0) {
                        sigma_clim = theta_vec[0];
                  }
                  if (std::isfinite(theta_vec[1]) && theta_vec[1] > 0.0) {
                        sigma_geog = theta_vec[1];
                  }
            }
            param1 = -1.0 / (2.0 * sigma_clim * sigma_clim);
            param2 = -1.0 / (2.0 * sigma_geog * sigma_geog);
            break;
      }

      case WeightCode::INVERSE_JOINT: {
            // param1 = eps_clim, param2 = eps_geog
            double eps_clim = 1e-12;
            double eps_geog = 1e-6;
            if (theta_vec.size() >= 2) {
                  if (std::isfinite(theta_vec[0]) && theta_vec[0] > 0.0) {
                        eps_clim = theta_vec[0];
                  }
                  if (std::isfinite(theta_vec[1]) && theta_vec[1] > 0.0) {
                        eps_geog = theta_vec[1];
                  }
            }
            param1 = eps_clim;
            param2 = eps_geog;
            break;
      }
      }

      return std::make_pair(param1, param2);
}

// Weight calculation from codes with optimized parameters
// For efficiency, pre-computed parameters are passed in
inline double weight_from_codes(
            WeightCode wc,
            double clim_dist,
            double geog_dist,
            double param1,
            double param2)
{
      switch (wc) {
      case WeightCode::UNIFORM:
            return 1.0;

      case WeightCode::INVERSE_CLIM: {
            // param1 = epsilon
            return 1.0 / (clim_dist + param1);
      }

      case WeightCode::INVERSE_GEOG: {
            // param1 = epsilon
            return 1.0 / (geog_dist + param1);
      }

      case WeightCode::GAUSSIAN_CLIM: {
            // param1 = -1/(2*sigma^2)
            const double d2 = clim_dist * clim_dist;
            return std::exp(param1 * d2);
      }

      case WeightCode::GAUSSIAN_GEOG: {
            // param1 = -1/(2*sigma^2)
            const double d2 = geog_dist * geog_dist;
            return std::exp(param1 * d2);
      }

      case WeightCode::GAUSSIAN_JOINT: {
            // param1 = -1/(2*sigma_clim^2), param2 = -1/(2*sigma_geog^2)
            const double d2_clim = clim_dist * clim_dist;
            const double d2_geog = geog_dist * geog_dist;
            return std::exp(param1 * d2_clim + param2 * d2_geog);
      }

      case WeightCode::INVERSE_JOINT: {
            // param1 = eps_clim, param2 = eps_geog
            const double d_clim_adj = clim_dist + param1;
            const double d_geog_adj = geog_dist + param2;
            const double joint_dist = std::sqrt(d_clim_adj * d_clim_adj +
                                                d_geog_adj * d_geog_adj);
            return 1.0 / joint_dist;
      }

      default: // NONE or unknown
            return 1.0;
      }
}

} // namespace analogs
