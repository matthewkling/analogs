#pragma once

#include <cmath>
#include <limits>
#include <vector>
#include <utility>

namespace analogs {

// Selection code enums (R->C++ mapping for select parameter)
enum class SelectCode : int {
      KNN_ENV = 0,
            KNN_GEOG = 1,
            ALL      = 2
};

// Aggregation code enums (R->C++ mapping for stat parameter)
enum class AggregateCode : int {
      NONE           = 0,
            COUNT          = 1,
            SUM_WEIGHTS    = 2,
            MEAN_WEIGHTS   = 3,
            SUM            = 4,
            MEAN           = 5,
            WEIGHTED_SUM   = 6,
            WEIGHTED_MEAN  = 7,
            ESS            = 8,
            REGRESSION     = 9,
            TABULATE       = 10
};

// Compute Euclidean environment distance (scalar) and/or per-var checks.
// Returns (ok, env_dist) where ok means thresholds satisfied.
inline std::pair<bool, double>
      env_ok_and_dist(const double* f_env_col,
                       const double* r_env_col,
                       int n_env,
                       int stride_f,
                       int stride_r,
                       bool use_pervar_env,
                       const std::vector<double>& max_env_pervar,
                       bool use_scalar_env,
                       double max_env_scalar)
      {
            double sumsq = 0.0;

            if (use_pervar_env) {
                  for (int k = 0; k < n_env; ++k) {
                        const double df = f_env_col[k * stride_f] - r_env_col[k * stride_r];
                        if (std::fabs(df) > max_env_pervar[k]) {
                              return std::make_pair(false, 0.0);
                        }
                        sumsq += df * df;
                  }
                  return std::make_pair(true, std::sqrt(sumsq));
            } else {
                  // scalar threshold or just distance
                  for (int k = 0; k < n_env; ++k) {
                        const double df = f_env_col[k * stride_f] - r_env_col[k * stride_r];
                        sumsq += df * df;
                  }
                  const double d = std::sqrt(sumsq);
                  if (use_scalar_env && d > max_env_scalar) {
                        return std::make_pair(false, d);
                  }
                  return std::make_pair(true, d);
            }
      }

} // namespace analogs
