#pragma once

#include <cmath>
#include <limits>
#include <vector>
#include <utility>

namespace analogs {

// Selection code enums (R->C++ mapping for select parameter)
enum class SelectCode : int {
      KNN_CLIM = 0,
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
            WEIGHTED_MEAN  = 7
};

// Compute Euclidean climate distance (scalar) and/or per-var checks.
// Returns (ok, clim_dist) where ok means thresholds satisfied.
inline std::pair<bool, double>
      clim_ok_and_dist(const double* f_clim_col,
                       const double* r_clim_col,
                       int n_clim,
                       int stride_f,
                       int stride_r,
                       bool use_pervar_clim,
                       const std::vector<double>& max_clim_pervar,
                       bool use_scalar_clim,
                       double max_clim_scalar)
      {
            double sumsq = 0.0;

            if (use_pervar_clim) {
                  for (int k = 0; k < n_clim; ++k) {
                        const double df = f_clim_col[k * stride_f] - r_clim_col[k * stride_r];
                        if (std::fabs(df) > max_clim_pervar[k]) {
                              return std::make_pair(false, 0.0);
                        }
                        sumsq += df * df;
                  }
                  return std::make_pair(true, std::sqrt(sumsq));
            } else {
                  // scalar threshold or just distance
                  for (int k = 0; k < n_clim; ++k) {
                        const double df = f_clim_col[k * stride_f] - r_clim_col[k * stride_r];
                        sumsq += df * df;
                  }
                  const double d = std::sqrt(sumsq);
                  if (use_scalar_clim && d > max_clim_scalar) {
                        return std::make_pair(false, d);
                  }
                  return std::make_pair(true, d);
            }
      }

} // namespace analogs
