#pragma once

#include <cmath>
#include <limits>
#include <vector>
#include <utility>

namespace analogs {

// Mode code enums (consistent with R wrapper mapping)
enum class ModeCode : int {
      KNN_CLIM = 0,
            KNN_GEOG = 1,
            COUNT    = 2,
            SUM      = 3,
            MEAN     = 4,
            ALL      = 5
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
