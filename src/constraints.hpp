#pragma once
#include "types.hpp"
#include <cmath>

namespace analogs {

// Climate constraint supports scalar Euclidean radius (preferred) or per-var bands.
struct ClimateConstraint {
      const double* max_abs_diff = nullptr; // length nvars when in vector mode
      double radius = INFINITY;             // scalar Euclidean radius
      size_tu nvars = 0;
      bool use_scalar = true;               // true => Euclidean(a,b) <= radius

      inline bool pass(const double* a, const double* b) const {
            if (use_scalar) {
                  double s = 0.0;
                  for (size_tu k=0; k<nvars; ++k) { double d = a[k]-b[k]; s += d*d; }
                  return std::sqrt(s) <= radius;
            } else {
                  for (size_tu k=0; k<nvars; ++k) {
                        if (std::fabs(a[k]-b[k]) > max_abs_diff[k]) return false;
                  }
                  return true;
            }
      }

      // Halfwidth used for bin-window enumeration; valid for both modes.
      inline double window_halfwidth(size_tu k) const {
            (void)k; // unused in scalar mode
            return use_scalar ? radius : max_abs_diff[k];
      }
};

struct GeoRadius { double max_km = INFINITY; };

} // namespace analogs
