#pragma once

#include <Rcpp.h>
#include <cmath>
#include <limits>

namespace analogs {

// Per-family kernel shape. The overall weight for a candidate is the PRODUCT
// of the climate-family weight and the geographic-family weight, each computed
// independently from its own distance and shape. This product model replaces
// the old fused WeightCode enum: it reproduces every separable kernel (a
// single-family kernel is just the other family = UNIFORM) and additionally
// supports mixed shapes (e.g. INVERSE climate x GAUSSIAN geography). The old
// non-separable INVERSE_JOINT (a coupled 1/||(c,g)|| norm) is intentionally
// dropped; INVERSE x INVERSE now gives 1/((1+c/tc)(1+g/tg)), a clean product.
enum class FamilyKernel : int {
      UNIFORM  = 0,   // constant weight 1 (also the "no weighting" case)
            GAUSSIAN = 1,   // exp(-d^2 / (2 theta^2))
            INVERSE  = 2    // 1 / (1 + d / theta)   [reparameterized: theta = half-weight scale]
};

// Pre-compute the single per-family weight parameter from theta, so the hot
// loop avoids repeated divides. Returns:
//   GAUSSIAN -> -1/(2 theta^2)   (multiply by d^2 inside exp)
//   INVERSE  ->  1/theta         (multiply by d, add 1, reciprocate)
//   UNIFORM  ->  0.0             (unused)
// theta defaults: 1.0 if not finite/positive. For INVERSE, theta is the
// half-weight distance (weight = 1/2 at d = theta), bounded at d=0 (weight 1)
// with no epsilon needed.
inline double precompute_family_param(FamilyKernel k, double theta) {
      switch (k) {
      case FamilyKernel::GAUSSIAN: {
            double sigma = (std::isfinite(theta) && theta > 0.0) ? theta : 1.0;
            return -1.0 / (2.0 * sigma * sigma);
      }
      case FamilyKernel::INVERSE: {
            double scale = (std::isfinite(theta) && theta > 0.0) ? theta : 1.0;
            return 1.0 / scale;
      }
      case FamilyKernel::UNIFORM:
      default:
            return 0.0;
      }
}

// Compute one family's weight from its shape, distance, and precomputed param.
// UNIFORM short-circuits to 1.0 (skipping exp/divide), so a query that only
// weights one family pays nothing for the other.
inline double per_family_weight(FamilyKernel k, double dist, double wparam) {
      switch (k) {
      case FamilyKernel::GAUSSIAN:
            // wparam = -1/(2 sigma^2)
            return std::exp(wparam * dist * dist);
      case FamilyKernel::INVERSE:
            // wparam = 1/theta;  weight = 1 / (1 + d/theta)
            return 1.0 / (1.0 + dist * wparam);
      case FamilyKernel::UNIFORM:
      default:
            return 1.0;
      }
}

// Combined weight = product of the two per-family weights. This is the single
// entry point used by the workers.
inline double weight_from_families(FamilyKernel clim_k, double clim_dist, double clim_wparam,
                                   FamilyKernel geog_k, double geog_dist, double geog_wparam) {
      return per_family_weight(clim_k, clim_dist, clim_wparam) *
            per_family_weight(geog_k, geog_dist, geog_wparam);
}

} // namespace analogs
