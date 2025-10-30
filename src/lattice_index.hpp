#pragma once
#include "types.hpp"
#include "constraints.hpp"
#include <vector>
#include <cstdint>
#include <unordered_map>
#include <cmath>

namespace analogs {

// ---- Bin coordinate key ---------------------------------------------------
struct BinCoord {
      std::vector<int32_t> idx;              // length = nvars
      bool operator==(const BinCoord& o) const noexcept { return idx == o.idx; }
};

struct BinCoordHash {
      std::size_t operator()(const BinCoord& k) const noexcept {
            std::size_t h = 1469598103934665603ull;
            for (size_t i = 0; i < k.idx.size(); ++i) {
                  std::size_t v = static_cast<std::size_t>(static_cast<uint32_t>(k.idx[i]));
                  h ^= v + 0x9e3779b97f4a7c15ull + (h<<6) + (h>>2);
            }
            return h;
      }
};

// ---- Lattice index with quantile-based binning ----------------------------
// OPTIMIZATION: Uses squared distances throughout to avoid expensive sqrt operations
struct LatticeIndex {
      MatrixView ref_clim; // [n_ref x n_vars]
      MatrixView ref_geo;  // [n_ref x 2]
      size_tu    nref = 0, nvars = 0;

      // Quantile-based binning parameters
      size_tu target_occupancy = 20;  // target points per bin
      std::vector<size_tu> n_bins_per_dim;  // computed adaptively
      std::vector<std::vector<double>> quantile_edges;  // [nvars][n_bins+1]

      // CSR payload
      std::vector<index_t> bin_ids_flat; // concatenated ids across all bins
      std::vector<size_tu> bin_offsets;  // size = n_bins+1
      std::vector<BinCoord> bin_coords;  // one coord per bin (same order as offsets)

      // Coordinate -> bin index for fast lookup
      std::unordered_map<BinCoord, size_tu, BinCoordHash> coord_to_bin;

      // Constructor
      LatticeIndex(MatrixView ref_clim_, MatrixView ref_geo_, size_tu nvars_,
                   size_tu target_occupancy_ = 20);

      // Build index with quantile binning
      void build();

      // Legacy rectangular window filter (unordered)
      void hard_filter(const double* focal_vars,
                       const ClimateConstraint& clim,
                       std::vector<index_t>& out_ids) const;

      // Best-first search for top-k or radius filtering using AABB lower bounds
      // OPTIMIZED: Uses squared distances internally
      void best_first_search(const double* focal_vars,
                             const ClimateConstraint& clim,
                             std::vector<index_t>& out_ids,
                             index_t k) const;

      // Back-compat API (delegates to best-first)
      void expanding_search(const double* focal_vars,
                            const ClimateConstraint& clim,
                            std::vector<index_t>& out_ids,
                            index_t k,
                            size_tu max_bins_to_visit = 0) const;

      // Diagnostic information
      void get_diagnostics(size_tu& total_bins, double& avg_occupancy,
                           double& min_occupancy, double& max_occupancy) const;

private:
      // Compute adaptive number of bins per dimension
      size_tu compute_n_bins() const;

      // Compute quantile edges for a single dimension
      void compute_quantile_edges(size_tu dim_idx, size_tu n_bins);

      // Assign a value to a bin using binary search on quantile edges
      inline int32_t bin_index(double value, size_tu dim_idx) const {
            const auto& edges = quantile_edges[dim_idx];
            if (edges.empty()) return 0;

            // Binary search to find the bin
            auto it = std::upper_bound(edges.begin(), edges.end(), value);
            if (it == edges.begin()) return 0;
            if (it == edges.end()) return static_cast<int32_t>(edges.size() - 2);

            int32_t bin = static_cast<int32_t>(std::distance(edges.begin(), it) - 1);
            return bin;
      }

      // AABB lower-bound SQUARED distance from focal to a bin
      // OPTIMIZATION: Returns squared distance to avoid sqrt
      inline double aabb_lower_bound_squared(const double* focal_vars,
                                             const BinCoord& coord) const {
            double sum_sq = 0.0;

            for (size_tu k = 0; k < nvars; ++k) {
                  int32_t bin_idx = coord.idx[k];

                  if (bin_idx < 0 || bin_idx >= static_cast<int32_t>(quantile_edges[k].size() - 1)) {
                        continue;
                  }

                  const double low = quantile_edges[k][bin_idx];
                  const double high = quantile_edges[k][bin_idx + 1];
                  const double fk = focal_vars[k];

                  double d;
                  if (fk < low) {
                        d = low - fk;
                  } else if (fk > high) {
                        d = fk - high;
                  } else {
                        d = 0.0;
                  }

                  sum_sq += d * d;
            }

            return sum_sq;  // Return squared distance (no sqrt!)
      }

      // Legacy window enumeration (for hard_filter)
      void bins_for_window(const double* focal_vars,
                           const ClimateConstraint& clim,
                           std::vector<size_tu>& out_bin_ids) const;
};

} // namespace analogs
