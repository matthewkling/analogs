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

// ---- Lattice index (uniform-width bins for now) ---------------------------
// NOTE: Supports best-first traversal using AABB lower bounds. When we move to
// quantile (uneven) bins, we'll store per-bin edges; traversal logic remains.
struct LatticeIndex {
      MatrixView ref_clim; // [n_ref x n_vars]
      MatrixView ref_geo;  // [n_ref x 2]
      size_tu    nref = 0, nvars = 0;

      // Binning parameters (uniform for now)
      std::vector<double> mins;         // per-var min of ref_clim
      std::vector<double> bin_widths;   // per-var bin widths (Inf collapses dimension)
      std::vector<double> inv_bw;       // cached 1/bin_widths (0 if Inf)

      // CSR payload
      std::vector<index_t> bin_ids_flat; // concatenated ids across all bins
      std::vector<size_tu> bin_offsets;  // size = n_bins+1
      std::vector<BinCoord> bin_coords;  // one coord per bin (same order as offsets)

      // Coordinate -> bin index (will be replaced by sort+binary-search later)
      std::unordered_map<BinCoord, size_tu, BinCoordHash> coord_to_bin;

      LatticeIndex(MatrixView ref_clim_, MatrixView ref_geo_, const double* bin_w, size_tu nvars_);
      void build();

      // Legacy rectangular window filter (unordered)
      void hard_filter(const double* focal_vars,
                       const ClimateConstraint& clim,
                       std::vector<index_t>& out_ids) const;

      // Best-first search for top-k or radius filtering using AABB lower bounds
      void best_first_search(const double* focal_vars,
                             const ClimateConstraint& clim,
                             std::vector<index_t>& out_ids,
                             index_t k) const;

      // Back-compat API (implemented in .cpp; delegates to best-first)
      void expanding_search(const double* focal_vars,
                            const ClimateConstraint& clim,
                            std::vector<index_t>& out_ids,
                            index_t k,
                            size_tu max_bins_to_visit = 0) const;

private:
      inline int32_t bin_index(double value, size_tu k) const {
            if (inv_bw[k] == 0.0) return 0; // collapsed dimension
            double t = (value - mins[k]) * inv_bw[k];
            long long bi = static_cast<long long>(std::floor(t));
            if (bi < INT32_MIN) bi = INT32_MIN;
            if (bi > INT32_MAX) bi = INT32_MAX;
            return static_cast<int32_t>(bi);
      }

      // AABB lower-bound distance from focal to a bin (uniform-width bins).
      double aabb_lower_bound(const double* focal_vars, const BinCoord& coord) const;

      // Legacy window enumeration (uniform bins only)
      void bins_for_window(const double* focal_vars,
                           const ClimateConstraint& clim,
                           std::vector<size_tu>& out_bin_ids) const;
};



} // namespace analogs
