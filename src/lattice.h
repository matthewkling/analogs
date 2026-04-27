#pragma once

#include "types.h"

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <array>
#include <cmath>
#include <limits>
#include <random>
#include <algorithm>


namespace analogs {

// Bin structure: stores point IDs and sampling weight
struct Bin {
      std::vector<index_t> point_ids;
      double sample_weight;  // 1.0 = no downsampling, >1.0 = each point represents multiple

      Bin() : sample_weight(1.0) {}
};

// Full-dimensional regular lattice over
//   dims = [x, y, clim1, ..., climP].
//
// - Uses 2 geographic dims (x, y) and all climate dims.
// - Bins are regular in each dimension.
// - Cells are stored sparsely in an unordered_map keyed by a flattened index.
// - Used only for *candidate generation*; all exact geo/clim tests
//   still happen in core.cpp.
class Lattice {
public:
      MetricType metric_type;

      size_tu n_points;     // number of reference points (before downsampling)
      size_tu n_geo_dims;   // number of geographic variables (2 for planar; 2 or 3 in general)
      size_tu n_clim_dims;  // number of climate variables
      size_tu n_dims;       // = n_geo_dims + n_clim_dims

      // Per-dimension metadata
      std::vector<double> mins;    // min value per dimension
      std::vector<double> maxs;    // max value per dimension
      std::vector<double> res;     // bin width per dimension
      std::vector<size_tu> n_bins; // number of bins per dimension
      std::vector<size_tu> strides; // flattening strides per dimension

      // Sparse bins: key (flattened bin index) -> Bin struct
      std::unordered_map<size_tu, Bin> bins;

      // Diagnostics
      size_tu total_bins;        // product of n_bins[d]
      size_tu min_cell_occ;      // min occupancy among occupied cells (before downsampling)
      size_tu max_cell_occ;      // max occupancy among occupied cells (before downsampling)
      size_tu n_cells_nonempty;  // number of occupied cells
      double downsample_actual;  // actual downsampling rate achieved (points kept / total points)

      Lattice()
            : metric_type(MetricType::Planar),
              n_points(0),
              n_geo_dims(2),
              n_clim_dims(0),
              n_dims(2),
              total_bins(1),
              min_cell_occ(std::numeric_limits<size_tu>::max()),
              max_cell_occ(0),
              n_cells_nonempty(0),
              downsample_actual(1.0) {}

      // Build lattice over all dims with optional downsampling
      void build(const double* ref_ptr,
                 size_tu n_ref,
                 size_tu n_geo,
                 size_tu n_clim,
                 size_tu stride_r,
                 MetricType metric,
                 double max_dist,
                 bool use_scalar_clim,
                 const std::vector<double>& max_clim_pervar,
                 double max_clim_scalar,
                 int bins_per_dim,
                 bool use_geo_lattice,
                 bool use_clim_lattice,
                 double downsample_rate = 1.0,
                 unsigned int seed = 0);

      // Get candidates within [lo, hi] bin ranges
      void get_candidates(const std::vector<size_tu>& lo,
                          const std::vector<size_tu>& hi,
                          std::vector<index_t>& out_indices,
                          std::vector<double>& out_weights) const;

      // Query candidate indices for a focal point
      void query(const double* focal_geo,
                 const double* focal_clim,
                 double max_dist,
                 bool use_scalar_clim,
                 const std::vector<double>& max_clim_pervar,
                 double max_clim_scalar,
                 std::vector<index_t>& out_indices,
                 std::vector<double>& out_weights) const;

      // KNN query with priority-queue based search
      // May contain fewer than k if not enough admissible refs.
      void knn_query(const double* focal_geo,
                     const double* focal_clim,
                     const double* ref_ptr,
                     size_tu stride_r,
                     size_tu n_clim,
                     bool rank_by_geog,
                     double max_geog,
                     bool use_scalar_clim,
                     const std::vector<double>& max_clim_pervar,
                     double max_clim_scalar,
                     int k,
                     std::vector<index_t>& out_indices,
                     std::vector<double>& out_weights) const;


private:
      void enumerate_cells(const std::vector<size_tu>& lo,
                           const std::vector<size_tu>& hi,
                           std::vector<size_tu>& idx,
                           size_tu dim,
                           std::vector<index_t>& out_indices,
                           std::vector<double>& out_weights) const;

private:
      inline size_tu flatten_idx(const std::array<size_tu, 8>& idx) const {
            size_tu key = 0;
            for (size_tu d = 0; d < n_dims; ++d) {
                  key += idx[d] * strides[d];
            }
            return key;
      }

};

// ---- Implementation -------------------------------------------------------

inline void Lattice::build(const double* ref_ptr,
                           size_tu n_ref,
                           size_tu n_geo,
                           size_tu n_clim,
                           size_tu stride_r,
                           MetricType metric,
                           double max_dist,
                           bool use_scalar_clim,
                           const std::vector<double>& max_clim_pervar,
                           double max_clim_scalar,
                           int    bins_per_dim,
                           bool   use_geo_lattice,
                           bool   use_clim_lattice,
                           double downsample_rate,
                           unsigned int seed) {
      metric_type    = metric;
      n_points       = n_ref;
      n_geo_dims     = n_geo;
      n_clim_dims    = n_clim;
      n_dims         = n_geo_dims + n_clim_dims;

      mins.assign(n_dims, 0.0);
      maxs.assign(n_dims, 0.0);
      res.assign(n_dims, 1.0);
      n_bins.assign(n_dims, 1);
      strides.assign(n_dims, 1);
      bins.clear();
      total_bins    = 1;
      min_cell_occ  = std::numeric_limits<size_tu>::max();
      max_cell_occ  = 0;
      n_cells_nonempty  = 0;
      downsample_actual = 1.0;

      if (n_ref == 0) {
            return;
      }

      // First pass: compute mins/maxs over all dims.
      for (size_tu j = 0; j < n_ref; ++j) {
            // geo dims
            for (size_tu g = 0; g < n_geo_dims; ++g) {
                  double v = ref_ptr[j + g * stride_r];
                  if (j == 0) {
                        mins[g] = maxs[g] = v;
                  } else {
                        if (v < mins[g]) mins[g] = v;
                        if (v > maxs[g]) maxs[g] = v;
                  }
            }

            // climate dims
            for (size_tu k = 0; k < n_clim_dims; ++k) {
                  size_tu d = n_geo_dims + k;
                  double v = ref_ptr[j + d * stride_r];
                  if (j == 0) {
                        mins[d] = maxs[d] = v;
                  } else {
                        if (v < mins[d]) mins[d] = v;
                        if (v > maxs[d]) maxs[d] = v;
                  }
            }
      }

      // Choose per-dimension bin widths and bin counts.
      //
      // - bins_per_dim > 0: requested resolution for *active* dims
      // - inactive dims get a single bin (no partitioning)
      size_tu B = (bins_per_dim > 0)
            ? static_cast<size_tu>(bins_per_dim)
                  : static_cast<size_tu>(10);

      for (size_tu d = 0; d < n_dims; ++d) {
            double span = maxs[d] - mins[d];

            bool active = false;
            if (d < n_geo_dims) {
                  active = use_geo_lattice;
            } else {
                  active = use_clim_lattice;
            }

            if (!active) {
                  // Dimension not used in lattice: single bin.
                  n_bins[d] = 1;
                  res[d]    = (span > 0.0) ? span : 1.0;
            } else {
                  // Active dimension: equal-width bins_per_dim.
                  n_bins[d] = (B > 0) ? B : 1;
                  res[d]    = (span > 0.0 && n_bins[d] > 0)
                        ? (span / static_cast<double>(n_bins[d]))
                        : 1.0;
            }
      }

      // Compute strides and total number of bins.
      strides.assign(n_dims, 0);
      if (n_dims > 0) {
            strides[n_dims - 1] = 1;
            for (int d = static_cast<int>(n_dims) - 2; d >= 0; --d) {
                  strides[d] = strides[d + 1] * n_bins[d + 1];
            }
            total_bins = strides[0] * n_bins[0];
      } else {
            total_bins = 0;
      }


      // Second pass: assign points to bins.
      std::vector<size_tu> idx(n_dims);

      for (size_tu j = 0; j < n_ref; ++j) {
            // geo dims
            for (size_tu g = 0; g < n_geo_dims; ++g) {
                  double v = ref_ptr[j + g * stride_r];
                  double pos = (v - mins[g]) / res[g];
                  if (pos < 0.0) pos = 0.0;
                  size_tu ib = static_cast<size_tu>(pos);
                  if (ib >= n_bins[g]) ib = n_bins[g] - 1;
                  idx[g] = ib;
            }

            // climate dims
            for (size_tu k = 0; k < n_clim_dims; ++k) {
                  size_tu d = n_geo_dims + k;
                  double v = ref_ptr[j + d * stride_r];
                  double pos = (v - mins[d]) / res[d];
                  if (pos < 0.0) pos = 0.0;
                  size_tu ib = static_cast<size_tu>(pos);
                  if (ib >= n_bins[d]) ib = n_bins[d] - 1;
                  idx[d] = ib;
            }

            // Flatten multi-index to key.
            size_tu key = 0;
            for (size_tu d = 0; d < n_dims; ++d) {
                  key += idx[d] * strides[d];
            }

            // Insert point into bin (store 0-based index like original code)
            bins[key].point_ids.push_back(static_cast<index_t>(j));
      }

      // Compute diagnostics before downsampling
      for (const auto& pair : bins) {
            size_tu occ = pair.second.point_ids.size();
            if (occ < min_cell_occ) min_cell_occ = occ;
            if (occ > max_cell_occ) max_cell_occ = occ;
      }
      n_cells_nonempty = bins.size();

      // Apply downsampling if requested
      if (downsample_rate < 1.0) {
            // Create RNG for randomization
            std::mt19937 rng(seed);

            // Target total points to keep
            size_tu target_total = static_cast<size_tu>(std::ceil(n_ref * downsample_rate));

            // Collect all bin sizes
            std::vector<size_tu> bin_sizes;
            bin_sizes.reserve(bins.size());
            for (const auto& pair : bins) {
                  bin_sizes.push_back(pair.second.point_ids.size());
            }

            // Sort descending for efficiency
            std::sort(bin_sizes.begin(), bin_sizes.end(), std::greater<size_tu>());

            // Binary search for optimal target occupancy
            // Note: if bins.size() > target_total, we'll keep 1 per bin (spatial coverage priority)
            size_tu lo = 1;
            size_tu hi = bin_sizes.empty() ? 1 : bin_sizes[0];
            size_tu target_occ = 1;

            while (lo <= hi) {
                  size_tu mid = lo + (hi - lo) / 2;

                  // Calculate how many points we'd keep with this target
                  size_tu kept = 0;
                  for (size_tu size : bin_sizes) {
                        kept += std::min(size, mid);
                  }

                  if (kept <= target_total) {
                        // Can increase target (this is valid, try to find larger)
                        target_occ = mid;
                        lo = mid + 1;
                  } else {
                        // Would keep too many, need to decrease target
                        hi = mid - 1;
                  }
            }

            // Ensure at least 1 (maintain spatial coverage - don't drop entire bins)
            target_occ = std::max(static_cast<size_tu>(1), target_occ);

            // Apply downsampling to bins
            size_tu total_kept = 0;

            for (auto& pair : bins) {
                  Bin& bin = pair.second;
                  size_tu original_size = bin.point_ids.size();

                  if (original_size > target_occ) {
                        // Shuffle and keep first target_occ points
                        std::shuffle(bin.point_ids.begin(), bin.point_ids.end(), rng);
                        bin.point_ids.resize(target_occ);
                        bin.sample_weight = static_cast<double>(original_size) / target_occ;
                  } else {
                        // Keep all points in this bin (no downsampling at this location)
                        bin.sample_weight = 1.0;
                  }

                  total_kept += bin.point_ids.size();
            }

            // Record actual downsampling rate achieved
            downsample_actual = static_cast<double>(total_kept) / n_ref;

            // Recompute occupancy statistics after downsampling
            min_cell_occ = std::numeric_limits<size_tu>::max();
            max_cell_occ = 0;
            for (const auto& pair : bins) {
                  size_tu occ = pair.second.point_ids.size();
                  if (occ < min_cell_occ) min_cell_occ = occ;
                  if (occ > max_cell_occ) max_cell_occ = occ;
            }
            if (bins.empty()) {
                  min_cell_occ = 0;
            }
      } else {
            // No downsampling
            downsample_actual = 1.0;
      }
}

inline void Lattice::get_candidates(const std::vector<size_tu>& lo,
                                    const std::vector<size_tu>& hi,
                                    std::vector<index_t>& out_indices,
                                    std::vector<double>& out_weights) const {
      std::vector<size_tu> idx(n_dims, 0);
      enumerate_cells(lo, hi, idx, 0, out_indices, out_weights);
}

inline void Lattice::query(const double* focal_geo,
                           const double* focal_clim,
                           double max_dist,
                           bool use_scalar_clim,
                           const std::vector<double>& max_clim_pervar,
                           double max_clim_scalar,
                           std::vector<index_t>& out_indices,
                           std::vector<double>& out_weights) const {
      out_indices.clear();
      out_weights.clear();

      if (n_points == 0 || n_dims == 0) return;

      std::vector<size_tu> lo(n_dims);
      std::vector<size_tu> hi(n_dims);

      // Per-dimension allowed bin index ranges
      for (size_tu d = 0; d < n_dims; ++d) {
            double minv = mins[d];
            double maxv = maxs[d];

            if (d < n_geo_dims) {
                  // Geo dims: bound by max_dist if finite
                  if (std::isfinite(max_dist) && max_dist > 0.0) {
                        double q = focal_geo[d];
                        minv = q - max_dist;
                        maxv = q + max_dist;
                  }
            } else {
                  // Climate dims: bound by per-var or scalar climate thresholds
                  size_tu clim_idx = d - n_geo_dims;
                  double q = focal_clim[clim_idx];

                  if (!max_clim_pervar.empty() && clim_idx < max_clim_pervar.size()) {
                        double thr = max_clim_pervar[clim_idx];
                        if (std::isfinite(thr) && thr > 0.0) {
                              minv = q - thr;
                              maxv = q + thr;
                        }
                  }

                  if (use_scalar_clim && std::isfinite(max_clim_scalar) && max_clim_scalar > 0.0) {
                        double thr2 = max_clim_scalar;
                        minv = (minv == mins[d]) ? (q - thr2) : std::max(minv, q - thr2);
                        maxv = (maxv == maxs[d]) ? (q + thr2) : std::min(maxv, q + thr2);
                  }
            }

            // Convert to bin indices
            if (minv < mins[d]) minv = mins[d];
            if (maxv > maxs[d]) maxv = maxs[d];

            double pos_lo = (minv - mins[d]) / res[d];
            double pos_hi = (maxv - mins[d]) / res[d];

            if (pos_lo < 0.0) pos_lo = 0.0;
            if (pos_hi < 0.0) pos_hi = 0.0;

            size_tu ib_lo = static_cast<size_tu>(pos_lo);
            size_tu ib_hi = static_cast<size_tu>(pos_hi);

            if (ib_lo >= n_bins[d]) ib_lo = n_bins[d] - 1;
            if (ib_hi >= n_bins[d]) ib_hi = n_bins[d] - 1;

            lo[d] = ib_lo;
            hi[d] = ib_hi;
      }

      get_candidates(lo, hi, out_indices, out_weights);
}

inline void Lattice::enumerate_cells(const std::vector<size_tu>& lo,
                                     const std::vector<size_tu>& hi,
                                     std::vector<size_tu>& idx,
                                     size_tu dim,
                                     std::vector<index_t>& out_indices,
                                     std::vector<double>& out_weights) const {
      if (dim == n_dims) {
            // Compute flattened key
            size_tu key = 0;
            for (size_tu d = 0; d < n_dims; ++d) {
                  key += idx[d] * strides[d];
            }

            auto it = bins.find(key);
            if (it != bins.end()) {
                  const Bin& bin = it->second;
                  double weight = bin.sample_weight;
                  for (index_t pid : bin.point_ids) {
                        out_indices.push_back(pid);
                        out_weights.push_back(weight);
                  }
            }
            return;
      }

      for (size_tu i = lo[dim]; i <= hi[dim]; ++i) {
            idx[dim] = i;
            enumerate_cells(lo, hi, idx, dim + 1, out_indices, out_weights);
      }
}

// Forward declarations for knn_query helper functions
inline double lb_geo_projected(const Lattice& lat,
                               const std::array<size_tu, 8>& idx,
                               const double* focal_geo);

inline double lb_geo_chord3d(const Lattice& lat,
                             const std::array<size_tu, 8>& idx,
                             const double* focal_geo);

inline double lb_clim(const Lattice& lat,
                      const std::array<size_tu, 8>& idx,
                      const double* focal_clim);

inline bool clim_ok_and_dist_knn(const double* focal_clim,
                                 const double* ref_clim_col,
                                 size_tu n_clim,
                                 size_tu stride_r,
                                 bool use_scalar_clim,
                                 const std::vector<double>& max_clim_pervar,
                                 double max_clim_scalar,
                                 double& clim_dist_out,
                                 bool compute_dist);

// KNN query implementation with weight tracking
inline void Lattice::knn_query(const double* focal_geo,
                               const double* focal_clim,
                               const double* ref_ptr,
                               size_tu stride_r,
                               size_tu n_clim,
                               bool rank_by_geog,
                               double max_geog,
                               bool use_scalar_clim,
                               const std::vector<double>& max_clim_pervar,
                               double max_clim_scalar,
                               int k,
                               std::vector<index_t>& out_indices,
                               std::vector<double>& out_weights) const {
      out_indices.clear();
      out_weights.clear();

      if (k <= 0 || n_points == 0 || n_dims == 0) return;

      // Only implement expanding search for Euclidean metrics
      if (rank_by_geog && !(metric_type == MetricType::Planar || metric_type == MetricType::Chord3D)) {
            return;
      }

      if (n_clim != n_clim_dims) {
            return;
      }

      const size_tu D = n_dims;
      const size_tu n_geo = n_geo_dims;

      // Compute threshold bin ranges
      std::vector<size_tu> lo(D), hi(D);
      for (size_tu d = 0; d < D; ++d) {
            double q = (d < n_geo) ? focal_geo[d] : focal_clim[d - n_geo];

            double minv = mins[d];
            double maxv = maxs[d];

            if (d < n_geo) {
                  if (std::isfinite(max_geog) && max_geog > 0.0) {
                        minv = q - max_geog;
                        maxv = q + max_geog;
                  }
            } else {
                  size_tu cidx = d - n_geo;
                  if (!max_clim_pervar.empty() && cidx < max_clim_pervar.size()) {
                        double t = max_clim_pervar[cidx];
                        if (std::isfinite(t) && t > 0.0) {
                              minv = q - t;
                              maxv = q + t;
                        }
                  }
                  if (use_scalar_clim && std::isfinite(max_clim_scalar) && max_clim_scalar > 0.0) {
                        minv = (minv == mins[d]) ? (q - max_clim_scalar) : std::max(minv, q - max_clim_scalar);
                        maxv = (maxv == maxs[d]) ? (q + max_clim_scalar) : std::min(maxv, q + max_clim_scalar);
                  }
            }

            if (minv < mins[d]) minv = mins[d];
            if (maxv > maxs[d]) maxv = maxs[d];

            double pos_lo = (minv - mins[d]) / res[d];
            double pos_hi = (maxv - mins[d]) / res[d];
            if (pos_lo < 0.0) pos_lo = 0.0;
            if (pos_hi < 0.0) pos_hi = 0.0;

            size_tu ib_lo = static_cast<size_tu>(pos_lo);
            size_tu ib_hi = static_cast<size_tu>(pos_hi);
            if (ib_lo >= n_bins[d]) ib_lo = n_bins[d] - 1;
            if (ib_hi >= n_bins[d]) ib_hi = n_bins[d] - 1;

            lo[d] = ib_lo;
            hi[d] = ib_hi;
      }

      // Find center cell
      std::array<size_tu, 8> center;
      for (size_tu d = 0; d < D; ++d) {
            double q = (d < n_geo) ? focal_geo[d] : focal_clim[d - n_geo];
            double pos = (q - mins[d]) / res[d];
            if (pos < 0.0) pos = 0.0;
            size_tu ic = static_cast<size_tu>(pos);
            if (ic >= n_bins[d]) ic = n_bins[d] - 1;
            if (ic < lo[d]) ic = lo[d];
            if (ic > hi[d]) ic = hi[d];
            center[d] = ic;
      }

      // Cell lower bound function
      auto cell_lb = [&](const std::array<size_tu, 8>& idx) -> double {
            if (rank_by_geog) {
                  if (metric_type == MetricType::Planar) {
                        return lb_geo_projected(*this, idx, focal_geo);
                  } else if (metric_type == MetricType::Chord3D) {
                        return lb_geo_chord3d(*this, idx, focal_geo);
                  }
                  return 0.0;
            } else {
                  return lb_clim(*this, idx, focal_clim);
            }
      };

      struct CellState {
            double lb;
            std::array<size_tu, 8> idx;
      };

      struct CellCmp {
            bool operator()(const CellState& a, const CellState& b) const {
                  return a.lb > b.lb;
            }
      };

      std::priority_queue<CellState, std::vector<CellState>, CellCmp> pq;
      std::unordered_set<size_tu> visited;
      visited.reserve(128);

      CellState start;
      start.idx = center;
      start.lb = cell_lb(center);
      size_tu start_key = flatten_idx(center);
      visited.insert(start_key);
      pq.push(start);

      // kNN heap with weights
      struct Neighbor {
            double dist;
            index_t idx;
            double weight;

            Neighbor(double d, index_t i, double w) : dist(d), idx(i), weight(w) {}

            bool operator<(const Neighbor& other) const {
                  return dist < other.dist;  // max-heap: top has largest
            }
      };

      std::priority_queue<Neighbor> knn;

      const bool use_geo_constraint = (std::isfinite(max_geog) && max_geog > 0.0);
      const double max_geog2 = use_geo_constraint ? max_geog * max_geog : std::numeric_limits<double>::infinity();
      double d_k = std::numeric_limits<double>::infinity();

      // Expanding search
      while (!pq.empty()) {
            CellState current = pq.top();
            pq.pop();

            if (!knn.empty() && current.lb > d_k) break;

            size_tu key = flatten_idx(current.idx);
            auto it_cell = bins.find(key);

            if (it_cell != bins.end()) {
                  const Bin& bin = it_cell->second;
                  double bin_weight = bin.sample_weight;

                  for (index_t j : bin.point_ids) {
                        // j is already 0-based, use directly

                        // Geo distance
                        double gdist2 = 0.0;
                        double gdist = 0.0;
                        if (rank_by_geog || use_geo_constraint) {
                              double sumsq = 0.0;
                              for (size_tu g = 0; g < n_geo; ++g) {
                                    double rg = ref_ptr[j + g * stride_r];
                                    double dg = focal_geo[g] - rg;
                                    sumsq += dg * dg;
                              }
                              gdist2 = sumsq;
                              gdist = std::sqrt(gdist2);
                        }

                        if (use_geo_constraint && gdist2 > max_geog2) continue;

                        // Climate constraints
                        double clim_dist = 0.0;
                        bool compute_clim = !rank_by_geog;
                        if (!clim_ok_and_dist_knn(focal_clim,
                                                  ref_ptr + j + n_geo * stride_r,
                                                  n_clim, stride_r,
                                                  use_scalar_clim, max_clim_pervar, max_clim_scalar,
                                                  clim_dist, compute_clim)) {
                              continue;
                        }

                        // Ranking distance
                        double key_dist = rank_by_geog ? gdist : clim_dist;

                        if (static_cast<int>(knn.size()) < k) {
                              knn.emplace(key_dist, j, bin_weight);
                              if (static_cast<int>(knn.size()) == k) {
                                    d_k = knn.top().dist;
                              }
                        } else if (key_dist < knn.top().dist) {
                              knn.pop();
                              knn.emplace(key_dist, j, bin_weight);
                              d_k = knn.top().dist;
                        }
                  }
            }

            // Expand to neighbors
            std::array<size_tu, 8> base = current.idx;
            for (size_tu d = 0; d < D; ++d) {
                  for (int offset = -1; offset <= 1; offset += 2) {
                        std::array<size_tu, 8> nb_idx = base;
                        size_tu id = nb_idx[d];

                        if (offset < 0) {
                              if (id == 0 || id <= lo[d]) continue;
                              --id;
                        } else {
                              if (id + 1 >= n_bins[d] || id >= hi[d]) continue;
                              ++id;
                        }
                        nb_idx[d] = id;

                        size_tu nb_key = flatten_idx(nb_idx);
                        if (visited.find(nb_key) != visited.end()) continue;
                        visited.insert(nb_key);

                        CellState nb;
                        nb.idx = nb_idx;
                        nb.lb = cell_lb(nb_idx);
                        if (nb.lb <= d_k) {
                              pq.push(nb);
                        }
                  }
            }
      }

      // Extract results in ascending distance order
      int m = static_cast<int>(knn.size());
      out_indices.resize(m);
      out_weights.resize(m);

      for (int pos = m - 1; pos >= 0; --pos) {
            const Neighbor& nb = knn.top();
            out_indices[pos] = nb.idx;
            out_weights[pos] = nb.weight;
            knn.pop();
      }
}

// Helper functions for knn_query
inline double lb_geo_projected(const Lattice& lat,
                               const std::array<size_tu, 8>& idx,
                               const double* focal_geo) {
      const double fx = focal_geo[0];
      const double fy = focal_geo[1];

      const double x_min = lat.mins[0] + static_cast<double>(idx[0]) * lat.res[0];
      const double x_max = x_min + lat.res[0];
      const double y_min = lat.mins[1] + static_cast<double>(idx[1]) * lat.res[1];
      const double y_max = y_min + lat.res[1];

      double dx = 0.0;
      if (fx < x_min) dx = x_min - fx;
      else if (fx > x_max) dx = fx - x_max;

      double dy = 0.0;
      if (fy < y_min) dy = y_min - fy;
      else if (fy > y_max) dy = fy - y_max;

      return std::sqrt(dx * dx + dy * dy);
}

inline double lb_geo_chord3d(const Lattice& lat,
                             const std::array<size_tu, 8>& idx,
                             const double* focal_geo) {
      double sumsq = 0.0;
      for (size_tu g = 0; g < lat.n_geo_dims; ++g) {
            const double fg = focal_geo[g];
            const double c_min = lat.mins[g] + static_cast<double>(idx[g]) * lat.res[g];
            const double c_max = c_min + lat.res[g];

            double d = 0.0;
            if (fg < c_min) d = c_min - fg;
            else if (fg > c_max) d = fg - c_max;

            sumsq += d * d;
      }
      return std::sqrt(sumsq);
}

inline double lb_clim(const Lattice& lat,
                      const std::array<size_tu, 8>& idx,
                      const double* focal_clim) {
      double sumsq = 0.0;
      for (size_tu c = 0; c < lat.n_clim_dims; ++c) {
            size_tu d = lat.n_geo_dims + c;
            const double fc = focal_clim[c];
            const double c_min = lat.mins[d] + static_cast<double>(idx[d]) * lat.res[d];
            const double c_max = c_min + lat.res[d];

            double diff = 0.0;
            if (fc < c_min) diff = c_min - fc;
            else if (fc > c_max) diff = fc - c_max;

            sumsq += diff * diff;
      }
      return std::sqrt(sumsq);
}

inline bool clim_ok_and_dist_knn(const double* focal_clim,
                                 const double* ref_clim_col,
                                 size_tu n_clim,
                                 size_tu stride_r,
                                 bool use_scalar_clim,
                                 const std::vector<double>& max_clim_pervar,
                                 double max_clim_scalar,
                                 double& clim_dist_out,
                                 bool compute_dist) {
      double sumsq = 0.0;

      for (size_tu c = 0; c < n_clim; ++c) {
            const double fc = focal_clim[c];
            const double rc = ref_clim_col[c * stride_r];
            const double diff = fc - rc;

            // Check per-variable threshold
            if (!max_clim_pervar.empty() && c < max_clim_pervar.size()) {
                  double thr = max_clim_pervar[c];
                  if (std::isfinite(thr) && thr > 0.0) {
                        if (std::fabs(diff) > thr) return false;
                  }
            }

            if (compute_dist || use_scalar_clim) {
                  sumsq += diff * diff;
            }
      }

      // Check scalar threshold
      if (use_scalar_clim && std::isfinite(max_clim_scalar) && max_clim_scalar > 0.0) {
            double dist = std::sqrt(sumsq);
            if (dist > max_clim_scalar) return false;
      }

      if (compute_dist) {
            clim_dist_out = std::sqrt(sumsq);
      }

      return true;
}

} // namespace analogs
