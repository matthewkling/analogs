#pragma once

#include "types.hpp"

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <array>
#include <cmath>
#include <limits>


namespace analogs {

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

      size_tu n_points;     // number of reference points
      size_tu n_geo_dims;   // always 2: x, y
      size_tu n_clim_dims;  // number of climate variables
      size_tu n_dims;       // = n_geo_dims + n_clim_dims

      // Per-dimension metadata
      std::vector<double> mins;    // min value per dimension
      std::vector<double> maxs;    // max value per dimension
      std::vector<double> res;     // bin width per dimension
      std::vector<size_tu> n_bins; // number of bins per dimension
      std::vector<size_tu> strides; // flattening strides per dimension

      // Sparse cells: key (flattened bin index) -> list of reference indices
      std::unordered_map<size_tu, std::vector<index_t> > cells;

      // Diagnostics
      size_tu total_bins;        // product of n_bins[d]
      size_tu min_cell_occ;      // min occupancy among occupied cells
      size_tu max_cell_occ;      // max occupancy among occupied cells
      size_tu n_cells_nonempty;  // number of occupied cells

      Lattice()
            : metric_type(MetricType::Planar),
              n_points(0),
              n_geo_dims(2),
              n_clim_dims(0),
              n_dims(2),
              total_bins(1),
              min_cell_occ(std::numeric_limits<size_tu>::max()),
              max_cell_occ(0),
              n_cells_nonempty(0) {}

      // Build lattice over all dims.
      //
      // ref_ptr: pointer to ref matrix (column-major R layout)
      // n_ref:   number of rows
      // n_clim:  number of climate columns (total columns = 2 + n_clim)
      // stride_r: n_ref (for column-major indexing)
      // metric:   geographic metric type (lon/lat vs planar)
      // max_dist: geographic threshold (km); may be Inf
      // use_scalar_clim: whether scalar max_clim is active
      // max_clim_pervar: per-variable climate thresholds (length n_clim)
      // max_clim_scalar: scalar Euclidean climate threshold (or Inf)
      void build(const double* ref_ptr,
                 size_tu n_ref,
                 size_tu n_clim,
                 size_tu stride_r,
                 MetricType metric,
                 double max_dist,
                 bool use_scalar_clim,
                 const std::vector<double>& max_clim_pervar,
                 double max_clim_scalar,
                 int    bins_per_dim,
                 bool   use_geo_lattice,
                 bool   use_clim_lattice);


      // Query candidate indices for a focal point.
      //
      // focal_geo: length-2 array [x, y]
      // focal_clim: length n_clim array [clim1, ..., climP]
      // max_dist, use_scalar_clim, max_clim_pervar, max_clim_scalar: same
      //   semantics as in build().
      // out_indices: will be filled with 0-based ref indices belonging to
      //   cells whose binned coordinates are consistent with the thresholds.
      void query(const double* focal_geo,
                 const double* focal_clim,
                 double max_dist,
                 bool use_scalar_clim,
                 const std::vector<double>& max_clim_pervar,
                 double max_clim_scalar,
                 std::vector<index_t>& out_indices) const;

      // KNN query using expanding search over lattice cells (Chebyshev shells).
      //
      // - focal_geo: length-2 [x, y]
      // - focal_clim: length n_clim (contiguous)
      // - ref_ptr: pointer to ref matrix (col-major: [x, y, clim1, ...])
      // - stride_r: n_ref (rows of ref matrix)
      // - n_clim: number of climate variables
      // - rank_by_geog: if true, kNN in geographic space; else in climate space
      // - max_geog: geographic constraint (km); may be Inf
      // - use_scalar_clim, max_clim_pervar, max_clim_scalar: climate constraints
      // - k: number of neighbors desired
      // - out_indices: filled with 0-based row indices (sorted by increasing
      //   ranking metric). May contain fewer than k if not enough admissible refs.
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
                     std::vector<index_t>& out_indices) const;


private:
      void enumerate_cells(const std::vector<size_tu>& lo,
                           const std::vector<size_tu>& hi,
                           std::vector<size_tu>& idx,
                           size_tu dim,
                           std::vector<index_t>& out_indices) const;

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
                           size_tu n_clim,
                           size_tu stride_r,
                           MetricType metric,
                           double max_dist,
                           bool use_scalar_clim,
                           const std::vector<double>& max_clim_pervar,
                           double max_clim_scalar,
                           int    bins_per_dim,
                           bool   use_geo_lattice,
                           bool   use_clim_lattice) {
      metric_type    = metric;
      n_points       = n_ref;
      n_geo_dims     = 2;
      n_clim_dims    = n_clim;
      n_dims         = n_geo_dims + n_clim_dims;

      mins.assign(n_dims, 0.0);
      maxs.assign(n_dims, 0.0);
      res.assign(n_dims, 1.0);
      n_bins.assign(n_dims, 1);
      strides.assign(n_dims, 1);
      cells.clear();
      total_bins    = 1;
      min_cell_occ  = std::numeric_limits<size_tu>::max();
      max_cell_occ  = 0;
      n_cells_nonempty  = 0;

      if (n_ref == 0) {
            return;
      }

      // First pass: compute mins/maxs over all dims.
      for (size_tu j = 0; j < n_ref; ++j) {
            // geo dims
            double x = ref_ptr[j];                // col 0
            double y = ref_ptr[j + stride_r];     // col 1

            if (j == 0) {
                  mins[0] = maxs[0] = x;
                  mins[1] = maxs[1] = y;
            } else {
                  if (x < mins[0]) mins[0] = x;
                  if (x > maxs[0]) maxs[0] = x;
                  if (y < mins[1]) mins[1] = y;
                  if (y > maxs[1]) maxs[1] = y;
            }

            // climate dims
            for (size_tu k = 0; k < n_clim; ++k) {
                  double v = ref_ptr[j + (2 + k) * stride_r];
                  size_tu d = 2 + k;
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


      // Second pass: assign points to cells.
      std::vector<size_tu> idx(n_dims);

      for (size_tu j = 0; j < n_ref; ++j) {
            // dim 0: x
            double v0 = ref_ptr[j];
            double pos0 = (v0 - mins[0]) / res[0];
            if (pos0 < 0.0) pos0 = 0.0;
            size_tu i0 = static_cast<size_tu>(pos0);
            if (i0 >= n_bins[0]) i0 = n_bins[0] - 1;
            idx[0] = i0;

            // dim 1: y
            double v1 = ref_ptr[j + stride_r];
            double pos1 = (v1 - mins[1]) / res[1];
            if (pos1 < 0.0) pos1 = 0.0;
            size_tu i1 = static_cast<size_tu>(pos1);
            if (i1 >= n_bins[1]) i1 = n_bins[1] - 1;
            idx[1] = i1;

            // climate dims
            for (size_tu k = 0; k < n_clim; ++k) {
                  size_tu d = 2 + k;
                  double v = ref_ptr[j + (2 + k) * stride_r];
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

            std::vector<index_t>& cell = cells[key];
            cell.push_back(static_cast<index_t>(j));

            size_tu occ = cell.size();
            if (occ > max_cell_occ) max_cell_occ = occ;
            if (occ < min_cell_occ) min_cell_occ = occ;
      }

      n_cells_nonempty = static_cast<size_tu>(cells.size());

      if (cells.empty()) {
            min_cell_occ = 0;
            max_cell_occ = 0;
      }
}

inline void Lattice::enumerate_cells(const std::vector<size_tu>& lo,
                                     const std::vector<size_tu>& hi,
                                     std::vector<size_tu>& idx,
                                     size_tu dim,
                                     std::vector<index_t>& out_indices) const {
      if (dim == n_dims) {
            size_tu key = 0;
            for (size_tu d = 0; d < n_dims; ++d) {
                  key += idx[d] * strides[d];
            }
            typename std::unordered_map<size_tu, std::vector<index_t> >::const_iterator it =
                  cells.find(key);
            if (it != cells.end()) {
                  const std::vector<index_t>& cell = it->second;
                  out_indices.insert(out_indices.end(), cell.begin(), cell.end());
            }
            return;
      }

      for (size_tu v = lo[dim]; v <= hi[dim]; ++v) {
            idx[dim] = v;
            enumerate_cells(lo, hi, idx, dim + 1, out_indices);
      }
}

inline void Lattice::query(const double* focal_geo,
                           const double* focal_clim,
                           double max_dist,
                           bool use_scalar_clim,
                           const std::vector<double>& max_clim_pervar,
                           double max_clim_scalar,
                           std::vector<index_t>& out_indices) const {
      out_indices.clear();
      if (n_points == 0 || n_dims == 0) return;

      std::vector<size_tu> lo(n_dims);
      std::vector<size_tu> hi(n_dims);
      std::vector<size_tu> idx(n_dims);

      // Per-dimension allowed bin index ranges.
      for (size_tu d = 0; d < n_dims; ++d) {
            double minv = mins[d];
            double maxv = maxs[d];

            if (d < n_geo_dims) {
                  // Geo dims: bound by max_dist if finite.
                  if (std::isfinite(max_dist) && max_dist > 0.0) {
                        double q = focal_geo[d];
                        minv = q - max_dist;
                        maxv = q + max_dist;
                  }
            } else {
                  // Climate dims: bound by per-var or scalar climate thresholds.
                  size_tu k = d - n_geo_dims;
                  double band = std::numeric_limits<double>::infinity();

                  if (k < max_clim_pervar.size() &&
                      std::isfinite(max_clim_pervar[k]) &&
                      max_clim_pervar[k] > 0.0) {
                        band = max_clim_pervar[k];
                  } else if (use_scalar_clim &&
                        std::isfinite(max_clim_scalar) &&
                        max_clim_scalar > 0.0) {
                        band = max_clim_scalar;
                  }

                  if (std::isfinite(band) && band > 0.0) {
                        double q = focal_clim[k];
                        minv = q - band;
                        maxv = q + band;
                  }
            }

            // Convert value range to bin index range.
            double pos_lo = (minv - mins[d]) / res[d];
            double pos_hi = (maxv - mins[d]) / res[d];
            if (pos_lo < 0.0) pos_lo = 0.0;
            if (pos_hi < 0.0) pos_hi = 0.0;

            size_tu ilo = static_cast<size_tu>(pos_lo);
            size_tu ihi = static_cast<size_tu>(pos_hi);

            if (ilo >= n_bins[d]) ilo = n_bins[d] - 1;
            if (ihi >= n_bins[d]) ihi = n_bins[d] - 1;
            if (ilo > ihi) {
                  size_tu tmp = ilo;
                  ilo = ihi;
                  ihi = tmp;
            }
            lo[d] = ilo;
            hi[d] = ihi;
      }

      // Enumerate all cells within the bin ranges and collect their points.
      enumerate_cells(lo, hi, idx, 0, out_indices);
}

// Forward declarations for helpers used by knn_query
inline double lb_geo_projected(const Lattice& lat,
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
                               std::vector<index_t>& out_indices) const {
      out_indices.clear();
      if (k <= 0 || n_points == 0 || n_dims == 0) return;

      // For now, only implement expanding search for projected geo (planar).
      if (metric_type != MetricType::Planar && rank_by_geog) {
            // Caller should fall back to old path for non-planar metrics.
            // We just bail here.
            return;
      }

      // Sanity: n_clim argument should match lattice's climate dims.
      if (n_clim != n_clim_dims) {
            // You could throw or assert here; for now we just bail.
            return;
      }

      const size_tu D = n_dims;
      const size_tu n_geo = n_geo_dims; // 2

      // 1. Compute per-dimension bin ranges [lo[d], hi[d]] based on thresholds
      std::vector<size_tu> lo(D), hi(D);

      for (size_tu d = 0; d < D; ++d) {
            double minv = mins[d];
            double maxv = maxs[d];

            if (d < n_geo) {
                  // Geo dims: bound by max_geog if finite
                  if (std::isfinite(max_geog) && max_geog > 0.0) {
                        const double q = focal_geo[d];
                        minv = q - max_geog;
                        maxv = q + max_geog;
                  }
            } else {
                  // Climate dims: bound by per-var or scalar climate thresholds
                  const size_tu kdim = d - n_geo;
                  double band = std::numeric_limits<double>::infinity();

                  if (kdim < max_clim_pervar.size() &&
                      std::isfinite(max_clim_pervar[kdim]) &&
                      max_clim_pervar[kdim] > 0.0) {
                        band = max_clim_pervar[kdim];
                  } else if (use_scalar_clim &&
                        std::isfinite(max_clim_scalar) &&
                        max_clim_scalar > 0.0) {
                        band = max_clim_scalar;
                  }

                  if (std::isfinite(band) && band > 0.0) {
                        const double q = focal_clim[kdim];
                        minv = q - band;
                        maxv = q + band;
                  }
            }

            // Convert to bin indices
            double pos_lo = (minv - mins[d]) / res[d];
            double pos_hi = (maxv - mins[d]) / res[d];
            if (pos_lo < 0.0) pos_lo = 0.0;
            if (pos_hi < 0.0) pos_hi = 0.0;

            size_tu ilo = static_cast<size_tu>(pos_lo);
            size_tu ihi = static_cast<size_tu>(pos_hi);

            if (ilo >= n_bins[d]) ilo = n_bins[d] - 1;
            if (ihi >= n_bins[d]) ihi = n_bins[d] - 1;
            if (ilo > ihi) {
                  size_tu tmp = ilo;
                  ilo = ihi;
                  ihi = tmp;
            }
            lo[d] = ilo;
            hi[d] = ihi;
      }

      // 2. Seed cell index (center of Chebyshev shells)
      std::array<size_tu, 8> center;
      for (size_tu d = 0; d < D; ++d) {
            double q;
            if (d < n_geo) {
                  q = focal_geo[d];
            } else {
                  q = focal_clim[d - n_geo];
            }
            double pos = (q - mins[d]) / res[d];
            if (pos < 0.0) pos = 0.0;
            size_tu ic = static_cast<size_tu>(pos);
            if (ic >= n_bins[d]) ic = n_bins[d] - 1;

            // Clamp to [lo[d], hi[d]] so we're inside the threshold box
            if (ic < lo[d]) ic = lo[d];
            if (ic > hi[d]) ic = hi[d];

            center[d] = ic;
      }

      auto cell_lb = [&](const std::array<size_tu, 8>& idx) -> double {
            if (rank_by_geog) {
                  return lb_geo_projected(*this, idx, focal_geo);
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
                  // min-heap by lb: top has smallest lb
                  return a.lb > b.lb;
            }
      };

      std::priority_queue<CellState, std::vector<CellState>, CellCmp> pq;
      std::unordered_set<size_tu> visited;
      visited.reserve(128);

      CellState start;
      start.idx = center;
      start.lb  = cell_lb(center);
      const size_tu start_key = flatten_idx(center);
      visited.insert(start_key);
      pq.push(start);

      // kNN heap: max-heap by distance (we keep smallest k)
      using Neighbor = std::pair<double, int>; // (distance, row_index_0based)
      struct NeighborCmp {
            bool operator()(const Neighbor& a, const Neighbor& b) const {
                  return a.first < b.first; // top has largest distance
            }
      };
      std::priority_queue<Neighbor, std::vector<Neighbor>, NeighborCmp> knn;

      const bool use_geo_constraint = (std::isfinite(max_geog) && max_geog > 0.0);
      const double max_geog2 = (use_geo_constraint ? max_geog * max_geog
                                      : std::numeric_limits<double>::infinity());

      double d_k = std::numeric_limits<double>::infinity();

      // 3. Expanding search over cells
      std::array<size_tu, 8> idx;

      while (!pq.empty()) {
            CellState cell = pq.top();
            pq.pop();

            // Early exit: no unvisited cell can improve on current k-th neighbor
            if (!knn.empty() && static_cast<int>(knn.size()) >= k && cell.lb > d_k) {
                  break;
            }

            idx = cell.idx;
            const size_tu key = flatten_idx(idx);

            // Look up refs in this cell
            auto it = cells.find(key);
            if (it != cells.end()) {
                  const std::vector<index_t>& refs = it->second;

                  for (size_tu t = 0; t < refs.size(); ++t) {
                        const int j = static_cast<int>(refs[t]); // 0-based row index

                        const double rx = ref_ptr[j];
                        const double ry = ref_ptr[j + stride_r];
                        const double* r_clim_col = ref_ptr + j + 2 * stride_r;

                        // Geographic distance (planar)
                        double gdist2 = 0.0;
                        double gdist = 0.0;
                        if (rank_by_geog || use_geo_constraint) {
                              const double dx = focal_geo[0] - rx;
                              const double dy = focal_geo[1] - ry;
                              gdist2 = dx * dx + dy * dy;
                              gdist  = std::sqrt(gdist2);
                        }

                        if (use_geo_constraint && gdist2 > max_geog2) {
                              continue;
                        }

                        // Climate constraints and (optionally) distance
                        double clim_dist = 0.0;
                        const bool compute_clim_dist = !rank_by_geog; // need it for knn_clim
                        if (!clim_ok_and_dist_knn(
                                    focal_clim, r_clim_col,
                                    n_clim, stride_r,
                                    use_scalar_clim, max_clim_pervar, max_clim_scalar,
                                    clim_dist, compute_clim_dist)) {
                              continue;
                        }

                        // Ranking distance: geog or clim
                        double key_dist = 0.0;
                        if (rank_by_geog) {
                              key_dist = gdist;
                        } else {
                              key_dist = clim_dist;
                        }

                        if (static_cast<int>(knn.size()) < k) {
                              knn.emplace(key_dist, j);
                        } else if (!knn.empty() && key_dist < knn.top().first) {
                              knn.pop();
                              knn.emplace(key_dist, j);
                        }

                        if (static_cast<int>(knn.size()) == k) {
                              d_k = knn.top().first;
                        }
                  }
            }

            // Expand neighbors in Chebyshev shells: +/-1 in each dimension
            for (size_tu d = 0; d < D; ++d) {
                  // -1
                  if (idx[d] > lo[d]) {
                        std::array<size_tu, 8> nb_idx = idx;
                        nb_idx[d] = idx[d] - 1;
                        const size_tu nb_key = flatten_idx(nb_idx);
                        if (visited.insert(nb_key).second) {
                              CellState nb;
                              nb.idx = nb_idx;
                              nb.lb = cell_lb(nb_idx);
                              pq.push(nb);
                        }
                  }
                  // +1
                  if (idx[d] + 1 <= hi[d]) {
                        std::array<size_tu, 8> nb_idx = idx;
                        nb_idx[d] = idx[d] + 1;
                        const size_tu nb_key = flatten_idx(nb_idx);
                        if (visited.insert(nb_key).second) {
                              CellState nb;
                              nb.idx = nb_idx;
                              nb.lb = cell_lb(nb_idx);
                              pq.push(nb);
                        }
                  }
            }
      }

      // 4. Extract neighbors in ascending distance order
      const int m = static_cast<int>(knn.size());
      out_indices.resize(m);
      // we want smallest distance first, so pop into reverse
      for (int pos = m - 1; pos >= 0; --pos) {
            const Neighbor nb = knn.top();
            out_indices[pos] = static_cast<index_t>(nb.second); // 0-based row index
            knn.pop();
      }
}


// Lower bound on geographic distance (projected) from focal to a cell.
inline double lb_geo_projected(const Lattice& lat,
                               const std::array<size_tu, 8>& idx,
                               const double* focal_geo) {
      // dim 0: x, dim 1: y
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

// Lower bound on climate distance (Euclidean) from focal to a cell.
inline double lb_clim(const Lattice& lat,
                      const std::array<size_tu, 8>& idx,
                      const double* focal_clim) {
      const size_tu n_geo = lat.n_geo_dims;   // 2
      const size_tu n_clim = lat.n_clim_dims;

      double sumsq = 0.0;
      for (size_tu k = 0; k < n_clim; ++k) {
            const size_tu d = n_geo + k;
            const double c_min = lat.mins[d] + static_cast<double>(idx[d]) * lat.res[d];
            const double c_max = c_min + lat.res[d];
            const double cf = focal_clim[k];

            double dc = 0.0;
            if (cf < c_min) dc = c_min - cf;
            else if (cf > c_max) dc = cf - c_max;

            sumsq += dc * dc;
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
      const bool use_pervar = !max_clim_pervar.empty();
      const double max_clim2 = (use_scalar_clim && std::isfinite(max_clim_scalar) && max_clim_scalar > 0.0)
            ? max_clim_scalar * max_clim_scalar
      : std::numeric_limits<double>::infinity();

      for (size_tu k = 0; k < n_clim; ++k) {
            const double fk = focal_clim[k];
            const double rk = ref_clim_col[k * stride_r];
            const double diff = fk - rk;
            const double adiff = std::fabs(diff);

            if (use_pervar && k < max_clim_pervar.size()) {
                  const double band = max_clim_pervar[k];
                  if (std::isfinite(band) && band > 0.0 && adiff > band) {
                        return false;
                  }
            }

            sumsq += diff * diff;
      }

      if (max_clim2 < std::numeric_limits<double>::infinity()) {
            if (sumsq > max_clim2) return false;
      }

      if (compute_dist) {
            clim_dist_out = std::sqrt(sumsq);
      }
      return true;
}


} // namespace analogs
