#pragma once

#include "types.hpp"

#include <vector>
#include <unordered_map>
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
    size_tu total_bins;   // product of n_bins[d]
    size_tu min_cell_occ; // min occupancy among occupied cells
    size_tu max_cell_occ; // max occupancy among occupied cells

    Lattice()
        : metric_type(MetricType::Planar),
          n_points(0),
          n_geo_dims(2),
          n_clim_dims(0),
          n_dims(2),
          total_bins(1),
          min_cell_occ(std::numeric_limits<size_tu>::max()),
          max_cell_occ(0) {}

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
               double max_clim_scalar);

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

private:
    void enumerate_cells(const std::vector<size_tu>& lo,
                         const std::vector<size_tu>& hi,
                         std::vector<size_tu>& idx,
                         size_tu dim,
                         std::vector<index_t>& out_indices) const;
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
                           double max_clim_scalar) {
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

    // Choose per-dimension bin widths.
    for (size_tu d = 0; d < n_dims; ++d) {
        double span = maxs[d] - mins[d];
        if (d < n_geo_dims) {
            // Geographic dims: use max_dist if finite, otherwise 10 bins.
            if (std::isfinite(max_dist) && max_dist > 0.0) {
                res[d] = max_dist;
            } else {
                res[d] = (span > 0.0) ? (span / 10.0) : 1.0;
            }
        } else {
            // Climate dims: use per-variable or scalar thresholds if finite,
            // otherwise 10 bins over the range.
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
                res[d] = band;
            } else {
                res[d] = (span > 0.0) ? (span / 10.0) : 1.0;
            }
        }
    }

    // Compute bin counts and strides.
    total_bins = 1;
    for (size_tu d = 0; d < n_dims; ++d) {
        double span = maxs[d] - mins[d];
        size_tu nb = 1;
        if (span > 0.0) {
            nb = static_cast<size_tu>(std::ceil(span / res[d])) + 1;
        }
        if (nb == 0) nb = 1;
        n_bins[d] = nb;
    }

    // strides: last dimension has stride 1
    strides[n_dims - 1] = 1;
    for (int d = static_cast<int>(n_dims) - 2; d >= 0; --d) {
        strides[d] = strides[d + 1] * n_bins[d + 1];
    }
    total_bins = strides[0] * n_bins[0];

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

} // namespace analogs
