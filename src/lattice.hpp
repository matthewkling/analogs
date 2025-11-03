#pragma once

#include "types.hpp"
#include "metrics.hpp"
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>

namespace analogs {

// ---- Coordinate System Handling ------------------------------------------
//
// We support two geographic cases controlled by MetricType:
//   * MetricType::Haversine : (geo1, geo2) = (lat, lon) in degrees;
//                              distances are great-circle.
//   * otherwise             : (geo1, geo2) = projected (x, y) coordinates;
//                              distances are Euclidean in the projected CRS.
//
// In both cases we store simple 2-D bounding boxes for each occupied
// (geo × clim) lattice cell and use them to cheaply lower-bound the
// distance from a query point to the cell.

// Axis-aligned bounding box in geographic (lat/lon) space.
struct AABB {
    double lat_min, lat_max;
    double lon_min, lon_max;

    AABB()
        : lat_min(std::numeric_limits<double>::infinity()),
          lat_max(-std::numeric_limits<double>::infinity()),
          lon_min(std::numeric_limits<double>::infinity()),
          lon_max(-std::numeric_limits<double>::infinity()) {}

    // Lower bound on great-circle distance from query point to any point
    // in this bin (assuming lat/lon in degrees).
    double lower_bound_distance(double qlat, double qlon) const {
        // Clamp to the closest point in the rectangle and compute
        // great-circle distance to that point.
        const double closest_lat = std::clamp(qlat, lat_min, lat_max);
        const double closest_lon = std::clamp(qlon, lon_min, lon_max);
        return haversine_distance(qlat, qlon, closest_lat, closest_lon);
    }
};

// Simple 2-D box for projected / planar coordinates.
struct PlanarBox {
    double x_min, x_max;
    double y_min, y_max;

    PlanarBox()
        : x_min(std::numeric_limits<double>::infinity()),
          x_max(-std::numeric_limits<double>::infinity()),
          y_min(std::numeric_limits<double>::infinity()),
          y_max(-std::numeric_limits<double>::infinity()) {}

    // Lower bound on Euclidean distance from query point to any point
    // in this bin (x/y units should match the distance units).
    double lower_bound_distance(double qx, double qy) const {
        const double x_dist =
            std::max(0.0, std::max(x_min - qx, qx - x_max));
        const double y_dist =
            std::max(0.0, std::max(y_min - qy, qy - y_max));
        return std::sqrt(x_dist * x_dist + y_dist * y_dist);
    }
};

// ---- Lattice Structure ---------------------------------------------------
//
// Sparse (geo × clim) lattice over reference points.  The discrete
// climate dimension is always treated as exact equality (i.e., we only
// search a single climate bin per query); geographic bins are pruned
// using the bounding boxes above.

class Lattice {
public:
    // CSR sparse structure over the climate dimension, with rows indexed
    // by geo_bin.  row_ptr has length n_geo_bins + 1, col_ind and
    // neighbors have length equal to the number of occupied cells (nnz).
    std::vector<size_tu> row_ptr;
    std::vector<size_tu> col_ind;                 // climate bin for each cell
    std::vector<std::vector<index_t>> neighbors;  // reference indices per cell

    // Bin metadata (only one of these is used, depending on metric_type).
    std::vector<AABB>      aabbs;   // for Haversine / lon-lat
    std::vector<PlanarBox> boxes;   // for projected / Euclidean

    // Resolution / domain info
    size_tu n_geo_bins;
    size_tu n_clim_bins;
    double geo_resolution;   // bin width in geo1 units
    double clim_resolution;  // bin width in climate units
    MetricType metric_type;

    double geo_min, geo_max;
    double clim_min, clim_max;

    Lattice()
        : n_geo_bins(0),
          n_clim_bins(0),
          geo_resolution(0),
          clim_resolution(0),
          metric_type(MetricType::Haversine),
          geo_min(0),
          geo_max(0),
          clim_min(0),
          clim_max(0) {}

    // Build lattice from reference data.
    //
    // ref_geo:  n × 2 matrix, columns = (geo1, geo2)
    // ref_clim: n × 1 vector/matrix
    //
    // geo_res:   bin width in geo1 units (e.g., degrees latitude or km)
    // clim_res:  bin width in climate units
    void build(const MatrixView& ref_geo,   // n × 2: geo1, geo2
               const MatrixView& ref_clim,  // n × 1: climate variable
               double geo_res,
               double clim_res,
               MetricType metric);

    // Query interface – return indices into `neighbors` (i.e., occupied
    // lattice cells) that belong to the same climate bin as `clim_val`
    // and whose geographic box could contain a point within
    // `geo_threshold` of (geo1, geo2) under the active metric.
    std::vector<size_tu> query_bins(double geo1,
                                    double geo2,
                                    double clim_val,
                                    double geo_threshold) const;

    // Collect all reference indices in the specified occupied cells.
    std::vector<index_t> get_neighbors_in_bins(
        const std::vector<size_tu>& bin_indices) const;

private:
    // Bin assignment in 1-D geo1 / climate space.
    size_tu geo_bin_index(double geo_val) const;
    size_tu clim_bin_index(double clim_val) const;
    size_tu joint_bin_index(size_tu geo_idx, size_tu clim_idx) const;

    // Update per-cell geographic bounds.
    void update_bounds(size_tu bin_idx, double geo1, double geo2);
};

// ---- Implementation -------------------------------------------------------

inline void Lattice::build(const MatrixView& ref_geo,
                           const MatrixView& ref_clim,
                           double geo_res,
                           double clim_res,
                           MetricType metric) {
    metric_type     = metric;
    geo_resolution  = geo_res;
    clim_resolution = clim_res;

    const size_tu n = ref_geo.nrow;
    if (n == 0) {
        // Fully reset state in case of reuse.
        n_geo_bins = n_clim_bins = 0;
        row_ptr.clear();
        col_ind.clear();
        neighbors.clear();
        aabbs.clear();
        boxes.clear();
        return;
    }

    // Compute 1-D bounds in geo1 and climate.
    geo_min  = geo_max  = ref_geo.data[0];
    clim_min = clim_max = ref_clim.data[0];

    for (size_tu i = 0; i < n; ++i) {
        const double g1 = ref_geo.data[i];
        const double c  = ref_clim.data[i];

        geo_min  = std::min(geo_min, g1);
        geo_max  = std::max(geo_max, g1);
        clim_min = std::min(clim_min, c);
        clim_max = std::max(clim_max, c);
    }

    // Determine number of bins.
    n_geo_bins  = static_cast<size_tu>(
        std::ceil((geo_max - geo_min) / geo_res)) + 1;
    n_clim_bins = static_cast<size_tu>(
        std::ceil((clim_max - clim_min) / clim_res)) + 1;

    // Temporary structure: map joint bin index to list of point indices.
    std::vector<std::vector<index_t>> temp_bins(n_geo_bins * n_clim_bins);

    // Assign points to bins.
    for (size_tu i = 0; i < n; ++i) {
        const double g1 = ref_geo.data[i];
        const double c  = ref_clim.data[i];

        const size_tu geo_idx  = geo_bin_index(g1);
        const size_tu clim_idx = clim_bin_index(c);
        const size_tu bin      = joint_bin_index(geo_idx, clim_idx);

        temp_bins[bin].push_back(static_cast<index_t>(i));
    }

    // Build CSR structure.
    row_ptr.clear();
    row_ptr.reserve(n_geo_bins + 1);
    col_ind.clear();
    neighbors.clear();

    aabbs.clear();
    boxes.clear();

    for (size_tu geo_idx = 0; geo_idx < n_geo_bins; ++geo_idx) {
        row_ptr.push_back(col_ind.size());

        for (size_tu clim_idx = 0; clim_idx < n_clim_bins; ++clim_idx) {
            const size_tu bin = joint_bin_index(geo_idx, clim_idx);

            if (!temp_bins[bin].empty()) {
                col_ind.push_back(clim_idx);
                neighbors.push_back(std::move(temp_bins[bin]));

                // Create and populate bounds for this occupied cell.
                const size_tu cell_idx = neighbors.size() - 1;

                if (metric_type == MetricType::Haversine) {
                    aabbs.emplace_back();
                } else {
                    boxes.emplace_back();
                }

                for (index_t pt_idx : neighbors[cell_idx]) {
                    const double g1 = ref_geo.data[pt_idx];
                    const double g2 = ref_geo.data[pt_idx + n];
                    update_bounds(cell_idx, g1, g2);
                }
            }
        }
    }

    row_ptr.push_back(col_ind.size());
}

inline std::vector<size_tu> Lattice::query_bins(double geo1,
                                                double geo2,
                                                double clim_val,
                                                double geo_threshold) const {
    std::vector<size_tu> result;

    if (n_geo_bins == 0 || n_clim_bins == 0) {
        return result;
    }

    const size_tu query_clim_idx = clim_bin_index(clim_val);

    // For now we scan all geo rows and prune only by climate bin and
    // bounding-box distance. If needed, we can later restrict the
    // geo_idx loop using geo_threshold and geo_resolution.
    for (size_tu geo_idx = 0; geo_idx < n_geo_bins; ++geo_idx) {
        const size_tu row_begin = row_ptr[geo_idx];
        const size_tu row_end   = row_ptr[geo_idx + 1];

        for (size_tu j = row_begin; j < row_end; ++j) {
            const size_tu clim_idx = col_ind[j];

            // Only consider the matching climate bin.
            if (clim_idx != query_clim_idx) {
                continue;
            }

            double lb_dist;
            if (metric_type == MetricType::Haversine) {
                lb_dist = aabbs[j].lower_bound_distance(geo1, geo2);
            } else {
                lb_dist = boxes[j].lower_bound_distance(geo1, geo2);
            }

            if (lb_dist <= geo_threshold) {
                // j indexes into neighbors / aabbs / boxes.
                result.push_back(j);
            }
        }
    }

    return result;
}

inline std::vector<index_t> Lattice::get_neighbors_in_bins(
    const std::vector<size_tu>& bin_indices) const {
    std::vector<index_t> result;

    for (size_tu bin_idx : bin_indices) {
        const auto& bin_neighbors = neighbors[bin_idx];
        result.insert(result.end(), bin_neighbors.begin(), bin_neighbors.end());
    }

    return result;
}

inline size_tu Lattice::geo_bin_index(double geo_val) const {
    const size_tu idx =
        static_cast<size_tu>((geo_val - geo_min) / geo_resolution);
    return std::min(idx, n_geo_bins - 1);
}

inline size_tu Lattice::clim_bin_index(double clim_val) const {
    const size_tu idx =
        static_cast<size_tu>((clim_val - clim_min) / clim_resolution);
    return std::min(idx, n_clim_bins - 1);
}

inline size_tu Lattice::joint_bin_index(size_tu geo_idx,
                                        size_tu clim_idx) const {
    return geo_idx * n_clim_bins + clim_idx;
}

inline void Lattice::update_bounds(size_tu bin_idx,
                                   double geo1,
                                   double geo2) {
    if (metric_type == MetricType::Haversine) {
        AABB& box = aabbs[bin_idx];
        box.lat_min = std::min(box.lat_min, geo1);
        box.lat_max = std::max(box.lat_max, geo1);
        box.lon_min = std::min(box.lon_min, geo2);
        box.lon_max = std::max(box.lon_max, geo2);
    } else {
        PlanarBox& box = boxes[bin_idx];
        box.x_min = std::min(box.x_min, geo1);
        box.x_max = std::max(box.x_max, geo1);
        box.y_min = std::min(box.y_min, geo2);
        box.y_max = std::max(box.y_max, geo2);
    }
}

}  // namespace analogs
w