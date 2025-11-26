#pragma once
#include <cstdint>
#include <cstddef>
#include <utility>

namespace analogs {
using size_tu = std::size_t;
using index_t = int32_t;

// Common pair type for kNN results: (distance, ref_index)
using Neighbor = std::pair<double, index_t>;

enum class MetricType {
      Haversine,  // lon/lat on sphere
      Planar,     // projected / Euclidean geometry
      Chord3D     // 3D Euclidean chord distance (for ECEF)
};

struct MatrixView {
      const double* data;
      size_tu nrow;
      size_tu ncol;
};
}
