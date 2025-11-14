#pragma once
#include <cstdint>
#include <cstddef>

namespace analogs {
    using size_tu = std::size_t;
    using index_t = int32_t;

    enum class MetricType {
        Haversine,  // lon/lat on sphere
        Planar      // projected / Euclidean geometry
    };
    
    struct MatrixView {
        const double* data;
        size_tu nrow;
        size_tu ncol;
    };
}