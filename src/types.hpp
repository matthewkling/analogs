#pragma once
#include <cstdint>
#include <cstddef>

namespace analogs {
    using size_tu = std::size_t;
    using index_t = int32_t;
    
    struct MatrixView {
        const double* data;
        size_tu nrow;
        size_tu ncol;
    };
}