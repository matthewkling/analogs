#pragma once
#include <cstdint>
#include <cstddef>

namespace analogs {

using index_t = int32_t;
using size_tu = std::size_t;

// Column-major non-owning matrix view over R's NumericMatrix memory.
struct MatrixView {
      const double* data = nullptr;
      size_tu nrow = 0, ncol = 0; // like R: column-major
      inline double operator()(size_tu i, size_tu j) const { return data[i + j*nrow]; }
};

// Simple non-owning int vector view
struct IntVectorView {
      const index_t* data = nullptr;
      size_tu n = 0;
      inline index_t operator[](size_tu i) const { return data[i]; }
};

// Geographic coordinate space selector
enum class GeoSpace { LonLat, Planar };

} // namespace analogs
