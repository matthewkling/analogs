#pragma once
#include "types.hpp"
#include <cmath>

namespace analogs {

struct Metric {
      virtual ~Metric() {}
      virtual double dist(const double* a, const double* b, size_tu len) const = 0;
};

struct Euclidean : Metric {
      double dist(const double* a, const double* b, size_tu len) const {
            double s = 0.0;
            for (size_tu k=0; k<len; ++k) { double d = a[k]-b[k]; s += d*d; }
            return std::sqrt(s);
      }
};

// Great-circle distance in km given lon/lat in degrees
struct Haversine2D : Metric {
      double dist(const double* a, const double* b, size_tu /*len*/) const {
            const double PI = 3.14159265358979323846;
            const double kRad = PI/180.0;
            const double R = 6371.0088; // mean Earth radius (km)
            double lon1 = a[0]*kRad, lat1 = a[1]*kRad;
            double lon2 = b[0]*kRad, lat2 = b[1]*kRad;
            double dlon = lon2 - lon1, dlat = lat2 - lat1;
            double sdlat = std::sin(dlat*0.5), sdlon = std::sin(dlon*0.5);
            double h = sdlat*sdlat + std::cos(lat1)*std::cos(lat2)*sdlon*sdlon;
            double ang = 2.0 * std::asin(std::sqrt(h));
            return R * ang;
      }
};

} // namespace analogs
