#pragma once
#include "types.hpp"
#include <cmath>
#include <algorithm>

namespace analogs {

// Simple abstract base metric (mainly for clarity / future extensibility)
struct Metric {
    virtual ~Metric() {}
    virtual double dist(const double* a,
                        const double* b,
                        size_tu len) const = 0;
};

// Euclidean metric in R^len
struct Euclidean : Metric {
    double dist(const double* a,
                const double* b,
                size_tu len) const override {
        double s = 0.0;
        for (size_tu k = 0; k < len; ++k) {
            double d = a[k] - b[k];
            s += d * d;
        }
        return std::sqrt(s);
    }
};

// Great-circle distance on a sphere using lon/lat in degrees.
// Expects a[0]=lon, a[1]=lat, same for b. Ignores len>2.
struct Haversine2D : Metric {
    double dist(const double* a,
                const double* b,
                size_tu /*len*/) const override {
        const double PI   = 3.14159265358979323846;
        const double kRad = PI / 180.0;
        const double R    = 6371.0088; // mean Earth radius (km)

        double lon1 = a[0] * kRad, lat1 = a[1] * kRad;
        double lon2 = b[0] * kRad, lat2 = b[1] * kRad;

        double dlon = lon2 - lon1;
        double dlat = lat2 - lat1;

        double sdlat = std::sin(0.5 * dlat);
        double sdlon = std::sin(0.5 * dlon);

        double h = sdlat * sdlat +
                   std::cos(lat1) * std::cos(lat2) * sdlon * sdlon;

        // Guard against tiny numerical issues
        double ang = 2.0 * std::asin(std::sqrt(std::min(1.0, h)));
        return R * ang;
    }
};

// Convenience helper for lon/lat points (deg) to km
inline double haversine_distance(double lon1, double lat1,
                                 double lon2, double lat2) {
    double a[2] = { lon1, lat1 };
    double b[2] = { lon2, lat2 };
    Haversine2D h;
    return h.dist(a, b, 2);
}

} // namespace analogs
