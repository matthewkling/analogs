#pragma once

#include "types.h"
#include "metrics.h"
#include <cmath>

namespace analogs {

// Helper: geographic distance (km), lon/lat vs projected.
inline double geo_distance_km(double x1, double y1,
                              double x2, double y2,
                              bool use_haversine) {
      if (use_haversine) {
            double a[2] = { x1, y1 };
            double b[2] = { x2, y2 };
            Haversine2D h;
            return h.dist(a, b, 2);
      } else {
            const double dx = x1 - x2;
            const double dy = y1 - y2;
            return std::sqrt(dx * dx + dy * dy);
      }
}

// Convert lon/lat (degrees) to ECEF (Earth-Centered Earth-Fixed) coordinates (km)
inline void lonlat_to_ecef(double lon_deg, double lat_deg, double R,
                           double &X, double &Y, double &Z) {
      const double lon = lon_deg * (M_PI / 180.0);
      const double lat = lat_deg * (M_PI / 180.0);

      const double cos_lat = std::cos(lat);

      X = R * cos_lat * std::cos(lon);
      Y = R * cos_lat * std::sin(lon);
      Z = R * std::sin(lat);
}

// Convert geographic distance threshold (km) to ECEF chord distance (km)
inline double km_to_chord(double km_dist, double R) {
      // arc_distance = km_dist
      // central_angle = arc_distance / R (radians)
      // chord_distance = 2 * R * sin(central_angle / 2)
      const double central_angle = km_dist / R;
      return 2.0 * R * std::sin(0.5 * central_angle);
}

} // namespace analogs
