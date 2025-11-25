// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>

#include "lattice.hpp"
#include "metrics.hpp"
#include "types.hpp"
#include "profiling.hpp"

#include <queue>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#include <chrono>

using namespace Rcpp;
using namespace analogs;
using namespace RcppParallel;

namespace analogs { thread_local ProfileTimer g_profiler; }

namespace {

// ----- Modes / enums consistent with R wrapper mapping -------------------
enum class ModeCode : int {
      KNN_CLIM = 0,
            KNN_GEOG = 1,
            COUNT    = 2,
            SUM      = 3,
            MEAN     = 4,
            ALL      = 5
};

enum class WeightCode : int {
      NONE            = 0,
            UNIFORM         = 1,
            INVERSE_CLIM    = 2,
            INVERSE_GEOG    = 3,
            GAUSSIAN_CLIM   = 4,
            GAUSSIAN_GEOG   = 5,
            GAUSSIAN_JOINT  = 6,
            INVERSE_JOINT   = 7
};

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

// Compute Euclidean climate distance (scalar) and/or per-var checks.
// Returns (ok, clim_dist) where ok means thresholds satisfied.
inline std::pair<bool,double>
      clim_ok_and_dist(const double* f_clim_col,
                       const double* r_clim_col,
                       int n_clim,
                       int stride_f,
                       int stride_r,
                       bool use_pervar_clim,
                       const std::vector<double>& max_clim_pervar,
                       bool use_scalar_clim,
                       double max_clim_scalar)
      {
            double sumsq = 0.0;

            if (use_pervar_clim) {
                  for (int k = 0; k < n_clim; ++k) {
                        const double df = f_clim_col[k * stride_f] - r_clim_col[k * stride_r];
                        if (std::fabs(df) > max_clim_pervar[k]) {
                              return std::make_pair(false, 0.0);
                        }
                        sumsq += df * df;
                  }
                  return std::make_pair(true, std::sqrt(sumsq));
            } else {
                  // scalar threshold or just distance
                  for (int k = 0; k < n_clim; ++k) {
                        const double df = f_clim_col[k * stride_f] - r_clim_col[k * stride_r];
                        sumsq += df * df;
                  }
                  const double d = std::sqrt(sumsq);
                  if (use_scalar_clim && d > max_clim_scalar) {
                        return std::make_pair(false, d);
                  }
                  return std::make_pair(true, d);
            }
      }

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

// Weight calculation from codes with optimized parameters
// For efficiency, pre-computed parameters are passed in
inline double weight_from_codes(WeightCode wc,
                                double clim_dist,
                                double geog_dist,
                                double param1,
                                double param2)
{
      switch (wc) {
      case WeightCode::UNIFORM:
            return 1.0;

      case WeightCode::INVERSE_CLIM: {
            // param1 = epsilon
            return 1.0 / (clim_dist + param1);
      }

      case WeightCode::INVERSE_GEOG: {
            // param1 = epsilon
            return 1.0 / (geog_dist + param1);
      }

      case WeightCode::GAUSSIAN_CLIM: {
            // param1 = -1/(2*sigma^2)
            const double d2 = clim_dist * clim_dist;
            return std::exp(param1 * d2);
      }

      case WeightCode::GAUSSIAN_GEOG: {
            // param1 = -1/(2*sigma^2)
            const double d2 = geog_dist * geog_dist;
            return std::exp(param1 * d2);
      }

      case WeightCode::GAUSSIAN_JOINT: {
            // param1 = -1/(2*sigma_clim^2), param2 = -1/(2*sigma_geog^2)
            const double d2_clim = clim_dist * clim_dist;
            const double d2_geog = geog_dist * geog_dist;
            return std::exp(param1 * d2_clim + param2 * d2_geog);
      }

      case WeightCode::INVERSE_JOINT: {
            // param1 = eps_clim, param2 = eps_geog
            const double d_clim_adj = clim_dist + param1;
            const double d_geog_adj = geog_dist + param2;
            const double joint_dist = std::sqrt(d_clim_adj * d_clim_adj +
                                                d_geog_adj * d_geog_adj);
            return 1.0 / joint_dist;
      }

      default: // NONE or unknown
            return 1.0;
      }
}

// Pre-compute weight parameters for efficiency
// Returns (param1, param2) for use in weight_from_codes
inline std::pair<double, double> precompute_weight_params(WeightCode wc,
                                                          const NumericVector& theta_vec) {
      double param1 = 0.0;
      double param2 = 0.0;

      switch (wc) {
      case WeightCode::UNIFORM:
      case WeightCode::NONE:
            break;

      case WeightCode::INVERSE_CLIM:
      case WeightCode::INVERSE_GEOG: {
            // param1 = epsilon
            if (theta_vec.size() > 0 && std::isfinite(theta_vec[0]) && theta_vec[0] > 0.0) {
            param1 = theta_vec[0];
      } else {
            // Default epsilon
            param1 = (wc == WeightCode::INVERSE_CLIM) ? 1e-12 : 1e-6;
      }
      break;
      }

      case WeightCode::GAUSSIAN_CLIM:
      case WeightCode::GAUSSIAN_GEOG: {
            // param1 = -1/(2*sigma^2)
            double sigma = 1.0;
            if (theta_vec.size() > 0 && std::isfinite(theta_vec[0]) && theta_vec[0] > 0.0) {
                  sigma = theta_vec[0];
            }
            param1 = -1.0 / (2.0 * sigma * sigma);
            break;
      }

      case WeightCode::GAUSSIAN_JOINT: {
            // param1 = -1/(2*sigma_clim^2), param2 = -1/(2*sigma_geog^2)
            double sigma_clim = 1.0;
            double sigma_geog = 1.0;
            if (theta_vec.size() >= 2) {
                  if (std::isfinite(theta_vec[0]) && theta_vec[0] > 0.0) {
                        sigma_clim = theta_vec[0];
                  }
                  if (std::isfinite(theta_vec[1]) && theta_vec[1] > 0.0) {
                        sigma_geog = theta_vec[1];
                  }
            }
            param1 = -1.0 / (2.0 * sigma_clim * sigma_clim);
            param2 = -1.0 / (2.0 * sigma_geog * sigma_geog);
            break;
      }

      case WeightCode::INVERSE_JOINT: {
            // param1 = eps_clim, param2 = eps_geog
            double eps_clim = 1e-12;
            double eps_geog = 1e-6;
            if (theta_vec.size() >= 2) {
                  if (std::isfinite(theta_vec[0]) && theta_vec[0] > 0.0) {
                        eps_clim = theta_vec[0];
                  }
                  if (std::isfinite(theta_vec[1]) && theta_vec[1] > 0.0) {
                        eps_geog = theta_vec[1];
                  }
            }
            param1 = eps_clim;
            param2 = eps_geog;
            break;
      }
      }

      return std::make_pair(param1, param2);
}


// -------------------------------------------------------------------------
// Worker for pair-returning modes: knn_clim, knn_geog, all
// Writes into out_indices[i] a vector<int> of 1-based ref indices.
// -------------------------------------------------------------------------
struct PairWorker : public Worker {
      const double* focal_ptr;      // focal matrix (original coords + clim)
      const double* ref_ptr;        // reference matrix (original coords + clim)
      const double* ref_latt_ptr;   // matrix used by lattice (may be ref_ptr or ECEF+clim)
      int n_focal;
      int n_ref;
      int n_clim;
      int stride_f;
      int stride_r;                 // stride for ref_ptr (n_ref)
      int stride_latt_r;            // stride for ref_latt_ptr

      bool use_lattice;
      bool use_geog_filter;
      bool use_haversine;
      bool use_scalar_clim;
      bool use_pervar_clim;
      bool use_ecef;                // NEW: true when lonlat + geog filtering → use ECEF

      double max_clim_scalar;
      double max_geog;
      double max_geog_chord;        // NEW: chord distance threshold for ECEF mode
      std::vector<double> max_clim_pervar;

      ModeCode mcode;
      int k;                        // k for kNN modes (>=1), ignored for ALL
      Lattice* lattice_ptr;         // may be nullptr if !use_lattice

      double R_earth;               // Earth radius (km)

      std::vector< std::vector<int> >& out_indices;

      PairWorker(const NumericMatrix& focal_mm,
                 const NumericMatrix& ref_mm,
                 const double* ref_latt_ptr_,
                 int stride_latt_r_,
                 bool use_lattice_,
                 bool use_geog_filter_,
                 bool use_haversine_,
                 bool use_scalar_clim_,
                 bool use_pervar_clim_,
                 double max_clim_scalar_,
                 double max_geog_,
                 double max_geog_chord_,
                 const std::vector<double>& max_clim_pervar_,
                 ModeCode mcode_,
                 int k_,
                 Lattice* lattice_ptr_,
                 bool use_ecef_,
                 double R_earth_,
                 std::vector< std::vector<int> >& out_indices_)
            : focal_ptr(REAL(focal_mm)),
              ref_ptr(REAL(ref_mm)),
              ref_latt_ptr(ref_latt_ptr_),
              n_focal(focal_mm.nrow()),
              n_ref(ref_mm.nrow()),
              n_clim(focal_mm.ncol() - 2),
              stride_f(focal_mm.nrow()),
              stride_r(ref_mm.nrow()),
              stride_latt_r(stride_latt_r_),
              use_lattice(use_lattice_),
              use_geog_filter(use_geog_filter_),
              use_haversine(use_haversine_),
              use_scalar_clim(use_scalar_clim_),
              use_pervar_clim(use_pervar_clim_),
              use_ecef(use_ecef_),
              max_clim_scalar(max_clim_scalar_),
              max_geog(max_geog_),
              max_geog_chord(max_geog_chord_),
              max_clim_pervar(max_clim_pervar_),
              mcode(mcode_),
              k(k_),
              lattice_ptr(lattice_ptr_),
              R_earth(R_earth_),
              out_indices(out_indices_)
      {}

      void operator()(std::size_t begin, std::size_t end) {
            const bool knn_mode      = (mcode == ModeCode::KNN_CLIM || mcode == ModeCode::KNN_GEOG);
            const bool rank_by_clim  = (mcode == ModeCode::KNN_CLIM);
            const bool rank_by_geog  = (mcode == ModeCode::KNN_GEOG);

            for (std::size_t i = begin; i < end; ++i) {
                  const double fx = focal_ptr[i];
                  const double fy = focal_ptr[i + stride_f];
                  const double* f_clim_col = focal_ptr + i + 2 * stride_f;

                  // -----------------------------------------------------------------
                  // Fast path: KNN modes + lattice + (planar OR ECEF)
                  // -----------------------------------------------------------------
                  if (knn_mode &&
                      use_lattice &&
                      lattice_ptr != nullptr &&
                      (!use_haversine || use_ecef)) {

                        const size_tu n_geo = lattice_ptr->n_geo_dims;
                        std::vector<double> fgeo_vec(n_geo);

                        if (use_ecef) {
                              // Convert lon/lat (degrees) to ECEF (km)
                              double X, Y, Z;
                              lonlat_to_ecef(fx, fy, R_earth, X, Y, Z);
                              fgeo_vec[0] = X;
                              if (n_geo > 1) fgeo_vec[1] = Y;
                              if (n_geo > 2) fgeo_vec[2] = Z;
                        } else {
                              // Planar projected
                              fgeo_vec[0] = fx;
                              if (n_geo > 1) fgeo_vec[1] = fy;
                        }

                        std::vector<double> fclim_vec(n_clim);
                        for (int kdim = 0; kdim < n_clim; ++kdim) {
                              fclim_vec[kdim] = f_clim_col[kdim * stride_f];
                        }

                        // Use chord distance threshold for ECEF, regular for planar
                        const double geog_thresh = use_ecef ? max_geog_chord : max_geog;

                        std::vector<index_t> cand0;
                        lattice_ptr->knn_query(
                                    fgeo_vec.data(),
                                    fclim_vec.data(),
                                    ref_latt_ptr,
                                    static_cast<size_tu>(stride_latt_r),
                                    static_cast<size_tu>(n_clim),
                                    /*rank_by_geog*/ rank_by_geog,
                                    geog_thresh,
                                    use_scalar_clim,
                                    max_clim_pervar,
                                    max_clim_scalar,
                                    k,
                                    cand0
                        );

                        // Convert 0-based indices to 1-based for R
                        std::vector<int> keep(cand0.size());
                        for (std::size_t t = 0; t < cand0.size(); ++t) {
                              keep[t] = static_cast<int>(cand0[t]) + 1;
                        }
                        out_indices[i] = std::move(keep);
                        continue;
                  }

                  // -----------------------------------------------------------------
                  // Otherwise: use existing candidate generation (lattice query or brute)
                  // -----------------------------------------------------------------
                  std::vector<index_t> cand;

                  if (use_lattice && lattice_ptr != nullptr) {
                        cand.reserve(128);

                        // Query lattice - use ECEF coords if applicable
                        double q_geo[3]; // max 3 for ECEF
                        if (use_ecef) {
                              double X, Y, Z;
                              lonlat_to_ecef(fx, fy, R_earth, X, Y, Z);
                              q_geo[0] = X;
                              q_geo[1] = Y;
                              q_geo[2] = Z;
                        } else {
                              q_geo[0] = fx;
                              q_geo[1] = fy;
                        }

                        std::vector<double> q_clim(n_clim);
                        for (int kdim = 0; kdim < n_clim; ++kdim) {
                              q_clim[kdim] = f_clim_col[kdim * stride_f];
                        }

                        const double geog_thresh = use_ecef ? max_geog_chord : max_geog;

                        lattice_ptr->query(q_geo,
                                           q_clim.data(),
                                           geog_thresh,
                                           use_scalar_clim,
                                           max_clim_pervar,
                                           use_scalar_clim ? max_clim_scalar
                                                 : std::numeric_limits<double>::infinity(),
                                                   cand);
                  } else {
                        cand.reserve(n_ref);
                        for (int j = 0; j < n_ref; ++j) {
                              cand.push_back(static_cast<index_t>(j));
                        }
                  }

                  // ALL mode: keep all matches passing filters
                  if (mcode == ModeCode::ALL) {
                        std::vector<int> keep;
                        keep.reserve(cand.size());

                        // ECEF focal coords (if needed)
                        double fx_ecef = 0, fy_ecef = 0, fz_ecef = 0;
                        if (use_ecef) {
                              lonlat_to_ecef(fx, fy, R_earth, fx_ecef, fy_ecef, fz_ecef);
                        }

                        for (size_t t = 0; t < cand.size(); ++t) {
                              const int j = static_cast<int>(cand[t]);
                              const double rx = ref_ptr[j];
                              const double ry = ref_ptr[j + stride_r];

                              // Geog filter
                              if (use_geog_filter) {
                                    double gdist;
                                    if (use_ecef) {
                                          // Use ECEF chord distance
                                          const double rx_ecef = ref_latt_ptr[j];
                                          const double ry_ecef = ref_latt_ptr[j + stride_latt_r];
                                          const double rz_ecef = ref_latt_ptr[j + 2 * stride_latt_r];
                                          const double dx = fx_ecef - rx_ecef;
                                          const double dy = fy_ecef - ry_ecef;
                                          const double dz = fz_ecef - rz_ecef;
                                          gdist = std::sqrt(dx*dx + dy*dy + dz*dz);
                                          if (gdist > max_geog_chord) continue;
                                    } else {
                                          gdist = geo_distance_km(fx, fy, rx, ry, use_haversine);
                                          if (gdist > max_geog) continue;
                                    }
                              }

                              // Climate checks
                              const double* r_clim_col = ref_ptr + j + 2 * stride_r;
                              const auto okd = clim_ok_and_dist(
                                    f_clim_col, r_clim_col,
                                    n_clim, stride_f, stride_r,
                                    use_pervar_clim, max_clim_pervar,
                                    use_scalar_clim, max_clim_scalar
                              );
                              if (!okd.first) continue;

                              keep.push_back(j + 1); // 1-based
                        }

                        out_indices[i] = std::move(keep);
                        continue;
                  }

                  // kNN modes (Climate or Geog) without expanding-search lattice path
                  using Neighbor = std::pair<double, int>; // (key_dist, ref_index_1based)
                  auto cmp = [](const Neighbor& a, const Neighbor& b) {
                        return a.first < b.first; // max-heap: top has largest distance
                  };
                  std::priority_queue<Neighbor, std::vector<Neighbor>, decltype(cmp)> pq(cmp);

                  // ECEF focal coords (if needed)
                  double fx_ecef = 0, fy_ecef = 0, fz_ecef = 0;
                  if (use_ecef) {
                        lonlat_to_ecef(fx, fy, R_earth, fx_ecef, fy_ecef, fz_ecef);
                  }

                  for (size_t t = 0; t < cand.size(); ++t) {
                        const int j = static_cast<int>(cand[t]);
                        const double rx = ref_ptr[j];
                        const double ry = ref_ptr[j + stride_r];

                        // Geog distance & filter
                        double gdist;
                        if (use_ecef) {
                              // Use ECEF chord distance
                              const double rx_ecef = ref_latt_ptr[j];
                              const double ry_ecef = ref_latt_ptr[j + stride_latt_r];
                              const double rz_ecef = ref_latt_ptr[j + 2 * stride_latt_r];
                              const double dx = fx_ecef - rx_ecef;
                              const double dy = fy_ecef - ry_ecef;
                              const double dz = fz_ecef - rz_ecef;
                              gdist = std::sqrt(dx*dx + dy*dy + dz*dz);
                              if (use_geog_filter && gdist > max_geog_chord) continue;
                        } else {
                              gdist = geo_distance_km(fx, fy, rx, ry, use_haversine);
                              if (use_geog_filter && gdist > max_geog) continue;
                        }

                        // Climate checks & distance
                        const double* r_clim_col = ref_ptr + j + 2 * stride_r;
                        const auto okd = clim_ok_and_dist(
                              f_clim_col, r_clim_col,
                              n_clim, stride_f, stride_r,
                              use_pervar_clim, max_clim_pervar,
                              use_scalar_clim, max_clim_scalar
                        );
                        if (!okd.first) continue;
                        const double clim_dist = okd.second;

                        const double key = rank_by_clim ? clim_dist : gdist;
                        const int ref_index_1based = j + 1;

                        if (static_cast<int>(pq.size()) < k) {
                              pq.emplace(key, ref_index_1based);
                        } else if (!pq.empty() && key < pq.top().first) {
                              pq.pop();
                              pq.emplace(key, ref_index_1based);
                        }
                  }

                  const int m = static_cast<int>(pq.size());
                  std::vector<int> idx_vec(m);
                  for (int pos = m - 1; pos >= 0; --pos) {
                        const Neighbor& nb = pq.top();
                        idx_vec[pos] = nb.second;
                        pq.pop();
                  }
                  out_indices[i] = std::move(idx_vec);
            }
      }
};


// -------------------------------------------------------------------------
// Worker for aggregate modes: COUNT / SUM / MEAN
// Writes into agg[i] the scalar aggregate for focal i.
// -------------------------------------------------------------------------
struct AggWorker : public Worker {
      const double* focal_ptr;
      const double* ref_ptr;
      const double* ref_latt_ptr;   // ECEF+clim matrix if use_ecef, else ref_ptr
      int n_focal;
      int n_ref;
      int n_clim;
      int stride_f;
      int stride_r;
      int stride_latt_r;

      bool use_lattice;
      bool use_geog_filter;
      bool use_haversine;
      bool use_scalar_clim;
      bool use_pervar_clim;
      bool use_ecef;                // NEW: use ECEF for geog filtering

      double max_clim_scalar;
      double max_geog;
      double max_geog_chord;        // NEW: chord threshold for ECEF
      std::vector<double> max_clim_pervar;

      ModeCode mcode;
      WeightCode wcode;
      double weight_param1;         // Pre-computed weight parameter 1
      double weight_param2;         // Pre-computed weight parameter 2

      Lattice* lattice_ptr;
      double R_earth;

      std::vector<double>& agg; // output

      AggWorker(const NumericMatrix& focal_mm,
                const NumericMatrix& ref_mm,
                const double* ref_latt_ptr_,
                int stride_latt_r_,
                bool use_lattice_,
                bool use_geog_filter_,
                bool use_haversine_,
                bool use_scalar_clim_,
                bool use_pervar_clim_,
                double max_clim_scalar_,
                double max_geog_,
                double max_geog_chord_,
                const std::vector<double>& max_clim_pervar_,
                ModeCode mcode_,
                WeightCode wcode_,
                double weight_param1_,
                double weight_param2_,
                Lattice* lattice_ptr_,
                bool use_ecef_,
                double R_earth_,
                std::vector<double>& agg_)
            : focal_ptr(REAL(focal_mm)),
              ref_ptr(REAL(ref_mm)),
              ref_latt_ptr(ref_latt_ptr_),
              n_focal(focal_mm.nrow()),
              n_ref(ref_mm.nrow()),
              n_clim(focal_mm.ncol() - 2),
              stride_f(focal_mm.nrow()),
              stride_r(ref_mm.nrow()),
              stride_latt_r(stride_latt_r_),
              use_lattice(use_lattice_),
              use_geog_filter(use_geog_filter_),
              use_haversine(use_haversine_),
              use_scalar_clim(use_scalar_clim_),
              use_pervar_clim(use_pervar_clim_),
              use_ecef(use_ecef_),
              max_clim_scalar(max_clim_scalar_),
              max_geog(max_geog_),
              max_geog_chord(max_geog_chord_),
              max_clim_pervar(max_clim_pervar_),
              mcode(mcode_),
              wcode(wcode_),
              weight_param1(weight_param1_),
              weight_param2(weight_param2_),
              lattice_ptr(lattice_ptr_),
              R_earth(R_earth_),
              agg(agg_)
      {}

      void operator()(std::size_t begin, std::size_t end) {
            for (std::size_t i = begin; i < end; ++i) {
                  const double fx = focal_ptr[i];
                  const double fy = focal_ptr[i + stride_f];
                  const double* f_clim_col = focal_ptr + i + 2 * stride_f;

                  std::vector<index_t> cand;

                  if (use_lattice && lattice_ptr != nullptr) {
                        cand.reserve(128);

                        const size_tu n_geo = lattice_ptr->n_geo_dims;
                        std::vector<double> q_geo(n_geo);

                        if (use_ecef && n_geo == 3) {
                              // Lattice was built in ECEF; query in ECEF.
                              double X, Y, Z;
                              lonlat_to_ecef(fx, fy, R_earth, X, Y, Z);
                              q_geo[0] = X;
                              q_geo[1] = Y;
                              q_geo[2] = Z;
                        } else {
                              // Planar / non-ECEF lattice: use raw x,y.
                              q_geo[0] = fx;
                              if (n_geo > 1) q_geo[1] = fy;
                        }

                        std::vector<double> q_clim(n_clim);
                        for (int kdim = 0; kdim < n_clim; ++kdim) {
                              q_clim[kdim] = f_clim_col[kdim * stride_f];
                        }

                        const double geog_thresh = use_ecef ? max_geog_chord : max_geog;

                        lattice_ptr->query(
                                    q_geo.data(),
                                    q_clim.data(),
                                    geog_thresh,
                                    use_scalar_clim,
                                    max_clim_pervar,
                                    use_scalar_clim
                                    ? max_clim_scalar
                        : std::numeric_limits<double>::infinity(),
                          cand
                        );
                  } else {
                        cand.reserve(n_ref);
                        for (int j = 0; j < n_ref; ++j) {
                              cand.push_back(static_cast<index_t>(j));
                        }
                  }

                  // Aggregate over candidates
                  double sum = 0.0;
                  double sum_w = 0.0;
                  int    count = 0;

                  // ECEF focal coords (if needed)
                  double fx_ecef = 0, fy_ecef = 0, fz_ecef = 0;
                  if (use_ecef) {
                        lonlat_to_ecef(fx, fy, R_earth, fx_ecef, fy_ecef, fz_ecef);
                  }

                  for (size_t t = 0; t < cand.size(); ++t) {
                        const int j = static_cast<int>(cand[t]);
                        const double rx = ref_ptr[j];
                        const double ry = ref_ptr[j + stride_r];

                        // Geog filter (in original space: planar, Haversine, or ECEF)
                        double gdist;
                        if (use_ecef) {
                              // Use ECEF chord distance
                              const double rx_ecef = ref_latt_ptr[j];
                              const double ry_ecef = ref_latt_ptr[j + stride_latt_r];
                              const double rz_ecef = ref_latt_ptr[j + 2 * stride_latt_r];
                              const double dx = fx_ecef - rx_ecef;
                              const double dy = fy_ecef - ry_ecef;
                              const double dz = fz_ecef - rz_ecef;
                              gdist = std::sqrt(dx*dx + dy*dy + dz*dz);
                              if (use_geog_filter && gdist > max_geog_chord) continue;
                        } else {
                              gdist = geo_distance_km(fx, fy, rx, ry, use_haversine);
                              if (use_geog_filter && gdist > max_geog) continue;
                        }

                        // Climate checks & distance
                        const double* r_clim_col = ref_ptr + j + 2 * stride_r;
                        const auto okd = clim_ok_and_dist(
                              f_clim_col, r_clim_col,
                              n_clim, stride_f, stride_r,
                              use_pervar_clim, max_clim_pervar,
                              use_scalar_clim, max_clim_scalar
                        );
                        if (!okd.first) continue;
                        const double clim_dist = okd.second;

                        ++count;

                        if (mcode == ModeCode::SUM || mcode == ModeCode::MEAN) {
                              double w = weight_from_codes(wcode, clim_dist, gdist,
                                                           weight_param1, weight_param2);
                              sum   += w * 1.0;
                              sum_w += w;
                        }

                  }

                  double val = NA_REAL;
                  if (mcode == ModeCode::COUNT) {
                        val = static_cast<double>(count);
                  } else if ((mcode == ModeCode::SUM || mcode == ModeCode::MEAN) &&
                        sum_w > 0.0) {
                        double mu = sum / sum_w;
                        if (mcode == ModeCode::SUM) {
                              val = sum;
                        } else { // MEAN
                              val = mu;
                        }
                  }

                  agg[i] = val;
            }
      }
};


} // anonymous namespace


// =========================================================================
// NEW: Build analog index (lattice) and return as Xptr
// =========================================================================

// [[Rcpp::export]]
SEXP build_analog_index_cpp(const NumericMatrix& ref_mm,
                            const std::string& coord_type,
                            int index_res)
{
      const int n_ref = ref_mm.nrow();
      const int ncol_ref = ref_mm.ncol();

      if (ncol_ref < 3) {
            stop("Reference data must have at least 2 coordinate columns and 1 climate variable");
      }

      const int n_clim = ncol_ref - 2;
      const bool use_haversine = (coord_type == "lonlat");
      const bool use_ecef = use_haversine; // Use ECEF for lonlat coords
      const double R_earth = 6371.0088;

      // Create lattice object
      Lattice* lattice = new Lattice();

      // Storage for coordinate ranges (for metadata)
      std::vector<double> coord_mins(2), coord_maxs(2);
      std::vector<double> clim_mins(n_clim), clim_maxs(n_clim);

      // Compute ranges
      const double* ref_ptr = REAL(ref_mm);
      const int stride = n_ref;

      for (int i = 0; i < n_ref; ++i) {
            // Coordinate ranges
            for (int d = 0; d < 2; ++d) {
                  double v = ref_ptr[i + d * stride];
                  if (i == 0) {
                        coord_mins[d] = coord_maxs[d] = v;
                  } else {
                        if (v < coord_mins[d]) coord_mins[d] = v;
                        if (v > coord_maxs[d]) coord_maxs[d] = v;
                  }
            }

            // Climate ranges
            for (int k = 0; k < n_clim; ++k) {
                  double v = ref_ptr[i + (2 + k) * stride];
                  if (i == 0) {
                        clim_mins[k] = clim_maxs[k] = v;
                  } else {
                        if (v < clim_mins[k]) clim_mins[k] = v;
                        if (v > clim_maxs[k]) clim_maxs[k] = v;
                  }
            }
      }

      // Prepare data for lattice building
      const double* lattice_data_ptr;
      std::vector<double> ecef_storage;
      int lattice_stride;
      size_t n_geo_dims;

      if (use_ecef) {
            // Build ECEF (X,Y,Z) + climate lattice
            n_geo_dims = 3;
            const size_t stride_e = static_cast<size_t>(n_ref);
            const size_t n_cols = n_geo_dims + static_cast<size_t>(n_clim);

            ecef_storage.assign(n_cols * n_ref, 0.0);

            for (size_t j = 0; j < static_cast<size_t>(n_ref); ++j) {
                  const double lon = ref_ptr[j];
                  const double lat = ref_ptr[j + stride_e];

                  double X, Y, Z;
                  lonlat_to_ecef(lon, lat, R_earth, X, Y, Z);

                  ecef_storage[j] = X;
                  ecef_storage[j + stride_e] = Y;
                  ecef_storage[j + 2 * stride_e] = Z;
            }

            // Copy climate variables
            for (size_t k = 0; k < static_cast<size_t>(n_clim); ++k) {
                  const double* src_col = ref_ptr + (2 + k) * stride_e;
                  double* dst_col = ecef_storage.data() + (3 + k) * stride_e;
                  std::copy(src_col, src_col + n_ref, dst_col);
            }

            lattice_data_ptr = ecef_storage.data();
            lattice_stride = n_ref;
      } else {
            // Use original projected coordinates
            n_geo_dims = 2;
            lattice_data_ptr = ref_ptr;
            lattice_stride = n_ref;
      }

      // Build lattice (always index both geo and climate dimensions)
      MetricType metric = use_ecef ? MetricType::Chord3D : MetricType::Planar;

      lattice->build(
                  lattice_data_ptr,
                  static_cast<size_tu>(n_ref),
                  static_cast<size_tu>(n_geo_dims),
                  static_cast<size_tu>(n_clim),
                  static_cast<size_tu>(lattice_stride),
                  metric,
                  std::numeric_limits<double>::infinity(),  // No geog constraint at build time
                  false,                                     // No scalar clim constraint at build time
                  std::vector<double>(),                     // No pervar clim constraint at build time
                  std::numeric_limits<double>::infinity(),
                  index_res,
                  true,  // use_geo_lattice
                  true   // use_clim_lattice
      );

      // Create Xptr
      XPtr<Lattice> lattice_xptr(lattice, true);

      // Store ECEF data if needed (for queries)
      SEXP ecef_xptr = R_NilValue;
      if (use_ecef) {
            std::vector<double>* ecef_copy = new std::vector<double>(ecef_storage);
            ecef_xptr = XPtr< std::vector<double> >(ecef_copy, true);
      }

      // Return list with lattice Xptr and metadata
      List result = List::create(
            Named("lattice_xptr") = lattice_xptr,
            Named("ecef_xptr") = ecef_xptr,
            Named("coord_type") = coord_type,
            Named("n_ref") = n_ref,
            Named("n_clim") = n_clim,
            Named("index_res") = index_res,
            Named("coord_mins") = coord_mins,
            Named("coord_maxs") = coord_maxs,
            Named("clim_mins") = clim_mins,
            Named("clim_maxs") = clim_maxs,
            Named("use_ecef") = use_ecef,
            Named("total_bins") = static_cast<double>(lattice->total_bins),
            Named("n_cells_nonempty") = static_cast<double>(lattice->n_cells_nonempty),
            Named("min_cell_occ") = static_cast<double>(lattice->min_cell_occ),
            Named("max_cell_occ") = static_cast<double>(lattice->max_cell_occ)
      );

      return result;
}


// =========================================================================
// NEW: Query analog index with focal points
// =========================================================================

// [[Rcpp::export]]
SEXP query_analog_index_cpp(SEXP index_list,
                            const NumericMatrix& focal_mm,
                            const NumericMatrix& ref_mm,
                            int k,
                            const NumericVector& max_clim,
                            double max_geog,
                            int mode_code,
                            int weight_code,
                            const NumericVector& theta)
{
      // Extract lattice and metadata from index
      List idx = as<List>(index_list);

      XPtr<Lattice> lattice_xptr = idx["lattice_xptr"];
      Lattice* lattice_ptr = lattice_xptr.get();

      if (lattice_ptr == nullptr) {
            stop("Invalid lattice pointer in index");
      }

      std::string coord_type = as<std::string>(idx["coord_type"]);
      int n_ref = as<int>(idx["n_ref"]);
      int n_clim = as<int>(idx["n_clim"]);
      bool use_ecef = as<bool>(idx["use_ecef"]);

      // Validate dimensions
      if (ref_mm.nrow() != n_ref) {
            stop("Reference data dimensions don't match index");
      }
      if (ref_mm.ncol() - 2 != n_clim) {
            stop("Number of climate variables doesn't match index");
      }
      if (focal_mm.ncol() != ref_mm.ncol()) {
            stop("Focal and reference data must have same number of columns");
      }

      const int n_focal = focal_mm.nrow();
      const bool use_haversine = (coord_type == "lonlat");
      const double R_earth = 6371.0088;

      // Parse climate constraints
      bool use_scalar_clim = false;
      bool use_pervar_clim = false;
      double max_clim_scalar = std::numeric_limits<double>::infinity();
      std::vector<double> max_clim_pervar_std(n_clim, std::numeric_limits<double>::infinity());

      if (max_clim.size() == 1) {
            double v = max_clim[0];
            if (std::isfinite(v)) {
                  use_scalar_clim = true;
                  max_clim_scalar = v;
            }
      } else if (max_clim.size() == n_clim) {
            for (int i = 0; i < n_clim; ++i) {
                  max_clim_pervar_std[i] = max_clim[i];
            }
            use_pervar_clim = true;
      } else if (max_clim.size() > 1) {
            stop("max_clim must be length 1 or equal to number of climate variables");
      }

      const bool use_geog_filter = std::isfinite(max_geog);
      const double max_geog_chord = (use_ecef && std::isfinite(max_geog))
            ? km_to_chord(max_geog, R_earth)
                  : max_geog;

      const ModeCode mcode = static_cast<ModeCode>(mode_code);
      const WeightCode wcode = static_cast<WeightCode>(weight_code);

      // Pre-compute weight parameters for efficiency
      auto weight_params = precompute_weight_params(wcode, theta);
      double weight_param1 = weight_params.first;
      double weight_param2 = weight_params.second;

      // Get ECEF data pointer if applicable
      const double* ref_latt_ptr;
      int stride_latt_r;

      if (use_ecef) {
            SEXP ecef_sexp = idx["ecef_xptr"];
            if (ecef_sexp == R_NilValue) {
                  stop("ECEF data missing from index");
            }
            XPtr< std::vector<double> > ecef_xptr(ecef_sexp);
            if (ecef_xptr.get() == nullptr) {
                  stop("ECEF data pointer is null");
            }
            ref_latt_ptr = ecef_xptr->data();
            stride_latt_r = n_ref;
      } else {
            ref_latt_ptr = REAL(ref_mm);
            stride_latt_r = n_ref;
      }

      // Execute query using workers
      const bool return_pairs = (mcode == ModeCode::KNN_CLIM ||
                                 mcode == ModeCode::KNN_GEOG ||
                                 mcode == ModeCode::ALL);

      if (return_pairs) {
            const int k_knn = (mcode == ModeCode::ALL ? 0 : k);
            std::vector< std::vector<int> > out_indices(n_focal);

            PairWorker worker(focal_mm,
                              ref_mm,
                              ref_latt_ptr,
                              stride_latt_r,
                              true,  // use_lattice
                              use_geog_filter,
                              use_haversine,
                              use_scalar_clim,
                              use_pervar_clim,
                              max_clim_scalar,
                              max_geog,
                              max_geog_chord,
                              max_clim_pervar_std,
                              mcode,
                              k_knn,
                              lattice_ptr,
                              use_ecef,
                              R_earth,
                              out_indices);

            parallelFor(0, static_cast<std::size_t>(n_focal), worker);

            // Convert to List<IntegerVector>
            List out(n_focal);
            for (int i = 0; i < n_focal; ++i) {
                  const std::vector<int>& v = out_indices[i];
                  IntegerVector idx_vec(v.size());
                  for (std::size_t j = 0; j < v.size(); ++j) {
                        idx_vec[j] = v[j];
                  }
                  out[i] = idx_vec;
            }

            // Attach diagnostics
            out.attr("n_focal") = n_focal;
            out.attr("n_ref") = n_ref;
            out.attr("n_clim") = n_clim;
            out.attr("max_dist") = max_geog;
            out.attr("max_clim") = max_clim;
            out.attr("geo_mode") = coord_type;
            out.attr("binning_method") = use_ecef ? "lattice_ecef" : "lattice";
            out.attr("total_bins") = static_cast<double>(lattice_ptr->total_bins);
            out.attr("n_bins_nonempty") = static_cast<double>(lattice_ptr->n_cells_nonempty);
            out.attr("min_bin_occupancy") = static_cast<double>(lattice_ptr->min_cell_occ);
            out.attr("max_bin_occupancy") = static_cast<double>(lattice_ptr->max_cell_occ);

            double avg_bin_occupancy = 0.0;
            if (lattice_ptr->total_bins > 0) {
                  avg_bin_occupancy = static_cast<double>(lattice_ptr->n_points) /
                        static_cast<double>(lattice_ptr->total_bins);
            }
            out.attr("avg_bin_occupancy") = avg_bin_occupancy;

            double avg_nonempty_bin_occupancy = 0.0;
            if (lattice_ptr->n_cells_nonempty > 0) {
                  avg_nonempty_bin_occupancy = static_cast<double>(lattice_ptr->n_points) /
                        static_cast<double>(lattice_ptr->n_cells_nonempty);
            }
            out.attr("avg_nonempty_bin_occupancy") = avg_nonempty_bin_occupancy;

            return out;
      }

      // Aggregate modes
      std::vector<double> agg_vals(n_focal, NA_REAL);

      AggWorker aworker(focal_mm,
                        ref_mm,
                        ref_latt_ptr,
                        stride_latt_r,
                        true,  // use_lattice
                        use_geog_filter,
                        use_haversine,
                        use_scalar_clim,
                        use_pervar_clim,
                        max_clim_scalar,
                        max_geog,
                        max_geog_chord,
                        max_clim_pervar_std,
                        mcode,
                        wcode,
                        weight_param1,
                        weight_param2,
                        lattice_ptr,
                        use_ecef,
                        R_earth,
                        agg_vals);

      parallelFor(0, static_cast<std::size_t>(n_focal), aworker);

      NumericVector agg(n_focal);
      for (int i = 0; i < n_focal; ++i) {
            agg[i] = agg_vals[i];
      }

      agg.attr("n_focal") = n_focal;
      agg.attr("n_ref") = n_ref;
      agg.attr("n_clim") = n_clim;
      agg.attr("max_dist") = max_geog;
      agg.attr("max_clim") = max_clim;
      agg.attr("geo_mode") = coord_type;
      agg.attr("binning_method") = use_ecef ? "lattice_ecef" : "lattice";
      agg.attr("total_bins") = static_cast<double>(lattice_ptr->total_bins);
      agg.attr("n_bins_nonempty") = static_cast<double>(lattice_ptr->n_cells_nonempty);
      agg.attr("min_bin_occupancy") = static_cast<double>(lattice_ptr->min_cell_occ);
      agg.attr("max_bin_occupancy") = static_cast<double>(lattice_ptr->max_cell_occ);

      double avg_bin_occupancy = 0.0;
      if (lattice_ptr->total_bins > 0) {
            avg_bin_occupancy = static_cast<double>(lattice_ptr->n_points) /
                  static_cast<double>(lattice_ptr->total_bins);
      }
      agg.attr("avg_bin_occupancy") = avg_bin_occupancy;

      double avg_nonempty_bin_occupancy = 0.0;
      if (lattice_ptr->n_cells_nonempty > 0) {
            avg_nonempty_bin_occupancy = static_cast<double>(lattice_ptr->n_points) /
                  static_cast<double>(lattice_ptr->n_cells_nonempty);
      }
      agg.attr("avg_nonempty_bin_occupancy") = avg_nonempty_bin_occupancy;

      return agg;
}


// =========================================================================
// EXISTING: Original find_analogs_core (updated with new weight functions)
// =========================================================================

// [[Rcpp::export]]
SEXP find_analogs_core(const NumericMatrix& focal_mm,
                       const NumericMatrix& ref_mm,
                       int k,
                       const NumericVector& max_clim,
                       double max_geog,
                       const std::string& geo_mode,
                       int mode_code,
                       int weight_code,
                       const NumericVector& theta,
                       int lattice_res)
{
      PROFILE_START("TOTAL");

      const int n_focal = focal_mm.nrow();
      const int n_ref   = ref_mm.nrow();

      const int ncol_focal = focal_mm.ncol();
      const int ncol_ref   = ref_mm.ncol();
      if (ncol_focal != ncol_ref) {
            stop("focal and ref must have the same number of columns");
      }
      if (ncol_focal < 3) {
            stop("Need at least 2 coordinate columns and 1 climate variable");
      }

      // Parse geometry mode
      const bool use_haversine = (geo_mode == "lonlat");

      // Climate dims
      const int n_clim = ncol_focal - 2;

      // Interpret max_clim
      bool   use_scalar_clim = false;
      bool   use_pervar_clim = false;
      double max_clim_scalar = std::numeric_limits<double>::infinity();
      std::vector<double> max_clim_pervar_std(n_clim, std::numeric_limits<double>::infinity());

      if (max_clim.size() == 1) {
            double v = max_clim[0];
            if (std::isfinite(v)) {
                  use_scalar_clim = true;
                  max_clim_scalar = v;
            }
      } else if (max_clim.size() == n_clim) {
            for (int i = 0; i < n_clim; ++i) {
                  max_clim_pervar_std[i] = max_clim[i];
            }
            use_pervar_clim = true;
      } else if (max_clim.size() > 1) {
            stop("max_clim must be length 1 or equal to the number of climate variables");
      }

      // Geographic threshold
      const bool use_geog_filter = std::isfinite(max_geog);

      // Decide whether to use lattice index
      const bool use_lattice = (use_geog_filter || use_scalar_clim || use_pervar_clim);
      const ModeCode   mcode = static_cast<ModeCode>(mode_code);
      const WeightCode wcode = static_cast<WeightCode>(weight_code);

      // Pre-compute weight parameters for efficiency
      auto weight_params = precompute_weight_params(wcode, theta);
      double weight_param1 = weight_params.first;
      double weight_param2 = weight_params.second;

      // Dimension roles
      const bool geo_filter  = use_geog_filter;
      const bool clim_filter = (use_scalar_clim || use_pervar_clim);
      const bool geo_sort    = (mcode == ModeCode::KNN_GEOG);
      const bool clim_sort   = (mcode == ModeCode::KNN_CLIM);

      // Which dims participate in the lattice
      const bool use_geo_lattice  = (geo_filter  || geo_sort);
      const bool use_clim_lattice = (clim_filter || clim_sort);

      // Lattice resolution (bins per active dimension)
      int bins_per_dim = lattice_res;
      if (bins_per_dim <= 0) {
            bins_per_dim = 10;
      }

      // Decide if we should use ECEF
      // Use ECEF whenever we have lonlat coords AND geographic filtering
      const bool use_ecef = (use_haversine && (geo_sort || use_geog_filter));

      // Lattice and diagnostics
      Lattice lattice;
      double  total_bins = 1.0;
      double  avg_bin_occupancy = static_cast<double>(n_ref);
      double  min_bin_occupancy = static_cast<double>(n_ref);
      double  max_bin_occupancy = static_cast<double>(n_ref);
      double  n_bins_nonempty = 1.0;
      double  avg_nonempty_bin_occupancy = static_cast<double>(n_ref);
      std::string binning_method = "none";

      // Pointer and storage for coordinates used by the lattice.
      const double* ref_latt_ptr = REAL(ref_mm);
      int stride_latt_r = n_ref;
      std::vector<double> ref_latt_storage;

      const double R_earth = 6371.0088; // km

      // Convert km threshold to chord distance
      const double max_geog_chord = (use_ecef && std::isfinite(max_geog))
            ? km_to_chord(max_geog, R_earth)
                  : max_geog;

      if (use_lattice) {

            MetricType metric = MetricType::Planar;

            if (use_ecef) {
                  PROFILE_START("ECEF_CONV");

                  // Build ECEF (X,Y,Z in km) + climate lattice
                  const size_tu n_geo = 3;
                  const size_tu stride_e = static_cast<size_tu>(n_ref);
                  const size_tu n_cols = n_geo + static_cast<size_tu>(n_clim);

                  ref_latt_storage.assign(static_cast<std::size_t>(n_cols * n_ref), 0.0);

                  const double* ref_ptr = REAL(ref_mm);

                  for (size_tu j = 0; j < static_cast<size_tu>(n_ref); ++j) {
                        const double lon = ref_ptr[j];
                        const double lat = ref_ptr[j + stride_e];

                        double X, Y, Z;
                        lonlat_to_ecef(lon, lat, R_earth, X, Y, Z);

                        ref_latt_storage[j]                = X;
                        ref_latt_storage[j + stride_e]     = Y;
                        ref_latt_storage[j + 2 * stride_e] = Z;
                  }

                  // Copy climate variables
                  for (size_tu k = 0; k < static_cast<size_tu>(n_clim); ++k) {
                        const double* src_col = REAL(ref_mm) + (2 + k) * stride_e;
                        double* dst_col = ref_latt_storage.data() + (3 + k) * stride_e;
                        std::copy(src_col, src_col + n_ref, dst_col);
                  }


                  PROFILE_STOP("ECEF_CONV");

                  metric = MetricType::Chord3D;

                  PROFILE_START("LATTICE_BUILD");

                  lattice.build(ref_latt_storage.data(),
                                static_cast<size_tu>(n_ref),
                                static_cast<size_tu>(3),       // n_geo = 3
                                static_cast<size_tu>(n_clim),
                                stride_e,
                                metric,
                                max_geog_chord,  // Use chord threshold
                                use_scalar_clim,
                                max_clim_pervar_std,
                                use_scalar_clim ? max_clim_scalar
                                      : std::numeric_limits<double>::infinity(),
                                        bins_per_dim,
                                        use_geo_lattice,
                                        use_clim_lattice);

                  PROFILE_STOP("LATTICE_BUILD");

                  PROFILE_COUNT("LATTICE_TOTAL_BINS", lattice.total_bins);
                  PROFILE_COUNT("LATTICE_NONEMPTY_BINS", lattice.n_cells_nonempty);
                  PROFILE_COUNT("LATTICE_MAX_OCCUPANCY", lattice.max_cell_occ);

                  ref_latt_ptr  = ref_latt_storage.data();
                  stride_latt_r = n_ref;
                  binning_method = "lattice_ecef";


            } else {
                  // Normal (2D projected) lattice path
                  metric = MetricType::Planar;

                  PROFILE_START("LATTICE_BUILD");

                  lattice.build(REAL(ref_mm),
                                static_cast<size_tu>(n_ref),
                                static_cast<size_tu>(2),
                                static_cast<size_tu>(n_clim),
                                static_cast<size_tu>(n_ref),
                                metric,
                                max_geog,
                                use_scalar_clim,
                                max_clim_pervar_std,
                                use_scalar_clim ? max_clim_scalar
                                      : std::numeric_limits<double>::infinity(),
                                        bins_per_dim,
                                        use_geo_lattice,
                                        use_clim_lattice);

                  PROFILE_STOP("LATTICE_BUILD");

                  PROFILE_COUNT("LATTICE_TOTAL_BINS", lattice.total_bins);
                  PROFILE_COUNT("LATTICE_NONEMPTY_BINS", lattice.n_cells_nonempty);
                  PROFILE_COUNT("LATTICE_MAX_OCCUPANCY", lattice.max_cell_occ);

                  ref_latt_ptr  = REAL(ref_mm);
                  stride_latt_r = n_ref;
                  binning_method = "lattice";
            }

            // Diagnostics from lattice
            total_bins = static_cast<double>(lattice.total_bins);
            if (lattice.total_bins > 0) {
                  avg_bin_occupancy =
                        static_cast<double>(lattice.n_points) /
                              static_cast<double>(lattice.total_bins);
            }
            n_bins_nonempty = static_cast<double>(lattice.n_cells_nonempty);
            if (lattice.n_cells_nonempty > 0) {
                  avg_nonempty_bin_occupancy =
                        static_cast<double>(lattice.n_points) /
                              static_cast<double>(lattice.n_cells_nonempty);
            } else {
                  avg_nonempty_bin_occupancy = 0.0;
            }

            if (lattice.min_cell_occ == std::numeric_limits<size_tu>::max()) {
                  min_bin_occupancy = 0.0;
                  max_bin_occupancy = 0.0;
            } else {
                  min_bin_occupancy = static_cast<double>(lattice.min_cell_occ);
                  max_bin_occupancy = static_cast<double>(lattice.max_cell_occ);
            }

            // Record problem dimensions
            PROFILE_COUNT("N_FOCAL", n_focal);
            PROFILE_COUNT("N_REF", n_ref);
            PROFILE_COUNT("N_CLIM", n_clim);

      } else {
            binning_method = "none";
      }


      const bool return_pairs =
            (mcode == ModeCode::KNN_CLIM ||
            mcode == ModeCode::KNN_GEOG ||
            mcode == ModeCode::ALL);

      if (return_pairs) {

            const int k_knn = (mcode == ModeCode::ALL ? 0 : k);

            std::vector< std::vector<int> > out_indices(n_focal);

            PROFILE_START("PAIR_WORKERS");

            PairWorker worker(focal_mm,
                              ref_mm,
                              ref_latt_ptr,
                              stride_latt_r,
                              use_lattice,
                              use_geog_filter,
                              use_haversine,
                              use_scalar_clim,
                              use_pervar_clim,
                              max_clim_scalar,
                              max_geog,
                              max_geog_chord,
                              max_clim_pervar_std,
                              mcode,
                              k_knn,
                              use_lattice ? &lattice : nullptr,
                              use_ecef,
                              R_earth,
                              out_indices);

            parallelFor(0, static_cast<std::size_t>(n_focal), worker);

            // Quick diagnostic: which code path are focals taking?
            size_t fast_path_count = 0;
            size_t regular_path_count = 0;
            size_t no_analog_count = 0;

            for (int i = 0; i < n_focal; ++i) {
                  if (out_indices[i].empty()) {
                        no_analog_count++;
                  } else {
                        // Heuristic: fast path typically finds exactly k neighbors
                        if (out_indices[i].size() == static_cast<size_t>(k)) {
                              fast_path_count++;
                        } else {
                              regular_path_count++;
                        }
                  }
            }

            PROFILE_COUNT("FAST_PATH_FOCALS", fast_path_count);
            PROFILE_COUNT("REGULAR_PATH_FOCALS", regular_path_count);
            PROFILE_COUNT("NO_ANALOG_FOCALS", no_analog_count);

            PROFILE_STOP("PAIR_WORKERS");

            // Count statistics from results
            size_t total_pairs = 0;
            for (int i = 0; i < n_focal; ++i) {
                  total_pairs += out_indices[i].size();
            }

            PROFILE_COUNT("TOTAL_PAIRS", total_pairs);
            if (n_focal > 0) {
                  PROFILE_COUNT("AVG_PAIRS_PER_FOCAL", total_pairs / n_focal);
            }

            PROFILE_START("EMIT_PAIRS");

            // Convert to List<IntegerVector>
            List out(n_focal);
            for (int i = 0; i < n_focal; ++i) {
                  const std::vector<int>& v = out_indices[i];
                  IntegerVector idx_vec(v.size());
                  for (std::size_t j = 0; j < v.size(); ++j) {
                        idx_vec[j] = v[j];
                  }
                  out[i] = idx_vec;
            }

            // Attach diagnostics
            out.attr("n_focal") = n_focal;
            out.attr("n_ref")   = n_ref;
            out.attr("n_clim")  = n_clim;
            out.attr("max_dist") = max_geog;
            out.attr("max_clim") = max_clim;
            out.attr("geo_mode") = geo_mode;
            out.attr("binning_method")    = binning_method;
            out.attr("total_bins")        = total_bins;
            out.attr("avg_bin_occupancy") = avg_bin_occupancy;
            out.attr("min_bin_occupancy") = min_bin_occupancy;
            out.attr("max_bin_occupancy") = max_bin_occupancy;
            out.attr("n_bins_nonempty")               = n_bins_nonempty;
            out.attr("avg_nonempty_bin_occupancy")    = avg_nonempty_bin_occupancy;

            PROFILE_STOP("EMIT_PAIRS");

            PROFILE_STOP("TOTAL");

            return out;
      }

      // Aggregate modes: COUNT / SUM / MEAN
      std::vector<double> agg_vals(n_focal, NA_REAL);

      PROFILE_START("AGG_WORKERS");

      AggWorker aworker(focal_mm,
                        ref_mm,
                        ref_latt_ptr,
                        stride_latt_r,
                        use_lattice,
                        use_geog_filter,
                        use_haversine,
                        use_scalar_clim,
                        use_pervar_clim,
                        max_clim_scalar,
                        max_geog,
                        max_geog_chord,
                        max_clim_pervar_std,
                        mcode,
                        wcode,
                        weight_param1,
                        weight_param2,
                        use_lattice ? &lattice : nullptr,
                        use_ecef,
                        R_earth,
                        agg_vals);

      parallelFor(0, static_cast<std::size_t>(n_focal), aworker);

      PROFILE_STOP("AGG_WORKERS");

      PROFILE_START("FORMAT_AGG");

      NumericVector agg(n_focal);
      for (int i = 0; i < n_focal; ++i) {
            agg[i] = agg_vals[i];
      }

      agg.attr("n_focal") = n_focal;
      agg.attr("n_ref")   = n_ref;
      agg.attr("n_clim")  = n_clim;
      agg.attr("max_dist") = max_geog;
      agg.attr("max_clim") = max_clim;
      agg.attr("geo_mode") = geo_mode;
      agg.attr("binning_method")    = binning_method;
      agg.attr("total_bins")        = total_bins;
      agg.attr("avg_bin_occupancy") = avg_bin_occupancy;
      agg.attr("min_bin_occupancy") = min_bin_occupancy;
      agg.attr("max_bin_occupancy") = max_bin_occupancy;
      agg.attr("n_bins_nonempty")               = n_bins_nonempty;
      agg.attr("avg_nonempty_bin_occupancy")    = avg_nonempty_bin_occupancy;

      PROFILE_STOP("FORMAT_AGG");

      PROFILE_STOP("TOTAL");

      return agg;
}


// [[Rcpp::export]]
Rcpp::List profile_find_analogs(const NumericMatrix& focal_mm,
                                const NumericMatrix& ref_mm,
                                int k,
                                const NumericVector& max_clim,
                                double max_geog,
                                const std::string& geo_mode,
                                int mode_code,
                                int weight_code,
                                const NumericVector& theta,
                                int lattice_res,
                                bool enable_profiling)
{
      if (enable_profiling) {
            analogs::g_profiler.enable();
      }

      SEXP result = find_analogs_core(focal_mm, ref_mm, k, max_clim, max_geog,
                                      geo_mode, mode_code, weight_code, theta, lattice_res);

      // Collect results
      auto events = analogs::g_profiler.get_events();
      auto counters = analogs::g_profiler.get_counters();

      std::vector<std::string> event_names;
      std::vector<double> event_times;
      std::vector<int> event_counts;

      for (const auto& kv : events) {
            event_names.push_back(kv.first);
            event_times.push_back(kv.second.duration_ms);
            event_counts.push_back(static_cast<int>(kv.second.count));
      }

      std::vector<std::string> counter_names;
      std::vector<double> counter_values;

      for (const auto& kv : counters) {
            counter_names.push_back(kv.first);
            counter_values.push_back(static_cast<double>(kv.second));
      }

      analogs::g_profiler.clear();
      analogs::g_profiler.disable();

      return Rcpp::List::create(
            Rcpp::Named("result") = result,
            Rcpp::Named("profile") = Rcpp::List::create(
                  Rcpp::Named("event_names") = event_names,
                  Rcpp::Named("event_times_ms") = event_times,
                  Rcpp::Named("event_counts") = event_counts,
                  Rcpp::Named("counter_names") = counter_names,
                  Rcpp::Named("counter_values") = counter_values
            )
      );
}
