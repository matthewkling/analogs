// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>

#include "lattice.hpp"
#include "metrics.hpp"
#include "types.hpp"

#include <queue>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#include <chrono>

using namespace Rcpp;
using namespace analogs;
using namespace RcppParallel;

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
      NONE          = 0, // not used (knn/all/count)
            UNIFORM       = 1,
            INVERSE_CLIM  = 2,
            INVERSE_GEOG  = 3
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

inline double weight_from_codes(WeightCode wc,
                                double clim_dist,
                                double geog_dist,
                                double theta)
{
      switch (wc) {
      case WeightCode::UNIFORM:
            return 1.0;
      case WeightCode::INVERSE_CLIM: {
            const double eps = (theta > 0.0 && std::isfinite(theta)) ? theta : 1e-12;
            return 1.0 / (clim_dist + eps);
      }
      case WeightCode::INVERSE_GEOG: {
            const double eps = (theta > 0.0 && std::isfinite(theta)) ? theta : 1e-6;
            return 1.0 / (geog_dist + eps);
      }
      default: // NONE or unknown
            return 1.0;
      }
}

// -------------------------------------------------------------------------
// Worker for pair-returning modes: knn_clim, knn_geog, all
// Writes into out_indices[i] a vector<int> of 1-based ref indices.
// -------------------------------------------------------------------------
struct PairWorker : public Worker {
      const double* focal_ptr;
      const double* ref_ptr;
      int n_focal;
      int n_ref;
      int n_clim;
      int stride_f;
      int stride_r;

      bool use_lattice;
      bool use_geog_filter;
      bool use_haversine;
      bool use_scalar_clim;
      bool use_pervar_clim;

      double max_clim_scalar;
      double max_geog;
      std::vector<double> max_clim_pervar;

      ModeCode mcode;
      int k;                    // k for kNN modes (>=1), ignored for ALL
      Lattice* lattice_ptr;     // may be nullptr if !use_lattice

      std::vector< std::vector<int> >& out_indices;

      PairWorker(const NumericMatrix& focal_mm,
                 const NumericMatrix& ref_mm,
                 bool use_lattice_,
                 bool use_geog_filter_,
                 bool use_haversine_,
                 bool use_scalar_clim_,
                 bool use_pervar_clim_,
                 double max_clim_scalar_,
                 double max_geog_,
                 const std::vector<double>& max_clim_pervar_,
                 ModeCode mcode_,
                 int k_,
                 Lattice* lattice_ptr_,
                 std::vector< std::vector<int> >& out_indices_)
            : focal_ptr(REAL(focal_mm)),
              ref_ptr(REAL(ref_mm)),
              n_focal(focal_mm.nrow()),
              n_ref(ref_mm.nrow()),
              n_clim(focal_mm.ncol() - 2),
              stride_f(focal_mm.nrow()),
              stride_r(ref_mm.nrow()),
              use_lattice(use_lattice_),
              use_geog_filter(use_geog_filter_),
              use_haversine(use_haversine_),
              use_scalar_clim(use_scalar_clim_),
              use_pervar_clim(use_pervar_clim_),
              max_clim_scalar(max_clim_scalar_),
              max_geog(max_geog_),
              max_clim_pervar(max_clim_pervar_),
              mcode(mcode_),
              k(k_),
              lattice_ptr(lattice_ptr_),
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
                  // Fast path: KNN modes + lattice + projected coordinates
                  // Use Lattice::knn_query with expanding search over cells.
                  // -----------------------------------------------------------------
                  if (knn_mode &&
                      use_lattice &&
                      lattice_ptr != nullptr &&
                      !use_haversine) {  // planar only for expanding-search v1

                      double fgeo[2] = { fx, fy };
                        std::vector<double> fclim_vec(n_clim);
                        for (int kdim = 0; kdim < n_clim; ++kdim) {
                              fclim_vec[kdim] = f_clim_col[kdim * stride_f];
                        }

                        std::vector<index_t> cand0;
                        lattice_ptr->knn_query(
                                    fgeo,
                                    fclim_vec.data(),
                                    ref_ptr,
                                    static_cast<size_tu>(stride_r),
                                    static_cast<size_tu>(n_clim),
                                    /*rank_by_geog*/ rank_by_geog,
                                    max_geog,
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
                        double q_geo[2] = { fx, fy };
                        std::vector<double> q_clim(n_clim);
                        for (int kdim = 0; kdim < n_clim; ++kdim) {
                              q_clim[kdim] = f_clim_col[kdim * stride_f];
                        }
                        lattice_ptr->query(q_geo,
                                           q_clim.data(),
                                           max_geog,
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

                        for (size_t t = 0; t < cand.size(); ++t) {
                              const int j = static_cast<int>(cand[t]);
                              const double rx = ref_ptr[j];
                              const double ry = ref_ptr[j + stride_r];

                              // Geog filter
                              if (use_geog_filter) {
                                    const double gdist = geo_distance_km(fx, fy, rx, ry, use_haversine);
                                    if (gdist > max_geog) continue;
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

                  for (size_t t = 0; t < cand.size(); ++t) {
                        const int j = static_cast<int>(cand[t]);
                        const double rx = ref_ptr[j];
                        const double ry = ref_ptr[j + stride_r];

                        // Geog distance & filter
                        const double gdist = geo_distance_km(fx, fy, rx, ry, use_haversine);
                        if (use_geog_filter && gdist > max_geog) continue;

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
      int n_focal;
      int n_ref;
      int n_clim;
      int stride_f;
      int stride_r;

      bool use_lattice;
      bool use_geog_filter;
      bool use_haversine;
      bool use_scalar_clim;
      bool use_pervar_clim;

      double max_clim_scalar;
      double max_geog;
      std::vector<double> max_clim_pervar;

      ModeCode mcode;
      WeightCode wcode;
      double theta;

      Lattice* lattice_ptr;    // may be nullptr if !use_lattice

      std::vector<double>& agg; // output

      AggWorker(const NumericMatrix& focal_mm,
                const NumericMatrix& ref_mm,
                bool use_lattice_,
                bool use_geog_filter_,
                bool use_haversine_,
                bool use_scalar_clim_,
                bool use_pervar_clim_,
                double max_clim_scalar_,
                double max_geog_,
                const std::vector<double>& max_clim_pervar_,
                ModeCode mcode_,
                WeightCode wcode_,
                double theta_,
                Lattice* lattice_ptr_,
                std::vector<double>& agg_)
            : focal_ptr(REAL(focal_mm)),
              ref_ptr(REAL(ref_mm)),
              n_focal(focal_mm.nrow()),
              n_ref(ref_mm.nrow()),
              n_clim(focal_mm.ncol() - 2),
              stride_f(focal_mm.nrow()),
              stride_r(ref_mm.nrow()),
              use_lattice(use_lattice_),
              use_geog_filter(use_geog_filter_),
              use_haversine(use_haversine_),
              use_scalar_clim(use_scalar_clim_),
              use_pervar_clim(use_pervar_clim_),
              max_clim_scalar(max_clim_scalar_),
              max_geog(max_geog_),
              max_clim_pervar(max_clim_pervar_),
              mcode(mcode_),
              wcode(wcode_),
              theta(theta_),
              lattice_ptr(lattice_ptr_),
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
                        double q_geo[2] = { fx, fy };
                        std::vector<double> q_clim(n_clim);
                        for (int kdim = 0; kdim < n_clim; ++kdim) {
                              q_clim[kdim] = f_clim_col[kdim * stride_f];
                        }
                        lattice_ptr->query(q_geo,
                                           q_clim.data(),
                                           max_geog,
                                           use_scalar_clim,
                                           max_clim_pervar,
                                           use_scalar_clim ? max_clim_scalar : std::numeric_limits<double>::infinity(),
                                           cand);
                  } else {
                        cand.reserve(n_ref);
                        for (int j = 0; j < n_ref; ++j) cand.push_back(static_cast<index_t>(j));
                  }

                  double count = 0.0;
                  double sum_w = 0.0;

                  for (size_t t = 0; t < cand.size(); ++t) {
                        const int j = static_cast<int>(cand[t]);
                        const double rx = ref_ptr[j];
                        const double ry = ref_ptr[j + stride_r];

                        // Geog filter
                        const double gdist = geo_distance_km(fx, fy, rx, ry, use_haversine);
                        if (use_geog_filter && gdist > max_geog) continue;

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
                              const double w = weight_from_codes(wcode, clim_dist, gdist, theta);
                              sum_w += w;
                        }
                  }

                  if (mcode == ModeCode::COUNT) {
                        agg[i] = count;
                  } else if (mcode == ModeCode::SUM) {
                        agg[i] = sum_w;
                  } else { // MEAN
                        agg[i] = (count > 0.0) ? (sum_w / count) : 0.0;
                  }
            }
      }
};

} // anonymous namespace


// -------------------------------------------------------------------------
// [[Rcpp::export]]
SEXP find_analogs_core(const NumericMatrix& focal_mm,
                       const NumericMatrix& ref_mm,
                       int k,
                       const NumericVector& max_clim,
                       double max_geog,
                       const std::string& geo_mode,
                       int mode_code,
                       int weight_code,
                       double theta,
                       int lattice_res)
{
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
            bins_per_dim = 10; // simple default; will later be replaced by choose_res
      }

      // Build lattice if useful
      Lattice lattice;
      double  total_bins = 1.0;
      double  avg_bin_occupancy = static_cast<double>(n_ref);
      double  min_bin_occupancy = static_cast<double>(n_ref);
      double  max_bin_occupancy = static_cast<double>(n_ref);
      double  n_bins_nonempty = 1.0;
      double  avg_nonempty_bin_occupancy = static_cast<double>(n_ref);
      std::string binning_method = "none";

      if (use_lattice) {
            MetricType metric = use_haversine ? MetricType::Haversine
            : MetricType::Planar;
            lattice.build(REAL(ref_mm),
                          static_cast<size_tu>(n_ref),
                          static_cast<size_tu>(n_clim),
                          static_cast<size_tu>(n_ref), // stride_r
                          metric,
                          max_geog,
                          use_scalar_clim,
                          max_clim_pervar_std,
                          use_scalar_clim ? max_clim_scalar
                                : std::numeric_limits<double>::infinity(),
                                  bins_per_dim,
                                  use_geo_lattice,
                                  use_clim_lattice);

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
            binning_method = "multi_dim_lattice";

      }

      const bool return_pairs =
            (mcode == ModeCode::KNN_CLIM ||
            mcode == ModeCode::KNN_GEOG ||
            mcode == ModeCode::ALL);

      if (return_pairs) {
            // Use k as provided for kNN, ignore for ALL (we treat ALL separately)
            const int k_knn = (mcode == ModeCode::ALL ? 0 : k);

            // Parallel: compute neighbor indices per focal into std::vector-of-vectors
            std::vector< std::vector<int> > out_indices(n_focal);

            PairWorker worker(focal_mm,
                              ref_mm,
                              use_lattice,
                              use_geog_filter,
                              use_haversine,
                              use_scalar_clim,
                              use_pervar_clim,
                              max_clim_scalar,
                              max_geog,
                              max_clim_pervar_std,
                              mcode,
                              k_knn,
                              use_lattice ? &lattice : nullptr,
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
            out.attr("n_ref")   = n_ref;
            out.attr("n_clim")  = n_clim;
            out.attr("max_dist") = max_geog;    // keep legacy attr name
            out.attr("max_clim") = max_clim;
            out.attr("geo_mode") = geo_mode;
            out.attr("binning_method")    = binning_method;
            out.attr("total_bins")        = total_bins;
            out.attr("avg_bin_occupancy") = avg_bin_occupancy;
            out.attr("min_bin_occupancy") = min_bin_occupancy;
            out.attr("max_bin_occupancy") = max_bin_occupancy;
            out.attr("n_bins_nonempty")               = n_bins_nonempty;
            out.attr("avg_nonempty_bin_occupancy")    = avg_nonempty_bin_occupancy;

            return out;
      }

      // Aggregate modes: COUNT / SUM / MEAN
      std::vector<double> agg_vals(n_focal, NA_REAL);

      AggWorker aworker(focal_mm,
                        ref_mm,
                        use_lattice,
                        use_geog_filter,
                        use_haversine,
                        use_scalar_clim,
                        use_pervar_clim,
                        max_clim_scalar,
                        max_geog,
                        max_clim_pervar_std,
                        mcode,
                        wcode,
                        theta,
                        use_lattice ? &lattice : nullptr,
                        agg_vals);

      parallelFor(0, static_cast<std::size_t>(n_focal), aworker);

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

      return agg;
}



// -------------------------------------------------------------------------
// KNN lattice benchmark helper
//
// - Builds Lattice with the same logic as find_analogs_core
// - Requires geo_mode = "projected" so !use_haversine and KNN lattice path
// - For each focal, calls Lattice::knn_query, measures total KNN time
// - Returns lattice diagnostics and neighbor count summaries
// -------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List bench_knn_core(const NumericMatrix& focal_mm,
                          const NumericMatrix& ref_mm,
                          int k,
                          const NumericVector& max_clim,
                          double max_geog,
                          const std::string& geo_mode,
                          int mode_code)
{
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

      if (k <= 0) {
            stop("k must be positive for bench_knn_core");
      }

      // Parse geometry mode
      const bool use_haversine = (geo_mode == "lonlat");
      if (use_haversine) {
            stop("bench_knn_core currently supports only geo_mode = 'projected' (planar) for KNN lattice benchmarking");
      }

      // Climate dims
      const int n_clim = ncol_focal - 2;

      // Interpret max_clim exactly as in find_analogs_core
      bool   use_scalar_clim = false;
      bool   use_pervar_clim = false;
      double max_clim_scalar = std::numeric_limits<double>::infinity();
      std::vector<double> max_clim_pervar_std(
                  n_clim, std::numeric_limits<double>::infinity()
      );

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

      const bool use_geog_filter = std::isfinite(max_geog);
      const bool use_lattice     = (use_geog_filter || use_scalar_clim || use_pervar_clim);

      if (!use_lattice) {
            stop("bench_knn_core: lattice will not be used (no finite max_geog or max_clim); nothing to benchmark");
      }

      // Mode (must be one of the KNN modes)
      ModeCode mcode = static_cast<ModeCode>(mode_code);
      const bool knn_mode     = (mcode == ModeCode::KNN_CLIM || mcode == ModeCode::KNN_GEOG);
      const bool rank_by_geog = (mcode == ModeCode::KNN_GEOG);

      if (!knn_mode) {
            stop("bench_knn_core: mode_code must correspond to KNN_CLIM or KNN_GEOG");
      }

      const bool geo_filter  = use_geog_filter;
      const bool clim_filter = (use_scalar_clim || use_pervar_clim);
      const bool geo_sort    = (mcode == ModeCode::KNN_GEOG);
      const bool clim_sort   = (mcode == ModeCode::KNN_CLIM);

      const bool use_geo_lattice  = (geo_filter  || geo_sort);
      const bool use_clim_lattice = (clim_filter || clim_sort);

      int bins_per_dim = 10; // hard-coded for now; can be tuned later


      // Build lattice (planar metric)
      Lattice lattice;
      double  total_bins         = 1.0;
      double  avg_bin_occupancy  = static_cast<double>(n_ref);
      double  min_bin_occupancy  = static_cast<double>(n_ref);
      double  max_bin_occupancy  = static_cast<double>(n_ref);
      double  build_time_ms      = 0.0;
      double  knn_time_ms        = 0.0;

      {
            MetricType metric = MetricType::Planar;

            auto t0 = std::chrono::high_resolution_clock::now();
            lattice.build(REAL(ref_mm),
                          static_cast<size_tu>(n_ref),
                          static_cast<size_tu>(n_clim),
                          static_cast<size_tu>(n_ref), // stride_r
                          metric,
                          max_geog,
                          use_scalar_clim,
                          max_clim_pervar_std,
                          use_scalar_clim ? max_clim_scalar
                                : std::numeric_limits<double>::infinity(),
                                  bins_per_dim,
                                  use_geo_lattice,
                                  use_clim_lattice);
            auto t1 = std::chrono::high_resolution_clock::now();
            build_time_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

            total_bins = static_cast<double>(lattice.total_bins);
            if (lattice.total_bins > 0) {
                  avg_bin_occupancy =
                        static_cast<double>(lattice.n_points) /
                              static_cast<double>(lattice.total_bins);
            }
            if (lattice.min_cell_occ == std::numeric_limits<size_tu>::max()) {
                  min_bin_occupancy = 0.0;
                  max_bin_occupancy = 0.0;
            } else {
                  min_bin_occupancy = static_cast<double>(lattice.min_cell_occ);
                  max_bin_occupancy = static_cast<double>(lattice.max_cell_occ);
            }
      }

      // KNN over focal points (mirrors PairWorker fast path, but sequential)
      const double* focal_ptr = REAL(focal_mm);
      const double* ref_ptr   = REAL(ref_mm);
      const int     stride_f  = n_focal;
      const int     stride_r  = n_ref;

      std::vector<double> neighbor_counts(n_focal, 0.0);

      auto tq0 = std::chrono::high_resolution_clock::now();

      std::vector<index_t> cand0;
      cand0.reserve(static_cast<std::size_t>(k));

      for (int i = 0; i < n_focal; ++i) {
            const double fx = focal_ptr[i];
            const double fy = focal_ptr[i + stride_f];
            const double* f_clim_col = focal_ptr + i + 2 * stride_f;

            double fgeo[2] = { fx, fy };
            std::vector<double> fclim_vec(n_clim);
            for (int kdim = 0; kdim < n_clim; ++kdim) {
                  fclim_vec[kdim] = f_clim_col[kdim * stride_f];
            }

            cand0.clear();
            lattice.knn_query(
                  fgeo,
                  fclim_vec.data(),
                  ref_ptr,
                  static_cast<size_tu>(stride_r),
                  static_cast<size_tu>(n_clim),
                  /*rank_by_geog*/ rank_by_geog,
                  max_geog,
                  use_scalar_clim,
                  max_clim_pervar_std,
                  max_clim_scalar,
                  k,
                  cand0
            );

            neighbor_counts[i] = static_cast<double>(cand0.size());
      }

      auto tq1 = std::chrono::high_resolution_clock::now();
      knn_time_ms = std::chrono::duration<double, std::milli>(tq1 - tq0).count();

      // Summaries of neighbor counts (normally <= k)
      double mean_neighbors = 0.0;
      double p95_neighbors  = 0.0;
      double min_neighbors  = 0.0;
      double max_neighbors  = 0.0;

      if (n_focal > 0) {
            std::vector<double> tmp = neighbor_counts;
            double sum = 0.0;
            min_neighbors = neighbor_counts[0];
            max_neighbors = neighbor_counts[0];

            for (int i = 0; i < n_focal; ++i) {
                  const double v = neighbor_counts[i];
                  sum += v;
                  if (v < min_neighbors) min_neighbors = v;
                  if (v > max_neighbors) max_neighbors = v;
            }
            mean_neighbors = sum / static_cast<double>(n_focal);

            std::sort(tmp.begin(), tmp.end());
            std::size_t idx95 = static_cast<std::size_t>(
                  std::floor(0.95 * (static_cast<double>(n_focal) - 1.0))
            );
            p95_neighbors = tmp[idx95];
      }

      return Rcpp::List::create(
            Rcpp::Named("n_focal")            = n_focal,
            Rcpp::Named("n_ref")              = n_ref,
            Rcpp::Named("n_clim")             = n_clim,
            Rcpp::Named("k")                  = k,
            Rcpp::Named("geo_mode")           = geo_mode,
            Rcpp::Named("mode_code")          = mode_code,
            Rcpp::Named("use_lattice")        = use_lattice,
            Rcpp::Named("build_time_ms")      = build_time_ms,
            Rcpp::Named("knn_time_ms")        = knn_time_ms,
            Rcpp::Named("total_bins")         = total_bins,
            Rcpp::Named("avg_bin_occupancy")  = avg_bin_occupancy,
            Rcpp::Named("min_bin_occupancy")  = min_bin_occupancy,
            Rcpp::Named("max_bin_occupancy")  = max_bin_occupancy,
            Rcpp::Named("neighbor_counts")    = neighbor_counts,
            Rcpp::Named("mean_neighbors")     = mean_neighbors,
            Rcpp::Named("p95_neighbors")      = p95_neighbors,
            Rcpp::Named("min_neighbors")      = min_neighbors,
            Rcpp::Named("max_neighbors")      = max_neighbors
      );
}


// -------------------------------------------------------------------------
// Internal lattice benchmark helper
//
// - Builds lattice with same logic as find_analogs_core.
// - For each focal, calls lattice.query() and records candidate count
//   *before* the expensive distance checks.
// - Returns timing and sparsity diagnostics.
// -------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List bench_lattice_core(const NumericMatrix& focal_mm,
                              const NumericMatrix& ref_mm,
                              const NumericVector& max_clim,
                              double max_geog,
                              const std::string& geo_mode)
{
      const int n_focal = focal_mm.nrow();
      const int n_ref   = ref_mm.nrow();

      const int ncol_focal = focal_mm.ncol();
      const int ncol_ref   = ref_mm.ncol();
      if (ncol_focal != ncol_ref) {
            Rcpp::stop("focal and ref must have the same number of columns");
      }
      if (ncol_focal < 3) {
            Rcpp::stop("Need at least 2 coordinate columns and 1 climate variable");
      }

      const bool use_haversine = (geo_mode == "lonlat");
      const int  n_clim        = ncol_focal - 2;

      // Interpret max_clim (same logic as find_analogs_core)
      bool   use_scalar_clim = false;
      bool   use_pervar_clim = false;
      double max_clim_scalar = std::numeric_limits<double>::infinity();
      std::vector<double> max_clim_pervar_std(
                  n_clim, std::numeric_limits<double>::infinity()
      );

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
            Rcpp::stop("max_clim must be length 1 or equal to the number of climate variables");
      }

      const bool use_geog_filter = std::isfinite(max_geog);
      const bool use_lattice     = (use_geog_filter || use_scalar_clim || use_pervar_clim);

      Lattice lattice;
      const bool geo_filter  = use_geog_filter;
      const bool clim_filter = (use_scalar_clim || use_pervar_clim);
      const bool use_geo_lattice  = geo_filter;
      const bool use_clim_lattice = clim_filter;

      int bins_per_dim = 10; // hard-coded dev default

      double  total_bins = 1.0;
      double  avg_bin_occupancy = static_cast<double>(n_ref);
      double  min_bin_occupancy = static_cast<double>(n_ref);
      double  max_bin_occupancy = static_cast<double>(n_ref);
      double  n_bins_nonempty   = 1.0;
      double  avg_nonempty_bin_occupancy = static_cast<double>(n_ref);

      double build_time_ms = 0.0;
      double query_time_ms = 0.0;

      if (use_lattice) {
            MetricType metric = use_haversine ? MetricType::Haversine
            : MetricType::Planar;

            auto t0 = std::chrono::high_resolution_clock::now();
            lattice.build(REAL(ref_mm),
                          static_cast<size_tu>(n_ref),
                          static_cast<size_tu>(n_clim),
                          static_cast<size_tu>(n_ref), // stride_r
                          metric,
                          max_geog,
                          use_scalar_clim,
                          max_clim_pervar_std,
                          use_scalar_clim ? max_clim_scalar
                                : std::numeric_limits<double>::infinity(),
                                  bins_per_dim,
                                  use_geo_lattice,
                                  use_clim_lattice);

            auto t1 = std::chrono::high_resolution_clock::now();
            build_time_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

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
      }

      // Per-focal candidate counts from lattice.query()
      std::vector<double> cand_counts(n_focal, 0.0);

      if (use_lattice) {
            auto tq0 = std::chrono::high_resolution_clock::now();

            const double* focal_ptr = REAL(focal_mm);
            const int     stride_f  = n_focal;

            std::vector<index_t> cand;
            cand.reserve(128);

            for (int i = 0; i < n_focal; ++i) {
                  const double fx = focal_ptr[i];
                  const double fy = focal_ptr[i + stride_f];
                  const double* f_clim_col = focal_ptr + i + 2 * stride_f;

                  double q_geo[2] = { fx, fy };
                  std::vector<double> q_clim(n_clim);
                  for (int kdim = 0; kdim < n_clim; ++kdim) {
                        q_clim[kdim] = f_clim_col[kdim * stride_f];
                  }

                  cand.clear();
                  lattice.query(q_geo,
                                q_clim.data(),
                                max_geog,
                                use_scalar_clim,
                                max_clim_pervar_std,
                                use_scalar_clim ? max_clim_scalar
                                      : std::numeric_limits<double>::infinity(),
                                        cand);

                  cand_counts[i] = static_cast<double>(cand.size());
            }

            auto tq1 = std::chrono::high_resolution_clock::now();
            query_time_ms = std::chrono::duration<double, std::milli>(tq1 - tq0).count();
      } else {
            // No lattice: all refs are candidates
            for (int i = 0; i < n_focal; ++i) {
                  cand_counts[i] = static_cast<double>(n_ref);
            }
      }

      // Summaries of candidate counts
      double mean_candidates = 0.0;
      double p95_candidates  = 0.0;
      double max_candidates  = 0.0;

      if (n_focal > 0) {
            double sum = 0.0;
            for (int i = 0; i < n_focal; ++i) {
                  sum += cand_counts[i];
                  if (cand_counts[i] > max_candidates) {
                        max_candidates = cand_counts[i];
                  }
            }
            mean_candidates = sum / static_cast<double>(n_focal);

            std::vector<double> tmp = cand_counts;
            std::sort(tmp.begin(), tmp.end());
            std::size_t idx95 = static_cast<std::size_t>(
                  std::floor(0.95 * (static_cast<double>(n_focal) - 1.0))
            );
            p95_candidates = tmp[idx95];
      }

      return Rcpp::List::create(
            Rcpp::Named("n_focal")                     = n_focal,
            Rcpp::Named("n_ref")                       = n_ref,
            Rcpp::Named("n_clim")                      = n_clim,
            Rcpp::Named("geo_mode")                    = geo_mode,
            Rcpp::Named("use_lattice")                 = use_lattice,
            Rcpp::Named("build_time_ms")               = build_time_ms,
            Rcpp::Named("query_time_ms")               = query_time_ms,
            Rcpp::Named("total_bins")                  = total_bins,
            Rcpp::Named("n_bins_nonempty")             = n_bins_nonempty,
            Rcpp::Named("avg_bin_occupancy")           = avg_bin_occupancy,
            Rcpp::Named("avg_nonempty_bin_occupancy")  = avg_nonempty_bin_occupancy,
            Rcpp::Named("min_bin_occupancy")           = min_bin_occupancy,
            Rcpp::Named("max_bin_occupancy")           = max_bin_occupancy,
            Rcpp::Named("cand_counts")                 = cand_counts,
            Rcpp::Named("mean_candidates")             = mean_candidates,
            Rcpp::Named("p95_candidates")              = p95_candidates,
            Rcpp::Named("max_candidates")              = max_candidates
      );
}
