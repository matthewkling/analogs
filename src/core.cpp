// [[Rcpp::depends(RcppParallel)]]

#include "lattice.hpp"
#include "metrics.hpp"
#include "types.hpp"

#include <Rcpp.h>
#include <queue>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

using namespace Rcpp;
using namespace analogs;

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

} // anonymous namespace


// -------------------------------------------------------------------------
// NOTE: signature updated to accept mode_code, weight_code, theta.
// mode_code:   0=knn_clim, 1=knn_geog, 2=count, 3=sum, 4=mean, 5=all
// weight_code: 0=none (unused), 1=uniform, 2=inverse_clim, 3=inverse_geog
// theta:       numeric hyperparameter for weights (epsilon), ignored if not used
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
                       double theta)
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

      // Pointers / strides (R column-major)
      const double* focal_ptr = REAL(focal_mm);
      const double* ref_ptr   = REAL(ref_mm);
      const int stride_f      = n_focal;
      const int stride_r      = n_ref;

      // Decide whether to use lattice index
      const bool use_lattice = (use_geog_filter || use_scalar_clim || use_pervar_clim);

      // Build lattice if useful
      Lattice lattice;
      double  total_bins = 1.0;
      double  avg_bin_occupancy = static_cast<double>(n_ref);
      double  min_bin_occupancy = static_cast<double>(n_ref);
      double  max_bin_occupancy = static_cast<double>(n_ref);
      std::string binning_method = "none";

      if (use_lattice) {
            MetricType metric = use_haversine ? MetricType::Haversine
            : MetricType::Planar;
            lattice.build(ref_ptr,
                          static_cast<size_tu>(n_ref),
                          static_cast<size_tu>(n_clim),
                          static_cast<size_tu>(stride_r),
                          metric,
                          max_geog,
                          use_scalar_clim,
                          max_clim_pervar_std,
                          use_scalar_clim ? max_clim_scalar : std::numeric_limits<double>::infinity());

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
            binning_method = "multi_dim_lattice";
      }

      const ModeCode   mcode = static_cast<ModeCode>(mode_code);
      const WeightCode wcode = static_cast<WeightCode>(weight_code);

      // ---------- Branch by output family: pairs vs aggregates ---------------

      const bool return_pairs =
            (mcode == ModeCode::KNN_CLIM ||
            mcode == ModeCode::KNN_GEOG ||
            mcode == ModeCode::ALL);

      if (return_pairs) {

            // Return List<IntegerVector> (indices per focal)
            List out(n_focal);

            for (int i = 0; i < n_focal; ++i) {
                  const double fx = focal_ptr[i];
                  const double fy = focal_ptr[i + stride_f];
                  const double* f_clim_col = focal_ptr + i + 2 * stride_f;

                  // Candidate set
                  std::vector<index_t> cand;
                  cand.reserve(use_lattice ? 128 : n_ref);

                  if (use_lattice) {
                        double q_geo[2] = { fx, fy };
                        std::vector<double> q_clim(n_clim);
                        for (int kdim = 0; kdim < n_clim; ++kdim) {
                              q_clim[kdim] = f_clim_col[kdim * stride_f];
                        }
                        lattice.query(q_geo,
                                      q_clim.data(),
                                      max_geog,
                                      use_scalar_clim,
                                      max_clim_pervar_std,
                                      use_scalar_clim ? max_clim_scalar : std::numeric_limits<double>::infinity(),
                                      cand);
                  } else {
                        cand.reserve(n_ref);
                        for (int j = 0; j < n_ref; ++j) cand.push_back(static_cast<index_t>(j));
                  }

                  if (mcode == ModeCode::ALL) {
                        // Collect all matches passing filters
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

                              // Climate checks (and distance)
                              const double* r_clim_col = ref_ptr + j + 2 * stride_r;
                              const auto okd = clim_ok_and_dist(f_clim_col, r_clim_col,
                                                                n_clim, stride_f, stride_r,
                                                                use_pervar_clim, max_clim_pervar_std,
                                                                use_scalar_clim, max_clim_scalar);
                              if (!okd.first) continue;

                              keep.push_back(j + 1); // 1-based
                        }

                        IntegerVector idx_vec(keep.size());
                        for (size_t u = 0; u < keep.size(); ++u) idx_vec[u] = keep[u];
                        out[i] = idx_vec;
                        continue;
                  }

                  // kNN modes (Climate or Geog)
                  const bool rank_by_clim = (mcode == ModeCode::KNN_CLIM);
                  using Neighbor = std::pair<double, int>; // (key_dist, ref_index_1based)

                  auto cmp = [](const Neighbor& a, const Neighbor& b) {
                        return a.first < b.first; // max-heap
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
                        const auto okd = clim_ok_and_dist(f_clim_col, r_clim_col,
                                                          n_clim, stride_f, stride_r,
                                                          use_pervar_clim, max_clim_pervar_std,
                                                          use_scalar_clim, max_clim_scalar);
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
                  IntegerVector idx_vec(m);
                  for (int pos = m - 1; pos >= 0; --pos) {
                        const Neighbor& nb = pq.top();
                        idx_vec[pos] = nb.second;
                        pq.pop();
                  }
                  out[i] = idx_vec;
            }

            // Attach diagnostics
            out.attr("n_focal") = n_focal;
            out.attr("n_ref")   = n_ref;
            out.attr("n_clim")  = n_clim;
            out.attr("max_dist") = max_geog;    // keep legacy attr name for now
            out.attr("max_clim") = max_clim;
            out.attr("geo_mode") = geo_mode;
            out.attr("binning_method")    = binning_method;
            out.attr("total_bins")        = total_bins;
            out.attr("avg_bin_occupancy") = avg_bin_occupancy;
            out.attr("min_bin_occupancy") = min_bin_occupancy;
            out.attr("max_bin_occupancy") = max_bin_occupancy;

            return out;
      }

      // ----------------- Aggregate modes: COUNT / SUM / MEAN -----------------

      NumericVector agg(n_focal, NA_REAL);

      for (int i = 0; i < n_focal; ++i) {
            const double fx = focal_ptr[i];
            const double fy = focal_ptr[i + stride_f];
            const double* f_clim_col = focal_ptr + i + 2 * stride_f;

            std::vector<index_t> cand;
            cand.reserve(use_lattice ? 128 : n_ref);

            if (use_lattice) {
                  double q_geo[2] = { fx, fy };
                  std::vector<double> q_clim(n_clim);
                  for (int kdim = 0; kdim < n_clim; ++kdim) {
                        q_clim[kdim] = f_clim_col[kdim * stride_f];
                  }
                  lattice.query(q_geo,
                                q_clim.data(),
                                max_geog,
                                use_scalar_clim,
                                max_clim_pervar_std,
                                use_scalar_clim ? max_clim_scalar : std::numeric_limits<double>::infinity(),
                                cand);
            } else {
                  cand.reserve(n_ref);
                  for (int j = 0; j < n_ref; ++j) cand.push_back(static_cast<index_t>(j));
            }

            // Streaming aggregation
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
                  const auto okd = clim_ok_and_dist(f_clim_col, r_clim_col,
                                                    n_clim, stride_f, stride_r,
                                                    use_pervar_clim, max_clim_pervar_std,
                                                    use_scalar_clim, max_clim_scalar);
                  if (!okd.first) continue;
                  const double clim_dist = okd.second;

                  // Update aggregates
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

      // Attach diagnostics to vector
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

      return agg;
}
