// [[Rcpp::depends(RcppParallel)]]

#include "metrics.hpp"
#include "types.hpp"
#include <Rcpp.h>
#include <queue>
#include <algorithm>
#include <cmath>

using namespace Rcpp;
using namespace analogs;

// ---- helper: geo distance -------------------------------------------------

static inline double geo_distance_km(double x1, double y1,
                                     double x2, double y2,
                                     bool use_haversine) {
    if (use_haversine) {
        double a[2] = { x1, y1 };
        double b[2] = { x2, y2 };
        Haversine2D h;
        return h.dist(a, b, 2);
    } else {
        double dx = x1 - x2;
        double dy = y1 - y2;
        return std::sqrt(dx * dx + dy * dy);
    }
}

// ---- main core function ---------------------------------------------------
//
// focal_mm, ref_mm: numeric matrices with columns [x, y, clim1, clim2, ...].
// k:   top-k to return per focal (k > 0). If k == 0, return all matches
//      that satisfy the constraints (used for "all", "sum", "mean", "count").
// max_clim: scalar (Euclidean radius) or vector (per-variable abs diff).
// max_dist: max geographic distance in km (Inf = no constraint).
// geo_mode: "lonlat" (Haversine) or "projected" (Euclidean).
//
// Return: list of length n_focal; element i is an integer vector of
// 1-based indices into ref_mm of matching analogs for focal i.
//
 // [[Rcpp::export]]
List find_analogs_core(const NumericMatrix& focal_mm,
                       const NumericMatrix& ref_mm,
                       int k,
                       const NumericVector& max_clim,
                       double max_dist,
                       const std::string& geo_mode) {

    const int n_focal = focal_mm.nrow();
    const int n_ref   = ref_mm.nrow();

    if (n_focal == 0 || n_ref == 0) {
        // Return empty list with diagnostics
        List out(n_focal);
        out.attr("n_focal") = n_focal;
        out.attr("n_ref")   = n_ref;
        out.attr("n_clim")  = 0;
        out.attr("binning_method") = "none";
        out.attr("total_bins") = 1;
        out.attr("avg_bin_occupancy") = n_ref;
        out.attr("min_bin_occupancy") = n_ref;
        out.attr("max_bin_occupancy") = n_ref;
        return out;
    }

    const int ncol_focal = focal_mm.ncol();
    const int ncol_ref   = ref_mm.ncol();

    if (ncol_focal != ncol_ref) {
        stop("focal and ref must have the same number of columns");
    }
    if (ncol_focal < 3) {
        stop("Need at least 2 coordinate columns and 1 climate variable");
    }

    const int n_clim = ncol_focal - 2;

    // Interpret geographic mode
    const bool use_haversine = (geo_mode == "lonlat");

    // Interpret max_dist (Inf => no filter)
    const bool use_dist_filter = std::isfinite(max_dist);

    // Interpret max_clim:
    //
    //   * length 1: scalar radius in Euclidean climate space
    //   * length p: per-variable abs difference thresholds
    //
    bool use_scalar_clim = false;
    bool use_pervar_clim = false;
    double max_clim_scalar = std::numeric_limits<double>::infinity();
    NumericVector max_clim_pervar;

    if (max_clim.size() == 1) {
        double v = max_clim[0];
        if (std::isfinite(v)) {
            use_scalar_clim = true;
            max_clim_scalar = v;
        }
    } else if (max_clim.size() == n_clim) {
        max_clim_pervar = clone(max_clim);
        use_pervar_clim = true;
    } else if (max_clim.size() > 1) {
        stop("max_clim must be length 1 or length equal to the number of climate variables");
    }

    // Climate distance metric (Euclidean in whitened space)
    Euclidean clim_metric;

    // Output list of index vectors
    List out(n_focal);

    // Data pointers for a slightly cheaper access pattern
    const double* focal_ptr = REAL(focal_mm);
    const double* ref_ptr   = REAL(ref_mm);
    const int stride_f      = n_focal; // column-major: nrow
    const int stride_r      = n_ref;

    // For each focal site
    for (int i = 0; i < n_focal; ++i) {
        // Coordinates of focal i
        const double fx = focal_ptr[i];             // col 0
        const double fy = focal_ptr[i + stride_f];  // col 1

        // Climate pointer for focal i (start of clim columns)
        const double* f_clim = focal_ptr + i + 2 * stride_f;

        // Storage for matches
        if (k > 0) {
            // Top-k by climate distance within constraints

            using Neighbor = std::pair<double, int>; // (clim_dist, ref_index_1based)

            auto cmp = [](const Neighbor& a, const Neighbor& b) {
                // max-heap by distance: largest distance at top
                return a.first < b.first;
            };
            std::priority_queue<Neighbor,
                                std::vector<Neighbor>,
                                decltype(cmp)> pq(cmp);

            for (int j = 0; j < n_ref; ++j) {
                // Geographic filter
                if (use_dist_filter) {
                    const double rx = ref_ptr[j];             // col 0
                    const double ry = ref_ptr[j + stride_r];  // col 1;
                    double gdist = geo_distance_km(fx, fy, rx, ry, use_haversine);
                    if (gdist > max_dist) continue;
                }

                // Climate pointer for ref j
                const double* r_clim = ref_ptr + j + 2 * stride_r;

                // Climate filters and distance
                double clim_dist = 0.0;

                if (use_pervar_clim) {
                    // Per-variable threshold AND Euclidean distance
                    double sumsq = 0.0;
                    bool ok = true;
                    for (int kdim = 0; kdim < n_clim; ++kdim) {
                        double df = f_clim[kdim * stride_f] - r_clim[kdim * stride_r];
                        if (std::fabs(df) > max_clim_pervar[kdim]) {
                            ok = false;
                            break;
                        }
                        sumsq += df * df;
                    }
                    if (!ok) continue;
                    clim_dist = std::sqrt(sumsq);
                } else {
                    // Only scalar threshold (if finite) and/or ordering
                    // Build temporary contiguous arrays for metric convenience
                    std::vector<double> tmp_f(n_clim);
                    std::vector<double> tmp_r(n_clim);
                    for (int kdim = 0; kdim < n_clim; ++kdim) {
                        tmp_f[kdim] = f_clim[kdim * stride_f];
                        tmp_r[kdim] = r_clim[kdim * stride_r];
                    }
                    clim_dist = clim_metric.dist(tmp_f.data(), tmp_r.data(),
                                                 static_cast<size_tu>(n_clim));
                    if (use_scalar_clim && clim_dist > max_clim_scalar) {
                        continue;
                    }
                }

                // Candidate survives all filters; add to top-k structure
                int ref_index_1based = j + 1;

                if (static_cast<int>(pq.size()) < k) {
                    pq.emplace(clim_dist, ref_index_1based);
                } else if (!pq.empty() && clim_dist < pq.top().first) {
                    pq.pop();
                    pq.emplace(clim_dist, ref_index_1based);
                }
            }

            // Extract and sort by increasing climate distance
            const int m = static_cast<int>(pq.size());
            IntegerVector idx(m);
            std::vector<double> dists(m);

            for (int pos = m - 1; pos >= 0; --pos) {
                const Neighbor& nb = pq.top();
                idx[pos]   = nb.second;
                dists[pos] = nb.first;
                pq.pop();
            }

            // (We don't return distances; R will recompute as needed.)
            out[i] = idx;
        } else {
            // k == 0: return all matches satisfying filters, no ordering
            std::vector<int> matches;
            matches.reserve(64); // small default; will grow if needed

            for (int j = 0; j < n_ref; ++j) {
                // Geographic filter
                if (use_dist_filter) {
                    const double rx = ref_ptr[j];
                    const double ry = ref_ptr[j + stride_r];
                    double gdist = geo_distance_km(fx, fy, rx, ry, use_haversine);
                    if (gdist > max_dist) continue;
                }

                // Climate pointer for ref j
                const double* r_clim = ref_ptr + j + 2 * stride_r;

                // Climate filters (no need for distance unless scalar threshold)
                if (use_pervar_clim) {
                    bool ok = true;
                    for (int kdim = 0; kdim < n_clim; ++kdim) {
                        double df = f_clim[kdim * stride_f] - r_clim[kdim * stride_r];
                        if (std::fabs(df) > max_clim_pervar[kdim]) {
                            ok = false;
                            break;
                        }
                    }
                    if (!ok) continue;
                } else if (use_scalar_clim) {
                    std::vector<double> tmp_f(n_clim);
                    std::vector<double> tmp_r(n_clim);
                    for (int kdim = 0; kdim < n_clim; ++kdim) {
                        tmp_f[kdim] = f_clim[kdim * stride_f];
                        tmp_r[kdim] = r_clim[kdim * stride_r];
                    }
                    double clim_dist = clim_metric.dist(tmp_f.data(), tmp_r.data(),
                                                        static_cast<size_tu>(n_clim));
                    if (clim_dist > max_clim_scalar) continue;
                }

                matches.push_back(j + 1); // 1-based for R
            }

            IntegerVector idx(matches.size());
            for (std::size_t m = 0; m < matches.size(); ++m) {
                idx[m] = matches[m];
            }
            out[i] = idx;
        }
    }

    // Attach basic diagnostics (index-related fields are placeholders for now)
    out.attr("n_focal") = n_focal;
    out.attr("n_ref")   = n_ref;
    out.attr("n_clim")  = n_clim;
    out.attr("max_dist") = max_dist;
    out.attr("max_clim") = max_clim;
    out.attr("geo_mode") = geo_mode;

    out.attr("binning_method")     = "none";   // lattice not yet used here
    out.attr("total_bins")         = 1;
    out.attr("avg_bin_occupancy")  = n_ref;
    out.attr("min_bin_occupancy")  = n_ref;
    out.attr("max_bin_occupancy")  = n_ref;

    return out;
}
