#include <RcppArmadillo.h>
#include <algorithm>
#include <numeric>
#include <set>
#include <unordered_set>
#include <memory>
#include <cmath>

#include "analog_utils.cpp"   // WeightType, parse_weight_type, calc_weight, find_top_n_weights
#include "lattice_index.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// ---------------------------------------------------------------------
// Aggregation type (matrix-only pipeline)
// ---------------------------------------------------------------------

enum AggType { NMAX, SUM, MEAN, COUNT, ALL };

inline AggType parse_agg_type(const std::string& fun_str) {
      if (fun_str == "nmax")  return NMAX;
      if (fun_str == "sum")   return SUM;
      if (fun_str == "mean")  return MEAN;
      if (fun_str == "count") return COUNT;
      if (fun_str == "all")   return ALL;
      return NMAX; // default
}

// ---------------------------------------------------------------------
// Climate lattice index (1–4 dims)
// ---------------------------------------------------------------------
struct ClimateIndexBank {
      bool available{false};
      int dims{0};
      std::unique_ptr<LatticeIndex<1>> idx1;
      std::unique_ptr<LatticeIndex<2>> idx2;
      std::unique_ptr<LatticeIndex<3>> idx3;
      std::unique_ptr<LatticeIndex<4>> idx4;

      void build_if_needed(const arma::mat& ref_climate, bool want_idx, int bins_per_dim = 15) {
            if (!want_idx) { available = false; dims = 0; return; }
            dims = (int)ref_climate.n_cols;
            available = (dims >= 1 && dims <= 4);
            if (!available) return;
            if (dims == 1) { idx1.reset(new LatticeIndex<1>()); idx1->build(ref_climate, bins_per_dim); }
            else if (dims == 2) { idx2.reset(new LatticeIndex<2>()); idx2->build(ref_climate, bins_per_dim); }
            else if (dims == 3) { idx3.reset(new LatticeIndex<3>()); idx3->build(ref_climate, bins_per_dim); }
            else if (dims == 4) { idx4.reset(new LatticeIndex<4>()); idx4->build(ref_climate, bins_per_dim); }
      }

      std::vector<int> query_box(const arma::rowvec& focal_clim, double threshold) const {
            if (!available) return {};
            if (dims == 1) return idx1->query_box(focal_clim, threshold);
            if (dims == 2) return idx2->query_box(focal_clim, threshold);
            if (dims == 3) return idx3->query_box(focal_clim, threshold);
            if (dims == 4) return idx4->query_box(focal_clim, threshold);
            return {};
      }
};

// ---------------------------------------------------------------------
// Small helpers
// ---------------------------------------------------------------------
static inline void sort_unique(std::vector<int>& v) {
      std::sort(v.begin(), v.end());
      v.erase(std::unique(v.begin(), v.end()), v.end());
}

static inline std::vector<int> intersect_sorted(const std::vector<int>& a, const std::vector<int>& b) {
      std::vector<int> out; out.reserve(std::min(a.size(), b.size()));
      size_t i = 0, j = 0;
      while (i < a.size() && j < b.size()) {
            if (a[i] < b[j]) ++i;
            else if (b[j] < a[i]) ++j;
            else { out.push_back(a[i]); ++i; ++j; }
      }
      return out;
}

// ---------------------------------------------------------------------
// Vectorized geo distance for a subset (builds from ref_coords rows(idx))
// ---------------------------------------------------------------------
static inline arma::vec geo_dist_subset(const arma::rowvec& fc,
                                        const arma::mat& ref_coords_rows) {
      // Inputs: fc = [x,y] in degrees; ref_coords_rows = kx2 in degrees
      const double lat1 = fc[1] * arma::datum::pi / 180.0;
      const double lon1 = fc[0] * arma::datum::pi / 180.0;
      const double cos_lat1 = std::cos(lat1);

      arma::vec lat2 = (arma::datum::pi / 180.0) * ref_coords_rows.col(1);
      arma::vec lon2 = (arma::datum::pi / 180.0) * ref_coords_rows.col(0);
      arma::vec dlat = lat2 - lat1;
      arma::vec dlon = lon2 - lon1;

      arma::vec sin_dlat2 = arma::sin(0.5 * dlat);
      arma::vec sin_dlon2 = arma::sin(0.5 * dlon);
      arma::vec a = arma::square(sin_dlat2) + cos_lat1 * arma::cos(lat2) % arma::square(sin_dlon2);
      a = arma::clamp(a, 0.0, 1.0);
      return 2.0 * 6371.0 * arma::asin(arma::sqrt(a));
}

// ---------------------------------------------------------------------
// Matrix–Matrix entry point (only exported routine now)
// ---------------------------------------------------------------------
// [[Rcpp::export]]
DataFrame find_analogs_matrix_matrix(
            arma::mat focal_coords,
            arma::mat focal_climate,
            arma::mat ref_coords,
            arma::mat ref_climate,
            Nullable<double> max_dist,
            Nullable<double> max_clim,
            std::string weight,
            std::string fun,
            Nullable<int> n
) {
      const int n_focal = focal_coords.n_rows;
      const int n_ref   = ref_coords.n_rows;
      const int n_vars  = focal_climate.n_cols;

      WeightType weight_type = parse_weight_type(weight);
      AggType    agg_type    = parse_agg_type(fun);
      int n_analogs = (agg_type == NMAX && !n.isNull()) ? as<int>(n) : 1;

      const double radius_km      = max_dist.isNotNull() ? as<double>(max_dist) : R_PosInf;
      const double clim_threshold = max_clim.isNotNull() ? as<double>(max_clim) : R_PosInf;

      // Build climate index once if useful
      const bool want_clim = (n_vars >= 1 && n_vars <= 4 && n_ref > 5000 && max_clim.isNotNull());
      ClimateIndexBank clim_idx; clim_idx.build_if_needed(ref_climate, want_clim, 15);

      // Build geo index once up-front if it's plausibly needed globally
      const bool need_geo_idx = (max_dist.isNotNull() || (agg_type == NMAX && weight_type == INVERSE_DIST)) && (n_ref > 5000);
      std::unique_ptr<LatticeIndex<2>> geo_idx;
      bool have_geo_idx = false;
      if (need_geo_idx) { geo_idx.reset(new LatticeIndex<2>()); geo_idx->build(ref_coords, 20); have_geo_idx = true; }

      // Expansion heuristics (adaptive to n)
      const int base_cutoff       = std::max(400, std::min(2000, 40 * n_analogs));
      const int target_candidates = std::max(500, 40 * n_analogs);

      // Outputs (preallocate some)
      std::vector<int> result_focal_idx, result_analog_idx; result_focal_idx.reserve(n_focal * n_analogs);
      std::vector<double> result_focal_x, result_focal_y, result_analog_x, result_analog_y;
      std::vector<double> result_clim_dist, result_geo_dist, result_weight, result_value;

      // Reusable buffers
      std::vector<int> pool; pool.reserve(n_ref);
      std::vector<int> vi;   vi.reserve(1024);
      std::vector<double> gdv, cdv, wv; gdv.reserve(1024); cdv.reserve(1024); wv.reserve(1024);

      for (int f = 0; f < n_focal; ++f) {
            const arma::rowvec fc = focal_coords.row(f);
            const arma::rowvec fv = focal_climate.row(f);
            if (fv.has_nan()) continue;

            pool.clear(); vi.clear(); gdv.clear(); cdv.clear(); wv.clear();

            // Hard climate filter first
            bool used_clim_filter = false;
            if (clim_threshold < R_PosInf && clim_idx.available) {
                  auto cc = clim_idx.query_box(fv, clim_threshold);
                  pool.assign(cc.begin(), cc.end());
                  used_clim_filter = true;
            } else {
                  pool.resize(n_ref); std::iota(pool.begin(), pool.end(), 0);
            }

            // Build geo index lazily only if we actually need it
            if (!have_geo_idx && (max_dist.isNotNull() || (agg_type == NMAX && weight_type == INVERSE_DIST && (int)pool.size() > base_cutoff))) {
                  if (n_ref > 5000) { geo_idx.reset(new LatticeIndex<2>()); geo_idx->build(ref_coords, 20); have_geo_idx = true; }
            }

            // If a hard geo filter exists and we have an index, intersect; otherwise rely on precise check later
            if (radius_km < R_PosInf && have_geo_idx) {
                  auto gc = geo_idx->query_box(fc, radius_km);
                  if (used_clim_filter) {
                        sort_unique(pool); sort_unique(gc); pool = intersect_sorted(pool, gc);
                  } else {
                        pool.swap(gc);
                  }
            }

            // Optional expansion when pool still large and we only need top-n
            if (agg_type == NMAX && (int)pool.size() > base_cutoff) {
                  double ring   = (weight_type == INVERSE_DIST) ? 50.0 : 0.05; // km or clim units
                  const double grow   = 1.7;
                  std::unordered_set<int> seen(pool.begin(), pool.end());
                  while ((int)pool.size() < target_candidates) {
                        if (weight_type == INVERSE_DIST && have_geo_idx) {
                              auto g2 = geo_idx->query_box(fc, ring);
                              for (int i : g2) if (!seen.count(i)) { seen.insert(i); pool.push_back(i); }
                        } else if (weight_type == INVERSE_CLIM && clim_idx.available) {
                              auto c2 = clim_idx.query_box(fv, ring);
                              for (int i : c2) if (!seen.count(i)) { seen.insert(i); pool.push_back(i); }
                        } else break;
                        ring *= grow;
                        if (ring > 1e6) break;
                  }
            }

            if (pool.empty()) {
                  if (agg_type == COUNT || agg_type == SUM || agg_type == MEAN) {
                        result_focal_idx.push_back(f + 1);
                        result_focal_x.push_back(fc[0]);
                        result_focal_y.push_back(fc[1]);
                        result_value.push_back( (agg_type == COUNT) ? 0.0 : NA_REAL );
                  }
                  continue;
            }

            // Slice ref rows for this pool once
            arma::uvec idx(pool.size());
            for (size_t i = 0; i < pool.size(); ++i) idx(i) = static_cast<unsigned int>(pool[i]);
            arma::mat Rv = ref_climate.rows(idx);
            arma::mat Rc = ref_coords.rows(idx);

            // Vectorized climate distances
            arma::vec cdists = arma::sqrt(arma::sum(arma::square(Rv.each_row() - fv), 1));

            // Vectorized geographic distances (computed only for pool)
            arma::vec gdist = geo_dist_subset(fc, Rc);

            // Precise thresholds and accumulation
            for (size_t i = 0; i < pool.size(); ++i) {
                  if (!std::isfinite(cdists[i]) || !std::isfinite(gdist[i])) continue;
                  if (cdists[i] > clim_threshold || gdist[i] > radius_km) continue;
                  const double w = calc_weight(cdists[i], gdist[i], weight_type);
                  vi.push_back(pool[i]);
                  gdv.push_back(gdist[i]);
                  cdv.push_back(cdists[i]);
                  wv.push_back(w);
            }

            // Aggregate
            if (agg_type == NMAX) {
                  arma::vec wv_arma(wv);
                  arma::uvec top = find_top_n_weights(wv_arma, n_analogs);
                  for (auto j : top) {
                        int ri = vi[j];
                        result_focal_idx.push_back(f + 1);
                        result_focal_x.push_back(fc[0]);
                        result_focal_y.push_back(fc[1]);
                        result_analog_idx.push_back(ri + 1);
                        result_analog_x.push_back(ref_coords(ri, 0));
                        result_analog_y.push_back(ref_coords(ri, 1));
                        result_clim_dist.push_back(cdv[j]);
                        result_geo_dist.push_back(gdv[j]);
                        result_weight.push_back(wv[j]);
                  }
            } else if (agg_type == COUNT) {
                  result_focal_idx.push_back(f + 1);
                  result_focal_x.push_back(fc[0]);
                  result_focal_y.push_back(fc[1]);
                  result_value.push_back((double)vi.size());
            } else if (agg_type == SUM || agg_type == MEAN) {
                  double s = std::accumulate(wv.begin(), wv.end(), 0.0);
                  result_focal_idx.push_back(f + 1);
                  result_focal_x.push_back(fc[0]);
                  result_focal_y.push_back(fc[1]);
                  result_value.push_back((agg_type == MEAN && !vi.empty()) ? (s / (double)vi.size()) : s);
            } else if (agg_type == ALL) {
                  for (size_t i = 0; i < vi.size(); ++i) {
                        int ri = vi[i];
                        result_focal_idx.push_back(f + 1);
                        result_focal_x.push_back(fc[0]);
                        result_focal_y.push_back(fc[1]);
                        result_analog_idx.push_back(ri + 1);
                        result_analog_x.push_back(ref_coords(ri, 0));
                        result_analog_y.push_back(ref_coords(ri, 1));
                        result_clim_dist.push_back(cdv[i]);
                        result_geo_dist.push_back(gdv[i]);
                        result_weight.push_back(wv[i]);
                  }
            }
      }

      if (agg_type == COUNT || agg_type == SUM || agg_type == MEAN) {
            return DataFrame::create(
                  _["focal_index"] = result_focal_idx,
                  _["focal_x"]     = result_focal_x,
                  _["focal_y"]     = result_focal_y,
                  _["value"]       = result_value
            );
      } else {
            return DataFrame::create(
                  _["focal_index"] = result_focal_idx,
                  _["focal_x"]     = result_focal_x,
                  _["focal_y"]     = result_focal_y,
                  _["analog_index"] = result_analog_idx,
                  _["analog_x"]     = result_analog_x,
                  _["analog_y"]     = result_analog_y,
                  _["clim_dist"]    = result_clim_dist,
                  _["geog_dist"]    = result_geo_dist,
                  _["weight"]       = result_weight
            );
      }
}
