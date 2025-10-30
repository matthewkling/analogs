#include "lattice_index.hpp"
#include <cmath>
#include <algorithm>
#include <limits>
#include <queue>

namespace analogs {

LatticeIndex::LatticeIndex(MatrixView ref_clim_, MatrixView ref_geo_, const double* bin_w, size_tu nvars_)
      : ref_clim(ref_clim_), ref_geo(ref_geo_), nref(ref_clim_.nrow), nvars(nvars_),
        mins(nvars_, 0.0), bin_widths(nvars_, 0.0), inv_bw(nvars_, 0.0) {
      for (size_tu k = 0; k < nvars; ++k) {
            bin_widths[k] = bin_w ? bin_w[k] : INFINITY;
            inv_bw[k] = (std::isfinite(bin_widths[k]) && bin_widths[k] > 0.0) ? 1.0 / bin_widths[k] : 0.0;
      }
}

void LatticeIndex::build() {
      // per-var mins
      for (size_tu k = 0; k < nvars; ++k) {
            double m = ref_clim(0, k);
            for (size_tu i = 1; i < nref; ++i) m = std::min(m, ref_clim(i, k));
            mins[k] = m;
      }

      // bucket references by integer bin coordinate
      std::unordered_map<BinCoord, std::vector<index_t>, BinCoordHash> tmp;
      tmp.reserve(static_cast<size_tu>(nref * 1.3));
      BinCoord c; c.idx.resize(nvars);
      for (size_tu i = 0; i < nref; ++i) {
            for (size_tu k = 0; k < nvars; ++k) c.idx[k] = bin_index(ref_clim(i, k), k);
            tmp[c].push_back(static_cast<index_t>(i));
      }

      const size_tu n_bins = static_cast<size_tu>(tmp.size());
      bin_offsets.clear(); bin_offsets.reserve(n_bins + 1); bin_offsets.push_back(0);
      bin_ids_flat.clear(); bin_ids_flat.reserve(nref);
      bin_coords.clear(); bin_coords.reserve(n_bins);
      coord_to_bin.clear(); coord_to_bin.reserve(n_bins);

      size_tu bid = 0;
      for (auto &kv : tmp) {
            const BinCoord &coord = kv.first; std::vector<index_t> &ids = kv.second;
            bin_coords.push_back(coord);
            coord_to_bin[coord] = bid;
            bin_ids_flat.insert(bin_ids_flat.end(), ids.begin(), ids.end());
            bin_offsets.push_back(bin_ids_flat.size());
            ++bid;
      }
}

// Legacy window enumeration (uniform bins)
void LatticeIndex::bins_for_window(const double* focal_vars,
                                   const ClimateConstraint& clim,
                                   std::vector<size_tu>& out_bin_ids) const {
      std::vector<int32_t> lo(nvars), hi(nvars);
      for (size_tu k = 0; k < nvars; ++k) {
            if (inv_bw[k] == 0.0) { int32_t c = bin_index(focal_vars[k], k); lo[k] = c; hi[k] = c; }
            else {
                  double half = clim.window_halfwidth(k);
                  double a = focal_vars[k] - half, b = focal_vars[k] + half;
                  lo[k] = bin_index(a, k); hi[k] = bin_index(b, k);
                  if (hi[k] < lo[k]) std::swap(hi[k], lo[k]);
            }
      }
      BinCoord cur; cur.idx.resize(nvars);
      std::vector<int32_t> idx(nvars); for (size_tu k = 0; k < nvars; ++k) idx[k] = lo[k];
      while (true) {
            for (size_tu k = 0; k < nvars; ++k) cur.idx[k] = idx[k];
            auto it = coord_to_bin.find(cur); if (it != coord_to_bin.end()) out_bin_ids.push_back(it->second);
            size_tu d = 0; for (; d < nvars; ++d) { if (idx[d] < hi[d]) { idx[d]++; break; } idx[d] = lo[d]; }
            if (d == nvars) break;
      }
}

// AABB lower bound for uniform bins (per-dim width). For collapsed dims (Inf), delta=0.
double LatticeIndex::aabb_lower_bound(const double* focal_vars, const BinCoord& coord) const {
      double s = 0.0;
      for (size_tu k = 0; k < nvars; ++k) {
            if (inv_bw[k] == 0.0) continue; // collapsed dimension contributes 0
            double low = mins[k] + static_cast<double>(coord.idx[k]) * bin_widths[k];
            double high = low + bin_widths[k];
            double fk = focal_vars[k];
            double d = 0.0;
            if (fk < low) d = low - fk; else if (fk > high) d = fk - high; else d = 0.0;
            s += d * d;
      }
      return std::sqrt(s);
}

// Best-first traversal using min-heap of bins by LB(AABB)
void LatticeIndex::best_first_search(const double* focal_vars,
                                     const ClimateConstraint& clim,
                                     std::vector<index_t>& out_ids,
                                     index_t k) const {
      out_ids.clear();
      if (bin_coords.empty()) return;

      struct QNode { double lb; size_tu bid; };
      struct Cmp { bool operator()(const QNode& a, const QNode& b) const { return a.lb > b.lb; } };
      std::priority_queue<QNode, std::vector<QNode>, Cmp> pq; // min-heap by lb

      std::vector<char> visited(bin_coords.size(), 0);

      // Seed: bin containing focal if present; else seed all bins (rare)
      BinCoord start; start.idx.resize(nvars);
      for (size_tu kdim = 0; kdim < nvars; ++kdim) start.idx[kdim] = bin_index(focal_vars[kdim], kdim);
      auto it0 = coord_to_bin.find(start);
      if (it0 != coord_to_bin.end()) {
            pq.push({0.0, it0->second});
      } else {
            for (size_tu b = 0; b < bin_coords.size(); ++b) pq.push({aabb_lower_bound(focal_vars, bin_coords[b]), b});
      }

      // Neighbor pusher (±1 in index space per dim)
      auto push_neighbors = [&](size_tu bid){
            const BinCoord& base = bin_coords[bid];
            BinCoord nb; nb.idx.resize(nvars);
            for (size_tu kdim = 0; kdim < nvars; ++kdim) nb.idx[kdim] = base.idx[kdim];
            for (size_tu kdim = 0; kdim < nvars; ++kdim) {
                  nb.idx[kdim] = base.idx[kdim] - 1; auto itL = coord_to_bin.find(nb);
                  if (itL != coord_to_bin.end() && !visited[itL->second]) pq.push({aabb_lower_bound(focal_vars, nb), itL->second});
                  nb.idx[kdim] = base.idx[kdim] + 1; auto itR = coord_to_bin.find(nb);
                  if (itR != coord_to_bin.end() && !visited[itR->second]) pq.push({aabb_lower_bound(focal_vars, nb), itR->second});
                  nb.idx[kdim] = base.idx[kdim];
            }
      };

      // k-NN maintenance (max-heap on exact distances)
      struct KNode { double d; index_t id; };
      struct KCmp { bool operator()(const KNode& a, const KNode& b) const { return a.d < b.d; } }; // max-heap
      std::priority_queue<KNode, std::vector<KNode>, KCmp> topk;
      const bool do_knn = (k > 0);
      double kth = std::numeric_limits<double>::infinity();

      while (!pq.empty()) {
            QNode q = pq.top(); pq.pop();
            if (visited[q.bid]) continue; visited[q.bid] = 1;

            // Radius early-stop (k<=0 and scalar radius)
            if (!do_knn && clim.use_scalar && q.lb > clim.radius) break;

            // Scan exact within this bin
            size_tu start_off = bin_offsets[q.bid]; size_tu stop_off = bin_offsets[q.bid + 1];
            for (size_tu p = start_off; p < stop_off; ++p) {
                  index_t rid = bin_ids_flat[p];
                  if (clim.use_scalar) {
                        double s = 0.0;
                        for (size_tu kk = 0; kk < nvars; ++kk) { double d = focal_vars[kk] - ref_clim(static_cast<size_tu>(rid), kk); s += d * d; }
                        double dist = std::sqrt(s);
                        if (do_knn) {
                              if ((int)topk.size() < k) { topk.push({dist, rid}); if ((int)topk.size() == k) kth = topk.top().d; }
                              else if (dist < topk.top().d) { topk.pop(); topk.push({dist, rid}); kth = topk.top().d; }
                        } else {
                              if (dist <= clim.radius) out_ids.push_back(rid);
                        }
                  } else {
                        // Per-var band check
                        bool ok = true;
                        for (size_tu kk = 0; kk < nvars; ++kk) {
                              if (std::fabs(focal_vars[kk] - ref_clim(static_cast<size_tu>(rid), kk)) > clim.max_abs_diff[kk]) { ok = false; break; }
                        }
                        if (ok) {
                              if (do_knn) {
                                    double s = 0.0;
                                    for (size_tu kk = 0; kk < nvars; ++kk) { double d = focal_vars[kk] - ref_clim(static_cast<size_tu>(rid), kk); s += d * d; }
                                    double dist = std::sqrt(s);
                                    if ((int)topk.size() < k) { topk.push({dist, rid}); if ((int)topk.size() == k) kth = topk.top().d; }
                                    else if (dist < topk.top().d) { topk.pop(); topk.push({dist, rid}); kth = topk.top().d; }
                              } else {
                                    out_ids.push_back(rid);
                              }
                        }
                  }
            }

            // k-NN stopping: next bin LB cannot beat current kth
            if (do_knn && !pq.empty() && std::isfinite(kth)) {
                  if (pq.top().lb >= kth && (int)topk.size() >= k) break;
            }

            // Expand frontier
            push_neighbors(q.bid);
      }

      if (do_knn) {
            std::vector<index_t> tmp; tmp.reserve((size_tu)k);
            while (!topk.empty()) { tmp.push_back(topk.top().id); topk.pop(); }
            std::reverse(tmp.begin(), tmp.end());
            out_ids.swap(tmp);
      }
}

void LatticeIndex::hard_filter(const double* focal_vars,
                               const ClimateConstraint& clim,
                               std::vector<index_t>& out_ids) const {
      out_ids.clear();
      std::vector<size_tu> bins_to_visit; bins_to_visit.reserve(64);
      bins_for_window(focal_vars, clim, bins_to_visit);
      for (size_tu b = 0; b < bins_to_visit.size(); ++b) {
            size_tu bid = bins_to_visit[b]; size_tu start = bin_offsets[bid]; size_tu stop = bin_offsets[bid + 1];
            for (size_tu p = start; p < stop; ++p) {
                  index_t rid = bin_ids_flat[p];
                  bool ok;
                  if (clim.use_scalar) {
                        double s = 0.0; for (size_tu kk = 0; kk < nvars; ++kk) { double d = focal_vars[kk] - ref_clim(static_cast<size_tu>(rid), kk); s += d * d; }
                        ok = (std::sqrt(s) <= clim.radius);
                  } else {
                        ok = true;
                        for (size_tu kk = 0; kk < nvars; ++kk) {
                              if (std::fabs(focal_vars[kk] - ref_clim(static_cast<size_tu>(rid), kk)) > clim.max_abs_diff[kk]) { ok = false; break; }
                        }
                  }
                  if (ok) out_ids.push_back(rid);
            }
      }
}

} // namespace analogs
