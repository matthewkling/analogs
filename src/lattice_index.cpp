#include "lattice_index.hpp"
#include <cmath>
#include <algorithm>
#include <limits>
#include <queue>

namespace analogs {

LatticeIndex::LatticeIndex(MatrixView ref_clim_, MatrixView ref_geo_, size_tu nvars_,
                           size_tu target_occupancy_)
      : ref_clim(ref_clim_), ref_geo(ref_geo_),
        nref(ref_clim_.nrow), nvars(nvars_),
        target_occupancy(target_occupancy_) {

      quantile_edges.resize(nvars);
      n_bins_per_dim.resize(nvars);
}

size_tu LatticeIndex::compute_n_bins() const {
      // Adaptive bin sizing: n_bins = max(10, floor(pow(n/target, 1/d)))
      double ratio = static_cast<double>(nref) / static_cast<double>(target_occupancy);
      double exponent = 1.0 / static_cast<double>(nvars);
      size_tu n_bins = static_cast<size_tu>(std::floor(std::pow(ratio, exponent)));

      // Bounds
      if (n_bins < 10) n_bins = 10;
      if (n_bins > 1000) n_bins = 1000;

      return n_bins;
}

void LatticeIndex::compute_quantile_edges(size_tu dim_idx, size_tu n_bins) {
      // Extract values for this dimension
      std::vector<double> values;
      values.reserve(nref);
      for (size_tu i = 0; i < nref; ++i) {
            values.push_back(ref_clim(i, dim_idx));
      }

      // Sort to compute quantiles
      std::sort(values.begin(), values.end());

      // Compute bin edges at quantiles
      quantile_edges[dim_idx].resize(n_bins + 1);
      quantile_edges[dim_idx][0] = values.front();
      quantile_edges[dim_idx][n_bins] = values.back();

      // Interior edges at quantiles
      for (size_tu b = 1; b < n_bins; ++b) {
            double quantile = static_cast<double>(b) / static_cast<double>(n_bins);
            double pos = quantile * static_cast<double>(nref - 1);
            size_tu idx = static_cast<size_tu>(pos);

            if (idx + 1 < nref) {
                  double frac = pos - static_cast<double>(idx);
                  quantile_edges[dim_idx][b] = values[idx] * (1.0 - frac) +
                        values[idx + 1] * frac;
            } else {
                  quantile_edges[dim_idx][b] = values[idx];
            }
      }

      // Ensure strict monotonicity
      for (size_tu b = 1; b <= n_bins; ++b) {
            if (quantile_edges[dim_idx][b] <= quantile_edges[dim_idx][b-1]) {
                  quantile_edges[dim_idx][b] = quantile_edges[dim_idx][b-1] +
                        std::numeric_limits<double>::epsilon();
            }
      }
}

// Helper: Encode bin coordinate as single integer (mixed-radix)
// This avoids heap-backed vectors and enables efficient sorting
inline uint64_t encode_bin_coord(const std::vector<int32_t>& indices,
                                 const std::vector<size_tu>& n_bins_per_dim) {
      uint64_t key = 0;
      uint64_t multiplier = 1;

      for (size_tu k = 0; k < indices.size(); ++k) {
            // Offset by a large constant to handle negative bin indices
            // (shouldn't happen in practice, but safe)
            uint64_t bin_val = static_cast<uint64_t>(indices[k] + 100000);
            key += bin_val * multiplier;
            multiplier *= (n_bins_per_dim[k] + 200000);  // space for offset range
      }

      return key;
}

// Helper: Decode integer key back to bin coordinate (for debugging/lookup)
inline void decode_bin_coord(uint64_t key,
                             const std::vector<size_tu>& n_bins_per_dim,
                             std::vector<int32_t>& indices) {
      indices.resize(n_bins_per_dim.size());
      uint64_t remaining = key;

      for (size_tu k = 0; k < n_bins_per_dim.size(); ++k) {
            uint64_t divisor = n_bins_per_dim[k] + 200000;
            uint64_t bin_val = remaining % divisor;
            indices[k] = static_cast<int32_t>(bin_val) - 100000;
            remaining /= divisor;
      }
}

void LatticeIndex::build() {
      // Compute number of bins
      size_tu n_bins = compute_n_bins();
      for (size_tu k = 0; k < nvars; ++k) {
            n_bins_per_dim[k] = n_bins;
            compute_quantile_edges(k, n_bins);
      }

      // OPTIMIZATION: Use flat vector + sort instead of hash map
      // This eliminates heap allocations and gives better cache locality
      std::vector<std::pair<uint64_t, index_t>> bin_assignments;
      bin_assignments.reserve(nref);

      std::vector<int32_t> temp_indices(nvars);

      // Assign each reference point to a bin
      for (size_tu i = 0; i < nref; ++i) {
            for (size_tu k = 0; k < nvars; ++k) {
                  temp_indices[k] = bin_index(ref_clim(i, k), k);
            }
            uint64_t key = encode_bin_coord(temp_indices, n_bins_per_dim);
            bin_assignments.push_back({key, static_cast<index_t>(i)});
      }

      // Sort by bin key - O(n log n) with excellent cache locality!
      std::sort(bin_assignments.begin(), bin_assignments.end());

      // Build CSR structure by linear scan of sorted data
      bin_ids_flat.clear();
      bin_ids_flat.reserve(nref);

      bin_offsets.clear();
      bin_offsets.reserve(nref / 20);  // estimate
      bin_offsets.push_back(0);

      bin_coords.clear();
      bin_coords.reserve(nref / 20);

      coord_to_bin.clear();
      coord_to_bin.reserve(nref / 20);

      if (bin_assignments.empty()) return;

      // Process first bin
      uint64_t current_key = bin_assignments[0].first;
      BinCoord current_coord;
      decode_bin_coord(current_key, n_bins_per_dim, current_coord.idx);
      bin_coords.push_back(current_coord);
      coord_to_bin[current_coord] = 0;

      size_tu current_bin_idx = 0;

      // Scan through sorted assignments
      for (size_tu i = 0; i < bin_assignments.size(); ++i) {
            uint64_t key = bin_assignments[i].first;
            index_t ref_id = bin_assignments[i].second;

            // New bin?
            if (key != current_key) {
                  // Finalize previous bin
                  bin_offsets.push_back(bin_ids_flat.size());

                  // Start new bin
                  current_key = key;
                  current_bin_idx++;

                  decode_bin_coord(key, n_bins_per_dim, current_coord.idx);
                  bin_coords.push_back(current_coord);
                  coord_to_bin[current_coord] = current_bin_idx;
            }

            // Add ref to current bin
            bin_ids_flat.push_back(ref_id);
      }

      // Finalize last bin
      bin_offsets.push_back(bin_ids_flat.size());
}

void LatticeIndex::get_diagnostics(size_tu& total_bins, double& avg_occupancy,
                                   double& min_occupancy, double& max_occupancy) const {
      total_bins = bin_coords.size();

      if (total_bins == 0) {
            avg_occupancy = 0.0;
            min_occupancy = 0.0;
            max_occupancy = 0.0;
            return;
      }

      min_occupancy = std::numeric_limits<double>::max();
      max_occupancy = 0.0;
      double sum = 0.0;

      for (size_tu b = 0; b < total_bins; ++b) {
            double occ = static_cast<double>(bin_offsets[b + 1] - bin_offsets[b]);
            sum += occ;
            if (occ < min_occupancy) min_occupancy = occ;
            if (occ > max_occupancy) max_occupancy = occ;
      }

      avg_occupancy = sum / static_cast<double>(total_bins);
}

void LatticeIndex::bins_for_window(const double* focal_vars,
                                   const ClimateConstraint& clim,
                                   std::vector<size_tu>& out_bin_ids) const {
      std::vector<int32_t> lo(nvars), hi(nvars);

      for (size_tu k = 0; k < nvars; ++k) {
            double half = clim.window_halfwidth(k);
            double min_val = focal_vars[k] - half;
            double max_val = focal_vars[k] + half;

            lo[k] = bin_index(min_val, k);
            hi[k] = bin_index(max_val, k);

            if (hi[k] < lo[k]) std::swap(hi[k], lo[k]);
      }

      // Enumerate all bins in the rectangular window
      BinCoord cur;
      cur.idx.resize(nvars);
      std::vector<int32_t> idx(nvars);
      for (size_tu k = 0; k < nvars; ++k) idx[k] = lo[k];

      while (true) {
            for (size_tu k = 0; k < nvars; ++k) cur.idx[k] = idx[k];

            auto it = coord_to_bin.find(cur);
            if (it != coord_to_bin.end()) {
                  out_bin_ids.push_back(it->second);
            }

            size_tu d = 0;
            for (; d < nvars; ++d) {
                  if (idx[d] < hi[d]) {
                        idx[d]++;
                        break;
                  }
                  idx[d] = lo[d];
            }
            if (d == nvars) break;
      }
}

// Best-first search using SQUARED distances (no sqrt!)
void LatticeIndex::best_first_search(const double* focal_vars,
                                     const ClimateConstraint& clim,
                                     std::vector<index_t>& out_ids,
                                     index_t k) const {
      out_ids.clear();
      if (bin_coords.empty()) return;

      // Priority queue nodes store SQUARED lower bounds
      struct QNode { double lb_sq; size_tu bid; };
      struct Cmp { bool operator()(const QNode& a, const QNode& b) const {
            return a.lb_sq > b.lb_sq;
      }};
      std::priority_queue<QNode, std::vector<QNode>, Cmp> pq;

      std::vector<char> visited(bin_coords.size(), 0);

      // Seed: bin containing focal
      BinCoord start;
      start.idx.resize(nvars);
      for (size_tu kdim = 0; kdim < nvars; ++kdim) {
            start.idx[kdim] = bin_index(focal_vars[kdim], kdim);
      }

      auto it0 = coord_to_bin.find(start);
      if (it0 != coord_to_bin.end()) {
            pq.push({0.0, it0->second});
      } else {
            // Seed all bins (focal not in any bin)
            for (size_tu b = 0; b < bin_coords.size(); ++b) {
                  pq.push({aabb_lower_bound_squared(focal_vars, bin_coords[b]), b});
            }
      }

      // Neighbor pusher
      auto push_neighbors = [&](size_tu bid){
            const BinCoord& base = bin_coords[bid];
            BinCoord nb;
            nb.idx.resize(nvars);
            for (size_tu kdim = 0; kdim < nvars; ++kdim) nb.idx[kdim] = base.idx[kdim];

            for (size_tu kdim = 0; kdim < nvars; ++kdim) {
                  nb.idx[kdim] = base.idx[kdim] - 1;
                  auto itL = coord_to_bin.find(nb);
                  if (itL != coord_to_bin.end() && !visited[itL->second]) {
                        pq.push({aabb_lower_bound_squared(focal_vars, nb), itL->second});
                  }

                  nb.idx[kdim] = base.idx[kdim] + 1;
                  auto itR = coord_to_bin.find(nb);
                  if (itR != coord_to_bin.end() && !visited[itR->second]) {
                        pq.push({aabb_lower_bound_squared(focal_vars, nb), itR->second});
                  }

                  nb.idx[kdim] = base.idx[kdim];
            }
      };

      // k-NN heap stores SQUARED distances
      struct KNode { double d_sq; index_t id; };
      struct KCmp { bool operator()(const KNode& a, const KNode& b) const {
            return a.d_sq < b.d_sq;
      }};
      std::priority_queue<KNode, std::vector<KNode>, KCmp> topk;

      const bool do_knn = (k > 0);
      double kth_sq = std::numeric_limits<double>::infinity();

      // Pre-compute SQUARED radius for comparisons
      const double radius_sq = clim.radius * clim.radius;

      while (!pq.empty()) {
            QNode q = pq.top();
            pq.pop();

            if (visited[q.bid]) continue;
            visited[q.bid] = 1;

            // Early-stop for radius queries (compare SQUARED distances)
            if (!do_knn && clim.use_scalar && q.lb_sq > radius_sq) break;

            // Scan bin
            size_tu start_off = bin_offsets[q.bid];
            size_tu stop_off = bin_offsets[q.bid + 1];

            for (size_tu p = start_off; p < stop_off; ++p) {
                  index_t rid = bin_ids_flat[p];

                  if (clim.use_scalar) {
                        // Compute SQUARED distance (no sqrt!)
                        double dist_sq = 0.0;
                        for (size_tu kk = 0; kk < nvars; ++kk) {
                              double d = focal_vars[kk] - ref_clim(static_cast<size_tu>(rid), kk);
                              dist_sq += d * d;
                        }

                        if (do_knn) {
                              // k-NN: compare squared distances
                              if ((int)topk.size() < k) {
                                    topk.push({dist_sq, rid});
                                    if ((int)topk.size() == k) kth_sq = topk.top().d_sq;
                              } else if (dist_sq < topk.top().d_sq) {
                                    topk.pop();
                                    topk.push({dist_sq, rid});
                                    kth_sq = topk.top().d_sq;
                              }
                        } else {
                              // Radius filter: compare squared distances
                              if (dist_sq <= radius_sq) {
                                    out_ids.push_back(rid);
                              }
                        }
                  } else {
                        // Per-var band check
                        bool ok = true;
                        for (size_tu kk = 0; kk < nvars; ++kk) {
                              if (std::fabs(focal_vars[kk] - ref_clim(static_cast<size_tu>(rid), kk)) >
                                        clim.max_abs_diff[kk]) {
                                    ok = false;
                                    break;
                              }
                        }

                        if (ok) {
                              if (do_knn) {
                                    double dist_sq = 0.0;
                                    for (size_tu kk = 0; kk < nvars; ++kk) {
                                          double d = focal_vars[kk] - ref_clim(static_cast<size_tu>(rid), kk);
                                          dist_sq += d * d;
                                    }

                                    if ((int)topk.size() < k) {
                                          topk.push({dist_sq, rid});
                                          if ((int)topk.size() == k) kth_sq = topk.top().d_sq;
                                    } else if (dist_sq < topk.top().d_sq) {
                                          topk.pop();
                                          topk.push({dist_sq, rid});
                                          kth_sq = topk.top().d_sq;
                                    }
                              } else {
                                    out_ids.push_back(rid);
                              }
                        }
                  }
            }

            // k-NN early stopping
            if (do_knn && !pq.empty() && std::isfinite(kth_sq)) {
                  if (pq.top().lb_sq >= kth_sq && (int)topk.size() >= k) break;
            }

            push_neighbors(q.bid);
      }

      // Extract results
      if (do_knn) {
            std::vector<index_t> tmp;
            tmp.reserve((size_tu)k);
            while (!topk.empty()) {
                  tmp.push_back(topk.top().id);
                  topk.pop();
            }
            std::reverse(tmp.begin(), tmp.end());
            out_ids.swap(tmp);
      }
}

void LatticeIndex::hard_filter(const double* focal_vars,
                               const ClimateConstraint& clim,
                               std::vector<index_t>& out_ids) const {
      out_ids.clear();

      std::vector<size_tu> bins_to_visit;
      bins_to_visit.reserve(64);
      bins_for_window(focal_vars, clim, bins_to_visit);

      const double radius_sq = clim.radius * clim.radius;

      for (size_tu b = 0; b < bins_to_visit.size(); ++b) {
            size_tu bid = bins_to_visit[b];
            size_tu start = bin_offsets[bid];
            size_tu stop = bin_offsets[bid + 1];

            for (size_tu p = start; p < stop; ++p) {
                  index_t rid = bin_ids_flat[p];
                  bool ok;

                  if (clim.use_scalar) {
                        double dist_sq = 0.0;
                        for (size_tu kk = 0; kk < nvars; ++kk) {
                              double d = focal_vars[kk] - ref_clim(static_cast<size_tu>(rid), kk);
                              dist_sq += d * d;
                        }
                        ok = (dist_sq <= radius_sq);
                  } else {
                        ok = true;
                        for (size_tu kk = 0; kk < nvars; ++kk) {
                              if (std::fabs(focal_vars[kk] - ref_clim(static_cast<size_tu>(rid), kk)) >
                                        clim.max_abs_diff[kk]) {
                                    ok = false;
                                    break;
                              }
                        }
                  }

                  if (ok) out_ids.push_back(rid);
            }
      }
}

void LatticeIndex::expanding_search(const double* focal_vars,
                                    const ClimateConstraint& clim,
                                    std::vector<index_t>& out_ids,
                                    index_t k,
                                    size_tu max_bins_to_visit) const {
      (void)max_bins_to_visit;
      best_first_search(focal_vars, clim, out_ids, k);
}

} // namespace analogs
