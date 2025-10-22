#ifndef LATTICE_INDEX_H
#define LATTICE_INDEX_H

#include <RcppArmadillo.h>
#include <map>
#include <vector>
#include <queue>
#include <cmath>

// N-dimensional bin ID (e.g., {3, 5} for 2D, {2, 4, 1} for 3D)
typedef std::vector<int> BinID;

// Helper to compare BinIDs for map
struct BinIDCompare {
      bool operator()(const BinID& a, const BinID& b) const {
            return a < b;
      }
};

// Generic N-dimensional lattice index
template<int N_DIMS>
class LatticeIndex {
private:
      bool built;
      int n_bins_per_dim;
      arma::vec bin_min;
      arma::vec bin_max;
      arma::vec bin_width;
      std::map<BinID, std::vector<int>, BinIDCompare> bins;

      // Get bin ID for a point
      BinID get_bin_id(const arma::rowvec& point) const {
            BinID bid(N_DIMS);
            for (int d = 0; d < N_DIMS; d++) {
                  int b = std::floor((point[d] - bin_min[d]) / bin_width[d]);
                  b = std::max(0, std::min(n_bins_per_dim - 1, b));
                  bid[d] = b;
            }
            return bid;
      }

      // Calculate distance from point to bin (minimum distance to any point in bin)
      double distance_to_bin(const arma::rowvec& point, const BinID& bid) const {
            double dist_sq = 0.0;

            for (int d = 0; d < N_DIMS; d++) {
                  double bin_start = bin_min[d] + bid[d] * bin_width[d];
                  double bin_end = bin_start + bin_width[d];

                  // Distance to bin in this dimension
                  double d_dim = 0.0;
                  if (point[d] < bin_start) {
                        d_dim = bin_start - point[d];
                  } else if (point[d] > bin_end) {
                        d_dim = point[d] - bin_end;
                  }
                  // else point is inside bin range in this dimension

                  dist_sq += d_dim * d_dim;
            }

            return std::sqrt(dist_sq);
      }

      // Get all bins within L-infinity (box) distance of center
      std::vector<BinID> get_bins_in_box(const arma::rowvec& center, double radius) const {
            std::vector<BinID> result;

            // Find bin range in each dimension
            std::vector<int> min_bins(N_DIMS), max_bins(N_DIMS);

            for (int d = 0; d < N_DIMS; d++) {
                  double low = center[d] - radius;
                  double high = center[d] + radius;

                  min_bins[d] = std::max(0,
                                         (int)std::floor((low - bin_min[d]) / bin_width[d]));
                  max_bins[d] = std::min(n_bins_per_dim - 1,
                                         (int)std::floor((high - bin_min[d]) / bin_width[d]));
            }

            // Generate all bin combinations in range (recursive)
            BinID current(N_DIMS);
            enumerate_bins_recursive(min_bins, max_bins, 0, current, result);

            return result;
      }

      // Recursively enumerate all bins in range
      void enumerate_bins_recursive(const std::vector<int>& min_bins,
                                    const std::vector<int>& max_bins,
                                    int dim,
                                    BinID& current,
                                    std::vector<BinID>& result) const {
            if (dim == N_DIMS) {
                  result.push_back(current);
                  return;
            }

            for (int b = min_bins[dim]; b <= max_bins[dim]; b++) {
                  current[dim] = b;
                  enumerate_bins_recursive(min_bins, max_bins, dim + 1, current, result);
            }
      }

public:
      LatticeIndex() : built(false), n_bins_per_dim(0) {
            bin_min.set_size(N_DIMS);
            bin_max.set_size(N_DIMS);
            bin_width.set_size(N_DIMS);
      }

      bool is_built() const { return built; }

      // Build index from points matrix [n_points × N_DIMS]
      void build(const arma::mat& points, int bins_per_dim = 20) {
            if (points.n_cols != N_DIMS) {
                  Rcpp::stop("Points matrix must have N_DIMS columns");
            }

            n_bins_per_dim = bins_per_dim;

            // Calculate bin boundaries for each dimension
            for (int d = 0; d < N_DIMS; d++) {
                  bin_min[d] = points.col(d).min();
                  bin_max[d] = points.col(d).max();

                  // Add small epsilon to max to ensure all points fit
                  double range = bin_max[d] - bin_min[d];
                  bin_max[d] += range * 0.001;

                  bin_width[d] = (bin_max[d] - bin_min[d]) / n_bins_per_dim;
            }

            // Assign points to bins
            for (size_t i = 0; i < points.n_rows; i++) {
                  BinID bid = get_bin_id(points.row(i));
                  bins[bid].push_back(i);
            }

            built = true;
      }

      // Query: return all points within radius (box filter)
      std::vector<int> query_box(const arma::rowvec& center, double radius) const {
            if (!built) return {};

            std::vector<int> result;
            auto bin_ids = get_bins_in_box(center, radius);

            for (const auto& bid : bin_ids) {
                  auto it = bins.find(bid);
                  if (it != bins.end()) {
                        const auto& indices = it->second;
                        result.insert(result.end(), indices.begin(), indices.end());
                  }
            }

            return result;
      }

      // Query: expanding search (returns points in order of bin proximity)
      // Stops when max_count found or max_radius exceeded
      std::vector<int> query_expanding(const arma::rowvec& center,
                                       double max_radius,
                                       int max_count = -1) const {
            if (!built) return {};

            std::vector<int> result;

            // Priority queue: (distance_to_bin, bin_id)
            typedef std::pair<double, BinID> QueueItem;
            auto compare = [](const QueueItem& a, const QueueItem& b) {
                  return a.first > b.first;  // Min heap
            };
            std::priority_queue<QueueItem, std::vector<QueueItem>, decltype(compare)> pq(compare);

            std::set<BinID, BinIDCompare> visited;

            // Start with bin containing center
            BinID start_bin = get_bin_id(center);
            pq.push({0.0, start_bin});

            while (!pq.empty()) {
                  QueueItem top = pq.top();
                  pq.pop();

                  double dist = top.first;
                  BinID bid = top.second;

                  // Early exit if beyond max radius
                  if (dist > max_radius) break;

                  // Skip if already visited
                  if (visited.count(bid)) continue;
                  visited.insert(bid);

                  // Add points from this bin
                  auto it = bins.find(bid);
                  if (it != bins.end()) {
                        const auto& indices = it->second;
                        result.insert(result.end(), indices.begin(), indices.end());

                        // Early exit if we have enough
                        if (max_count > 0 && (int)result.size() >= max_count) {
                              break;
                        }
                  }

                  // Add neighboring bins to queue
                  for (int d = 0; d < N_DIMS; d++) {
                        // Check bins before and after in this dimension
                        for (int delta : {-1, 1}) {
                              BinID neighbor = bid;
                              neighbor[d] += delta;

                              // Check bounds
                              if (neighbor[d] >= 0 && neighbor[d] < n_bins_per_dim) {
                                    if (!visited.count(neighbor)) {
                                          double neighbor_dist = distance_to_bin(center, neighbor);
                                          pq.push({neighbor_dist, neighbor});
                                    }
                              }
                        }
                  }
            }

            return result;
      }

      // Get statistics (for debugging/tuning)
      void print_stats() const {
            if (!built) {
                  Rcpp::Rcout << "Index not built" << std::endl;
                  return;
            }

            int total_points = 0;
            int non_empty_bins = 0;
            int max_bin_size = 0;

            for (const auto& kv : bins) {
                  int size = kv.second.size();
                  if (size > 0) {
                        non_empty_bins++;
                        total_points += size;
                        max_bin_size = std::max(max_bin_size, size);
                  }
            }

            double avg_bin_size = non_empty_bins > 0 ?
            (double)total_points / non_empty_bins : 0;

            Rcpp::Rcout << "LatticeIndex<" << N_DIMS << "> stats:" << std::endl;
            Rcpp::Rcout << "  Bins per dim: " << n_bins_per_dim << std::endl;
            Rcpp::Rcout << "  Total bins: " << std::pow(n_bins_per_dim, N_DIMS) << std::endl;
            Rcpp::Rcout << "  Non-empty bins: " << non_empty_bins << std::endl;
            Rcpp::Rcout << "  Total points: " << total_points << std::endl;
            Rcpp::Rcout << "  Avg bin size: " << avg_bin_size << std::endl;
            Rcpp::Rcout << "  Max bin size: " << max_bin_size << std::endl;
      }
};

#endif // LATTICE_INDEX_H
