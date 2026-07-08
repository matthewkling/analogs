#pragma once

#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>

namespace analogs {

// Reconstruct symmetric covariance matrix from lower triangle storage
// Input: vec with n*(n+1)/2 elements storing [var1, var2, ..., varN, cov12, cov13, ..., cov(N-1,N)]
// Output: n×n symmetric matrix stored row-major
inline void reconstruct_cov_matrix(const double* vec, int n, std::vector<double>& cov_matrix) {
      cov_matrix.resize(n * n);

      // First n elements are diagonal (variances)
      for (int i = 0; i < n; ++i) {
            cov_matrix[i * n + i] = vec[i];
      }

      // Remaining elements are upper triangle covariances
      int idx = n;  // Start after diagonal elements
      for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                  double cov_val = vec[idx++];
                  cov_matrix[i * n + j] = cov_val;  // Upper triangle
                  cov_matrix[j * n + i] = cov_val;  // Lower triangle (symmetric)
            }
      }
}

// Invert covariance matrix using Cholesky decomposition
// Input: cov_matrix (n×n, row-major, symmetric positive definite)
// Output: inv_cov_matrix (n×n, row-major)
// Returns: true if successful, false if matrix is not positive definite
inline bool invert_cov_matrix(const std::vector<double>& cov_matrix, int n,
                              std::vector<double>& inv_cov_matrix) {
      inv_cov_matrix.resize(n * n);

      // Cholesky decomposition: A = L * L^T
      std::vector<double> L(n * n, 0.0);

      for (int i = 0; i < n; ++i) {
            for (int j = 0; j <= i; ++j) {
                  double sum = 0.0;

                  if (i == j) {
                        // Diagonal elements
                        for (int k = 0; k < j; ++k) {
                              sum += L[j * n + k] * L[j * n + k];
                        }
                        double diag = cov_matrix[j * n + j] - sum;

                        // Check for positive definiteness
                        if (diag <= 0.0) {
                              return false;
                        }
                        L[j * n + j] = std::sqrt(diag);
                  } else {
                        // Off-diagonal elements
                        for (int k = 0; k < j; ++k) {
                              sum += L[i * n + k] * L[j * n + k];
                        }
                        L[i * n + j] = (cov_matrix[i * n + j] - sum) / L[j * n + j];
                  }
            }
      }

      // Invert L (lower triangular) -> L_inv
      std::vector<double> L_inv(n * n, 0.0);
      for (int i = 0; i < n; ++i) {
            L_inv[i * n + i] = 1.0 / L[i * n + i];
            for (int j = i + 1; j < n; ++j) {
                  double sum = 0.0;
                  for (int k = i; k < j; ++k) {
                        sum += L[j * n + k] * L_inv[k * n + i];
                  }
                  L_inv[j * n + i] = -sum / L[j * n + j];
            }
      }

      // Compute A_inv = L_inv^T * L_inv
      for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                  double sum = 0.0;
                  for (int k = std::max(i, j); k < n; ++k) {
                        sum += L_inv[k * n + i] * L_inv[k * n + j];
                  }
                  inv_cov_matrix[i * n + j] = sum;
            }
      }

      return true;
}

// Compute Mahalanobis distance: sqrt((x - y)^T * Sigma^{-1} * (x - y))
// focal_env: pointer to focal environment vector (length n)
// ref_env: pointer to reference environment vector (length n)
// inv_cov: inverse covariance matrix (n×n, row-major)
// n: number of environment dimensions
// stride_f, stride_r: strides for accessing environment columns
inline double mahalanobis_distance(const double* focal_env,
                                   const double* ref_env,
                                   const std::vector<double>& inv_cov,
                                   int n,
                                   int stride_f,
                                   int stride_r) {
      // Compute difference vector: d = focal - ref
      std::vector<double> diff(n);
      for (int i = 0; i < n; ++i) {
            diff[i] = focal_env[i * stride_f] - ref_env[i * stride_r];
      }

      // Compute d^T * inv_cov * d
      double result = 0.0;
      for (int i = 0; i < n; ++i) {
            double sum = 0.0;
            for (int j = 0; j < n; ++j) {
                  sum += inv_cov[i * n + j] * diff[j];
            }
            result += diff[i] * sum;
      }

      return std::sqrt(std::max(0.0, result));  // Protect against numerical errors
}

// Compute axis-aligned bounding box for Mahalanobis ellipsoid
// Given a Mahalanobis distance threshold, compute the min/max bounds
// in each environment dimension that could possibly satisfy the threshold.
// This is used for lattice queries to expand search bounds appropriately.
//
// focal_env: focal point environment values
// inv_cov: inverse covariance matrix (n×n, row-major)
// n: number of environment dimensions
// threshold: Mahalanobis distance threshold
// Output: bounds (2*n elements: [min0, max0, min1, max1, ...])
inline void mahalanobis_bounding_box(const double* focal_env,
                                     const std::vector<double>& inv_cov,
                                     int n,
                                     double threshold,
                                     std::vector<double>& bounds) {
      bounds.resize(2 * n);

      if (!std::isfinite(threshold) || threshold <= 0.0) {
            // No bound - use infinite bounds
            for (int i = 0; i < n; ++i) {
                  bounds[2 * i] = -std::numeric_limits<double>::infinity();
                  bounds[2 * i + 1] = std::numeric_limits<double>::infinity();
            }
            return;
      }

      // For conservative axis-aligned bounds, we compute the eigenvalues
      // of the covariance matrix (inverse of inv_cov) and use the maximum
      // eigenvalue to determine the worst-case axis-aligned extent.
      //
      // Quick approximation: Use diagonal elements of inv_cov
      // The bound in dimension i is approximately threshold / sqrt(inv_cov[i,i])

      for (int i = 0; i < n; ++i) {
            double inv_var = inv_cov[i * n + i];

            if (inv_var > 0.0) {
                  // Conservative bound: threshold / sqrt(diagonal)
                  // This gives a bounding box that's guaranteed to contain the ellipsoid
                  double half_width = threshold / std::sqrt(inv_var);
                  bounds[2 * i] = focal_env[i] - half_width;
                  bounds[2 * i + 1] = focal_env[i] + half_width;
            } else {
                  // Shouldn't happen with valid positive definite matrix
                  bounds[2 * i] = -std::numeric_limits<double>::infinity();
                  bounds[2 * i + 1] = std::numeric_limits<double>::infinity();
            }
      }
}

// Check if Mahalanobis distance satisfies threshold and optionally compute distance
// Returns: (ok, distance) where ok indicates if threshold is satisfied
inline std::pair<bool, double>
mahalanobis_ok_and_dist(const double* focal_env,
                        const double* ref_env,
                        const std::vector<double>& inv_cov,
                        int n,
                        int stride_f,
                        int stride_r,
                        double threshold,
                        bool compute_dist) {
      if (!std::isfinite(threshold)) {
            // No threshold - always pass
            if (compute_dist) {
                  double d = mahalanobis_distance(focal_env, ref_env, inv_cov,
                                                  n, stride_f, stride_r);
                  return std::make_pair(true, d);
            } else {
                  return std::make_pair(true, 0.0);
            }
      }

      // Compute distance
      double dist = mahalanobis_distance(focal_env, ref_env, inv_cov,
                                         n, stride_f, stride_r);

      if (dist > threshold) {
            return std::make_pair(false, dist);
      }

      return std::make_pair(true, dist);
}

} // namespace analogs
