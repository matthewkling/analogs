#pragma once

// Small dense Cholesky solver for weighted ridge regression.
// Designed for p+1 <= ~10 (intercept + a few covariates).
// No external linear algebra dependencies.

#include <vector>
#include <cmath>
#include <limits>

namespace analogs {

// In-place Cholesky factorization of symmetric positive-definite matrix A.
// A is stored as a flat vector in row-major order, dimension dim x dim.
// On success, lower triangle of A contains L such that A = L L'.
// Returns false if matrix is not positive definite (diagonal element <= 0).
inline bool cholesky_factor(std::vector<double>& A, int dim) {
      for (int j = 0; j < dim; ++j) {
            double sum = 0.0;
            for (int k = 0; k < j; ++k) {
                  sum += A[j * dim + k] * A[j * dim + k];
            }
            double diag = A[j * dim + j] - sum;
            if (diag <= 0.0) return false;
            A[j * dim + j] = std::sqrt(diag);

            for (int i = j + 1; i < dim; ++i) {
                  double sum2 = 0.0;
                  for (int k = 0; k < j; ++k) {
                        sum2 += A[i * dim + k] * A[j * dim + k];
                  }
                  A[i * dim + j] = (A[i * dim + j] - sum2) / A[j * dim + j];
            }
      }
      return true;
}

// Solve L L' x = b given Cholesky factor L (lower triangle of A, row-major).
// Overwrites b with the solution x.
inline void cholesky_solve(const std::vector<double>& A, int dim, std::vector<double>& b) {
      // Forward substitution: L y = b
      for (int i = 0; i < dim; ++i) {
            for (int k = 0; k < i; ++k) {
                  b[i] -= A[i * dim + k] * b[k];
            }
            b[i] /= A[i * dim + i];
      }
      // Back substitution: L' x = y
      for (int i = dim - 1; i >= 0; --i) {
            for (int k = i + 1; k < dim; ++k) {
                  b[i] -= A[k * dim + i] * b[k];
            }
            b[i] /= A[i * dim + i];
      }
}

// Solve weighted ridge regression for one focal cell.
//
// Given n_analogs observations with:
//   - response values y[v] for each of n_vars value variables
//   - covariates cov_ptr[j + p * cov_stride] for analog j, covariate p
//   - weights w[j] (combined distance × sample weight)
//   - ridge penalty lambda (applied to covariate coefficients only, not intercept)
//
// Writes (n_covs + 1) coefficients per value variable into out_coeffs,
// laid out as: [intercept_v0, beta1_v0, ..., betap_v0, intercept_v1, ...].
//
// Returns false if the system is singular (only possible when lambda == 0).
inline bool solve_ridge(
            const std::vector<std::size_t>& analog_indices,  // 0-based indices into pool
            const std::vector<double>& weights,              // combined weights for each analog
            const double* values_ptr,                        // pool values matrix (column-major)
            int values_stride,                               // stride for values (n_ref)
            int n_vars,                                      // number of value variables
            const double* cov_ptr,                           // pool covariates matrix (column-major)
            int cov_stride,                                  // stride for covariates (n_ref)
            int n_covs,                                      // number of covariate columns
            double lambda,                                   // ridge penalty
            double* out_coeffs                               // output: n_vars * (n_covs + 1)
) {
      const int dim = n_covs + 1;  // intercept + covariates
      const int n = static_cast<int>(analog_indices.size());

      // Build X'WX (symmetric, stored as flat row-major dim x dim)
      std::vector<double> XtWX(dim * dim, 0.0);

      // Build X'Wy for each value variable
      // Layout: [XtWy_var0[0..dim-1], XtWy_var1[0..dim-1], ...]
      std::vector<double> XtWy(n_vars * dim, 0.0);

      for (int t = 0; t < n; ++t) {
            const std::size_t j = analog_indices[t];
            const double w = weights[t];

            // Build row of X: [1, cov1, cov2, ..., covp]
            // We don't need to store the full design matrix — just accumulate
            // the outer product w * x_t * x_t' into XtWX

            // Intercept × intercept
            XtWX[0] += w;  // 1 * w * 1

            // Intercept × covariates and covariates × intercept
            for (int p = 0; p < n_covs; ++p) {
                  double cp = cov_ptr[j + p * cov_stride];
                  double wcp = w * cp;
                  XtWX[0 * dim + (p + 1)] += wcp;        // row 0, col p+1
                  XtWX[(p + 1) * dim + 0] += wcp;        // row p+1, col 0
            }

            // Covariate × covariate block
            for (int p1 = 0; p1 < n_covs; ++p1) {
                  double cp1 = cov_ptr[j + p1 * cov_stride];
                  double wcp1 = w * cp1;
                  for (int p2 = p1; p2 < n_covs; ++p2) {
                        double cp2 = cov_ptr[j + p2 * cov_stride];
                        double val = wcp1 * cp2;
                        XtWX[(p1 + 1) * dim + (p2 + 1)] += val;
                        if (p2 != p1) {
                              XtWX[(p2 + 1) * dim + (p1 + 1)] += val;
                        }
                  }
            }

            // X'Wy for each value variable
            for (int v = 0; v < n_vars; ++v) {
                  double yv = values_ptr[j + v * values_stride];
                  // Skip NA values
                  if (ISNAN(yv)) continue;

                  double wyv = w * yv;
                  XtWy[v * dim + 0] += wyv;  // intercept component

                  for (int p = 0; p < n_covs; ++p) {
                        double cp = cov_ptr[j + p * cov_stride];
                        XtWy[v * dim + (p + 1)] += wyv * cp;
                  }
            }
      }

      // Add ridge penalty to covariate diagonal (not intercept)
      for (int p = 0; p < n_covs; ++p) {
            XtWX[(p + 1) * dim + (p + 1)] += lambda;
      }

      // Cholesky factorization of XtWX
      if (!cholesky_factor(XtWX, dim)) {
            // Singular system — write NAs for all coefficients
            for (int v = 0; v < n_vars; ++v) {
                  for (int d = 0; d < dim; ++d) {
                        out_coeffs[v * dim + d] = NA_REAL;
                  }
            }
            return false;
      }

      // Solve for each value variable using the same factorization
      for (int v = 0; v < n_vars; ++v) {
            // Copy XtWy for this variable (cholesky_solve overwrites it)
            std::vector<double> rhs(XtWy.begin() + v * dim,
                                    XtWy.begin() + (v + 1) * dim);
            cholesky_solve(XtWX, dim, rhs);

            // Write coefficients to output
            for (int d = 0; d < dim; ++d) {
                  out_coeffs[v * dim + d] = rhs[d];
            }
      }

      return true;
}

} // namespace analogs
