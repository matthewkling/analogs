#pragma once

// Small dense Cholesky solver for weighted ridge regression.
// Designed for p+1 <= ~10 (intercept + a few covariates).
// No external linear algebra dependencies.
// Returns coefficients and standard errors.

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

// Compute diagonal of (X'WX + lambda*I)^{-1} from its Cholesky factor L.
// Uses the identity: if A = L L', then A^{-1} = L^{-T} L^{-1},
// and diag(A^{-1})_i = sum_k (L^{-1})_{ki}^2.
// We compute L^{-1} column by column via forward substitution.
inline void cholesky_inv_diag(const std::vector<double>& L, int dim,
                              std::vector<double>& inv_diag) {
      inv_diag.resize(dim);

      // For each column j of L^{-1}, solve L * col_j = e_j
      // Then inv_diag[j] = sum of squared elements in col_j
      for (int j = 0; j < dim; ++j) {
            std::vector<double> col(dim, 0.0);
            col[j] = 1.0;
            for (int i = j; i < dim; ++i) {
                  for (int k = j; k < i; ++k) {
                        col[i] -= L[i * dim + k] * col[k];
                  }
                  col[i] /= L[i * dim + i];
            }

            double sumsq = 0.0;
            for (int i = j; i < dim; ++i) {
                  sumsq += col[i] * col[i];
            }
            inv_diag[j] = sumsq;
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
// Writes output into out_coeffs and out_se, each of size n_vars * (n_covs + 1).
// Layout: [intercept_v0, beta1_v0, ..., betap_v0, intercept_v1, ...].
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
            double* out_coeffs,                              // output: n_vars * (n_covs + 1)
            double* out_se                                   // output: n_vars * (n_covs + 1)
) {
      const int dim = n_covs + 1;  // intercept + covariates
      const int n = static_cast<int>(analog_indices.size());

      // Build X'WX (symmetric, stored as flat row-major dim x dim)
      std::vector<double> XtWX(dim * dim, 0.0);

      // Build X'Wy for each value variable
      std::vector<double> XtWy(n_vars * dim, 0.0);

      // Track sum of weights for sigma^2 denominator
      double sum_w = 0.0;

      for (int t = 0; t < n; ++t) {
            const std::size_t j = analog_indices[t];
            const double w = weights[t];
            sum_w += w;

            // Intercept × intercept
            XtWX[0] += w;

            // Intercept × covariates and covariates × intercept
            for (int p = 0; p < n_covs; ++p) {
                  double cp = cov_ptr[j + p * cov_stride];
                  double wcp = w * cp;
                  XtWX[0 * dim + (p + 1)] += wcp;
                  XtWX[(p + 1) * dim + 0] += wcp;
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
                  if (ISNAN(yv)) continue;

                  double wyv = w * yv;
                  XtWy[v * dim + 0] += wyv;

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

      // Cholesky factorization
      std::vector<double> L = XtWX;
      if (!cholesky_factor(L, dim)) {
            // Singular system — write NAs for all coefficients and SEs
            for (int v = 0; v < n_vars; ++v) {
                  for (int d = 0; d < dim; ++d) {
                        out_coeffs[v * dim + d] = NA_REAL;
                        out_se[v * dim + d] = NA_REAL;
                  }
            }
            return false;
      }

      // Compute diagonal of (X'WX + lambda*I)^{-1} — shared across all value variables
      std::vector<double> inv_diag;
      cholesky_inv_diag(L, dim, inv_diag);

      // Solve for each value variable and compute SEs
      for (int v = 0; v < n_vars; ++v) {
            // Solve for coefficients
            std::vector<double> rhs(XtWy.begin() + v * dim,
                                    XtWy.begin() + (v + 1) * dim);
            cholesky_solve(L, dim, rhs);

            // Write coefficients
            for (int d = 0; d < dim; ++d) {
                  out_coeffs[v * dim + d] = rhs[d];
            }

            // Compute weighted residual variance: sigma^2 = sum(w * e^2) / (sum(w) - dim)
            double wss = 0.0;
            for (int t = 0; t < n; ++t) {
                  const std::size_t j = analog_indices[t];
                  const double w = weights[t];
                  double yv = values_ptr[j + v * values_stride];
                  if (ISNAN(yv)) continue;

                  double fitted = rhs[0];
                  for (int p = 0; p < n_covs; ++p) {
                        fitted += rhs[p + 1] * cov_ptr[j + p * cov_stride];
                  }
                  double resid = yv - fitted;
                  wss += w * resid * resid;
            }

            double denom = sum_w - static_cast<double>(dim);
            double sigma2;
            if (denom > 0.0) {
                  sigma2 = wss / denom;
            } else {
                  sigma2 = NA_REAL;
            }

            // Write standard errors: se_j = sqrt(sigma^2 * inv_diag_j)
            for (int d = 0; d < dim; ++d) {
                  if (ISNAN(sigma2) || sigma2 < 0.0) {
                        out_se[v * dim + d] = NA_REAL;
                  } else {
                        out_se[v * dim + d] = std::sqrt(sigma2 * inv_diag[d]);
                  }
            }
      }

      return true;
}

} // namespace analogs
