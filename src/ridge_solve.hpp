#pragma once

// Small dense Cholesky solver for weighted ridge regression.
// Designed for p+1 <= ~10 (intercept + a few covariates).
// No external linear algebra dependencies.
//
// Supports coefficient-only, ESS-based, and design-based (sandwich) SE
// computation, selected via SeCode.

#include <vector>
#include <cmath>
#include <limits>
#include "se_code.hpp"

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

// Compute diagonal of A^{-1} from Cholesky factor L (A = L L').
// Uses: A^{-1} = L^{-T} L^{-1}, so diag(A^{-1})_j = sum_k (L^{-1})_{kj}^2.
// We compute column j of L^{-1} by forward-solving L x = e_j.
inline void cholesky_inv_diag(const std::vector<double>& L, int dim,
                              std::vector<double>& inv_diag) {
      inv_diag.resize(dim);

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

// Compute full inverse A^{-1} from Cholesky factor L (A = L L'), stored
// row-major in out_inv (dim x dim). Only the symmetric inverse is produced.
// Used for the sandwich SE: needs full A^{-1}, not just diag.
inline void cholesky_inverse(const std::vector<double>& L, int dim,
                             std::vector<double>& out_inv) {
      out_inv.assign(dim * dim, 0.0);

      // Solve A x_col = e_j for each column j.
      // This is a triangular solve L y = e_j (forward) then L' x = y (back).
      std::vector<double> rhs(dim);
      for (int j = 0; j < dim; ++j) {
            std::fill(rhs.begin(), rhs.end(), 0.0);
            rhs[j] = 1.0;

            // Forward: L y = e_j
            for (int i = 0; i < dim; ++i) {
                  for (int k = 0; k < i; ++k) {
                        rhs[i] -= L[i * dim + k] * rhs[k];
                  }
                  rhs[i] /= L[i * dim + i];
            }
            // Back: L' x = y
            for (int i = dim - 1; i >= 0; --i) {
                  for (int k = i + 1; k < dim; ++k) {
                        rhs[i] -= L[k * dim + i] * rhs[k];
                  }
                  rhs[i] /= L[i * dim + i];
            }

            // Write column j of inverse (it's symmetric; fill whole column)
            for (int i = 0; i < dim; ++i) {
                  out_inv[i * dim + j] = rhs[i];
            }
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
// Writes output into out_coeffs (always) and out_se (only if se_code != NONE).
// Layout: [intercept_v0, beta1_v0, ..., betap_v0, intercept_v1, ...].
//
// SE formulas:
//   ESS:    Var(beta) = sigma2_ess * (X'WX + lambda I)^{-1}
//           sigma2_ess = (RSS_w / sum_w) * n_eff / (n_eff - p)
//           n_eff = (sum_w)^2 / sum_w_sq
//   DESIGN: Var(beta) = A^{-1} M A^{-1}, A = X'WX + lambda I
//           M_{ab} = sum_t w_t^2 r_t^2 x_{ta} x_{tb}
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
            SeCode se_code,                                  // which SE variant (or NONE)
            double* out_coeffs,                              // output: n_vars * (n_covs + 1)
            double* out_se                                   // output: n_vars * (n_covs + 1); ignored if se_code == NONE
) {
      const int dim = n_covs + 1;  // intercept + covariates
      const int n = static_cast<int>(analog_indices.size());
      const bool want_se = (se_code != SeCode::NONE);

      // Build X'WX (symmetric, stored as flat row-major dim x dim)
      std::vector<double> XtWX(dim * dim, 0.0);

      // Build X'Wy for each value variable
      std::vector<double> XtWy(n_vars * dim, 0.0);

      // Running sums for SE support
      double sum_w = 0.0;
      double sum_w_sq = 0.0;  // only used for SE

      for (int t = 0; t < n; ++t) {
            const std::size_t j = analog_indices[t];
            const double w = weights[t];
            sum_w += w;
            if (want_se) sum_w_sq += w * w;

            // Intercept × intercept
            XtWX[0] += w;

            // Intercept × covariates (symmetric)
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

      // Cholesky factorization of A = X'WX + lambda I_covs
      std::vector<double> L = XtWX;
      if (!cholesky_factor(L, dim)) {
            // Singular: write NAs for coefficients (and SEs if requested)
            for (int v = 0; v < n_vars; ++v) {
                  for (int d = 0; d < dim; ++d) {
                        out_coeffs[v * dim + d] = NA_REAL;
                        if (want_se) out_se[v * dim + d] = NA_REAL;
                  }
            }
            return false;
      }

      // ESS-based SE: one inv_diag suffices (shared across vars); sigma2 per var.
      std::vector<double> inv_diag;
      // DESIGN-based SE: full inverse of A needed to form A^{-1} M A^{-1}.
      std::vector<double> A_inv;
      if (want_se) {
            if (se_code == SeCode::ESS) {
                  cholesky_inv_diag(L, dim, inv_diag);
            } else {  // DESIGN
                  cholesky_inverse(L, dim, A_inv);
            }
      }

      // n_eff for ESS variant
      double n_eff = 0.0;
      if (want_se && se_code == SeCode::ESS && sum_w_sq > 0.0) {
            n_eff = (sum_w * sum_w) / sum_w_sq;
      }

      // Solve per variable, write coefficients, compute SEs if requested.
      for (int v = 0; v < n_vars; ++v) {
            // Solve for coefficients
            std::vector<double> rhs(XtWy.begin() + v * dim,
                                    XtWy.begin() + (v + 1) * dim);
            cholesky_solve(L, dim, rhs);

            for (int d = 0; d < dim; ++d) {
                  out_coeffs[v * dim + d] = rhs[d];
            }

            if (!want_se) continue;

            // ---- SE computation for this variable ----
            // Both variants need a pass over analogs to compute residuals.
            // ESS needs RSS_w = sum w * r^2.
            // DESIGN needs M = sum w^2 * r^2 * x x' (symmetric dim×dim).
            double rss_w = 0.0;
            std::vector<double> M;  // only used for design
            if (se_code == SeCode::DESIGN) {
                  M.assign(dim * dim, 0.0);
            }

            bool any_valid = false;

            for (int t = 0; t < n; ++t) {
                  const std::size_t j = analog_indices[t];
                  const double w = weights[t];
                  double yv = values_ptr[j + v * values_stride];
                  if (ISNAN(yv)) continue;
                  any_valid = true;

                  double fitted = rhs[0];
                  for (int p = 0; p < n_covs; ++p) {
                        fitted += rhs[p + 1] * cov_ptr[j + p * cov_stride];
                  }
                  const double resid = yv - fitted;

                  if (se_code == SeCode::ESS) {
                        rss_w += w * resid * resid;
                  } else {
                        // DESIGN: accumulate outer product w^2 r^2 x x'
                        const double w2r2 = w * w * resid * resid;
                        // Build x row implicitly: x = [1, cov1, cov2, ...]
                        // M[0][0] += w2r2
                        M[0] += w2r2;
                        for (int p = 0; p < n_covs; ++p) {
                              double cp = cov_ptr[j + p * cov_stride];
                              double w2r2_cp = w2r2 * cp;
                              M[0 * dim + (p + 1)] += w2r2_cp;
                              M[(p + 1) * dim + 0] += w2r2_cp;
                        }
                        for (int p1 = 0; p1 < n_covs; ++p1) {
                              double cp1 = cov_ptr[j + p1 * cov_stride];
                              double w2r2_cp1 = w2r2 * cp1;
                              for (int p2 = p1; p2 < n_covs; ++p2) {
                                    double cp2 = cov_ptr[j + p2 * cov_stride];
                                    double val = w2r2_cp1 * cp2;
                                    M[(p1 + 1) * dim + (p2 + 1)] += val;
                                    if (p2 != p1) {
                                          M[(p2 + 1) * dim + (p1 + 1)] += val;
                                    }
                              }
                        }
                  }
            }

            if (!any_valid) {
                  for (int d = 0; d < dim; ++d) out_se[v * dim + d] = NA_REAL;
                  continue;
            }

            if (se_code == SeCode::ESS) {
                  // sigma2_ess = (RSS_w / sum_w) * n_eff / (n_eff - p)
                  const double df = n_eff - static_cast<double>(dim);
                  if (sum_w <= 0.0 || n_eff <= 0.0 || df <= 0.0) {
                        for (int d = 0; d < dim; ++d) out_se[v * dim + d] = NA_REAL;
                        continue;
                  }
                  const double sigma2 = (rss_w / sum_w) * (n_eff / df);
                  if (!(sigma2 >= 0.0) || !std::isfinite(sigma2)) {
                        for (int d = 0; d < dim; ++d) out_se[v * dim + d] = NA_REAL;
                        continue;
                  }
                  for (int d = 0; d < dim; ++d) {
                        out_se[v * dim + d] = std::sqrt(sigma2 * inv_diag[d]);
                  }
            } else {
                  // DESIGN: diag(A^{-1} M A^{-1})_d = row_d(A^{-1} M) dot col_d(A^{-1})
                  // Since A^{-1} is symmetric, row_d(A^{-1}) = col_d(A^{-1}).
                  // Compute diag directly:
                  //   diag_d = sum_{a,b} A^{-1}_{da} * M_{ab} * A^{-1}_{bd}
                  for (int d = 0; d < dim; ++d) {
                        double acc = 0.0;
                        for (int a = 0; a < dim; ++a) {
                              double Aid_da = A_inv[d * dim + a];
                              double inner = 0.0;
                              for (int b = 0; b < dim; ++b) {
                                    inner += M[a * dim + b] * A_inv[b * dim + d];
                              }
                              acc += Aid_da * inner;
                        }
                        if (!(acc >= 0.0) || !std::isfinite(acc)) {
                              out_se[v * dim + d] = NA_REAL;
                        } else {
                              out_se[v * dim + d] = std::sqrt(acc);
                        }
                  }
            }
      }

      return true;
}

} // namespace analogs
