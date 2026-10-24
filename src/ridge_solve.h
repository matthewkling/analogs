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
#include "se_code.h"

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

// -----------------------------------------------------------------------------
// Weighted ridge regression for one neighborhood, split into reusable stages.
//
// The solve is decomposed so that every consumer shares everything downstream
// of the sufficient statistics. The exact engine sees individual analogs and
// fills the moments by accumulation; a binned (FFT) engine sees only
// kernel-weighted moment sums and fills them directly.
//
//   RidgeMoments        X'WX, X'Wy, sum_w, sum_w_sq (sufficient statistics),
//                       shared across response variables plus per-variable
//                       copies for variables with missing responses
//   accumulate_moments  fill RidgeMoments from an explicit analog list
//   factor_ridge        add the information-scaled penalty, Cholesky-factor
//   solve_coefs         coefficients for every response variable
//   kish_n_eff          Kish effective sample size from sum_w, sum_w_sq
//   ess_se_from_rss     ESS-framed SEs for one variable, given its RSS_w
//   sandwich_se_diag    design-framed SEs for one variable, given its meat M
//   solve_ridge         exact-engine entry point composing the above, plus the
//                       residual pass that ESS and DESIGN SEs need
//
// Model: n analogs with response values y[v] for each of n_vars variables,
// covariates x (intercept prepended), and combined weights w (distance x
// sample x area x user). The ridge penalty applies to covariate coefficients
// only, not the intercept.
//
// Penalty: the system solved is A = X'WX + (lambda * sum_w / info) I_covs.
// Since X'WX ~ sum_w * V (V = weighted covariate covariance; the unpenalized
// intercept absorbs covariate means), each slope is shrunk by roughly
// info*v / (info*v + lambda): lambda acts as lambda pseudo-observations.
// `info` is the neighborhood's information content in units of
// "full-relevance observations", computed by the caller as
// sum(kernel * area * user / mean_user), i.e. the combined weights with the
// downsampling weight excluded and user weights rescaled to mean 1. With no
// user weights and no downsampling, info == sum_w and the penalty is plain
// lambda. Rescaling all user weights by a constant leaves the fit unchanged.
// The penalty is computed once per neighborhood from the full analog set and
// shared by all response variables (see "Missing responses" below).
//
// Missing responses: an analog whose response v is NA is dropped from every
// term of variable v's fit (X'WX, X'Wy, sum_w, sum_w_sq, residuals), so each
// variable is fit on exactly the analogs that observe it. Variables with no
// missing responses in the neighborhood share one X'WX and one factorization.
// The penalty uses the full-neighborhood ratio sum_w / info for every
// variable; this equals the per-variable ratio exactly when no downsampling
// is active, and is a close approximation otherwise.
//
// SE formulas:
//   ESS:    Var(beta) = sigma2_ess * A^{-1}
//           sigma2_ess = (RSS_w / sum_w) * n_eff / (n_eff - p)
//           n_eff = (sum_w)^2 / sum_w_sq
//   DESIGN: Var(beta) = A^{-1} M A^{-1}
//           M_{ab} = sum_t w_t^2 r_t^2 x_{ta} x_{tb}
//
// Output layout (coefficients and SEs alike):
//   [intercept_v0, beta1_v0, ..., betap_v0, intercept_v1, ...]
// -----------------------------------------------------------------------------

// Sufficient statistics for one neighborhood's weighted regression.
//
// The "shared" system (XtWX, sum_w, sum_w_sq) covers every analog. A variable
// with at least one missing response among the analogs gets its own system,
// accumulated over only the analogs that observe it (own[v] >= 0 indexes it).
// Accumulating the subset directly, rather than subtracting the missing rows
// from the shared system, avoids cancellation error. XtWX is stored
// UNPENALIZED; factor_ridge() adds the penalty to a copy.
struct RidgeMoments {
      int dim = 0;                       // n_covs + 1 (intercept + covariates)
      int n_vars = 0;                    // number of response variables
      std::vector<double> XtWX;          // shared: dim x dim, row-major, symmetric
      double sum_w = 0.0;                // shared
      double sum_w_sq = 0.0;             // shared; only needed for ESS SEs
      std::vector<double> XtWy;          // n_vars x dim, row-major (NA rows skipped)
      std::vector<int> own;              // per variable: -1 = shared, else own-system slot
      std::vector<double> own_XtWX;      // n_own x dim x dim
      std::vector<double> own_sum_w;     // n_own
      std::vector<double> own_sum_w_sq;  // n_own

      void reset(int n_covs, int n_vars_in) {
            dim = n_covs + 1;
            n_vars = n_vars_in;
            XtWX.assign(dim * dim, 0.0);
            sum_w = 0.0;
            sum_w_sq = 0.0;
            XtWy.assign(n_vars * dim, 0.0);
            own.assign(n_vars, -1);
            own_XtWX.clear();
            own_sum_w.clear();
            own_sum_w_sq.clear();
      }

      // Give variable v its own (zeroed) system; returns its slot.
      int add_own(int v) {
            const int slot = static_cast<int>(own_sum_w.size());
            own[v] = slot;
            own_XtWX.resize(own_XtWX.size() + dim * dim, 0.0);
            own_sum_w.push_back(0.0);
            own_sum_w_sq.push_back(0.0);
            return slot;
      }

      // Variable v's (unpenalized) system.
      const double* xtwx(int v) const {
            return (own[v] < 0) ? XtWX.data() : own_XtWX.data() + own[v] * dim * dim;
      }
      double sw(int v) const   { return (own[v] < 0) ? sum_w    : own_sum_w[own[v]]; }
      double sw2(int v) const  { return (own[v] < 0) ? sum_w_sq : own_sum_w_sq[own[v]]; }
};

inline void fill_na(double* out, int n) {
      for (int i = 0; i < n; ++i) out[i] = NA_REAL;
}

// Add weight w times the outer product x x' (x = [1, covariates of row j])
// into a dim x dim row-major symmetric matrix.
inline void add_weighted_outer(double* A, int dim, double w,
                               const double* cov_ptr, int cov_stride,
                               std::size_t j) {
      const int n_covs = dim - 1;

      // Intercept x intercept
      A[0] += w;

      // Intercept x covariates (symmetric)
      for (int p = 0; p < n_covs; ++p) {
            double cp = cov_ptr[j + p * cov_stride];
            double wcp = w * cp;
            A[0 * dim + (p + 1)] += wcp;
            A[(p + 1) * dim + 0] += wcp;
      }

      // Covariate x covariate block
      for (int p1 = 0; p1 < n_covs; ++p1) {
            double cp1 = cov_ptr[j + p1 * cov_stride];
            double wcp1 = w * cp1;
            for (int p2 = p1; p2 < n_covs; ++p2) {
                  double cp2 = cov_ptr[j + p2 * cov_stride];
                  double val = wcp1 * cp2;
                  A[(p1 + 1) * dim + (p2 + 1)] += val;
                  if (p2 != p1) {
                        A[(p2 + 1) * dim + (p1 + 1)] += val;
                  }
            }
      }
}

// Accumulate moments over an explicit analog list (see RidgeMoments for the
// treatment of missing responses).
inline void accumulate_moments(
            const std::vector<std::size_t>& analog_indices,  // 0-based indices into pool
            const std::vector<double>& weights,              // combined weight per analog
            const double* values_ptr,                        // pool values (column-major)
            int values_stride,                               // stride for values (n_ref)
            const double* cov_ptr,                           // pool covariates (column-major)
            int cov_stride,                                  // stride for covariates (n_ref)
            bool want_w_sq,                                  // accumulate sum_w_sq?
            RidgeMoments& m                                  // reset() by caller
) {
      const int dim = m.dim;
      const int n_covs = dim - 1;
      const int n_vars = m.n_vars;
      const int n = static_cast<int>(analog_indices.size());

      // Pre-scan: which variables have a missing response in this neighborhood?
      for (int v = 0; v < n_vars; ++v) {
            for (int t = 0; t < n; ++t) {
                  if (ISNAN(values_ptr[analog_indices[t] + v * values_stride])) {
                        m.add_own(v);
                        break;
                  }
            }
      }
      const bool any_own = !m.own_sum_w.empty();

      for (int t = 0; t < n; ++t) {
            const std::size_t j = analog_indices[t];
            const double w = weights[t];

            // Shared system (all analogs)
            m.sum_w += w;
            if (want_w_sq) m.sum_w_sq += w * w;
            add_weighted_outer(m.XtWX.data(), dim, w, cov_ptr, cov_stride, j);

            for (int v = 0; v < n_vars; ++v) {
                  double yv = values_ptr[j + v * values_stride];
                  if (ISNAN(yv)) continue;

                  // Own system for variables with missing responses
                  if (any_own && m.own[v] >= 0) {
                        const int s = m.own[v];
                        m.own_sum_w[s] += w;
                        if (want_w_sq) m.own_sum_w_sq[s] += w * w;
                        add_weighted_outer(m.own_XtWX.data() + s * dim * dim, dim, w,
                                           cov_ptr, cov_stride, j);
                  }

                  // X'Wy
                  double wyv = w * yv;
                  m.XtWy[v * dim + 0] += wyv;
                  for (int p = 0; p < n_covs; ++p) {
                        double cp = cov_ptr[j + p * cov_stride];
                        m.XtWy[v * dim + (p + 1)] += wyv * cp;
                  }
            }
      }
}

// Ridge penalty added to the covariate diagonal (see header comment).
inline double ridge_penalty(const RidgeMoments& m, double lambda, double info) {
      return (info > 0.0) ? lambda * m.sum_w / info : 0.0;
}

// Form A = X'WX + penalty * I_covs from a dim x dim system and Cholesky-factor
// it into L. Returns false if A is not positive definite.
inline bool factor_ridge(const double* XtWX, int dim, double penalty,
                         std::vector<double>& L) {
      L.assign(XtWX, XtWX + dim * dim);
      for (int p = 1; p < dim; ++p) {
            L[p * dim + p] += penalty;
      }
      return cholesky_factor(L, dim);
}

// Solve A beta = X'Wy for variable v, given the factor L of v's system.
inline void solve_coefs(const std::vector<double>& L, const RidgeMoments& m,
                        int v, double* out_beta) {
      const int dim = m.dim;
      std::vector<double> rhs(m.XtWy.begin() + v * dim,
                              m.XtWy.begin() + (v + 1) * dim);
      cholesky_solve(L, dim, rhs);
      for (int d = 0; d < dim; ++d) out_beta[d] = rhs[d];
}

// Kish effective sample size; 0 when undefined.
inline double kish_n_eff(double sum_w, double sum_w_sq) {
      return (sum_w_sq > 0.0) ? (sum_w * sum_w) / sum_w_sq : 0.0;
}

// ESS-framed SEs for one variable (dim values written to out_se_v).
// inv_diag = diag(A^{-1}).
inline void ess_se_from_rss(double rss_w, double sum_w, double n_eff, int dim,
                            const std::vector<double>& inv_diag,
                            double* out_se_v) {
      const double df = n_eff - static_cast<double>(dim);
      if (sum_w <= 0.0 || n_eff <= 0.0 || df <= 0.0) {
            fill_na(out_se_v, dim);
            return;
      }
      const double sigma2 = (rss_w / sum_w) * (n_eff / df);
      if (!(sigma2 >= 0.0) || !std::isfinite(sigma2)) {
            fill_na(out_se_v, dim);
            return;
      }
      for (int d = 0; d < dim; ++d) {
            out_se_v[d] = std::sqrt(sigma2 * inv_diag[d]);
      }
}

// Design-framed (sandwich) SEs for one variable: sqrt(diag(A^{-1} M A^{-1})).
// A^{-1} is symmetric, so diag_d = sum_{a,b} A^{-1}_{da} M_{ab} A^{-1}_{bd}.
inline void sandwich_se_diag(const std::vector<double>& A_inv,
                             const std::vector<double>& M, int dim,
                             double* out_se_v) {
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
            out_se_v[d] = (!(acc >= 0.0) || !std::isfinite(acc))
                  ? NA_REAL : std::sqrt(acc);
      }
}

// Solve every response variable's system and write coefficients; delegate SEs
// to se_fn(v, beta, inv_diag, A_inv, se_v), which receives diag(A^{-1}) when
// se_code == ESS and the full A^{-1} when se_code == DESIGN. Variables without
// missing responses share one (lazily computed) factorization. Variables whose
// system is singular get NA outputs. Returns false if any variable was singular.
template <class SeFn>
inline bool solve_all_vars(const RidgeMoments& m, double penalty, SeCode se_code,
                           double* out_coeffs, double* out_se, SeFn se_fn) {
      const int dim = m.dim;
      const bool want_se = (se_code != SeCode::NONE);

      std::vector<double> L_shared, inv_diag_shared, A_inv_shared;
      bool shared_done = false, shared_ok = false;
      std::vector<double> L_own, inv_diag_own, A_inv_own;
      bool all_ok = true;

      for (int v = 0; v < m.n_vars; ++v) {
            double* beta = out_coeffs + v * dim;
            double* se_v = want_se ? out_se + v * dim : nullptr;

            // Factor (and, for SEs, invert) variable v's system.
            std::vector<double>* L;
            std::vector<double>* inv_diag;
            std::vector<double>* A_inv;
            bool ok;
            bool fresh;
            if (m.own[v] < 0) {
                  fresh = !shared_done;
                  if (fresh) {
                        shared_ok = factor_ridge(m.xtwx(v), dim, penalty, L_shared);
                        shared_done = true;
                  }
                  ok = shared_ok;
                  L = &L_shared; inv_diag = &inv_diag_shared; A_inv = &A_inv_shared;
            } else {
                  fresh = true;
                  ok = factor_ridge(m.xtwx(v), dim, penalty, L_own);
                  L = &L_own; inv_diag = &inv_diag_own; A_inv = &A_inv_own;
            }

            if (!ok) {
                  fill_na(beta, dim);
                  if (want_se) fill_na(se_v, dim);
                  all_ok = false;
                  continue;
            }

            solve_coefs(*L, m, v, beta);
            if (!want_se) continue;

            if (fresh) {
                  if (se_code == SeCode::ESS) {
                        cholesky_inv_diag(*L, dim, *inv_diag);
                  } else {
                        cholesky_inverse(*L, dim, *A_inv);
                  }
            }
            se_fn(v, beta, *inv_diag, *A_inv, se_v);
      }

      return all_ok;
}

// Exact-engine entry point: solve weighted ridge regression for one focal
// from its explicit analog list. Writes out_coeffs (always) and out_se (only
// if se_code != NONE). SEs use a residual pass over the analogs. Variables
// whose system is singular (e.g. no analog observes them) get NA outputs.
// Returns false if any variable was singular.
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
            double info,                                     // information content (see above)
            SeCode se_code,                                  // which SE variant (or NONE)
            double* out_coeffs,                              // output: n_vars * (n_covs + 1)
            double* out_se                                   // output: n_vars * (n_covs + 1); ignored if se_code == NONE
) {
      const int dim = n_covs + 1;
      const int n = static_cast<int>(analog_indices.size());

      RidgeMoments m;
      m.reset(n_covs, n_vars);
      accumulate_moments(analog_indices, weights, values_ptr, values_stride,
                         cov_ptr, cov_stride, se_code != SeCode::NONE, m);

      // Residual pass. ESS needs RSS_w = sum w r^2; DESIGN needs
      // M = sum w^2 r^2 x x'. Analogs missing response v are skipped.
      std::vector<double> M;
      auto residual_se = [&](int v, const double* beta,
                             const std::vector<double>& inv_diag,
                             const std::vector<double>& A_inv, double* se_v) {
            double rss_w = 0.0;
            if (se_code == SeCode::DESIGN) M.assign(dim * dim, 0.0);
            bool any_valid = false;

            for (int t = 0; t < n; ++t) {
                  const std::size_t j = analog_indices[t];
                  const double w = weights[t];
                  double yv = values_ptr[j + v * values_stride];
                  if (ISNAN(yv)) continue;
                  any_valid = true;

                  double fitted = beta[0];
                  for (int p = 0; p < n_covs; ++p) {
                        fitted += beta[p + 1] * cov_ptr[j + p * cov_stride];
                  }
                  const double resid = yv - fitted;

                  if (se_code == SeCode::ESS) {
                        rss_w += w * resid * resid;
                  } else {
                        add_weighted_outer(M.data(), dim, w * w * resid * resid,
                                           cov_ptr, cov_stride, j);
                  }
            }

            if (!any_valid) {
                  fill_na(se_v, dim);
            } else if (se_code == SeCode::ESS) {
                  ess_se_from_rss(rss_w, m.sw(v), kish_n_eff(m.sw(v), m.sw2(v)),
                                  dim, inv_diag, se_v);
            } else {
                  sandwich_se_diag(A_inv, M, dim, se_v);
            }
      };

      return solve_all_vars(m, ridge_penalty(m, lambda, info), se_code,
                            out_coeffs, out_se, residual_se);
}

// Moments-based entry point (binned / FFT engine): solve from moments filled
// directly by the caller. Supports SE codes NONE and ESS only (the design SE
// needs per-analog residuals). For ESS, the weighted RSS is recovered from
// moments as
//   RSS_w = sum_wy2[v] - 2 beta' X'Wy + beta' X'WX beta   (X'WX unpenalized)
// where sum_wy2[v] is sum w y^2 over the analogs observing variable v.
inline bool solve_ridge_from_moments(const RidgeMoments& m, const double* sum_wy2,
                                     double lambda, double info, SeCode se_code,
                                     double* out_coeffs, double* out_se) {
      const int dim = m.dim;
      auto moment_se = [&](int v, const double* beta,
                           const std::vector<double>& inv_diag,
                           const std::vector<double>& /*A_inv*/, double* se_v) {
            const double* A = m.xtwx(v);
            const double* b = m.XtWy.data() + v * dim;
            double bb = 0.0, bAb = 0.0;
            for (int a = 0; a < dim; ++a) {
                  bb += beta[a] * b[a];
                  double row = 0.0;
                  for (int c = 0; c < dim; ++c) row += A[a * dim + c] * beta[c];
                  bAb += beta[a] * row;
            }
            double rss_w = sum_wy2[v] - 2.0 * bb + bAb;
            if (rss_w < 0.0) rss_w = 0.0;  // fp / binning guard
            ess_se_from_rss(rss_w, m.sw(v), kish_n_eff(m.sw(v), m.sw2(v)),
                            dim, inv_diag, se_v);
      };
      return solve_all_vars(m, ridge_penalty(m, lambda, info), se_code,
                            out_coeffs, out_se, moment_se);
}

} // namespace analogs
