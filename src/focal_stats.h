#pragma once

// Per-focal aggregate statistics: output column layout, accumulators, and
// finalization. Shared by the exact engine (AggWorker, which accumulates over
// individual analogs) and the FFT engine (which fills the same accumulators
// from kernel-convolved fields), so both engines emit identical columns from
// identical formulas.
//
// Column layout of one focal's output row:
//   [regular stats]   COUNT, SUM_WEIGHTS, MEAN_WEIGHTS, ESS, in request order
//   [value stats]     for each y variable: SUM, MEAN, WEIGHTED_SUM,
//                     WEIGHTED_MEAN in request order, with the weighted_mean
//                     SE column directly after WEIGHTED_MEAN when requested
//   [tabulate]        for each y variable: one column per class
//   [regression]      for each y variable: intercept + slopes; then, when SEs
//                     are requested, the same block of SEs

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include "codes.h"
#include "se_code.h"

namespace analogs {

struct StatLayout {
      std::vector<AggregateCode> acodes;
      SeCode se_code;
      int n_vars;
      int n_covs;

      std::vector<int> regular_pos;  // positions in acodes of regular stats
      std::vector<int> value_pos;    // positions in acodes of value stats
      int wm_value_idx = -1;         // index into value_pos of WEIGHTED_MEAN, or -1
      bool has_regression = false;
      bool has_tabulate = false;
      bool need_weights = false;     // does any stat use distance weights?
      int n_regular = 0;
      int n_value = 0;

      bool want_se;
      bool wm_se;                    // weighted_mean SE accumulators needed
      bool wm_se_design;

      std::vector<int> n_classes_per_var;
      std::vector<int> tab_offsets;  // start of variable v's classes in the tabulate block
      int tab_total = 0;

      int reg_dim;                   // n_covs + 1 when regression is requested, else 0
      bool reg_se;

      StatLayout(const std::vector<AggregateCode>& acodes_in, SeCode se_code_in,
                 bool has_values, int n_vars_in, int n_covs_in,
                 const std::vector<int>& n_classes_per_var_in)
            : acodes(acodes_in), se_code(se_code_in), n_vars(n_vars_in),
              n_covs(n_covs_in), n_classes_per_var(n_classes_per_var_in) {

            for (int s = 0; s < static_cast<int>(acodes.size()); ++s) {
                  const AggregateCode a = acodes[s];
                  if (a == AggregateCode::SUM || a == AggregateCode::MEAN ||
                      a == AggregateCode::WEIGHTED_SUM ||
                      a == AggregateCode::WEIGHTED_MEAN) {
                        if (a == AggregateCode::WEIGHTED_MEAN) {
                              wm_value_idx = static_cast<int>(value_pos.size());
                        }
                        value_pos.push_back(s);
                  } else if (a == AggregateCode::REGRESSION) {
                        has_regression = true;
                  } else if (a == AggregateCode::TABULATE) {
                        has_tabulate = true;
                  } else {
                        regular_pos.push_back(s);  // COUNT, SUM_WEIGHTS, MEAN_WEIGHTS, ESS
                  }
                  if (a == AggregateCode::SUM_WEIGHTS ||
                      a == AggregateCode::MEAN_WEIGHTS ||
                      a == AggregateCode::WEIGHTED_SUM ||
                      a == AggregateCode::WEIGHTED_MEAN ||
                      a == AggregateCode::ESS ||
                      a == AggregateCode::REGRESSION ||
                      a == AggregateCode::TABULATE) {
                        need_weights = true;
                  }
            }
            n_regular = static_cast<int>(regular_pos.size());
            n_value = static_cast<int>(value_pos.size());

            want_se = (se_code != SeCode::NONE);
            wm_se = want_se && (wm_value_idx >= 0) && has_values;
            wm_se_design = wm_se && (se_code == SeCode::DESIGN);

            if (has_tabulate) {
                  tab_offsets.assign(n_vars, 0);
                  for (int v = 0; v < n_vars; ++v) {
                        tab_offsets[v] = tab_total;
                        if (v < static_cast<int>(n_classes_per_var.size())) {
                              tab_total += n_classes_per_var[v];
                        }
                  }
            }

            reg_dim = has_regression ? (n_covs + 1) : 0;
            reg_se = want_se && has_regression;
      }

      // Total output columns per focal.
      int n_cols() const {
            const int wm_se_cols = (want_se && wm_value_idx >= 0) ? n_vars : 0;
            const int reg_cols = has_regression
                  ? n_vars * (reg_dim + (want_se ? reg_dim : 0)) : 0;
            return n_regular + n_value * n_vars + wm_se_cols + tab_total + reg_cols;
      }
};

// Per-focal sufficient statistics for every non-regression stat.
struct FocalAccum {
      double count = 0.0;                 // integer-valued in the exact engine
      double sum_weights = 0.0;
      double sum_weights_squared = 0.0;
      std::vector<double> regular_accum;  // n_regular (SUM_WEIGHTS / MEAN_WEIGHTS slots)
      std::vector<double> value_accum;    // n_vars x n_value, index v * n_value + idx
      std::vector<double> tabulate_accum; // tab_total
      std::vector<double> wm_sum_wy2;     // n_vars (ESS and DESIGN SE)
      std::vector<double> wm_sum_w2y;     // n_vars (DESIGN SE)
      std::vector<double> wm_sum_w2y2;    // n_vars (DESIGN SE)
      double wm_sum_w2 = 0.0;             // DESIGN SE

      void reset(const StatLayout& L, bool has_values) {
            count = 0.0;
            sum_weights = 0.0;
            sum_weights_squared = 0.0;
            regular_accum.assign(L.n_regular, 0.0);
            value_accum.clear();
            if (has_values && L.n_value > 0) {
                  value_accum.assign(L.n_vars * L.n_value, 0.0);
            }
            tabulate_accum.clear();
            if (L.has_tabulate && L.tab_total > 0) {
                  tabulate_accum.assign(L.tab_total, 0.0);
            }
            wm_sum_wy2.clear();
            wm_sum_w2y.clear();
            wm_sum_w2y2.clear();
            wm_sum_w2 = 0.0;
            if (L.wm_se) {
                  wm_sum_wy2.assign(L.n_vars, 0.0);
                  if (L.wm_se_design) {
                        wm_sum_w2y.assign(L.n_vars, 0.0);
                        wm_sum_w2y2.assign(L.n_vars, 0.0);
                  }
            }
      }
};

// Write the regular, value, and tabulate columns for one focal into `row`.
// Returns the index of the next column (where regression output begins).
inline int write_nonregression_stats(const StatLayout& L, const FocalAccum& a,
                                     bool has_values, double* row) {
      int col_idx = 0;

      // Regular stats
      for (int idx = 0; idx < L.n_regular; ++idx) {
            const AggregateCode code = L.acodes[L.regular_pos[idx]];
            double result = NA_REAL;

            if (code == AggregateCode::COUNT) {
                  result = a.count;
            } else if (code == AggregateCode::SUM_WEIGHTS) {
                  result = a.regular_accum[idx];
            } else if (code == AggregateCode::MEAN_WEIGHTS) {
                  result = (a.count > 0) ? (a.regular_accum[idx] / a.count) : NA_REAL;
            } else if (code == AggregateCode::ESS) {
                  if (a.sum_weights_squared > 0.0 && a.sum_weights > 0.0) {
                        result = (a.sum_weights * a.sum_weights) / a.sum_weights_squared;
                  } else {
                        result = 0.0;
                  }
            }

            row[col_idx++] = result;
      }

      // Value stats, with the weighted_mean SE column directly after
      // WEIGHTED_MEAN when requested.
      if (has_values && L.n_value > 0) {
            for (int v = 0; v < L.n_vars; ++v) {
                  for (int idx = 0; idx < L.n_value; ++idx) {
                        const AggregateCode code = L.acodes[L.value_pos[idx]];
                        const int accum_idx = v * L.n_value + idx;
                        double result = NA_REAL;

                        if (code == AggregateCode::SUM) {
                              result = a.value_accum[accum_idx];
                        } else if (code == AggregateCode::MEAN) {
                              result = (a.count > 0) ? (a.value_accum[accum_idx] / a.count) : NA_REAL;
                        } else if (code == AggregateCode::WEIGHTED_SUM) {
                              result = a.value_accum[accum_idx];
                        } else if (code == AggregateCode::WEIGHTED_MEAN) {
                              result = (a.sum_weights > 0) ? (a.value_accum[accum_idx] / a.sum_weights) : NA_REAL;
                        }

                        row[col_idx++] = result;

                        if (L.wm_se && idx == L.wm_value_idx) {
                              double se_val = NA_REAL;
                              if (a.sum_weights > 0.0 && !ISNAN(result)) {
                                    const double ybar = result;
                                    const double sum_wy2 = a.wm_sum_wy2[v];

                                    if (L.se_code == SeCode::ESS) {
                                          // n_eff = (Σw)^2 / Σw^2;  var_w(y) = Σwy^2/Σw - ybar^2
                                          // SE = sqrt(var_w / n_eff)
                                          if (a.sum_weights_squared > 0.0) {
                                                const double n_eff =
                                                      (a.sum_weights * a.sum_weights) / a.sum_weights_squared;
                                                double var_w = (sum_wy2 / a.sum_weights) - (ybar * ybar);
                                                if (var_w < 0.0) var_w = 0.0;  // fp guard
                                                if (n_eff > 0.0 && std::isfinite(var_w)) {
                                                      se_val = std::sqrt(var_w / n_eff);
                                                }
                                          }
                                    } else {
                                          // DESIGN: SE = sqrt(Σ w^2 (y - ybar)^2) / Σw, with
                                          //   numerator = Σw^2 y^2 - 2 ybar Σw^2 y + ybar^2 Σw^2
                                          double num = a.wm_sum_w2y2[v]
                                                - 2.0 * ybar * a.wm_sum_w2y[v]
                                                + ybar * ybar * a.wm_sum_w2;
                                          if (num < 0.0) num = 0.0;  // fp guard
                                          if (std::isfinite(num)) {
                                                se_val = std::sqrt(num) / a.sum_weights;
                                          }
                                    }
                              }
                              row[col_idx++] = se_val;
                        }
                  }
            }
      }

      // Tabulate: per-class weighted vote totals, variable by variable.
      if (L.has_tabulate && L.tab_total > 0) {
            for (int b = 0; b < L.tab_total; ++b) {
                  row[col_idx++] = a.tabulate_accum[b];
            }
      }

      return col_idx;
}

// Write regression output starting at row[col_idx]: coefficients for every
// variable, then SEs when requested. coeffs/se are n_vars x reg_dim; pass
// nullptr for coeffs to write all-NA output (e.g. a focal with no analogs).
inline void write_regression_stats(const StatLayout& L, const double* coeffs,
                                   const double* se, double* row, int col_idx) {
      const int n = L.n_vars * L.reg_dim;
      for (int t = 0; t < n; ++t) {
            row[col_idx++] = coeffs ? coeffs[t] : NA_REAL;
      }
      if (L.reg_se) {
            for (int t = 0; t < n; ++t) {
                  row[col_idx++] = (coeffs && se) ? se[t] : NA_REAL;
            }
      }
}

} // namespace analogs
