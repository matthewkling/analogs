# Performance metrics for cross-validation output

Computes per-variable prediction-performance metrics from the output of
[`analog_cv()`](https://matthewkling.github.io/analogs/reference/analog_cv.md).
Handles both tabular (data.frame) and raster (SpatRaster) CV output,
both single-y and multi-y configurations, and continuous, binary, and
categorical outcomes.

## Usage

``` r
cv_performance(x, outcome_type = "auto", weights = NULL)
```

## Arguments

- x:

  Output from
  [`analog_cv()`](https://matthewkling.github.io/analogs/reference/analog_cv.md).
  Either a data.frame or a SpatRaster produced with
  `include_residuals = TRUE`.

- outcome_type:

  Controls outcome-type classification (continuous / binary only):

  - `"auto"` (default): auto-detect as described above.

  - `"continuous"` or `"binary"`: force all variables to this type.

  - A named character vector with one entry per variable, giving the
    type for each (e.g.,
    `c(biomass = "continuous", presence = "binary")`).

  Categorical CV results (from `stat = "tabulate"`) ignore this argument
  unless explicitly set, in which case they error.

- weights:

  Optional numeric vector of per-location weights (one value per
  row/cell of `x`). Used for weighted versions of all metrics. Default
  `NULL` gives unweighted metrics.

## Value

A data.frame in long format with columns `variable`, `type`, `metric`,
and `value`. One row per (variable, metric) pair.

## Output format

Results are returned in long format with one row per (variable, metric)
pair and columns:

- `variable`: name of the response variable.

- `type`: outcome type (`"continuous"`, `"binary"`, or `"categorical"`).

- `metric`: metric name (see below).

- `value`: numeric metric value.

This format is designed to accommodate additional outcome types and
metrics in the future without schema changes.

## Continuous metrics

For variables classified as continuous:

- `n`: number of locations with finite observed and predicted values.

- `rmse`: root mean squared error of held-out predictions.

- `mae`: mean absolute error.

- `bias`: mean signed residual (`mean(obs - pred)`). Positive values
  indicate systematic under-prediction.

- `r2`: out-of-sample R² (a.k.a. "predicted R²"), computed as
  `1 - SS_res / SS_tot` using the held-out residuals. Can be negative
  when predictions are worse than simply predicting the overall mean.

## Binary metrics

For variables classified as binary (observed values all in `[0, 1]` with
both classes present):

- `n`: number of locations with finite observed and predicted values.

- `auc`: area under the ROC curve (via Mann-Whitney U; handles ties). A
  threshold-independent measure of rank-based discrimination.

- `tss`: true skill statistic (`sensitivity + specificity - 1`) at the
  threshold that maximizes it over unique predicted values.

- `tss_threshold`: the threshold used for `tss`.

- `brier`: Brier score (mean squared error between prediction and 0/1
  outcome). A proper scoring rule; most interpretable when predictions
  are in `[0, 1]`.

## Categorical metrics

For results from
[`analog_cv()`](https://matthewkling.github.io/analogs/reference/analog_cv.md)
with `stat = "tabulate"`:

- `n`: number of locations with non-NA observed class and a non-empty
  analog neighborhood (i.e., at least one analog with a non-NA class).

- `accuracy`: proportion of locations where the predicted (primary)
  class matches the observed class.

- `brier`: mean per-focal multiclass Brier score, computed on
  row-normalized vote shares (`Σ_k (p_k - I[obs == k])^2`). Range
  `[0, 2]`; lower is better.

- `n_classes`: number of distinct classes (`K`).

- `confusion[<obs>|<pred>]`: count of locations with observed class
  `<obs>` predicted as `<pred>`. One row per `K * K` cell, with zero
  counts included. Filter via `startsWith(metric, "confusion[")` to pull
  just the confusion matrix entries.

## Outcome-type detection

The CV result's column structure determines its overall type:

- Categorical results (from `stat = "tabulate"`) have `obs_*`,
  `primary_*`, and `brier_*` columns. All variables are categorical.

- Continuous/binary results have `obs_*` and `residual_*` columns. When
  `outcome_type = "auto"` (the default), each variable is classified
  per-variable as `"binary"` if observed values are all in `[0, 1]`
  after removing NAs (with both classes present), or `"continuous"`
  otherwise.

Users can override classification for continuous/binary cases by passing
a scalar type name (applies to all variables) or a named character
vector (one entry per variable). `outcome_type` does not apply to
categorical CV output and passing anything other than `"auto"` for a
categorical result will error.

## See also

[`analog_cv()`](https://matthewkling.github.io/analogs/reference/analog_cv.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Continuous outcome
cv <- analog_cv(
  fun      = analog_impact,
  pool     = sites,
  y        = sites$biomass,
  env     = kernel("gaussian", theta = 0.2, max = 0.5),
  geog     = kernel(max = 100)
)
cv_performance(cv)

# Binary outcome (presence/absence)
cv_bin <- analog_cv(
  fun      = analog_impact,
  pool     = sites,
  y        = sites$presence,   # 0/1 values
  env     = kernel("gaussian", theta = 0.2, max = 0.5),
  geog     = kernel(max = 100)
)
cv_performance(cv_bin)

# Categorical outcome (e.g., vegetation type)
cv_cat <- analog_cv(
  fun      = analog_impact,
  pool     = sites,
  y        = factor(sites$vegetation),
  stat     = c("count", "sum_weights", "tabulate"),
  env     = kernel("gaussian", theta = 0.2, max = 0.5),
  geog     = kernel(max = 100)
)
perf <- cv_performance(cv_cat)

# Pull just the headline scalars
perf[!startsWith(perf$metric, "confusion["), ]

# Pull the confusion matrix in long form
perf[startsWith(perf$metric, "confusion["), ]

# Parameter tuning via AUC
thetas <- c(0.1, 0.2, 0.3, 0.5)
auc <- sapply(thetas, function(th) {
  cv <- analog_cv(
    fun = analog_impact, pool = sites, y = sites$presence,
    env = kernel("gaussian", theta = th, max = 0.5), geog = kernel(max = 100)
  )
  perf <- cv_performance(cv)
  perf$value[perf$metric == "auc"]
})
thetas[which.max(auc)]
} # }
```
