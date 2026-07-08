# Cross-validate an analog function

Runs an analog impact or regression analysis in cross-validation mode,
generating held-out predictions and residuals for each location in
`pool`. Each location is predicted using only neighbors that exclude
itself, providing an honest assessment of how well the specified
configuration predicts observed `y` values. Supports leave-one-out (LOO)
and k-fold cross-validation methods.

## Usage

``` r
analog_cv(
  fun,
  pool,
  y,
  covariates = NULL,
  cv_method = c("loo", "kfold"),
  n_folds = NULL,
  fold_id = NULL,
  include_residuals = TRUE,
  ...
)
```

## Arguments

- fun:

  An analog function to cross-validate. Must be one of
  [`analog_impact()`](https://matthewkling.github.io/analogs/reference/analog_impact.md),
  [`analog_regression()`](https://matthewkling.github.io/analogs/reference/analog_regression.md),
  or
  [`analog_search()`](https://matthewkling.github.io/analogs/reference/analog_search.md)
  (passed as a function object, not a string).

- pool:

  The reference dataset. Matrix/data.frame with columns x, y, and
  environmental variables, or a SpatRaster with environmental variable
  layers. Pre-built `analog_index` objects are not supported;
  `analog_cv()` builds indices internally per fold (for k-fold) or once
  (for LOO).

- y:

  Response variable(s). For continuous prediction targets
  (`weighted_mean`, `regression`): numeric vector, matrix, data.frame,
  or SpatRaster. For categorical (`tabulate`): factor or coercible-to-
  factor vector / character / integer / matrix / data.frame /
  SpatRaster. Must have exactly the same number of rows/cells as `pool`.

- covariates:

  Predictor variables (required for regression; must be supplied
  whenever `fun` will fit local regressions). Matrix, data.frame, or
  SpatRaster. Must have exactly the same number of rows/cells as `pool`.

- cv_method:

  One of `"loo"` (default) or `"kfold"`.

- n_folds:

  Integer number of folds for k-fold CV. Pool rows are randomly assigned
  to folds. Ignored when `cv_method = "loo"` or when `fold_id` is
  supplied.

- fold_id:

  Optional integer vector of length `nrow(pool)` giving a fold
  assignment for each pool row. Overrides `n_folds`. Can be used to
  manually specify nonrandom folds, such as for spatial block
  cross-validation.

- include_residuals:

  Logical; if `TRUE` (default), the output includes per-focal
  residual-equivalent columns (see `@return`).

- ...:

  Additional arguments passed to `fun` (e.g., `env`, `geog`, `k`,
  `lambda`, `select`, `se`, `weight`). Note: `fun` must accept
  `exclude_self` (directly or via `...`);
  [`analog_search()`](https://matthewkling.github.io/analogs/reference/analog_search.md)
  accepts it as a named parameter, and the wrapper helpers forward it
  via their own `...`.

## Value

A data.frame or SpatRaster (matching the format of `pool`) with one row
per pool location, containing all variables that `fun` would return,
plus the following residual-equivalent columns when
`include_residuals = TRUE` and a prediction target is identified.

For continuous prediction targets (`weighted_mean`, `regression`):

- `obs` / `obs_{yname}`: observed y value at this location.

- `residual` / `residual_{yname}`: observed minus held-out prediction.

For categorical prediction target (`tabulate`):

- `obs` / `obs_{yname}`: observed class label (character).

- `primary` / `primary_{yname}`: predicted (modal) class label
  (character; argmax across the per-class vote columns).

- `brier` / `brier_{yname}`: per-focal Brier score, computed on
  row-normalized vote shares (range 0-2).

The underlying `n_<level>` vote columns from the analog search are also
retained, so users can compute additional metrics (entropy, top-k
accuracy, custom losses) in postprocessing.

Always present:

- `fold`: fold assignment (k-fold only).

Rows are ordered to match `pool`'s input row order. For SpatRaster
output, character columns (`obs`, `primary`) are dropped with a message;
pass `pool` as a data.frame to retain them.

## Details

`analog_cv()` supports two CV methods:

- `"loo"` (leave-one-out): Each focal location excludes its own pool row
  from its neighborhood. Implemented as a single call with
  self-exclusion. Fast and the most granular form of CV.

- `"kfold"`: Pool is partitioned into `n_folds` folds (or user- supplied
  `fold_id`). Each fold's locations are predicted using the remaining
  folds as the pool. Implemented as `k` separate calls with the index
  rebuilt per fold. Reduces optimism from spatial autocorrelation by
  holding out larger contiguous sets of locations (if folds are
  spatially blocked).

Supported functions:
[`analog_impact()`](https://matthewkling.github.io/analogs/reference/analog_impact.md),
[`analog_regression()`](https://matthewkling.github.io/analogs/reference/analog_regression.md),
and
[`analog_search()`](https://matthewkling.github.io/analogs/reference/analog_search.md).
Other `analog_*()` functions have no `y` input and thus no prediction to
validate.

When `fun = analog_search`, residuals are computed against:

- the `weighted_mean` column if `stat` includes `"weighted_mean"` but
  not `"regression"` or `"tabulate"`;

- fitted values from regression coefficients if `stat` includes
  `"regression"` but not `"weighted_mean"` or `"tabulate"`;

- per-class weighted votes if `stat` includes `"tabulate"` (categorical
  y; produces Brier score and primary-class label rather than a numeric
  residual).

If `stat` includes more than one of these, the prediction target is
ambiguous and `analog_cv()` will error. If it includes none, residuals
are skipped and only the underlying search columns are returned.

For categorical CV (`stat = "tabulate"`), `y` must be a factor or
coercible-to-factor input (character, integer codes, single-layer
categorical SpatRaster). The output uses different residual-equivalent
columns: see `@return` below.

## See also

[`analog_search()`](https://matthewkling.github.io/analogs/reference/analog_search.md),
[`analog_impact()`](https://matthewkling.github.io/analogs/reference/analog_impact.md),
[`analog_regression()`](https://matthewkling.github.io/analogs/reference/analog_regression.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# LOO for AIM
cv <- analog_cv(
  fun      = analog_impact,
  pool     = sites,
  y        = sites$biomass,
  env      = kernel("gaussian", theta = 0.2, max = 0.5),
  geog     = kernel(max = 100)
)
rmse <- sqrt(mean(cv$residual^2, na.rm = TRUE))

# 10-fold CV for local regression
cv_reg <- analog_cv(
  fun         = analog_regression,
  pool        = sites,
  y           = sites$income,
  covariates  = data.frame(education = sites$edu),
  select      = "knn_geog",
  k           = 50,
  geog        = kernel("gaussian", theta = 20),
  cv_method   = "kfold",
  n_folds     = 10
)

# LOO categorical CV (vegetation projection)
cv_veg <- analog_cv(
  fun      = analog_impact,
  pool     = sites,
  y        = factor(sites$vegetation_type),
  stat     = "tabulate",
  env      = kernel("gaussian", theta = 0.2, max = 0.5),
  geog     = kernel(max = 100)
)
# Per-focal Brier score and primary class are in cv_veg$brier and
# cv_veg$primary; the n_<level> vote columns are also retained.
accuracy <- mean(cv_veg$primary == cv_veg$obs, na.rm = TRUE)
mean_brier <- mean(cv_veg$brier, na.rm = TRUE)
} # }
```
