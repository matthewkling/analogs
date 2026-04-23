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
  climate variables, or a SpatRaster with climate variable layers.
  Pre-built `analog_index` objects are not supported; `analog_cv()`
  builds indices internally per fold (for k-fold) or once (for LOO).

- y:

  Response variable(s). Numeric vector, matrix, data.frame, or
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

  Logical; if `TRUE` (default), the output includes `obs[_{yname}]` and
  `residual[_{yname}]` columns for each `y` variable (when a prediction
  target can be identified).

- ...:

  Additional arguments passed to `fun` (e.g., `max_clim`, `max_geog`,
  `kernel`, `theta`, `k`, `lambda`, `select`, `se`). Note: `fun` must
  accept `exclude_self` (directly or via `...`);
  [`analog_search()`](https://matthewkling.github.io/analogs/reference/analog_search.md)
  accepts it as a named parameter, and the wrapper helpers forward it
  via their own `...`.

## Value

A data.frame or SpatRaster (matching the format of `pool`) with one row
per pool location, containing all variables that `fun` would return,
plus:

- `obs` / `obs_{yname}`: observed y value at this location (when
  `include_residuals = TRUE` and a prediction target is identified).

- `residual` / `residual_{yname}`: observed minus held-out prediction.

- `fold`: fold assignment (k-fold only).

Rows are ordered to match `pool`'s input row order.

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
  not `"regression"`;

- fitted values from regression coefficients if `stat` includes
  `"regression"` but not `"weighted_mean"`;

If `stat` includes both, the prediction target is ambiguous and
`analog_cv()` will error. If it includes neither, residuals are skipped
and only the underlying search columns are returned.

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
  max_clim = 0.5,
  max_geog = 100,
  kernel   = "gaussian_clim",
  theta    = 0.2
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
  kernel      = "gaussian_geog",
  theta       = 20,
  cv_method   = "kfold",
  n_folds     = 10
)

# Power-user CV via analog_search with a custom stat
cv_custom <- analog_cv(
  fun      = analog_search,
  pool     = sites,
  y        = sites$biomass,
  select   = "knn_clim",
  k        = 10,
  stat     = c("count", "weighted_mean"),
  kernel   = "gaussian_clim",
  theta    = 0.3
)
} # }
```
