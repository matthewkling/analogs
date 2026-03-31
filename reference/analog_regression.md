# Local weighted regression across analog neighborhoods

Fits a weighted local regression of `values` on `covariates` within each
focal location's analog neighborhood. Analog neighborhoods are defined
by climatic similarity, geographic proximity, or both, while covariates
capture additional predictors that influence outcomes within each
neighborhood but are not captured by to the search dimensions. This
generalizes the weighted mean — which averages over all
within-neighborhood variation — by resolving variation driven by these
auxiliary predictors. Supports ordinary and ridge-penalized weighted
least squares. With purely geographic neighborhoods, this is equivalent
to geographically weighted regression (GWR); with climate-based
neighborhoods, it extends the analog impact model (AIM) framework to
incorporate local covariate effects.

## Usage

``` r
analog_regression(
  x,
  pool,
  values,
  covariates,
  max_geog = NULL,
  max_clim = NULL,
  select = "all",
  k = NULL,
  weight = c("gaussian_clim", "inverse_clim", "gaussian_geog", "inverse_geog",
    "gaussian_joint", "inverse_joint", "uniform"),
  theta = NULL,
  lambda = 0,
  stat = c("count", "ess", "regression"),
  x_cov = NULL,
  coord_type = "auto",
  index_res = "auto",
  n_threads = NULL,
  progress = FALSE
)
```

## Arguments

- x:

  Focal locations for which regressions will be fit. Should be a
  matrix/data.frame with columns x, y, and climate variables, or a
  SpatRaster with climate variable layers.

- pool:

  The reference dataset to search for analogs. Either:

  - Matrix/data.frame with columns x, y, and climate variables, or
    SpatRaster with climate variable layers, OR

  - An `analog_index` object created by
    [`build_analog_index()`](https://matthewkling.github.io/analogs/reference/build_analog_index.md)
    (for repeated queries).

- values:

  Response variable(s) to model via local regression. Can be a numeric
  vector, matrix, data.frame, or SpatRaster. Must have exactly the same
  number of rows/cells as `pool`. A separate regression is fit for each
  variable.

- covariates:

  Predictor variables for local regression. Can be a numeric vector
  (single covariate), matrix, data.frame, or SpatRaster. Must have
  exactly the same number of rows/cells as `pool`. Column/layer names
  carry through to output. These variables are NOT used for the analog
  search itself – only for regression within each neighborhood.

- max_geog:

  Maximum geographic distance constraint (default: NULL = no geographic
  constraint). When specified, only reference locations within this
  distance are considered. Radius units should be specified in
  kilometers if `coord_type = "lonlat"`, or in projected coordinate
  units if `coord_type = "projected"`.

- max_clim:

  Maximum climate distance constraint (default: NULL = no climate
  constraint). Can be either:

  - A scalar: Euclidean radius in climate space (e.g., 0.5)

  - A vector: Per-variable absolute differences (length must equal
    number of climate variables)

  Only reference locations within this climate distance are considered.
  When `x_cov` is provided, scalar thresholds are interpreted in
  Mahalanobis distance units.

- select:

  Character string specifying the analog selection strategy. One of:

  - `"all"` (default): Select all analogs that satisfy the `max_clim`
    and `max_geog` constraints.

  - `"knn_clim"`: For each focal, select up to `k` analogs with smallest
    climate distance, subject to filters.

  - `"knn_geog"`: For each focal, select up to `k` analogs with smallest
    geographic distance, subject to filters.

- k:

  Number of nearest analogs to return per focal location for kNN
  selection modes. Required when `select` is `"knn_geog"` or
  `"knn_clim"`; must be `NULL` for `select = "all"`.

- weight:

  Weighting function for matches, used only when `stat` includes
  `"sum_weights"` or `"mean_weights"`. One of:

  - `"uniform"`: All matches weighted equally (weight = 1.0).

  - `"inverse_clim"`: Inverse climate distance, weight = 1 /
    (climate_distance + eps), with epsilon given by `theta`.

  - `"inverse_geog"`: Inverse geographic distance, weight = 1 /
    (geographic_distance + eps), with epsilon given by `theta`.

  - `"gaussian_clim"`: Gaussian kernel on climate distance, weight =
    exp(-climate_distance^2 / (2*sigma^2)), with sigma given by `theta`.
    `"gaussian_geog"`: Gaussian kernel on geographic distance, weight =
    exp(-geographic_distance^2 / (2*sigma^2)), with sigma given by
    `theta`.

  - `"gaussian_joint"`: Gaussian kernel on combined distance, weight =
    exp(-(clim_dist^2 / (2*sigma_clim^2) + geog_dist^2 /
    (2*sigma_geog^2))), with sigmas given by `theta`.

  - `"inverse_joint"`: Inverse joint distance, weight = 1 /
    (sqrt(clim_dist^2 + geog_dist^2) + eps), with epsilon given by
    `theta`.

- theta:

  Optional numeric parameter used by weighting functions when `stat`
  includes `"sum_weights"` or `"mean_weights"` and `weight` is not
  `"uniform"`. Interpretation depends on `weight`:

  - For `"inverse_clim"` or `"inverse_geog"`: epsilon value added to
    distances (scalar; default: 1e-12 for climate, 1e-6 for geography).

  - For `"gaussian_clim"` or `"gaussian_geog"`: sigma bandwidth
    parameter (scalar; larger values = slower decay with distance).

  - For `"gaussian_joint"` or `"inverse_joint"`: 2-element vector
    `c(theta_clim, theta_geog)` (defaults: 1 for climate, 1 for
    geography).

- lambda:

  Ridge penalty parameter (default: 0, giving ordinary weighted least
  squares). Higher values shrink covariate coefficients toward zero,
  with the intercept approaching the weighted mean as `lambda -> Inf`.
  Useful when some neighborhoods have few analogs relative to the number
  of covariates.

- stat:

  Statistic(s) to compute. `"regression"` is always included. Additional
  stats like `"count"`, `"ess"`, and `"weighted_mean"` can be requested
  alongside regression coefficients. Default includes `"count"` and
  `"ess"` for diagnostics.

- x_cov:

  Optional focal-specific covariance matrices for Mahalanobis distance
  calculations. Should be a matrix or data.frame with one row per focal
  location and one column per unique covariance component, or a
  SpatRaster with a layer for each component. For n climate variables,
  there are n\*(n+1)/2 unique components, ordered as: variances first
  (diagonals), then covariances (upper triangle by row).

- coord_type:

  Coordinate system type:

  - `"auto"` (default): Automatically detect from coordinate ranges.

  - `"lonlat"`: Unprojected lon/lat coordinates (uses great-circle
    distance; assumes `max_geog` is in km).

  - `"projected"`: Projected XY coordinates (uses planar distance;
    assumes `max_geog` is in projection units).

- index_res:

  Tuning parameter giving the number of bins per dimension of the
  internally-used lattice search index. Either:

  - A positive integer.

  - `"auto"` (the default): Automatically tune the index resolution by
    optimizing compute time on a subsample of focal points. If focal has
    relatively few rows, auto-tuning is skipped and a default resolution
    of 16 is used.

  Ignored if `pool` is an `analog_index` (uses index's resolution).

- n_threads:

  Optional integer number of threads to use for the computation. If
  `NULL` (default), the global RcppParallel setting is used (see
  [`RcppParallel::setThreadOptions`](https://rdrr.io/pkg/RcppParallel/man/setThreadOptions.html)).

- progress:

  Logical; if `TRUE`, display a progress bar during computation.
  Progress tracking works by splitting the focal dataset into chunks and
  processing them sequentially. Useful for large datasets. Default is
  `FALSE`.

## Value

A data.frame (or SpatRaster if `x` is a SpatRaster) with one row per
focal location containing:

- `index`, `x`, `y`: Focal location identifiers

- Columns for any additional stats requested (e.g., `count`, `ess`)

- `intercept`: Regression intercept (predicted value when all covariates
  equal zero)

- One column per covariate: regression slope coefficients

- With multiple `values` variables: columns are named
  `{coeff}_{varname}` (e.g., `intercept_biomass`, `slope_biomass`)

## Details

This function is a wrapper that calls
[`analog_search()`](https://matthewkling.github.io/analogs/reference/analog_search.md)
with `"regression"` included in `stat`.

### Method

For each focal location, the function:

1.  Selects analog pool locations based on `select`, `max_clim`,
    `max_geog`, and `k`

2.  Computes distance-based weights for each analog (via `weight` and
    `theta`)

3.  Fits a weighted least squares regression of `values` on `covariates`
    using these weights, with optional ridge penalty `lambda`

4.  Returns the regression coefficients (intercept + slopes)

The math: `beta = (X'WX + lambda * I_p)^{-1} X'Wy`, where `W` is
diagonal weights, `X` is the design matrix (intercept + covariates), and
`I_p` penalizes covariate coefficients only (not the intercept).

### Relationship to Weighted Mean

With `lambda -> Inf`, covariate coefficients shrink to zero and the
intercept converges to the weighted mean. With `lambda = 0` (the
default), the full local regression is used. If covariates are centered
(weighted mean zero within each neighborhood), the intercept equals the
weighted mean at any lambda. The `lambda` parameter thus provides smooth
interpolation between a simple weighted average and a full local
regression.

### Common Configurations

- **Geographically weighted regression** (`select = "all"` or
  `"knn_geog"`, with `max_geog` and geographic weights): Local spatial
  regression using geographic proximity to define and weight
  neighborhoods, equivalent to GWR.

- **Climate-nearest regression** (`select = "knn_clim"`, with `max_geog`
  and `k`): Fixed-size neighborhoods of the k most similar climates
  within geographic range.

### Prediction

The function returns coefficients only. Prediction at new covariate
values (e.g., a fine-resolution topography grid) can be done via
regression algebra, e.g.:

    prediction <- intercept + beta_1 * covariate_1 + beta_2 * covariate_2

## See also

[`analog_search()`](https://matthewkling.github.io/analogs/reference/analog_search.md)
for the underlying flexible analog search function;
[`analog_impact()`](https://matthewkling.github.io/analogs/reference/analog_impact.md)
for the standard AIM workflow.

## Examples

``` r
if (FALSE) { # \dontrun{
# GWR-style spatial regression
gwr_result <- analog_regression(
  x = sites,
  pool = sites,
  values = sites$income,
  covariates = data.frame(education = sites$edu, access = sites$access),
  select = "knn_geog",
  k = 50,
  max_clim = NULL,
  weight = "gaussian_geog",
  theta = 20
)
} # }
```
