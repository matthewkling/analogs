# Tune Index Resolution

Automatically finds the optimal lattice index resolution for your data
and query pattern using adaptive bracketing search. Runs test queries
with different resolutions and recommends the one with the fastest
compute speed.

## Usage

``` r
tune_index_res(
  x,
  pool,
  downsample = 1,
  seed = NULL,
  select = "all",
  stat = NULL,
  max_clim = NULL,
  max_geog = NULL,
  k = NULL,
  kernel = NULL,
  theta = NULL,
  x_cov = NULL,
  y = NULL,
  covariates = NULL,
  lambda = 0,
  se = c("none", "ess", "design"),
  coord_type = c("auto", "lonlat", "projected"),
  n_threads = NULL,
  default_res = 16L,
  verbose = FALSE
)
```

## Arguments

- x:

  Focal locations for which analogs will be found. Should be a
  matrix/data.frame with columns x, y, and climate variables, or a
  SpatRaster with climate variable layers.

- pool:

  The reference dataset to search for analogs. Either:

  - Matrix/data.frame with columns x, y, and climate variables, or
    SpatRaster with climate variable layers, OR

  - An `analog_index` object created by
    [`build_analog_index()`](https://matthewkling.github.io/analogs/reference/build_analog_index.md)
    (for repeated queries).

- downsample:

  Optional downsampling rate (0-1) for the reference pool, indicating
  the proportion of points to retain. Values \< 1 reduce memory and
  improve speed at some cost to precision. Default is 1.0 (no
  downsampling). Ignored if `pool` is a pre-built index. When
  `downsample < 1`, `index_res` must be set explicitly (auto-tuning is
  not supported in this case; see the `index_res` parameter for
  details).

- seed:

  Optional random seed for reproducible downsampling. If `NULL`
  (default), uses current R random state. Ignored if `pool` is a
  pre-built index or `downsample = 1`.

- select:

  Character string specifying the analog selection strategy. One of:

  - `"all"` (default): Select all analogs that satisfy the `max_clim`
    and `max_geog` constraints.

  - `"knn_clim"`: For each focal, select up to `k` analogs with smallest
    climate distance, subject to filters.

  - `"knn_geog"`: For each focal, select up to `k` analogs with smallest
    geographic distance, subject to filters.

- stat:

  Statistic(s) used to aggregate selected analogs. Either:

  - `NULL` or `"none"`: Return all selected analog pairs as a
    data.frame.

  - `"count"`: For each focal, count the number of selected analogs.

  - `"sum_weights"`: For each focal, sum the weights of selected analogs
    (see `kernel` and `theta`). When `normalize = TRUE`, the reported
    value is the normalized density `D / D_max`, on roughly `[0, 1]`;
    otherwise it is the raw kernel-weight sum.

  - `"mean_weights"`: For each focal, mean of weights of selected
    analogs.

  - `"sum"`: Sum of `y` values across analogs (requires `y`).

  - `"mean"`: Mean of `y` values across analogs (requires `y`).

  - `"weighted_sum"`: Sum of (`y` × kernel weight) across analogs
    (requires `y` and `kernel`).

  - `"weighted_mean"`: Weighted mean of `y` values across analogs
    (requires `y` and `kernel`).

  - `"ess"`: Kish's effective sample size (ESS), computed as the squared
    sum of weights divided by the sum of squared weights (requires
    `kernel`).

  - `"regression"`: Weighted least squares (or ridge) regression of `y`
    on `covariates` within each analog neighborhood. Returns intercept
    and slope coefficients. Requires `y`, `covariates`, and `kernel`.
    See `lambda` for regularization.

  - `"tabulate"`: if `y` is categorical, separately sum the kernel
    weights of analogs matching each level of `y`. With
    `kernel = "uniform"` this reduces to a per-class vote count; with a
    distance-decay kernel it gives similarity-weighted support per
    class. Requires `y` (factor or coercible-to-factor) and `kernel`.
    Output has one column per class. `"tabulate"` is mutually exclusive
    with `"sum"`, `"mean"`, `"weighted_sum"`, `"weighted_mean"`, and
    `"regression"` (different `y` semantics); it can be combined with
    `"count"`, `"sum_weights"`, `"mean_weights"`, and `"ess"`.

  - A character vector combining multiple stats (e.g.,
    `c("count", "weighted_mean", "regression")`). Note: `"none"` cannot
    be combined with other stats.

- max_clim:

  Maximum climate distance constraint (default: NULL = no climate
  constraint). Can be either:

  - A scalar: Euclidean radius in climate space (e.g., 0.5)

  - A vector: Per-variable absolute differences (length must equal
    number of climate variables)

  Only reference locations within this climate distance are considered.
  When `x_cov` is provided, scalar thresholds are interpreted in
  Mahalanobis distance units.

- max_geog:

  Maximum geographic distance constraint (default: NULL = no geographic
  constraint). When specified, only reference locations within this
  distance are considered. Radius units should be specified in
  kilometers if `coord_type = "lonlat"`, or in projected coordinate
  units if `coord_type = "projected"`.

- k:

  Number of nearest analogs to return per focal location for kNN
  selection modes. Required when `select` is `"knn_geog"` or
  `"knn_clim"`; must be `NULL` for `select = "all"`.

- kernel:

  Kernel decay function for weighting matches, used only when `stat`
  includes a weighted aggregation (`"sum_weights"`, `"mean_weights"`,
  `"weighted_sum"`, `"weighted_mean"`, `"ess"`, `"regression"`, or
  `"tabulate"`). One of:

  - `"uniform"`: All matches weighted equally (kernel weight = 1.0).

  - `"inverse_clim"`: Inverse climate distance, kernel weight = 1 /
    (climate_distance + eps), with epsilon given by `theta`.

  - `"inverse_geog"`: Inverse geographic distance, kernel weight = 1 /
    (geographic_distance + eps), with epsilon given by `theta`.

  - `"gaussian_clim"`: Gaussian kernel on climate distance, kernel
    weight = exp(-climate_distance^2 / (2 sigma^2)), with sigma given by
    `theta`.

  - `"gaussian_geog"`: Gaussian kernel on geographic distance, kernel
    weight = exp(-geographic_distance^2 / (2 sigma^2)), with sigma given
    by `theta`.

  - `"gaussian_joint"`: Gaussian kernel on combined distance, kernel
    weight = exp(-(clim_dist^2 / (2 sigma_clim^2) + geog_dist^2 / (2
    sigma_geog^2))), with sigmas given by `theta`.

  - `"inverse_joint"`: Inverse joint distance, kernel weight = 1 /
    (sqrt(clim_dist^2 + geog_dist^2) + eps), with epsilon given by
    `theta`.

- theta:

  Optional numeric parameter controlling the shape of the weighting
  `kernel`, used whenever `kernel` is active (i.e. whenever `stat`
  includes a weighted aggregation) and `kernel` is not `"uniform"`.
  Interpretation depends on `kernel`:

  - For `"inverse_clim"` or `"inverse_geog"`: epsilon value added to
    distances (scalar; default: 1e-12 for climate, 1e-6 for geography).

  - For `"gaussian_clim"` or `"gaussian_geog"`: sigma bandwidth
    parameter (scalar; larger values = slower decay with distance).

  - For `"gaussian_joint"` or `"inverse_joint"`: 2-element vector
    `c(theta_clim, theta_geog)` (defaults: 1 for climate, 1 for
    geography).

  See
  [`kernel_params()`](https://matthewkling.github.io/analogs/reference/kernel_params.md)
  for help choosing `theta` and `max_clim` / `max_geog` values that work
  well together.

- x_cov:

  Optional focal-specific covariance matrices for Mahalanobis distance
  calculations. Should be a matrix or data.frame with one row per focal
  location and one column per unique covariance component, or a
  SpatRaster with a layer for each component. For n climate variables,
  there are n\*(n+1)/2 unique components, ordered as: variances first
  (diagonals), then covariances (upper triangle by row).

- y:

  Optional vector, factor, matrix/data.frame, or SpatRaster giving
  values for each reference location (must have same number of
  rows/cells as `pool`). Required for stats `"sum"`, `"mean"`,
  `"weighted_sum"`, `"weighted_mean"`, `"regression"`, and `"tabulate"`.
  Numeric for continuous stats; factor or coercible-to-factor
  (character, integer, logical) for `stat = "tabulate"`.

- covariates:

  Optional matrix/data.frame or SpatRaster giving covariate values for
  each reference location (must have same number of rows/cells as
  `pool`). Required when `stat` includes `"regression"`.

- lambda:

  Ridge penalty parameter for `stat = "regression"` (default: 0, giving
  ordinary weighted least squares). Higher values shrink covariate
  coefficients toward zero, with the intercept approaching the weighted
  mean as `lambda -> Inf`. Ignored when `"regression"` is not in `stat`.

- se:

  Standard-error framing to apply to SE-supporting stats
  (`"weighted_mean"` and `"regression"`). One of:

  - `"none"` (default): no SE columns are returned.

  - `"ess"`: effective-sample-size framing. For `weighted_mean`,
    `SE = sqrt(var_w(y) / n_eff)`, where `n_eff = (Σw)² / Σw²` is Kish's
    effective sample size and `var_w(y) = Σwy²/Σw - ȳ_w²`. For
    regression, `Var(β̂) = σ²_ess · (X'WX + λI)⁻¹`, with residual
    variance corrected using `n_eff - p` degrees of freedom.

  - `"design"`: design-based framing (no assumption that weights are
    precisions). For `weighted_mean`, `SE = sqrt(Σ w²(y - ȳ_w)²) / Σw`.

- coord_type:

  Coordinate system type:

  - `"auto"` (default): Automatically detect from coordinate ranges.

  - `"lonlat"`: Unprojected lon/lat coordinates (uses great-circle
    distance; assumes `max_geog` is in km).

  - `"projected"`: Projected XY coordinates (uses planar distance;
    assumes `max_geog` is in projection units).

- n_threads:

  Optional integer number of threads to use for the computation. If
  `NULL` (default), the global RcppParallel setting is used (see
  [`RcppParallel::setThreadOptions`](https://rdrr.io/pkg/RcppParallel/man/setThreadOptions.html)).

- default_res:

  Default resolution to use as starting point for search. Default is 16.

- verbose:

  Logical; if TRUE, print the selected resolution. Default is FALSE.

## Value

An integer giving the recommended index resolution (bins per dimension).

## Details

The function uses an adaptive bracketing algorithm:

1.  Starts with three resolutions: default/2, default, default\*2

2.  Evaluates elapsed time for each

3.  If minimum is at an edge, expands search in that direction

4.  Returns resolution with lowest elapsed time

This typically requires only 3-5 query evaluations total, making it much
faster than exhaustive grid search.

The function only performs tuning for non-trivial problem sizes (\>2000
focal points). For smaller datasets, it returns the default resolution.

A subsample of focal points is used for benchmarking to keep tuning fast
while still being representative of actual query performance.

## Examples

``` r
if (FALSE) { # \dontrun{
# Find optimal resolution for velocity queries
optimal_res <- tune_index_res(
  x = sample_sites,
  pool = climate_data,
  select = "knn_geog",
  stat = NULL,
  max_clim = 0.5,
  k = 1
)

# Use the optimized resolution
index <- build_analog_index(climate_data, index_res = optimal_res)
} # }
```
