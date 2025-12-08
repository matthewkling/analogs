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
  select = "all",
  stat = NULL,
  max_clim = NULL,
  max_geog = NULL,
  k = NULL,
  weight = NULL,
  theta = NULL,
  x_cov = NULL,
  values = NULL,
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

  - An `analog_index()` object created by
    [`build_analog_index()`](https://matthewkling.github.io/analogs/reference/build_analog_index.md)
    (for repeated queries).

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
    (see `weight` and `theta`).

  - `"mean_weights"`: For each focal, mean of weights of selected
    analogs.

  - `"sum"`: Sum of values across analogs (requires `values`).

  - `"mean"`: Mean of values across analogs (requires `values`).

  - `"weighted_sum"`: Sum of (value × weight) across analogs (requires
    `values` and `weight`).

  - `"weighted_mean"`: Weighted mean of values across analogs (requires
    `values` and `weight`).

  - A character vector combining multiple stats (e.g.,
    `c("count", "sum", "mean")`). Note: `"none"` cannot be combined with
    other stats.

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

  - `"gaussian_joint"`: Joint Gaussian kernel, weight =
    exp(-climate_distance^2/(2*sigma_c^2) -
    geographic_distance^2/(2*sigma_g^2)), with sigma values given by
    `theta` as a 2-element vector c(sigma_clim, sigma_geog).

  - `"inverse_joint"`: Joint inverse distance, weight = 1 /
    sqrt((climate_distance + eps_clim)^2 + (geographic_distance +
    eps_geog)^2), with epsilon values given by `theta` as a 2-element
    vector c(eps_clim, eps_geog).

  For `stat` not including weighted aggregations, `weight` must be
  `NULL`.

- theta:

  Optional numeric parameter used by weighting functions when `stat`
  includes `"sum_weights"` or `"mean_weights"` and `weight` is not
  `"uniform"`. Interpretation depends on `weight`:

  - For `"inverse_clim"` or `"inverse_geog"`: epsilon value added to
    distances (scalar; default: 1e-12 for climate, 1e-6 for geography).

  - For `"gaussian_clim"` or `"gaussian_geog"`: sigma bandwidth
    parameter (scalar; larger values = slower decay with distance).

  - For `"gaussian_joint"`: 2-element vector c(sigma_clim, sigma_geog)
    giving bandwidths for climate and geographic dimensions.

  - For `"inverse_joint"`: 2-element vector c(eps_clim, eps_geog) giving
    epsilon values for climate and geographic dimensions.

  If `theta` is `NULL`, sensible defaults are used for single-parameter
  weights. For `weight = "uniform"` or for non-weighted aggregations,
  `theta` must be `NULL`.

- x_cov:

  Optional focal-specific covariance matrices for Mahalanobis distance
  calculations. Should be a matrix or data.frame with one row per focal
  location and one column per unique covariance component. For n climate
  variables, there are n\*(n+1)/2 unique components, ordered as:
  variances first (diagonals), then covariances (upper triangle by row).

- values:

  Optional user-defined variables for each reference location to
  aggregate across selected analogs. Can be:

  - A numeric vector (single variable)

  - A matrix or data.frame with numeric columns (multiple variables)

  Must have exactly `nrow(pool)` rows (or number of reference locations
  if pool is an index). Each row corresponds to a reference location.

  When provided, enables value-based aggregation stats:

  - `"sum"`: Sum of values across analogs

  - `"mean"`: Mean of values across analogs

  - `"weighted_sum"`: Sum of (value × weight) - requires `weight`

  - `"weighted_mean"`: Sum of (value × weight) / sum of weights -
    requires `weight`

  For stat = NULL/"none" (pairs mode), value columns are included in
  output for each analog pair.

- coord_type:

  Coordinate system type (default: "auto"):

  - `"auto"`: Automatically detect from coordinate ranges.

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
