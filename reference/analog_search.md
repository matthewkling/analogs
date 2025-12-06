# Find Climate Analogs

Identifies locations in a reference dataset that are climatically
similar to focal locations, with optional constraints on climate
distance and geographic distance. This function uses a two-stage
approach: first selecting analogs based on specified criteria, then
optionally aggregating the results.

## Usage

``` r
analog_search(
  x,
  pool,
  select = "all",
  stat = NULL,
  max_clim = NULL,
  max_geog = NULL,
  x_cov = NULL,
  k = NULL,
  weight = NULL,
  theta = NULL,
  report_dist = TRUE,
  coord_type = c("auto", "lonlat", "projected"),
  index_res = "auto",
  n_threads = NULL
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
    [`build_analog_index`](https://matthewkling.github.io/analogs/reference/build_analog_index.md)
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

  - A character vector combining multiple stats (e.g.,
    `c("count", "sum_weights")`). Note: `"none"` cannot be combined with
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

- x_cov:

  Optional focal-specific covariance matrices for Mahalanobis distance
  calculations. Should be a matrix or data.frame with one row per focal
  location and one column per unique covariance component. For n climate
  variables, there are n\*(n+1)/2 unique components, ordered as:
  variances first (diagonals), then covariances (upper triangle by row).

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
    exp(-climate_distance^2 / (2\*sigma^2)), with sigma given by
    `theta`.

  - `"gaussian_geog"`: Gaussian kernel on geographic distance, weight =
    exp(-geographic_distance^2 / (2\*sigma^2)), with sigma given by
    `theta`.

  - `"gaussian_joint"`: Joint Gaussian kernel, weight =
    exp(-climate_distance^2/(2\*sigma_c^2) -
    geographic_distance^2/(2\*sigma_g^2)), with sigma values given by
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

- report_dist:

  Logical; if TRUE (default), include distance columns in output when
  `stat` is `NULL` or `"none"`. Set to FALSE for more compact output.

- coord_type:

  Coordinate system type (default: "auto"):

  - `"auto"`: Automatically detect from coordinate ranges.

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

## Value

The return value depends on the `stat` parameter:

\*\*For stat = NULL or "none"\*\*: A data.frame with one row per
focal-analog pair:

- `index`: Index of focal location (1-based).

- `x, y`: Coordinates of focal location.

- `analog_index`: Index of analog location in reference dataset
  (1-based).

- `analog_x, analog_y`: Coordinates of analog location.

- `clim_dist`: Climate distance (if `report_dist = TRUE`).

- `geog_dist`: Geographic distance in km (if `report_dist = TRUE`).

\*\*For stat = single aggregation or vector of aggregations\*\*: A
data.frame with one row per focal location:

- `index`: Index of focal location (1-based).

- `x, y`: Coordinates of focal location.

- One column per requested stat: `count`, `sum_weights`, and/or
  `mean_weights`.

All outputs include diagnostic attributes propagated from the C++ core.

## Details

The function uses a spatial indexing structure (lattice-based) to
quickly search through large reference datasets. Climate similarity is
measured using Euclidean distance in climate space (ideally
pre-whitened; see Details). Geographic distance can be computed for
lon/lat coordinates (great-circle distance) or projected coordinates
(planar distance).

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic pair queries
analog_search(x = focal, pool = ref, select = "all", max_clim = 0.5)
analog_search(x = focal, pool = ref, select = "knn_geog", max_clim = 0.5, k = 1)

# Single aggregation
analog_search(x = focal, pool = ref, select = "all", stat = "count", max_clim = 0.5)

# Multiple aggregations in one pass
analog_search(x = focal, pool = ref, select = "all",
              stat = c("count", "sum_weights", "mean_weights"),
              max_clim = 0.5, weight = "gaussian_clim", theta = 0.1)

# With pre-built index (for repeated queries)
index <- build_analog_index(ref)
analog_search(x = focal1, pool = index, select = "knn_geog", max_clim = 0.5, k = 1)
analog_search(x = focal2, pool = index, select = "all", stat = "count", max_clim = 0.3)
} # }
```
