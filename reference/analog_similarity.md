# Analog similarity: best climate analogs within a geographic envelope

Finds, for each focal location, the climate–nearest neighbor(s) in a
reference dataset that satisfy a specified geographic distance
threshold. Among other uses, this operation is often the first step in a
traditional analog impact modeling (AIM) analysis.

## Usage

``` r
analog_similarity(
  x,
  pool,
  x_cov = NULL,
  y = NULL,
  coord_type = "auto",
  max_geog,
  max_clim = NULL,
  k = 20,
  index_res = "auto",
  n_threads = NULL,
  downsample = 1,
  seed = NULL,
  progress = FALSE
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

- x_cov:

  Optional focal-specific covariance matrices for Mahalanobis distance
  calculations. Should be a matrix or data.frame with one row per focal
  location and one column per unique covariance component, or a
  SpatRaster with a layer for each component. For n climate variables,
  there are n\*(n+1)/2 unique components, ordered as: variances first
  (diagonals), then covariances (upper triangle by row).

- y:

  Optional user-defined variables for each reference location in `pool`
  to aggregate across selected analogs. Can be a numeric vector (single
  variable), matrix or data.frame with numeric columns (multiple
  variables), or a SpatRaster with one or more numeric layers. Must have
  exactly the same number of reference locations as `pool`.

  When provided, enables value-based aggregation stats `"sum"`,
  `"mean"`, `"weighted_sum"`, `"weighted_mean"`, and `"regression"`. For
  stat = NULL/"none" (pairs mode), y value columns are included in
  output for each analog pair.

- coord_type:

  Coordinate system type:

  - `"auto"` (default): Automatically detect from coordinate ranges.

  - `"lonlat"`: Unprojected lon/lat coordinates (uses great-circle
    distance; assumes `max_geog` is in km).

  - `"projected"`: Projected XY coordinates (uses planar distance;
    assumes `max_geog` is in projection units).

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

- k:

  Number of nearest analogs to return per focal location for kNN
  selection modes. Required when `select` is `"knn_geog"` or
  `"knn_clim"`; must be `NULL` for `select = "all"`.

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

- downsample:

  Optional downsampling rate (0-1) for the reference pool, indicating
  the proportion of points to retain. Values \< 1 reduce memory and
  improve speed at some cost to precision. Default is 1.0 (no
  downsampling). Ignored if `pool` is a pre-built index.

- seed:

  Optional random seed for reproducible downsampling. If `NULL`
  (default), uses current R random state. Ignored if `pool` is a
  pre-built index or `downsample = 1`.

- progress:

  Logical; if `TRUE`, display a progress bar during computation.
  Progress tracking works by splitting the focal dataset into chunks and
  processing them sequentially. Useful for large datasets. Default is
  `FALSE`.

## Value

Return type depends on input format and query mode.

Returns a data.frame, unless `x` is a SpatRaster and results have
exactly one record per input cell (aggregation mode, or pairwise with
`k = 1`), in which case returns a SpatRaster with one layer per output
variable.

Pairwise mode (`stat = NULL` or `"none"`) returns one row per
focal-analog pair, with the following variables:

- `index`, `x`, `y`: Focal location (1-based index and coordinates)
  corresponding to input `x`

- `analog_index`, `analog_x`, `analog_y`: Analog location corresponding
  to input `pool`

- `clim_dist`: Climate distance (Euclidean or Mahalanobis)

- `geog_dist`: Geographic distance (km for lonlat, projection units
  otherwise)

- Value columns (if `y` provided): one per variable

Aggregation mode (one or more `stat` values) returns one row per focal
location, with the following variables:

- `index`, `x`, `y`: Focal location

- One column per requested statistic. For `stat` with single `y`
  variable: column named by stat (e.g., `sum`, `mean`). For `stat` with
  multiple `y` variables: columns named `{stat}_{varname}` (e.g.,
  `sum_biomass`, `mean_richness`)

- For `stat = "regression"`: columns for `intercept` and each covariate
  name, or `intercept_{varname}` and `{covariate}_{varname}` with
  multiple `y` variables.

All results include metadata attributes (`select`, `stat`, `weight`,
etc.). Use
[`analog_summary()`](https://matthewkling.github.io/analogs/reference/analog_summary.md)
to view a formatted summary.

## Details

This function is a wrapper that calls
[`analog_search()`](https://matthewkling.github.io/analogs/reference/analog_search.md)
using `select = "knn_clim"`.

For each focal location, `analog_similarity()`:

1.  Identifies all reference points within `max_geog` km (and optional
    climate filter).

2.  Selects the `k` closest in *climate* distance.

This is the natural "inverse" of
[`analog_velocity`](https://matthewkling.github.io/analogs/reference/analog_velocity.md):
instead of finding where the focal climate moves geographically, it
finds the closest climatically similar conditions that are
geographically reachable.

## See also

[`analog_search()`](https://matthewkling.github.io/analogs/reference/analog_search.md)
for the underlying flexible analog search function;
[`tiled_analog_search()`](https://matthewkling.github.io/analogs/reference/tiled_analog_search.md)
for memory-safe searches on large raster datasets.

## Examples

``` r
if (FALSE) { # \dontrun{
# One-shot query
im <- analog_similarity(
  x = clim$clim1,
  pool = clim$clim2,
  max_geog = 100,
  k = 20
)

# With pre-built index (for repeated queries)
index <- build_analog_index(clim$clim2)
i1 <- analog_similarity(x = sites1, pool = index, max_geog = 100, k = 20)
i2 <- analog_similarity(x = sites2, pool = index, max_geog = 50, k = 10)
} # }
```
