# Analog similarity: best climate analogs within a geographic envelope

Finds, for each focal location, the climate–nearest neighbor(s) in a
reference dataset that satisfy a specified geographic distance
threshold. This function is a wrapper that calls
[`analog_search()`](https://matthewkling.github.io/analogs/reference/analog_search.md)
using `select = "knn_clim"`.

## Usage

``` r
analog_similarity(
  x,
  pool,
  x_cov = NULL,
  y = NULL,
  weight = NULL,
  coord_type = "auto",
  max_geog,
  max_clim = NULL,
  k = 20,
  index_res = "auto",
  cell_area_weight = "auto",
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

  Optional vector, factor, matrix/data.frame, or SpatRaster giving
  values for each reference location (must have same number of
  rows/cells as `pool`). Required for stats `"sum"`, `"mean"`,
  `"weighted_sum"`, `"weighted_mean"`, `"regression"`, and `"tabulate"`.
  Numeric for continuous stats; factor or coercible-to-factor
  (character, integer, logical) for `stat = "tabulate"`.

- weight:

  Optional pool site weights for use in aggregation. Numeric vector,
  single-column matrix/data.frame, or single-layer SpatRaster, with one
  value per row/cell of `pool`. For aggregation stats like
  `"weighted_mean"`, `"regression"`, etc., weights multiply through the
  weighted aggregation alongside any kernel weighting and cell-area
  weighting; they do not influence which analogs are selected by `knn_*`
  modes (selection remains distance-only). They are reported in pair
  mode as a `user_weight` column. Values must be non-negative; `NA` is
  allowed and treated as 0 (the point is excluded from aggregation).
  Default `NULL` means no user-supplied weights.

  If you want to exclude a static subset of pool sites entirely, masking
  `pool` (and any associated `y` / `covariates`) upfront is more
  efficient than passing `weight = 0` for those sites, since the lattice
  index will not have to scan or distance-compute against them. Use
  `weight = 0` for cases where the mask varies per query against a
  shared index, or where some sites have a continuous weight and others
  should be excluded.

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
    of 16 is used. Auto-tuning is **not supported** when
    `downsample < 1`, because the speed-optimal resolution can sometimes
    result in higher uncertainty of stat results under downsampling. In
    that case set `index_res` explicitly; finer values (e.g. 32)
    generally give better accuracy at the possible cost of query speed.

  Ignored if `pool` is an `analog_index` (uses index's resolution).

- cell_area_weight:

  Controls cell-area weighting when `pool` is a raster. One of `"auto"`
  (default; on for raster pools, off otherwise), `TRUE` (force on;
  errors if `pool` is not a SpatRaster), or `FALSE` (force off).
  Cell-area weights correct aggregation statistics for non-uniform cell
  areas (e.g. lonlat grids near the poles, or projected grids on
  non-equal-area projections); they are computed via
  [`terra::cellSize()`](https://rspatial.github.io/terra/reference/cellSize.html)
  and normalized to mean 1. When `pool` is a pre-built `analog_index`,
  this argument must agree with the index's stored configuration:
  `cell_area_weight = FALSE` errors if the index was built with
  cell-area weighting on (rebuild the index instead).

- n_threads:

  Optional integer number of threads to use for the computation. If
  `NULL` (default), the global RcppParallel setting is used (see
  [`RcppParallel::setThreadOptions`](https://rdrr.io/pkg/RcppParallel/man/setThreadOptions.html)).

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

- progress:

  Logical; if `TRUE`, display a progress bar during computation.
  Progress tracking works by splitting the focal dataset into chunks and
  processing them sequentially. Useful for large datasets. Default is
  `FALSE`.

## Value

A data.frame, or a SpatRaster when `x` is one and `k = 1`. Contains one
row per focal-analog pair with `index`, `x`, `y`, `analog_index`,
`analog_x`, `analog_y`, `clim_dist`, and `geog_dist`. See
[`analog_search()`](https://matthewkling.github.io/analogs/reference/analog_search.md)
for full column conventions and
[`metadata()`](https://matthewkling.github.io/analogs/reference/metadata.md)
for attached metadata attributes.

## Details

For each focal location, `analog_similarity()`:

1.  Identifies all reference points within `max_geog` km (and optional
    climate filter).

2.  Selects the `k` closest in *climate* distance.

This is the natural "inverse" of
[`analog_velocity`](https://matthewkling.github.io/analogs/reference/analog_velocity.md):
instead of finding where the focal climate moves geographically, it
finds the closest climatically similar conditions that are
geographically reachable.

Among other uses, this operation is often the first step in a
traditional analog impact modeling (AIM) analysis – though see
[`analog_impact()`](https://matthewkling.github.io/analogs/reference/analog_impact.md)
for a more complete AIM implementation.

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
