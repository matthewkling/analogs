# Climate impact: nearest climate analogs within a geographic envelope

Computes, for each focal location, the climate–nearest neighbor(s) in a
reference dataset that satisfy a specified geographic distance
threshold. This helper wraps
[`analog_search`](https://matthewkling.github.io/analogs/reference/analog_search.md)
using `mode = "knn_clim"`.

## Usage

``` r
analog_impact(
  x,
  pool,
  max_geog,
  x_cov = NULL,
  k = 20,
  max_clim = NULL,
  coord_type = "auto",
  report_dist = TRUE,
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
  For example:

  - 2 variables: c(var1, var2, cov12)

  - 3 variables: c(var1, var2, var3, cov12, cov13, cov23)

  When provided, all climate distances are computed as Mahalanobis
  distances using each focal's covariance structure. For focal-specific
  variances only (no covariance), set off-diagonal covariances to zero.
  Default is NULL (Euclidean climate distance).

- k:

  Number of nearest analogs to return per focal location for kNN modes.
  Required when `mode` is `"knn_geog"` or `"knn_clim"`; must be `NULL`
  for other modes.

- max_clim:

  Maximum climate distance constraint (default: NULL = no climate
  constraint). Can be either:

  - A scalar: Euclidean radius in climate space (e.g., 0.5)

  - A vector: Per-variable absolute differences (length must equal
    number of climate variables)

  Only reference locations within this climate distance are considered.
  When `x_cov` is provided, scalar thresholds are interpreted in
  Mahalanobis distance units.

- coord_type:

  Coordinate system type (default: "auto"):

  - `"auto"`: Automatically detect from coordinate ranges.

  - `"lonlat"`: Unprojected lon/lat coordinates (uses great-circle
    distance; assumes `max_geog` is in km).

  - `"projected"`: Projected XY coordinates (uses planar distance;
    assumes `max_geog` is in projection units).

- report_dist:

  Logical; if TRUE (default), include distance columns in output when
  `mode` is `"knn_geog"`, `"knn_clim"` or `"all"`. Set to FALSE for more
  compact output.

- index_res:

  Tuning parameter giving the number of bins per dimension of the
  internally-used lattice search index. Either:

  - A positive integer.

  - `"auto"` (the default): Automatically tune the intex resolution by
    optimizing compute time on a subsample of focal points. If focal has
    relatively few rows, auto-tuning is skipped and a default resolution
    of 16 is used.

  Ignored if `pool` is an `analog_index` (uses index's resolution).

- n_threads:

  Optional integer number of threads to use for the computation. If
  `NULL` (default), the global RcppParallel setting is used (see
  [`RcppParallel::setThreadOptions`](https://rdrr.io/pkg/RcppParallel/man/setThreadOptions.html)).

## Value

A data.frame with one row per focal–analog pair, including:

- `focal_index`, `analog_index`

- `focal_x`, `focal_y`, `analog_x`, `analog_y`

- `clim_dist`, `geog_dist` (if `report_dist = TRUE`)

Diagnostic attributes from the underlying spatial index are preserved.

## Details

It is useful for estimating the potential ecological impact of local
climate change: e.g., how climate conditions at a site compare to those
available within a species' dispersal range.

For each focal location, `analog_impact()`:

1.  Identifies all reference points within `max_geog` km (and optional
    climate filter).

2.  Selects the `k` closest in *climate* distance.

This is the natural "inverse" of
[`analog_velocity`](https://matthewkling.github.io/analogs/reference/analog_velocity.md):
instead of finding where the focal climate moves geographically, it
finds the closest climatically similar conditions that are
geographically reachable.

## Examples

``` r
if (FALSE) { # \dontrun{
# One-shot query
im <- analog_impact(
  x = clim$clim1,
  pool = clim$clim2,
  max_geog = 100,
  k = 20
)

# With pre-built index (for repeated queries)
index <- build_analog_index(clim$clim2)
i1 <- analog_impact(x = sites1, pool = index, max_geog = 100, k = 20)
i2 <- analog_impact(x = sites2, pool = index, max_geog = 50, k = 10)
} # }
```
