# Climate velocity: geographically nearest climate analogs

Computes, for each focal location, the geographic nearest neighbor(s) in
a reference dataset that satisfy a specified maximum climate distance
threshold. This helper wraps
[`analog_search`](https://matthewkling.github.io/analogs/reference/analog_search.md)
using `select = "knn_geog"` and is used for estimating analog-based
climate velocity.

## Usage

``` r
analog_velocity(
  x,
  pool,
  x_cov = NULL,
  values = NULL,
  coord_type = "auto",
  max_clim,
  max_geog = NULL,
  k = 1,
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

  - An `analog_index()` object created by
    [`build_analog_index()`](https://matthewkling.github.io/analogs/reference/build_analog_index.md)
    (for repeated queries).

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

  Coordinate system type:

  - `"auto"` (default): Automatically detect from coordinate ranges.

  - `"lonlat"`: Unprojected lon/lat coordinates (uses great-circle
    distance; assumes `max_geog` is in km).

  - `"projected"`: Projected XY coordinates (uses planar distance;
    assumes `max_geog` is in projection units).

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

- Value columns (if `values` provided): one per variable

Aggregation mode (one or more `stat` values) returns one row per focal
location, with the following variables:

- `index`, `x`, `y`: Focal location

- One column per requested statistic. For `stat` with single `values`
  variable: column named by stat (e.g., `sum`, `mean`). For `stat` with
  multiple `values` variables: columns named `{stat}_{varname}` (e.g.,
  `sum_biomass`, `mean_richness`)

All results include metadata attributes (`select`, `stat`, `weight`,
etc.). Use
[`analog_summary()`](https://matthewkling.github.io/analogs/reference/analog_summary.md)
to view a formatted summary.

## References

Hamann A, Roberts DR, Barber QE, Carroll C, Nielsen SE (2015). "Velocity
of climate change algorithms for guiding conservation and management."
*Global Change Biology*, **21**(2), 997-1004.
[doi:10.1111/gcb.12736](https://doi.org/10.1111/gcb.12736)

Dobrowski SZ, Parks SA (2016). "Climate change velocity underestimates
climate change exposure in mountainous regions." *Nature
Communications*, **7**, 12349.
[doi:10.1038/ncomms12349](https://doi.org/10.1038/ncomms12349)

## Examples

``` r
if (FALSE) { # \dontrun{
# One-shot query
v <- analog_velocity(
  x = clim$clim1,
  pool = clim$clim2,
  max_clim = 0.5,
  k = 1
)

# With pre-built index (for repeated queries)
index <- build_analog_index(clim$clim2)
v1 <- analog_velocity(x = sites1, pool = index, max_clim = 0.5, k = 1)
v2 <- analog_velocity(x = sites2, pool = index, max_clim = 0.3, k = 1)

# With focal-specific covariance matrices
v_mahal <- analog_velocity(
  x = clim$clim1,
  pool = clim$clim2,
  x_cov = baseline_covariances,
  max_clim = 2,  # In Mahalanobis distance units
  k = 1
)
} # }
```
