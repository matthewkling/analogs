# Climate velocity: nearest geographic analogs within a climate envelope

Computes, for each focal location, the geographic nearest neighbor(s) in
a reference dataset that satisfy a specified climate distance threshold.
This helper wraps
[`find_analogs`](https://matthewkling.github.io/analogs/reference/find_analogs.md)
using `mode = "knn_geog"` and is most commonly used for estimating
climate velocity (the rate and direction at which organisms would have
to move to track constant climate conditions).

## Usage

``` r
analog_velocity(
  x,
  pool,
  max_clim,
  k = 1,
  max_geog = NULL,
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

- max_clim:

  Maximum climate distance constraint (default: NULL = no climate
  constraint). Can be either:

  - A scalar: Euclidean radius in climate space (e.g., 0.5)

  - A vector: Per-variable absolute differences (length must equal
    number of climate variables)

  Only reference locations within this climate distance are considered.

- k:

  Number of nearest analogs to return per focal location for kNN modes.
  Required when `mode` is `"knn_geog"` or `"knn_clim"`; must be `NULL`
  for other modes.

- max_geog:

  Maximum geographic distance constraint (default: NULL = no geographic
  constraint). When specified, only reference locations within this
  distance are considered. Radius units should be specified in
  kilometers if `coord_type = "lonlat"`, or in projected coordinate
  units if `coord_type = "projected"`.

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

Diagnostic attributes (e.g., binning statistics) from the underlying
spatial index are preserved.

## Details

For each focal point, this function:

1.  Identifies all reference points satisfying the climate (and optional
    geographic) threshold(s).

2.  Among those, selects the `k` nearest in *geographic* distance.

This is the classical operation needed for estimating *climate
velocity*: the minimum relocation distance needed to maintain similar
climatic conditions under temporal change.

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
} # }
```
