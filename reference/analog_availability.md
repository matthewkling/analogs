# Analog availability: count of all analogs within climate/geographic limits

Computes, for each focal location, how many reference locations satisfy
the supplied climate and geographic constraints. This is useful for
mapping analog "availability" or environmental similarity density.

## Usage

``` r
analog_availability(
  x,
  pool,
  max_clim = NULL,
  max_geog = NULL,
  coord_type = "auto",
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

A data.frame with columns: - focal_index - focal_x, focal_y - value (the
count of analogs)

## Examples

``` r
if (FALSE) { # \dontrun{
# One-shot query
avail <- analog_availability(
  x = sites,
  pool = climate_data,
  max_clim = 0.5,
  max_geog = 100
)

# With pre-built index (for repeated queries)
index <- build_analog_index(climate_data)
a1 <- analog_availability(x = sites1, pool = index, max_clim = 0.5, max_geog = 100)
a2 <- analog_availability(x = sites2, pool = index, max_clim = 0.3, max_geog = 50)
} # }
```
