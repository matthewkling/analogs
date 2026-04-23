# Analog availability: count of all analogs within climate/geographic limits

Computes, for each focal location, how many reference locations satisfy
the supplied climate and geographic constraints. This is useful for
mapping analog "availability" or environmental similarity density.

## Usage

``` r
analog_availability(
  x,
  pool,
  x_cov = NULL,
  coord_type = "auto",
  max_clim = NULL,
  max_geog = NULL,
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

- For `stat = "regression"`: columns for `coef_intercept` and
  `coef_{covariate}`, or `coef_intercept_{varname}` and
  `coef_{covariate}_{varname}` with multiple `y` variables.

- When `se != "none"`: matching SE columns (`se_weighted_mean`,
  `se_intercept`, etc.) for each SE-supporting stat.

All results include metadata attributes (`select`, `stat`, `kernel`,
etc.). Use
[`analog_summary()`](https://matthewkling.github.io/analogs/reference/analog_summary.md)
to view a formatted summary.

## Details

This function is a wrapper that calls
[`analog_search()`](https://matthewkling.github.io/analogs/reference/analog_search.md)
using `select = "all"` and `stat = "count"`.

## References

Stralberg D, Carroll C, Pedlar JH, Wilsey CB, McKenney DW, Nielsen SE
(2018). "Macrorefugia for North American trees and songbirds: Climatic
limiting factors and multi-scale topographic influences." *Global
Ecology and Biogeography*, **27**(6), 690-703.
[doi:10.1111/geb.12731](https://doi.org/10.1111/geb.12731)

Carroll C, Lawler JJ, Roberts DR, Hamann A (2015). "Biotic and climatic
velocity identify contrasting areas of vulnerability to climate change."
*PLOS ONE*, **10**(10), e0140486.
[doi:10.1371/journal.pone.0140486](https://doi.org/10.1371/journal.pone.0140486)

Parks SA, Holsinger LM, Abatzoglou JT, Littlefield CE, Zeller KA (2023).
"Protected areas not likely to serve as steppingstones for species
undergoing climate-induced range shifts." *Global Change Biology*,
**29**(10), 2681-2696.
[doi:10.1111/gcb.16629](https://doi.org/10.1111/gcb.16629)

## See also

[`analog_search()`](https://matthewkling.github.io/analogs/reference/analog_search.md)
for the underlying flexible analog search function;
[`tiled_analog_search()`](https://matthewkling.github.io/analogs/reference/tiled_analog_search.md)
for memory-safe searches on large raster datasets.

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
