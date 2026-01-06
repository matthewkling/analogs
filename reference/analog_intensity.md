# Analog intensity: weighted sum of analogs within climate/geographic limits

Computes, for each focal location, the sum of weights of all reference
locations that satisfy the supplied climate and geographic constraints.
The weights are controlled by the `weight` and `theta` arguments and are
applied after filtering.

## Usage

``` r
analog_intensity(
  x,
  pool,
  x_cov = NULL,
  coord_type = "auto",
  max_clim = NULL,
  max_geog = NULL,
  weight = c("uniform", "inverse_clim", "inverse_geog", "gaussian_clim", "gaussian_geog",
    "gaussian_joint", "inverse_joint"),
  theta = NULL,
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

  - `"gaussian_joint"`: Gaussian kernel on combined distance, weight =
    exp(-(clim_dist^2 / (2*sigma_clim^2) + geog_dist^2 /
    (2*sigma_geog^2))), with sigmas given by `theta`.

  - `"inverse_joint"`: Inverse joint distance, weight = 1 /
    (sqrt(clim_dist^2 + geog_dist^2) + eps), with epsilon given by
    `theta`.

- theta:

  Optional numeric parameter used by weighting functions when `stat`
  includes `"sum_weights"` or `"mean_weights"` and `weight` is not
  `"uniform"`. Interpretation depends on `weight`:

  - For `"inverse_clim"` or `"inverse_geog"`: epsilon value added to
    distances (scalar; default: 1e-12 for climate, 1e-6 for geography).

  - For `"gaussian_clim"` or `"gaussian_geog"`: sigma bandwidth
    parameter (scalar; larger values = slower decay with distance).

  - For `"gaussian_joint"` or `"inverse_joint"`: 2-element vector
    `c(theta_clim, theta_geog)` (defaults: 1 for climate, 1 for
    geography).

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

## Details

This function is a wrapper that calls
[`analog_search()`](https://matthewkling.github.io/analogs/reference/analog_search.md)
using `select = "all"` and `stat = "sum_weights"`.

## References

Mahony CR, Cannon AJ, Wang T, Aitken SN (2017). "A closer look at novel
climates: New methods and insights at continental to landscape scales."
*Global Change Biology*, **23**(9), 3934-3955.
[doi:10.1111/gcb.13645](https://doi.org/10.1111/gcb.13645)

Abatzoglou JT, Dobrowski SZ, Parks SA (2020). "Multivariate climate
departures have outpaced univariate changes across global lands."
*Scientific Reports*, **10**(1), 3891.
[doi:10.1038/s41598-020-60270-5](https://doi.org/10.1038/s41598-020-60270-5)

Williams JW, Jackson ST, Kutzbach JE (2007). "Projected distributions of
novel and disappearing climates by 2100 AD." *Proceedings of the
National Academy of Sciences*, **104**(14), 5738-5742.
[doi:10.1073/pnas.0606292104](https://doi.org/10.1073/pnas.0606292104)

## See also

[`analog_search()`](https://matthewkling.github.io/analogs/reference/analog_search.md)
for the underlying flexible analog search function;
[`tiled_analog_search()`](https://matthewkling.github.io/analogs/reference/tiled_analog_search.md)
for memory-safe searches on large raster datasets.

## Examples

``` r
if (FALSE) { # \dontrun{
# One-shot query with inverse weighting
intens <- analog_intensity(
  x = sites,
  pool = climate_data,
  max_clim = 0.5,
  max_geog = 100,
  weight = "inverse_clim"
)

# Gaussian weighting by climate distance
intens_gauss <- analog_intensity(
  x = sites,
  pool = climate_data,
  max_clim = 0.5,
  max_geog = 100,
  weight = "gaussian_clim",
  theta = 0.2  # bandwidth parameter
)

# Joint Gaussian weighting (both climate and geography)
intens_joint <- analog_intensity(
  x = sites,
  pool = climate_data,
  max_clim = 0.5,
  max_geog = 100,
  weight = "gaussian_joint",
  theta = c(0.2, 50)  # c(clim_bandwidth, geog_bandwidth)
)

# With pre-built index (for repeated queries)
index <- build_analog_index(climate_data)
i1 <- analog_intensity(x = sites1, pool = index, max_clim = 0.5,
                       weight = "inverse_clim")
i2 <- analog_intensity(x = sites2, pool = index, max_geog = 100,
                       weight = "inverse_geog")
} # }
```
