# Analog intensity: weighted sum of analogs within climate/geographic limits

Computes, for each focal location, the weighted sum of all reference
locations that satisfy the supplied climate and geographic constraints.
The weights are controlled by the `weight` and `theta` arguments and are
applied after filtering.

## Usage

``` r
analog_intensity(
  x,
  pool,
  x_cov = NULL,
  values = NULL,
  coord_type = "auto",
  max_clim = NULL,
  max_geog = NULL,
  weight = c("uniform", "inverse_clim", "inverse_geog", "gaussian_clim", "gaussian_geog",
    "gaussian_joint", "inverse_joint"),
  theta = NULL,
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
