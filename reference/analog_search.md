# Find Climate Analogs

Identifies locations in a reference dataset that are climatically
similar to focal locations, with optional constraints on climate
distance and geographic distance. Analog searches use a two-stage
approach: first selecting analogs based on specified criteria, then
optionally aggregating the results.

## Usage

``` r
analog_search(
  x,
  pool,
  x_cov = NULL,
  values = NULL,
  coord_type = c("auto", "lonlat", "projected"),
  max_clim = NULL,
  max_geog = NULL,
  select = "all",
  k = NULL,
  stat = NULL,
  weight = NULL,
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

- select:

  Character string specifying the analog selection strategy. One of:

  - `"all"` (default): Select all analogs that satisfy the `max_clim`
    and `max_geog` constraints.

  - `"knn_clim"`: For each focal, select up to `k` analogs with smallest
    climate distance, subject to filters.

  - `"knn_geog"`: For each focal, select up to `k` analogs with smallest
    geographic distance, subject to filters.

- k:

  Number of nearest analogs to return per focal location for kNN
  selection modes. Required when `select` is `"knn_geog"` or
  `"knn_clim"`; must be `NULL` for `select = "all"`.

- stat:

  Statistic(s) used to aggregate selected analogs. Either:

  - `NULL` or `"none"`: Return all selected analog pairs as a
    data.frame.

  - `"count"`: For each focal, count the number of selected analogs.

  - `"sum_weights"`: For each focal, sum the weights of selected analogs
    (see `weight` and `theta`).

  - `"mean_weights"`: For each focal, mean of weights of selected
    analogs.

  - `"sum"`: Sum of values across analogs (requires `values`).

  - `"mean"`: Mean of values across analogs (requires `values`).

  - `"weighted_sum"`: Sum of (value × weight) across analogs (requires
    `values` and `weight`).

  - `"weighted_mean"`: Weighted mean of values across analogs (requires
    `values` and `weight`).

  - A character vector combining multiple stats (e.g.,
    `c("count", "sum", "mean")`). Note: `"none"` cannot be combined with
    other stats.

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

## Details

Parameters fall into four categories:

- *Data parameters* (`x`, `pool`, `x_cov`, `values`, `coord_type`) give
  attributes of the data on which to operate.

- *Selection parameters* (`select`, `max_clim`, `max_geog`, `k`) define
  which analogs to `select` from the `pool` for each `x`.

- *Aggregation parameters* (`stat`, `weight`, `theta`) control how
  selected analogs are summarized.

- *Computation parameters* (`index_res`, `n_threads`) specify internal
  behavior affecting compute performance.

The function uses a spatial indexing structure (lattice-based) to
quickly search through large reference datasets. Climate similarity is
measured using Euclidean distance in climate space (ideally
pre-whitened). Geographic distance can be computed for lon/lat
coordinates (great-circle distance) or projected coordinates (planar
distance).

## References

Mahony CR, Cannon AJ, Wang T, Aitken SN (2017). "A closer look at novel
climates: New methods and insights at continental to landscape scales."
*Global Change Biology*, **23**(9), 3934-3955.
[doi:10.1111/gcb.13645](https://doi.org/10.1111/gcb.13645)

Hamann A, Roberts DR, Barber QE, Carroll C, Nielsen SE (2015). "Velocity
of climate change algorithms for guiding conservation and management."
*Global Change Biology*, **21**(2), 997-1004.
[doi:10.1111/gcb.12736](https://doi.org/10.1111/gcb.12736)

Grenier P, Parent A-C, Huard D, Anctil F, Chaumont D (2013). "An
assessment of six dissimilarity metrics for climate analogs." *Journal
of Applied Meteorology and Climatology*, **52**(4), 733-752.
[doi:10.1175/JAMC-D-12-0170.1](https://doi.org/10.1175/JAMC-D-12-0170.1)

## See also

[`analog_summary()`](https://matthewkling.github.io/analogs/reference/analog_summary.md)
to print a summary of search result metadata.

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

# With user-defined values
values_df <- data.frame(
  biomass = runif(nrow(ref), 0, 100),
  richness = rpois(nrow(ref), 20)
)
#> Error: object 'ref' not found

analog_search(x = focal, pool = ref, values = values_df,
              stat = c("count", "sum", "mean"),
              max_clim = 0.5)
#> Error: object 'ref' not found
# Returns: index, x, y, count, sum_biomass, mean_biomass, sum_richness, mean_richness

# Weighted aggregation of values
analog_search(x = focal, pool = ref, values = values_df$biomass,
              stat = c("weighted_sum", "weighted_mean"),
              weight = "gaussian_clim", theta = 0.2,
              max_clim = 0.5)
#> Error: object 'ref' not found
# Returns: index, x, y, weighted_sum, weighted_mean
```
