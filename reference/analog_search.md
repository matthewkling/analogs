# Find Climate Analogs

Identifies locations in a reference dataset that are climatically
similar and/or geographically proximal to focal locations. Analog
searches use a two-stage approach: first selecting analogs based on
specified criteria, then optionally aggregating the results.

## Usage

``` r
analog_search(
  x,
  pool,
  x_cov = NULL,
  values = NULL,
  max_clim = NULL,
  max_geog = NULL,
  select = "all",
  k = NULL,
  stat = NULL,
  weight = NULL,
  theta = NULL,
  coord_type = c("auto", "lonlat", "projected"),
  downsample = 1,
  seed = NULL,
  index_res = "auto",
  n_threads = NULL,
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
  location and one column per unique covariance component, or a
  SpatRaster with a layer for each component. For n climate variables,
  there are n\*(n+1)/2 unique components, ordered as: variances first
  (diagonals), then covariances (upper triangle by row).

- values:

  Optional user-defined variables for each reference location in `pool`
  to aggregate across selected analogs. Can be a numeric vector (single
  variable), matrix or data.frame with numeric columns (multiple
  variables), or a SpatRaster with one or more numeic layers. Must have
  exactly the same number of reference locations as `pool`.

  When provided, enables value-based aggregation stats `"sum"`,
  `"mean"`, `"weighted_sum"`, and `"weighted_mean"`. For stat =
  NULL/"none" (pairs mode), value columns are included in output for
  each analog pair.

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

  - `"ess"`: Kish's effective sample size (ESS), computed as the squared
    sum of weights divided by the sum of squared weights (requires
    `weight`).

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

- coord_type:

  Coordinate system type:

  - `"auto"` (default): Automatically detect from coordinate ranges.

  - `"lonlat"`: Unprojected lon/lat coordinates (uses great-circle
    distance; assumes `max_geog` is in km).

  - `"projected"`: Projected XY coordinates (uses planar distance;
    assumes `max_geog` is in projection units).

- downsample:

  Optional downsampling rate (0-1) for the reference pool, indicating
  the proportion of points to retain. Values \< 1 reduce memory and
  improve speed at some cost to precision. Default is 1.0 (no
  downsampling). Ignored if `pool` is a pre-built index.

- seed:

  Optional random seed for reproducible downsampling. If `NULL`
  (default), uses current R random state. Ignored if `pool` is a
  pre-built index or `downsample = 1`.

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

### Parameter categories

- *Data parameters* (`x`, `pool`, `x_cov`, `values`, `coord_type`) give
  attributes of the data on which to operate.

- *Selection parameters* (`select`, `max_clim`, `max_geog`, `k`) define
  which analogs to `select` from the `pool` for each `x`.

- *Aggregation parameters* (`stat`, `weight`, `theta`) control how
  selected analogs are summarized.

- *Computation parameters* (`n_threads`, `index_res`, `downsample`,
  `seed`, `progress`) specify behavior for controlling compute
  performance.

### Distance metrics

Geographic distance can be computed for lon/lat coordinates
(great-circle distance) or projected coordinates (planar distance).

Climate similarity is measured using Euclidean or Mahalanobis distance
in climate space. In general, when multiple climate variables are used,
it is recommended to use pre-whitened (scaled) climate data, to avoid
major artifacts from climate variables with different units.
Pre-whitening can be done using
[`scale()`](https://rdrr.io/r/base/scale.html) for dataset-wide
Euclidean distances, or
[`mahalanobis_transform()`](https://matthewkling.github.io/analogs/reference/mahalanobis_transform.md)
for dataset-wide Mahalanobis distances.

The function also supports climate distance calculations based on *local
temporal* covariance structure at focal locations, via the `x_cov`
parameter. These local covariance values need to be pre-calculated.

### Computational optimization

The analog search architecture is designed with compute performance in
mind:

- All internal computations are done in C++.

- Searches use a lattice-based indexing structure to efficiently search
  through large reference datasets. By default, the lattice resolution
  is tuned for optimal performance.

- Parallel processing is available via the `threads` parameter.

- You can `downsample` prohibitively large reference pool datasets to
  improve speed and memory, using a stratified sampling scheme that
  reduces loss of precision relative to random sampling.

- For large datasets, enable `progress = TRUE` to display a progress bar
  during computation.

- For raster datasets that are too large to fit in memory,
  [`tiled_analog_search()`](https://matthewkling.github.io/analogs/reference/tiled_analog_search.md)
  offers a memory-safe option.

## References

Mahony CR, Cannon AJ, Wang T, Aitken SN (2017). "A closer look at novel
climates: new methods and insights at continental to landscape scales."
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

[`tiled_analog_search()`](https://matthewkling.github.io/analogs/reference/tiled_analog_search.md)
offers memory-safe searches on large raster datasets. Helper functions
such as
[`analog_impact()`](https://matthewkling.github.io/analogs/reference/analog_impact.md),
[`analog_velocity()`](https://matthewkling.github.io/analogs/reference/analog_velocity.md),
and
[`analog_intensity()`](https://matthewkling.github.io/analogs/reference/analog_intensity.md)
offer simpler interfaces for common search types.
