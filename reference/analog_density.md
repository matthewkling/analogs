# Analog density: kernel-weighted analog count within climate/geographic limits

Computes, for each focal location, the kernel-weighted sum of all
reference locations that satisfy the supplied climate and geographic
constraints.

## Usage

``` r
analog_density(
  x,
  pool,
  x_cov = NULL,
  weight = NULL,
  coord_type = "auto",
  stat = c("sum_weights"),
  max_clim = NULL,
  max_geog = NULL,
  kernel = c("uniform", "inverse_clim", "inverse_geog", "gaussian_clim", "gaussian_geog",
    "gaussian_joint", "inverse_joint"),
  theta = NULL,
  normalize = "auto",
  index_res = "auto",
  cell_area_weight = "auto",
  n_threads = NULL,
  downsample = 1,
  seed = NULL,
  progress = FALSE,
  ...
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

- stat:

  Character vector of one or more density-style summary statistics.
  Allowed options: `"sum_weights"` (default), `"mean_weights"`,
  `"count"`, `"ess"`. See
  [`analog_search()`](https://matthewkling.github.io/analogs/reference/analog_search.md)
  for detailed stat descriptions.

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

- kernel:

  Kernel decay function for weighting matches, used only when `stat`
  includes a weighted aggregation (`"sum_weights"`, `"mean_weights"`,
  `"weighted_sum"`, `"weighted_mean"`, `"ess"`, `"regression"`, or
  `"tabulate"`). One of:

  - `"uniform"`: All matches weighted equally (kernel weight = 1.0).

  - `"inverse_clim"`: Inverse climate distance, kernel weight = 1 /
    (climate_distance + eps), with epsilon given by `theta`.

  - `"inverse_geog"`: Inverse geographic distance, kernel weight = 1 /
    (geographic_distance + eps), with epsilon given by `theta`.

  - `"gaussian_clim"`: Gaussian kernel on climate distance, kernel
    weight = exp(-climate_distance^2 / (2 sigma^2)), with sigma given by
    `theta`.

  - `"gaussian_geog"`: Gaussian kernel on geographic distance, kernel
    weight = exp(-geographic_distance^2 / (2 sigma^2)), with sigma given
    by `theta`.

  - `"gaussian_joint"`: Gaussian kernel on combined distance, kernel
    weight = exp(-(clim_dist^2 / (2 sigma_clim^2) + geog_dist^2 / (2
    sigma_geog^2))), with sigmas given by `theta`.

  - `"inverse_joint"`: Inverse joint distance, kernel weight = 1 /
    (sqrt(clim_dist^2 + geog_dist^2) + eps), with epsilon given by
    `theta`.

- theta:

  Optional numeric parameter controlling the shape of the weighting
  `kernel`, used whenever `kernel` is active (i.e. whenever `stat`
  includes a weighted aggregation) and `kernel` is not `"uniform"`.
  Interpretation depends on `kernel`:

  - For `"inverse_clim"` or `"inverse_geog"`: epsilon value added to
    distances (scalar; default: 1e-12 for climate, 1e-6 for geography).

  - For `"gaussian_clim"` or `"gaussian_geog"`: sigma bandwidth
    parameter (scalar; larger values = slower decay with distance).

  - For `"gaussian_joint"` or `"inverse_joint"`: 2-element vector
    `c(theta_clim, theta_geog)` (defaults: 1 for climate, 1 for
    geography).

  See
  [`kernel_params()`](https://matthewkling.github.io/analogs/reference/kernel_params.md)
  for help choosing `theta` and `max_clim` / `max_geog` values that work
  well together.

- normalize:

  One of `TRUE`, `FALSE`, or `"auto"` (default). Only used if `stat`
  includes `"sum_weights"` or `"tabulate"` and `pool` is a raster. When
  active, results for these stats are divided by a global scalar so that
  they represent a fraction of a theoretically "perfect" scenario where
  the full search area within `max_geog` is occupied wall-to-wall by
  cells whose climate exactly matches `x`. See details under
  [`analog_search()`](https://matthewkling.github.io/analogs/reference/analog_search.md)
  for more info.

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

- ...:

  Additional arguments forwarded to
  [`analog_search()`](https://matthewkling.github.io/analogs/reference/analog_search.md)
  (e.g. `mean_cell_area` for tiled workflows).

## Value

A data.frame, or a SpatRaster when `x` is one. Contains `index`, `x`,
`y`, and `sum_weights`. When `normalize = TRUE`, the value of `D_max` is
attached as a result attribute. See
[`metadata()`](https://matthewkling.github.io/analogs/reference/metadata.md)
for attached metadata attributes.

## Details

This function is a wrapper that calls
[`analog_search()`](https://matthewkling.github.io/analogs/reference/analog_search.md)
using `select = "all"` and `stat = "sum_weights"`. It also supports
`"count"`, `"mean_weights"`, and `"ess"` as additional stats.

By default (`normalize = "auto"`), the returned `sum_weights` column is
divided by a global scalar `D_max` whenever `pool` meets the
prerequisites (raster `pool` with cell-area weighting active). This
gives values on roughly `[0, 1]`, interpretable as the fraction of
theoretical maximum analog density achieved at each focal. See
[`analog_search()`](https://matthewkling.github.io/analogs/reference/analog_search.md)
for more details.

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
# Normalized density (default): returns sum_weights on roughly [0, 1]
dens <- analog_density(
  x = sites,
  pool = climate_data,
  max_clim = 0.5,
  max_geog = 100,
  kernel = "gaussian_clim",
  theta = 0.2
)
attr(dens, "D_max")  # the global denominator used

# Raw (unnormalized) kernel-weighted sums
dens_raw <- analog_density(
  x = sites,
  pool = climate_data,
  max_clim = 0.5,
  max_geog = 100,
  kernel = "gaussian_clim",
  theta = 0.2,
  normalize = FALSE
)

# Joint Gaussian weighting (both climate and geography)
intens_joint <- analog_density(
  x = sites,
  pool = climate_data,
  max_clim = 0.5,
  max_geog = 100,
  kernel = "gaussian_joint",
  theta = c(0.2, 50)  # c(clim_bandwidth, geog_bandwidth)
)

# With pre-built index (for repeated queries)
index <- build_analog_index(climate_data)
i1 <- analog_density(x = sites1, pool = index, max_clim = 0.5,
                     max_geog = 100, kernel = "inverse_clim")
i2 <- analog_density(x = sites2, pool = index, max_geog = 100,
                     kernel = "gaussian_clim", theta = 0.3)
} # }
```
