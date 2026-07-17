# Analog density: kernel-weighted analog count within environmental/geographic limits

Computes, for each focal location, the kernel-weighted sum of all
reference locations that satisfy the supplied environmental and
geographic constraints.

## Usage

``` r
analog_density(
  x,
  pool,
  x_cov = NULL,
  weight = NULL,
  coord_type = "auto",
  stat = c("sum_weights"),
  env = NULL,
  geog = NULL,
  normalize = "auto",
  env_res_adj = "auto",
  geog_res_adj = "auto",
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
  matrix/data.frame with columns x, y, and environmental variables, or a
  SpatRaster with environmental variable layers.

- pool:

  The reference dataset to search for analogs. Either:

  - Matrix/data.frame with columns x, y, and environmental variables, or
    SpatRaster with environmental variable layers, OR

  - An `analog_index` object created by
    [`build_analog_index()`](https://matthewkling.github.io/analogs/reference/build_analog_index.md)
    (for repeated queries).

- x_cov:

  Optional focal-specific covariance matrices for Mahalanobis distance
  calculations. Should be a matrix or data.frame with one row per focal
  location and one column per unique covariance component, or a
  SpatRaster with a layer for each component. For n environmental
  variables, there are n\*(n+1)/2 unique components, ordered as:
  variances first (diagonals), then covariances (upper triangle by row).

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

- env, geog:

  Per-family distance treatment, each a
  [`kernel()`](https://matthewkling.github.io/analogs/reference/kernel.md)
  object (or `NULL`). A kernel bundles the hard distance threshold, the
  weighting kernel shape, and the kernel's scale for one family:
  environmental (`env`) or geography (`geog`).
  `kernel(weight, theta, max, min)` where:

  - `max`: hard upper distance threshold — candidates beyond it (in that
    family's distance) are excluded. For `env`, `max` may be a single
    Euclidean radius or a per-variable vector of absolute-difference
    thresholds (length equal to the number of environmental variables);
    scalar environmental thresholds are in Mahalanobis units when
    `x_cov` is supplied. For `geog`, `max` is a single radius
    (kilometers when `coord_type = "lonlat"`, projected units
    otherwise).

  - `min`: hard lower distance threshold — candidates closer than it are
    excluded, so retained candidates form an annulus `min <= d <= max`.
    Supported only for `geog` (a single radius, same units as the
    geographic `max`); setting `min` on `env` is an error. Mainly used
    to impose a spatial buffer for cross-validation (see
    [`analog_cv()`](https://matthewkling.github.io/analogs/reference/analog_cv.md)).

  - `weight`: kernel shape for weighted aggregations — `"uniform"` (no
    distance weighting), `"gaussian"` (`exp(-d^2 / (2 theta^2))`), or
    `"inverse"` (`1 / (1 + d / theta)`). The overall kernel weight is
    the product of the two families' weights, so shapes may be mixed
    (e.g. an inverse environmental kernel with a Gaussian geographic
    kernel).

  - `theta`: the kernel's scale (Gaussian bandwidth, or inverse
    half-weight distance). See
    [`kernel_params()`](https://matthewkling.github.io/analogs/reference/kernel_params.md)
    for calibrated values.

  A `NULL` kernel (the default for both) applies no threshold and no
  weighting for that family. See
  [`kernel()`](https://matthewkling.github.io/analogs/reference/kernel.md)
  for details.

- normalize:

  One of `TRUE`, `FALSE`, or `"auto"` (default). Only used if `stat`
  includes `"sum_weights"` or `"tabulate"` and `pool` is a raster. When
  active, results for these stats are divided by a global scalar so that
  they represent a fraction of a theoretically "perfect" scenario where
  the full search area within `max_geog` is occupied wall-to-wall by
  cells whose environment exactly matches `x`. See details under
  [`analog_search()`](https://matthewkling.github.io/analogs/reference/analog_search.md)
  for more info.

- env_res_adj, geog_res_adj:

  Control the lattice search-index resolution of the environmental and
  geographic families, each a multiplier on a data-dependent default
  (targeting ~50 pool points per occupied bin, split between families by
  effective dimensionality, so it scales with pool size). Each is
  either:

  - A non-negative number: `1` uses the default for that family, larger
    values are finer, smaller are coarser, and `0` deactivates it.

  - `"auto"` (the default for both): tune a single overall resolution
    scale by optimizing compute time on a subsample of focal points. If
    focal has relatively few rows, tuning is skipped. Not supported when
    `downsample < 1` (set explicit numeric values instead).

  A family that the query does not constrain (no corresponding `max_*`
  and not the knn sort key) is **automatically deactivated**, overriding
  any explicit value (with a message), since binning an unconstrained
  family only costs time. Ignored if `pool` is an `analog_index` (uses
  the index's resolution).

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
  `downsample < 1`, resolution must be set explicitly via `geog_res_adj`
  / `env_res_adj` (auto-tuning is not supported in this case; see those
  parameters for details).

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
  env = kernel("gaussian", theta = 0.2, max = 0.5),
  geog = kernel(max = 100)
)
attr(dens, "D_max")  # the global denominator used

# Raw (unnormalized) kernel-weighted sums
dens_raw <- analog_density(
  x = sites,
  pool = climate_data,
  env = kernel("gaussian", theta = 0.2, max = 0.5),
  geog = kernel(max = 100),
  normalize = FALSE
)

# Joint Gaussian weighting (both climate and geography)
intens_joint <- analog_density(
  x = sites,
  pool = climate_data,
  env = kernel("gaussian", theta = 0.2, max = 0.5),
  geog = kernel("gaussian", theta = 50, max = 100)
)

# With pre-built index (for repeated queries)
index <- build_analog_index(climate_data)
i1 <- analog_density(x = sites1, pool = index,
                     env = kernel("inverse", max = 0.5),
                     geog = kernel(max = 100))
i2 <- analog_density(x = sites2, pool = index,
                     env = kernel("gaussian", theta = 0.3),
                     geog = kernel(max = 100))
} # }
```
