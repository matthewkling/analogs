# Climate velocity: geographically nearest climate analogs

Finds, for each focal location, the geographic nearest neighbor(s) in a
reference dataset that satisfy a specified maximum environmental
distance threshold. Distances to these analogs, divided by time elapsed,
give analog-based climate velocity (Hamann et al. 2015; Dobrowski and
Parks 2016).

## Usage

``` r
analog_velocity(
  x,
  pool,
  x_cov = NULL,
  y = NULL,
  weight = NULL,
  coord_type = "auto",
  env,
  geog = NULL,
  k = 1,
  env_res_adj = "auto",
  geog_res_adj = "auto",
  cell_area_weight = "auto",
  n_threads = NULL,
  downsample = 1,
  seed = NULL,
  progress = FALSE
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

- y:

  Optional vector, factor, matrix/data.frame, or SpatRaster giving
  values for each reference location (must have same number of
  rows/cells as `pool`). Required for stats `"sum"`, `"mean"`,
  `"weighted_sum"`, `"weighted_mean"`, `"regression"`, and `"tabulate"`.
  Numeric for continuous stats; factor or coercible-to-factor
  (character, integer, logical) for `stat = "tabulate"`.

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

- k:

  Number of nearest analogs to return per focal location for kNN
  selection modes. Required when `select` is `"knn_geog"` or
  `"knn_env"`; must be `NULL` for `select = "all"`.

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

## Value

A data.frame, or a SpatRaster when `x` is one and `k = 1`. Contains one
row per focal-analog pair with `index`, `x`, `y`, `analog_index`,
`analog_x`, `analog_y`, `env_dist`, and `geog_dist`. See
[`analog_search()`](https://matthewkling.github.io/analogs/reference/analog_search.md)
for full column conventions and
[`metadata()`](https://matthewkling.github.io/analogs/reference/metadata.md)
for attached metadata attributes.

## Details

This function is a wrapper that calls
[`analog_search()`](https://matthewkling.github.io/analogs/reference/analog_search.md)
using `select = "knn_geog"` and `stat = "none"`. Note that it **does not
return velocity per se**—it returns geographic and environmental
distances to each focal site's nearest analog(s); to compute velocity,
you can divide these geographic distances by the length of time elapsed
between your `x` and `pool` datasets.

## References

Hamann A, Roberts DR, Barber QE, Carroll C, Nielsen SE (2015). "Velocity
of climate change algorithms for guiding conservation and management."
*Global Change Biology*, **21**(2), 997-1004.
[doi:10.1111/gcb.12736](https://doi.org/10.1111/gcb.12736)

Dobrowski SZ, Parks SA (2016). "Climate change velocity underestimates
climate change exposure in mountainous regions." *Nature
Communications*, **7**, 12349.
[doi:10.1038/ncomms12349](https://doi.org/10.1038/ncomms12349)

## See also

[`analog_search()`](https://matthewkling.github.io/analogs/reference/analog_search.md)
for the underlying flexible analog search function;
[`tiled_analog_search()`](https://matthewkling.github.io/analogs/reference/tiled_analog_search.md)
for memory-safe searches on large raster datasets.

## Examples

``` r
if (FALSE) { # \dontrun{
# One-shot query
v <- analog_velocity(
  x = clim$clim1,
  pool = clim$clim2,
  env = kernel(max = 0.5),
  k = 1
)

# With pre-built index (for repeated queries)
index <- build_analog_index(clim$clim2)
v1 <- analog_velocity(x = sites1, pool = index, env = kernel(max = 0.5), k = 1)
v2 <- analog_velocity(x = sites2, pool = index, env = kernel(max = 0.3), k = 1)

# With focal-specific covariance matrices
v_mahal <- analog_velocity(
  x = clim$clim1,
  pool = clim$clim2,
  x_cov = baseline_covariances,
  env = kernel(max = 2),  # In Mahalanobis distance units
  k = 1
)
} # }
```
