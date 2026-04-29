# Generalized analog search

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
  y = NULL,
  covariates = NULL,
  max_clim = NULL,
  max_geog = NULL,
  select = "all",
  k = NULL,
  stat = NULL,
  kernel = NULL,
  theta = NULL,
  lambda = 0,
  se = c("none", "ess", "design"),
  exclude_self = FALSE,
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

- y:

  Optional vector, factor, matrix/data.frame, or SpatRaster giving
  values for each reference location (must have same number of
  rows/cells as `pool`). Required for stats `"sum"`, `"mean"`,
  `"weighted_sum"`, `"weighted_mean"`, `"regression"`, and `"tabulate"`.
  Numeric for continuous stats; factor or coercible-to-factor
  (character, integer, logical) for `stat = "tabulate"`.

- covariates:

  Optional matrix/data.frame or SpatRaster giving covariate values for
  each reference location (must have same number of rows/cells as
  `pool`). Required when `stat` includes `"regression"`.

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
    (see `kernel` and `theta`).

  - `"mean_weights"`: For each focal, mean of weights of selected
    analogs.

  - `"sum"`: Sum of `y` values across analogs (requires `y`).

  - `"mean"`: Mean of `y` values across analogs (requires `y`).

  - `"weighted_sum"`: Sum of (`y` × kernel weight) across analogs
    (requires `y` and `kernel`).

  - `"weighted_mean"`: Weighted mean of `y` values across analogs
    (requires `y` and `kernel`).

  - `"ess"`: Kish's effective sample size (ESS), computed as the squared
    sum of weights divided by the sum of squared weights (requires
    `kernel`).

  - `"regression"`: Weighted least squares (or ridge) regression of `y`
    on `covariates` within each analog neighborhood. Returns intercept
    and slope coefficients. Requires `y`, `covariates`, and `kernel`.
    See `lambda` for regularization.

  - `"tabulate"`: For each focal and each level of categorical `y`, sum
    the kernel weights of analogs whose `y` falls in that level. With
    `kernel = "uniform"` this reduces to a per-class vote count; with a
    distance-decay kernel it gives similarity-weighted support per
    class. Requires `y` (factor or coercible-to-factor) and `kernel`.
    Output has one column per class. `"tabulate"` is mutually exclusive
    with `"sum"`, `"mean"`, `"weighted_sum"`, `"weighted_mean"`, and
    `"regression"` (different `y` semantics); it can be combined with
    `"count"`, `"sum_weights"`, `"mean_weights"`, and `"ess"`.

  - A character vector combining multiple stats (e.g.,
    `c("count", "weighted_mean", "regression")`). Note: `"none"` cannot
    be combined with other stats.

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

- lambda:

  Ridge penalty parameter for `stat = "regression"` (default: 0, giving
  ordinary weighted least squares). Higher values shrink covariate
  coefficients toward zero, with the intercept approaching the weighted
  mean as `lambda -> Inf`. Ignored when `"regression"` is not in `stat`.

- se:

  Standard-error framing to apply to SE-supporting stats
  (`"weighted_mean"` and `"regression"`). One of:

  - `"none"` (default): no SE columns are returned.

  - `"ess"`: effective-sample-size framing. For `weighted_mean`,
    `SE = sqrt(var_w(y) / n_eff)`, where `n_eff = (Σw)² / Σw²` is Kish's
    effective sample size and `var_w(y) = Σwy²/Σw - ȳ_w²`. For
    regression, `Var(β̂) = σ²_ess · (X'WX + λI)⁻¹`, with residual
    variance corrected using `n_eff - p` degrees of freedom.

  - `"design"`: design-based framing (no assumption that weights are
    precisions). For `weighted_mean`, `SE = sqrt(Σ w²(y - ȳ_w)²) / Σw`.

- exclude_self:

  Logical, default `FALSE`. `TRUE` is typically used for
  cross-validation, such as via
  [`analog_cv()`](https://matthewkling.github.io/analogs/reference/analog_cv.md),
  in which case each focal excludes the pool row at the same index from
  its analog neighborhood. This requires `x` and `pool` to be the same R
  object (checked via
  [`identical()`](https://rdrr.io/r/base/identical.html)), and is
  incompatible with pre-built indices, `downsample != 1`, and
  `progress = TRUE`.

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

- For `stat = "tabulate"`: one column per level of `y`, named
  `n_{level}` for a single unnamed `y`, or `{varname}_n_{level}` when
  `y` is named or has multiple columns.

- When `se != "none"`: matching SE columns (`se_weighted_mean`,
  `se_intercept`, etc.) for each SE-supporting stat.

All results include metadata attributes (`select`, `stat`, `kernel`,
etc.). Use
[`metadata()`](https://matthewkling.github.io/analogs/reference/metadata.md)
to retrieve them as a named list, or see
[`?metadata`](https://matthewkling.github.io/analogs/reference/metadata.md)
for a full reference.

## Details

### Parameter categories

- *Data parameters* (`x`, `pool`, `x_cov`, `y`, `covariates`,
  `coord_type`) give attributes of the data on which to operate.

- *Selection parameters* (`select`, `max_clim`, `max_geog`, `k`) define
  which analogs to `select` from the `pool` for each `x`.

- *Aggregation parameters* (`stat`, `kernel`, `theta`, `lambda`, `se`)
  control how selected analogs are summarized.

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

- Parallel processing is available via the `n_threads` parameter.

- You can `downsample` prohibitively large reference pool datasets to
  improve speed and memory, using a stratified sampling scheme that
  reduces loss of precision relative to random sampling.

- For large datasets, enable `progress = TRUE` to display a progress bar
  during computation.

- For raster datasets that are too large to fit in memory,
  [`tiled_analog_search()`](https://matthewkling.github.io/analogs/reference/tiled_analog_search.md)
  offers a memory-safe option.

### Cross-validation

For honest prediction error when `x` and `pool` are the same dataset,
use
[`analog_cv()`](https://matthewkling.github.io/analogs/reference/analog_cv.md)
or set `exclude_self = TRUE` to exclude each focal's own row from its
analog neighborhood.

## References

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
[`analog_cv()`](https://matthewkling.github.io/analogs/reference/analog_cv.md)
provides cross-validation workflows.
