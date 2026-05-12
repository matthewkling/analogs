# Climate impact assessment via analog impact model

Predicts ecological outcome variables based on climatic similarity and
geographic proximity. When focal and reference data are from the same
time period, performs a climate-informed spatial interpolation. When
they are from different eras, implements an analog impact model (AIM)
projecting the potential ecological state under climate change. For each
focal location's future climate, identifies locations with similar
baseline climates within a specified geographic range, then aggregates
their ecological characteristics weighted by climate similarity.
Aggregation can be a weighted mean for continuous outcomes or a weighted
class count for categorical outcomes.

## Usage

``` r
analog_impact(
  x,
  pool,
  y,
  weight = NULL,
  covariates = NULL,
  max_geog = NULL,
  max_clim = 1,
  kernel = c("gaussian_clim", "inverse_clim", "gaussian_joint", "inverse_joint"),
  theta = 0.25,
  stat = c("weighted_mean", "count", "sum_weights"),
  lambda = 0,
  se = c("none", "ess", "design"),
  normalize = "auto",
  x_cov = NULL,
  coord_type = "auto",
  index_res = "auto",
  cell_area_weight = "auto",
  n_threads = NULL,
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

- covariates:

  Optional matrix/data.frame or SpatRaster giving covariate values for
  each reference location (must have same number of rows/cells as
  `pool`). Required when `stat` includes `"regression"`.

- max_geog:

  Maximum geographic distance constraint (default: NULL = no geographic
  constraint). When specified, only reference locations within this
  distance are considered. Radius units should be specified in
  kilometers if `coord_type = "lonlat"`, or in projected coordinate
  units if `coord_type = "projected"`.

- max_clim:

  Maximum climate distance constraint (default: NULL = no climate
  constraint). Can be either:

  - A scalar: Euclidean radius in climate space (e.g., 0.5)

  - A vector: Per-variable absolute differences (length must equal
    number of climate variables)

  Only reference locations within this climate distance are considered.
  When `x_cov` is provided, scalar thresholds are interpreted in
  Mahalanobis distance units.

- kernel:

  Kernel decay function for weighting analogs during aggregation. Only
  weighting options that are based on *climate* are allowed:
  `"inverse_clim"` (default), `"gaussian_clim"`, `"inverse_joint"`,
  `"gaussian_joint"`. See
  [`analog_search()`](https://matthewkling.github.io/analogs/reference/analog_search.md)
  for details.

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

- stat:

  Statistic(s) to compute across analogs (default: c("count",
  "sum_weights", "weighted_mean")). See
  [`analog_search()`](https://matthewkling.github.io/analogs/reference/analog_search.md)
  for the full list of options. Common choices:

  - `"weighted_mean"`: Expected ecological state under continuous y

  - `"tabulate"`: Per-class similarity-weighted support under
    categorical y. Output has one column per class, named `n_<level>`
    (single-y) or `<varname>_n_<level>` (multi-y). `"tabulate"` is
    mutually exclusive with `weighted_mean`/`sum`/
    `mean`/`weighted_sum`/`regression`, but can be combined with
    `count`, `sum_weights`, `mean_weights`, and `ess`.

  - `"count"`: Analog availability (number of analogs)

  - `"sum_weights"`: Analog density (sum of climate-similarity weights
    across analogs)

  The default `c("count", "sum_weights", "weighted_mean")` is
  appropriate for continuous `y`; for categorical `y`, swap
  `weighted_mean` for `tabulate`.

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

- normalize:

  One of `TRUE`, `FALSE`, or `"auto"` (default). Only used if `stat`
  includes `"sum_weights"` or `"tabulate"` and `pool` is a raster. When
  active, results for these stats are divided by a global scalar so that
  they represent a fraction of a theoretically "perfect" scenario where
  the full search area within `max_geog` is occupied wall-to-wall by
  cells whose climate exactly matches `x`. See details under
  [`analog_search()`](https://matthewkling.github.io/analogs/reference/analog_search.md)
  for more info.

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

- progress:

  Logical; if `TRUE`, display a progress bar during computation.
  Progress tracking works by splitting the focal dataset into chunks and
  processing them sequentially. Useful for large datasets. Default is
  `FALSE`.

- ...:

  Additional arguments passed to
  [`analog_search()`](https://matthewkling.github.io/analogs/reference/analog_search.md).

## Value

A data.frame, or a SpatRaster when `x` is one. Contains `index`, `x`,
`y` plus one or more columns determined by `stat`. See
[`analog_search()`](https://matthewkling.github.io/analogs/reference/analog_search.md)
for column-naming conventions across stats and
[`metadata()`](https://matthewkling.github.io/analogs/reference/metadata.md)
for attached metadata attributes.

## Details

### The Analog Impact Model (AIM) Framework

When `x` represents future climate and `pool` represent baseline
climate, this function implements the "reverse analog" approach from the
climate change ecology literature. It addresses the question, "For a
location's future climate, what ecological conditions exist in current
locations with similar climates that are within dispersal range?"

The methodology:

1.  For each focal location's future climate conditions

2.  Find all current locations with similar climates (within `max_clim`)

3.  Constrain to dispersal-reachable distance (within `max_geog`)

4.  Weight each analog by climate similarity (via `kernel` function)

5.  Aggregate ecosystem characteristics across these weighted analogs

Unlike traditional AIM implementations that select k nearest climate
neighbors, this function uses all analogs within thresholds combined
with climate-distance-based kernel weighting. This approach eliminates
arbitrary choice of k, provides smoother, more continuous results, and
lets the kernel function (via `theta`) naturally control influence.
(Note that the traditional version can be implemented via
`analog_search(select = "knn_clim", stat = "mean", ...))`.)

### Choosing Parameters

- `max_geog`: Set based on species dispersal ability (e.g., 5-500 km)

- `max_clim`: Defines what counts as an "analog"

- `theta`: Controls kernel decay. The weight should decay to a small
  value at the `max_clim`/`max_geog` boundary. If `theta` is too large
  relative to thresholds, the hard cutoffs do most of the filtering and
  weighting becomes nearly uniform. For Gaussian kernels with three or
  fewer climate variables, a reasonable rule of thumb is to set `theta`
  to `max_* / 3`.

### Interpreting Results

- `weighted_mean`: Expected ecosystem state if colonized by species from
  analog locations.

- `tabulate`: For each class in `y`, the sum of analog weights that fall
  in that class. Each class's column gives the total
  climatic-similarity-weighted support among analogs (or, with
  `analog_search(stat = "tabulate", kernel = "uniform")`, a plain vote
  count). The "primary" projection at a focal location is
  [`which.max()`](https://rdrr.io/r/base/which.min.html) of these
  columns; "agreement" is the largest column value divided by the row
  sum (or, equivalently, divided by `sum_weights` when both are
  requested). Bray-Curtis dissimilarity between two such weighted-vote
  matrices (e.g., contemporary vs. future analogs) gives a measure of
  compositional vulnerability to ecological transformation (cf. Hoecker
  et al. 2026).

- `count`: How many analogs exist within max_clim and max_geog? Low
  counts indicate limited analog availability, while zero counts
  indicate climates that are novel within the geographic search radius.

- `sum_weights`: Total analog density. Low values indicate sparse or
  distant climate matches. This metric captures both the number and
  quality of analogs. Interpretation details vary based on the `kernel`
  parameter. Pass `normalize = TRUE` to express `sum_weights` (and any
  `tabulate` columns) as a fraction of theoretical maximum analog
  availability, on roughly `[0, 1]`. See
  [`analog_search()`](https://matthewkling.github.io/analogs/reference/analog_search.md)
  for preconditions.

## See also

[`analog_search()`](https://matthewkling.github.io/analogs/reference/analog_search.md)
for the underlying flexible analog search function;
[`analog_cv()`](https://matthewkling.github.io/analogs/reference/analog_cv.md)
for cross-validation of AIM fits.

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic climate impact assessment with continuous response
impact <- analog_impact(
  x = future_climate,
  pool = current_climate,
  y = current$habitat,
  max_geog = 100,    # 100 km dispersal range
  max_clim = 0.5     # Climate analog threshold
)

# With uncertainty quantification on weighted_mean
impact_se <- analog_impact(
  x = future_climate,
  pool = current_climate,
  y = current$habitat,
  max_geog = 100,
  max_clim = 0.5,
  se = "ess"
)

# Categorical response (e.g. vegetation type):
# compute weighted-vote counts per class within climatic neighborhoods
veg_votes <- analog_impact(
  x    = future_climate,
  pool = current_climate,
  y    = factor(current$vegetation_type),
  stat = c("count", "sum_weights", "tabulate"),
  kernel = "gaussian_clim",
  max_geog = 500,
  max_clim = 1
)
# Output has one `n_<level>` column per vegetation type.
# Primary class per focal:  apply(veg_votes[, grep("^n_", names(veg_votes))],
#                                   1, function(r) names(r)[which.max(r)])
# Agreement (max share):    apply(votes_mat, 1, max) / rowSums(votes_mat)
} # }
```
