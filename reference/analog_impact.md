# Climate impact assessment via analog impact model

Assesses potential climate change impacts using the analog impact model
(AIM) methodology. For each focal location's future climate, identifies
locations with similar baseline climates within a specified geographic
range, then aggregates their ecological characteristics weighted by
climate similarity. This quantifies what ecosystem conditions are likely
accessible via dispersal as climate changes.

## Usage

``` r
analog_impact(
  x,
  pool,
  y,
  covariates = NULL,
  max_geog = NULL,
  max_clim = 1,
  kernel = c("gaussian_clim", "inverse_clim", "gaussian_joint", "inverse_joint"),
  theta = 0.25,
  stat = c("count", "sum_weights", "weighted_mean"),
  lambda = 0,
  se = c("none", "ess", "design"),
  x_cov = NULL,
  coord_type = "auto",
  index_res = "auto",
  n_threads = NULL,
  progress = FALSE
)
```

## Arguments

- x:

  Focal locations (generally with future climate conditions). Should be
  a matrix/data.frame with columns x, y, and climate variables, or a
  SpatRaster with climate variable layers.

- pool:

  The reference dataset (generally representing baseline climate
  conditions). Either:

  - Matrix/data.frame with columns x, y, and climate variables, or
    SpatRaster with climate variable layers, OR

  - An `analog_index()` object created by
    [`build_analog_index()`](https://matthewkling.github.io/analogs/reference/build_analog_index.md)
    (for repeated queries).

- y:

  Ecological or environmental variable(s) for the same era as `pool`, to
  aggregate across climate analogs. Examples include occupancy of focal
  species, species richness, biomass, or any other ecological state
  variable. Can be a numeric vector (single variable), matrix or
  data.frame with numeric columns (multiple variables), or a SpatRaster
  with one or more numeric layers. Must have exactly the same number of
  reference locations as `pool`.

- covariates:

  Optional auxiliary predictor variables for each reference location in
  `pool`, used with `stat = "regression"`. Can be a numeric vector
  (single covariate), matrix or data.frame with numeric columns
  (multiple covariates), or a SpatRaster with one or more numeric
  layers. Must have exactly the same number of rows/cells as `pool`.
  Column names are used for output layer naming. These variables are NOT
  used for the analog search itself – only for local regression within
  each analog neighborhood.

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

- stat:

  Statistic(s) to compute across analogs (default: c("count",
  "sum_weights", "weighted_mean")). See
  [`analog_search()`](https://matthewkling.github.io/analogs/reference/analog_search.md)
  for options. The default statistics provide a complete picture:

  - `"count"`: Analog availability (number of analogs)

  - `"sum_weights"`: Analog intensity

  - `"weighted_mean"`: Expected ecological state

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
    For regression, the sandwich estimator
    `(X'WX + λI)⁻¹ M (X'WX + λI)⁻¹` with `M = Σ w² r² x xᵀ` (Huber-White
    / heteroskedasticity-consistent).

  Both variants are scale-invariant in the weights when `lambda = 0`.
  When `se != "none"` but no requested stat supports SE, a warning is
  issued and no SE columns are produced. SE-supporting stat columns are
  named `se_weighted_mean` / `se_weighted_mean_{varname}` and
  `se_{intercept|covname}` / `se_{intercept|covname}_{varname}`.

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

A data.frame with one row per focal location containing:

- `index`: Row number from input `x` data

- `x, y`: Coordinates of focal location

- One column per requested statistic

- For value statistics with multiple variables: `{stat}_{varname}`
  (e.g., `weighted_mean_habitat_quality`)

- For `"regression"`: `coef_intercept` and `coef_{covariate}`
  coefficient columns

- When `se != "none"`: SE columns for SE-supporting stats (e.g.,
  `se_weighted_mean`)

## Details

### The Analog Impact Model (AIM) Framework

This function implements the "reverse analog" approach from the climate
change ecology literature. It addresses the question, "For a location's
future climate, what ecological conditions exist in current locations
with similar climates that are within dispersal range?"

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

- `count`: How many analogs exist within max_clim and max_geog? Low
  counts indicate limited analog availability, while zero counts
  indicate climates that are novel within the geographic search radius.

- `sum_weights`: Total analog intensity. Low values indicate sparse or
  distant climate matches. This metric captures both the number and
  quality of analogs. Interpretation details vary based on the `kernel`
  parameter.

- `weighted_mean`: Expected ecosystem state if colonized by species from
  analog locations.

## See also

[`analog_search()`](https://matthewkling.github.io/analogs/reference/analog_search.md)
for the underlying flexible analog search function.

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic climate impact assessment
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
} # }
```
