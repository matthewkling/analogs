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
  covariates = NULL,
  max_geog = NULL,
  max_clim = 1,
  kernel = c("gaussian_clim", "inverse_clim", "gaussian_joint", "inverse_joint"),
  theta = 0.25,
  stat = c("weighted_mean", "count", "sum_weights"),
  lambda = 0,
  se = c("none", "ess", "design"),
  x_cov = NULL,
  coord_type = "auto",
  index_res = "auto",
  n_threads = NULL,
  progress = FALSE,
  ...
)
```

## Arguments

- x:

  Focal locations (usually future climate conditions for AIM, or
  baseline conditions for static interpolation). Should be a
  matrix/data.frame with columns x, y, and climate variables, or a
  SpatRaster with climate variable layers.

- pool:

  The reference dataset (generally representing baseline climate
  conditions). Either:

  - Matrix/data.frame with columns x, y, and climate variables, or
    SpatRaster with climate variable layers, OR

  - An `analog_index` object created by
    [`build_analog_index()`](https://matthewkling.github.io/analogs/reference/build_analog_index.md)
    (for repeated queries).

- y:

  Ecological or environmental variable(s) for the same era as `pool`, to
  aggregate across climate analogs. For continuous stats
  (`weighted_mean`, `sum`, etc.), `y` should be numeric: occupancy of
  focal species, species richness, biomass, or any other numeric
  ecological state. For `stat = "tabulate"` (categorical aggregation),
  `y` should be a factor or a vector that can be coerced to a factor
  (e.g., character vector of vegetation type names, integer class codes,
  or a single-layer categorical SpatRaster). Can be a vector / factor
  (single variable), matrix or data.frame (multiple variables), or a
  SpatRaster with one or more layers. Must have exactly the same number
  of reference locations as `pool`.

- covariates:

  Covariates to use with `stat = "regression"` (NULL otherwise; see
  [`analog_regression()`](https://matthewkling.github.io/analogs/reference/analog_regression.md)
  for the regression-focused wrapper).

- max_geog:

  Maximum geographic distance for analogs (km if lonlat; projection
  units otherwise). Acts as a hard dispersal constraint.

- max_clim:

  Maximum climate distance threshold defining what counts as an analog
  (default: 1.0).

- kernel:

  Kernel decay function for weighting analogs during aggregation. Only
  weighting options that are based on *climate* are allowed:
  `"inverse_clim"` (default), `"gaussian_clim"`, `"inverse_joint"`,
  `"gaussian_joint"`. See
  [`analog_search()`](https://matthewkling.github.io/analogs/reference/analog_search.md)
  for details.

- theta:

  Kernel bandwidth (default: 0.25). See
  [`analog_search()`](https://matthewkling.github.io/analogs/reference/analog_search.md).

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

  - `"sum_weights"`: Analog intensity (sum of climate-similarity weights
    across analogs)

  The default `c("count", "sum_weights", "weighted_mean")` is
  appropriate for continuous `y`; for categorical `y`, swap
  `weighted_mean` for `tabulate`.

- lambda:

  Ridge penalty for regression (default 0).

- se:

  One of "none" (default), "ess", or "design"; SE variant for
  `weighted_mean` and `regression`.

- x_cov:

  Optional covariance matrices for Mahalanobis distance (one per focal
  point), as expected by
  [`analog_search()`](https://matthewkling.github.io/analogs/reference/analog_search.md).

- coord_type:

  One of "auto", "lonlat", or "projected".

- index_res:

  Resolution parameter for the lattice index (or "auto").

- n_threads:

  Optional integer number of threads.

- progress:

  Logical; show a progress bar.

- ...:

  Additional arguments passed to
  [`analog_search()`](https://matthewkling.github.io/analogs/reference/analog_search.md).

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

- `sum_weights`: Total analog intensity. Low values indicate sparse or
  distant climate matches. This metric captures both the number and
  quality of analogs. Interpretation details vary based on the `kernel`
  parameter.

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
