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
  values,
  max_geog = NULL,
  max_clim = 1,
  weight = c("gaussian_clim", "inverse_clim", "gaussian_joint", "inverse_joint"),
  theta = 0.25,
  stat = c("count", "sum_weights", "weighted_mean"),
  x_cov = NULL,
  coord_type = "auto",
  index_res = "auto",
  n_threads = NULL
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

- values:

  Ecological or environmental variable(s) for the same era as `pool`, to
  aggregate across climate analogs. Must have exactly `nrow(pool)` rows.
  Examples include occupancy of focal species, species richness,
  biomass, or any other ecological state variable.

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

- weight:

  Function for weighting analogs during aggregation. Only weight options
  that are based on *climate* are allowed: `"inverse_clim"` (default),
  `"gaussian_clim"`, `"inverse_joint"`, `"gaussian_joint"`. See
  [`analog_search()`](https://matthewkling.github.io/analogs/reference/analog_search.md)
  for details.

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

- stat:

  Statistic(s) to compute across analogs (default: c("count",
  "sum_weights", "weighted_mean")). See
  [`analog_search()`](https://matthewkling.github.io/analogs/reference/analog_search.md)
  for options. The default statistics provide a complete picture:

  - `"count"`: Analog availability (number of analogs)

  - `"sum_weights"`: Analog intensity

  - `"weighted_mean"`: Expected ecological state

- x_cov:

  Optional focal-specific covariance matrices for Mahalanobis distance
  calculations. Should be a matrix or data.frame with one row per focal
  location and one column per unique covariance component. For n climate
  variables, there are n\*(n+1)/2 unique components, ordered as:
  variances first (diagonals), then covariances (upper triangle by row).

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

## Value

A data.frame with one row per focal location containing:

- `index`: Row number from input `x` data

- `x, y`: Coordinates of focal location

- One column per requested statistic

- For value statistics with multiple variables: `{stat}_{varname}`
  (e.g., `weighted_mean_habitat_quality`)

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

4.  Weight each analog by climate similarity (via `weight` function)

5.  Aggregate ecosystem characteristics across these weighted analogs

Unlike traditional AIM implementations that select k nearest climate
neighbors, this function uses all analogs within thresholds combined
with climate-distance-based weighting. This approach eliminates
arbitrary choice of k, provides smoother, more continuous results, and
lets the weight function (via `theta`) naturally control influence.
(Note that the traditional version can be implemented via
`analog_search(select = "knn_clim", stat = "mean", ...))`.)

### Choosing Parameters

- `max_geog`: Set based on species dispersal ability (e.g., 5-500 km)

- `max_clim`: Defines what counts as an "analog"

- `theta`: Controls weight decay. The weight should decay to a small
  value at the `max_clim`/`max_geog` boundary. If `theta` is too large
  relative to thresholds, the hard cutoffs do most of the filtering and
  weighting becomes nearly uniform. For Gaussian weights with three or
  fewer climate variables, a reasonable rule of thumb is to set `theta`
  to `max_* / 3`.

### Interpreting Results

- `count`: How many analogs exist within max_clim and max_geog? Low
  counts indicate limited analog availability, while zero counts
  indicate climates that are novel within the geographic search radius.

- `sum_weights`: Total analog intensity. Low values indicate sparse or
  distant climate matches. This metric captures both the number and
  quality of analogs. Interpretation details vary based on the `weight`
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
  values = current$habitat,
  max_geog = 100,    # 100 km dispersal range
  max_clim = 0.5     # Climate analog threshold
)

# Interpret results:
# - count: how many analogs?
# - sum_weights: analog availability (higher = more/better matches)
# - weighted_mean: expected habitat quality (higher = lower impact)

# Multiple ecosystem variables
values_df <- data.frame(
  habitat_quality = current$habitat,
  species_richness = current$richness,
  forest_cover = current$forest
)
impact_multi <- analog_impact(
  x = future_climate,
  pool = current_climate,
  values = values_df,
  max_geog = 150,
  max_clim = 0.4,
  theta = 0.25
)

# Custom statistics and weighting
impact_custom <- analog_impact(
  x = future_climate,
  pool = current_climate,
  values = current$biomass,
  stat = c("count", "weighted_mean", "weighted_sum"),
  max_clim = 0.6,
  max_geog = 200,
  weight = "gaussian_joint",    # Weight by both climate and geography
  theta = c(0.2, 50)           # Climate and geographic decay
)

# With pre-built index for multiple scenarios
current_index <- build_analog_index(current_climate)

impact_current <- analog_impact(current_climate, current_index,
                                 values = current$quality,
                                 max_geog = 100)
impact_ssp126 <- analog_impact(future_ssp126, current_index,
                                 values = current$quality,
                                 max_geog = 100)
impact_ssp585 <- analog_impact(future_ssp585, current_index,
                                 values = current$quality,
                                 max_geog = 100)
} # }
```
