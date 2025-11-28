# Find Climate Analogs

Identifies locations in a reference dataset that are climatically
similar to focal locations, with optional constraints on climate
distance and geographic distance. This function supports multiple use
cases including climate velocity analysis, analog availability mapping,
and climate impact assessment.

## Usage

``` r
analog_search(
  x,
  pool,
  mode,
  max_clim = NULL,
  max_geog = NULL,
  x_cov = NULL,
  k = NULL,
  weight = NULL,
  theta = NULL,
  report_dist = TRUE,
  coord_type = c("auto", "lonlat", "projected"),
  index_res = "auto",
  resolutions = NULL,
  n_sample = NULL,
  n_reps = NULL,
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

  - An `analog_index` object created by
    [`build_analog_index`](https://matthewkling.github.io/analogs/reference/build_analog_index.md)
    (for repeated queries).

- mode:

  Character string specifying the analog search mode. One of:

  - `"knn_clim"`: For each focal, return up to `k` analogs with smallest
    climate distance, subject to `max_clim` and `max_geog` filters.

  - `"knn_geog"`: For each focal, return up to `k` analogs with smallest
    geographic distance, subject to `max_clim` and `max_geog` filters.

  - `"all"`: Return all analogs that satisfy the filters.

  - `"count"`: For each focal, count how many analogs satisfy the
    filters.

  - `"sum"`: For each focal, sum weights of all analogs that satisfy the
    filters (see `weight` and `theta`).

  - `"mean"`: For each focal, mean of weights of all analogs that
    satisfy the filters.

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

- x_cov:

  Optional focal-specific covariance matrices for Mahalanobis distance
  calculations. Should be a matrix or data.frame with one row per focal
  location and one column per unique covariance component. For n climate
  variables, there are n\*(n+1)/2 unique components, ordered as:
  variances first (diagonals), then covariances (upper triangle by row).
  For example:

  - 2 variables: c(var1, var2, cov12)

  - 3 variables: c(var1, var2, var3, cov12, cov13, cov23)

  When provided, all climate distances are computed as Mahalanobis
  distances using each focal's covariance structure. For focal-specific
  variances only (no covariance), set off-diagonal covariances to zero.
  Default is NULL (Euclidean climate distance).

- k:

  Number of nearest analogs to return per focal location for kNN modes.
  Required when `mode` is `"knn_geog"` or `"knn_clim"`; must be `NULL`
  for other modes.

- weight:

  Weighting function for matches, used only when `mode` is `"sum"` or
  `"mean"`. One of:

  - `"uniform"`: All matches weighted equally (weight = 1.0).

  - `"inverse_clim"`: Weight = 1 / (climate_distance + epsilon), with
    epsilon given by `theta` (or a small default if `theta` is `NULL`).

  - `"inverse_geog"`: Weight = 1 / (geographic_distance + epsilon), with
    epsilon given by `theta` (or a small default if `theta` is `NULL`).

  - `"gaussian_clim"`: Gaussian kernel on climate distance, weight =
    exp(-climate_distance^2 / (2\*sigma^2)), with sigma (bandwidth)
    given by `theta`.

  - `"gaussian_geog"`: Gaussian kernel on geographic distance, weight =
    exp(-geographic_distance^2 / (2\*sigma^2)), with sigma (bandwidth)
    given by `theta`.

  - `"gaussian_joint"`: Bivariate Gaussian kernel (product of
    independent Gaussians over climate and geographic distances), with
    bandwidths given by `theta` as a 2-element vector c(sigma_clim,
    sigma_geog).

  - `"inverse_joint"`: Joint inverse distance, weight = 1 /
    sqrt((climate_distance + eps_clim)^2 + (geographic_distance +
    eps_geog)^2), with epsilon values given by `theta` as a 2-element
    vector c(eps_clim, eps_geog).

  For `mode` in `"knn_geog"`, `"knn_clim"`, `"count"`, or `"all"`,
  `weight` must be `NULL`.

- theta:

  Optional numeric parameter used by weighting functions when `mode` is
  `"sum"` or `"mean"` and `weight` is not `"uniform"`. Interpretation
  depends on `weight`:

  - For `"inverse_clim"` or `"inverse_geog"`: epsilon value added to
    distances (scalar; default: 1e-12 for climate, 1e-6 for geography).

  - For `"gaussian_clim"` or `"gaussian_geog"`: sigma bandwidth
    parameter (scalar; larger values = slower decay with distance).

  - For `"gaussian_joint"`: 2-element vector c(sigma_clim, sigma_geog)
    giving bandwidths for climate and geographic dimensions.

  - For `"inverse_joint"`: 2-element vector c(eps_clim, eps_geog) giving
    epsilon values for climate and geographic dimensions.

  If `theta` is `NULL`, sensible defaults are used for single-parameter
  weights. For `weight = "uniform"` or for non-aggregating modes,
  `theta` must be `NULL`.

- report_dist:

  Logical; if TRUE (default), include distance columns in output when
  `mode` is `"knn_geog"`, `"knn_clim"` or `"all"`. Set to FALSE for more
  compact output.

- coord_type:

  Coordinate system type (default: "auto"):

  - `"auto"`: Automatically detect from coordinate ranges.

  - `"lonlat"`: Unprojected lon/lat coordinates (uses great-circle
    distance; assumes `max_geog` is in km).

  - `"projected"`: Projected XY coordinates (uses planar distance;
    assumes `max_geog` is in projection units).

- index_res:

  Tuning parameter giving the number of bins per dimension of the
  internally-used lattice search index. Either:

  - A positive integer.

  - `"auto"` (the default): Automatically tune the intex resolution by
    optimizing compute time on a subsample of focal points. If focal has
    relatively few rows, auto-tuning is skipped and a default resolution
    of 16 is used.

  Ignored if `pool` is an `analog_index` (uses index's resolution).

- n_threads:

  Optional integer number of threads to use for the computation. If
  `NULL` (default), the global RcppParallel setting is used (see
  [`RcppParallel::setThreadOptions`](https://rdrr.io/pkg/RcppParallel/man/setThreadOptions.html)).

## Value

The return value depends on the `mode` parameter:

\*\*For mode = "knn_geog", "knn_clim" or "all"\*\*: A data.frame with
one row per focal-analog pair:

- `focal_index`: Index of focal location (1-based).

- `focal_x, focal_y`: Coordinates of focal location.

- `analog_index`: Index of analog location in reference dataset
  (1-based).

- `analog_x, analog_y`: Coordinates of analog location.

- `clim_dist`: Climate distance (if `report_dist = TRUE`).

- `geog_dist`: Geographic distance in km (if `report_dist = TRUE`).

\*\*For mode = "sum", "mean", or "count"\*\*: A data.frame with one row
per focal location:

- `focal_index`: Index of focal location (1-based).

- `focal_x, focal_y`: Coordinates of focal location.

- `value`: Aggregated value (count, sum of weights, or mean of weights).

All outputs include diagnostic attributes propagated from the C++ core,
including:

- `total_bins`: Number of spatial bins in the lattice index.

- `avg_bin_occupancy`: Average points per bin.

- `min_bin_occupancy, max_bin_occupancy`: Range of bin occupancy.

- `binning_method`: Method used ("multi_dim_lattice" or "none").

- `n_ref, n_clim`: Size of reference dataset and number of climate
  variables.

## Details

The function uses a spatial indexing structure (lattice-based) to
quickly search through large reference datasets. Climate similarity is
measured using Euclidean distance in climate space (ideally
pre-whitened; see Details). Geographic distance can be computed for
lon/lat coordinates (great-circle distance) or projected coordinates
(planar distance).

\*\*Common Use Cases:\*\*

**Climate Velocity** (nearest geographic neighbor with similar climate):

    analog_search(
      x        = clim$clim1,
      pool     = clim$clim2,
      mode     = "knn_geog",
      max_clim = 0.5,
      max_geog = NULL,
      k        = 1
    )

**Climate Impact** (climatically similar locations within dispersal
range):

    analog_search(
      x        = clim$clim1,
      pool     = clim$clim2,
      mode     = "knn_clim",
      max_clim = 0.5,
      max_geog = 100,
      k        = 20
    )

**Analog Availability** (count of suitable locations):

    analog_search(
      x        = clim$clim1,
      pool     = clim$clim1,
      mode     = "count",
      max_clim = 0.5,
      max_geog = 100
    )

**Weighted Analog Intensity** (e.g., distance-weighted availability):

    analog_search(
      x        = clim$clim1,
      pool     = clim$clim1,
      mode     = "sum",
      max_clim = 0.5,
      max_geog = 100,
      weight   = "inverse_geog",
      theta    = 1e-6
    )

**Using a Pre-built Index** (for repeated queries):

    # Build index once
    index <- build_analog_index(clim$clim2, index_res = 16)

    # Query multiple times
    v1 <- analog_search(x = sites1, pool = index, mode = "knn_geog", max_clim = 0.5, k = 1)
    v2 <- analog_search(x = sites2, pool = index, mode = "knn_geog", max_clim = 0.3, k = 1)

**Focal-specific Mahalanobis Distance**:

    # With focal-specific covariance matrices
    analog_search(
      x        = clim$clim1,
      pool     = clim$clim2,
      x_cov    = focal_covariances,  # n_focal x 3 matrix for 2 climate vars
      mode     = "knn_geog",
      max_clim = 2,  # In Mahalanobis distance units
      k        = 1
    )
