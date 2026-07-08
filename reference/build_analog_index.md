# Build Analog Index

Pre-builds a reusable lattice index from reference environmental data.
The index can be queried multiple times with different focal points and
parameters, avoiding the need to rebuild the lattice for each query.

## Usage

``` r
build_analog_index(
  pool,
  coord_type = c("auto", "lonlat", "projected"),
  geog_res_adj = 1,
  env_res_adj = 1,
  downsample = 1,
  seed = NULL,
  cell_area_weight = "auto",
  mean_cell_area = NULL
)
```

## Arguments

- pool:

  The reference dataset to search for analogs. Should be a
  matrix/data.frame with columns x, y, and environmental variables, or a
  SpatRaster with environmental variable layers.

- coord_type:

  Coordinate system type:

  - `"auto"` (default): Automatically detect from coordinate ranges.

  - `"lonlat"`: Unprojected lon/lat coordinates (uses great-circle
    distance; assumes `max_geog` is in km).

  - `"projected"`: Projected XY coordinates (uses planar distance;
    assumes `max_geog` is in projection units).

- env_res_adj, geog_res_adj:

  Non-negative scalars adjusting the lattice resolution of the
  environmental and geographic variable families, each as a multiplier
  on a data-dependent default. The default allocation targets an average
  of roughly 50 pool points per occupied bin and splits that budget
  between the two families by their effective dimensionality;
  `env_res_adj = geog_res_adj = 1` (the default) uses that allocation.
  Larger values give a finer lattice for that family (more, smaller
  bins), smaller values a coarser one, and `0` deactivates the family
  entirely (one bin per axis), which is appropriate when a query does
  not constrain that family. Because the default scales with pool size,
  these adjustments travel across datasets of different sizes. Within
  each family, bins are allocated in proportion to each axis' data range
  so a search radius spans a comparable number of bins on every axis.
  [`analog_search()`](https://matthewkling.github.io/analogs/reference/analog_search.md)
  sets these automatically from the query's constraints (deactivating an
  unconstrained family). Ignored if `pool` is an `analog_index`.

- downsample:

  Optional downsampling rate (0-1) indicating the proportion of points
  in `pool` to retain. Downsampling reduces memory use and improves
  query speed at the cost of some precision; adaptive stratified
  sampling is used to minimize loss of precision. The default is 1.0 (no
  downsampling). See Details for more info.

- seed:

  Optional random seed for reproducible downsampling. If `NULL`
  (default), uses current R random state.

- cell_area_weight:

  Controls cell-area weighting for raster pools. One of:

  - `"auto"` (default): Compute cell-area weights when `pool` is a
    SpatRaster, and skip them otherwise. This corrects aggregation
    statistics for non-uniform cell areas (e.g. lonlat grids where cell
    area shrinks toward the poles, or projected grids on non-equal-area
    projections).

  - `TRUE`: Force cell-area weighting on. Errors if `pool` is not a
    SpatRaster.

  - `FALSE`: Force cell-area weighting off; treat all pool points as
    having equal weight.

  - A numeric vector of length `nrow(pool)`: Use these caller-supplied
    weights as-is, without any further normalization. This is intended
    for advanced workflows like
    [`tiled_analog_search()`](https://matthewkling.github.io/analogs/reference/tiled_analog_search.md)
    that need to maintain a globally consistent normalization across
    multiple per-tile index builds; most users should use one of the
    three options above.

  When `"auto"` or `TRUE` triggers computation, weights are computed via
  [`terra::cellSize()`](https://rspatial.github.io/terra/reference/cellSize.html)
  and normalized to mean 1 over finite values, so absolute magnitudes of
  stats like `sum_weights` remain comparable to the unweighted case. The
  weights are stored on the returned index and used during all
  subsequent queries.

- mean_cell_area:

  Optional scalar mean cell area (in km^2) to attach to the index,
  overriding any value auto-computed from the raster pool. Intended for
  internal use by
  [`tiled_analog_search()`](https://matthewkling.github.io/analogs/reference/tiled_analog_search.md)
  to propagate a globally-consistent mean area across per-tile index
  builds (so that `analog_density(normalize = TRUE)` produces consistent
  values across tiles). Most users should leave this `NULL`.

## Value

An S3 object of class `"analog_index"` containing:

- The compiled lattice index (internal C++ structure)

- Reference data

- Metadata: coordinate type, dimensions, ranges, resolution

- Diagnostics: bin counts, occupancy statistics, and downsampling info

## Details

The lattice index is a multidimensional grid of bins, built over both
geographic and environmental dimensions. This structure enables
efficient analog searches by first filtering and sorting bins of similar
points before computing exact results. For lon/lat coordinates, the
index uses ECEF (Earth-Centered Earth-Fixed) space internally for
optimal performance.

Index resolution is controlled per family by `geog_res_adj` and
`env_res_adj`, each a multiplier on a pool-size-dependent default
(targeting ~50 points per bin, split between the families by effective
dimensionality). `1` uses the default for that family, larger values are
finer, smaller are coarser, and `0` deactivates a family (one bin per
axis). The optimal values depend on your data and query patterns;
[`analog_search()`](https://matthewkling.github.io/analogs/reference/analog_search.md)
sets them automatically from the query's constraints, or accept the
defaults of `1`, which work well for many applications.

### Downsampling

For very large datasets, downsampling can significantly improve memory
usage and query speed, at the cost of some precision. The `downsample`
parameter controls the target fraction of the data points in `pool` that
are retained in the index. Downsampling uses an adaptive stratified
approach: densely-packed bins are thinned more aggressively while sparse
bins are preserved, which helps reduce imprecision in sparse regions
compared to fully random sampling. Note: The actual rate may be higher
than requested if maintaining at least one point per occupied bin
requires it (common with sparse data or fine-grained binning); check
`index$downsample_actual`.

Each remaining analog in the downsampled pool gets a `sample_weight`
indicating the number of points it represents in the original pool; this
weight is the inverse of the sampling rate in the analog's index bin.
For pair queries (`stat = "none"`), results include each analog's
`sample_weight`. For aggregation stats (count, sum, mean, etc.),
sampling weights are used internally to automatically correct for the
downsampling bias.

## Examples

``` r
if (FALSE) { # \dontrun{
# Build index with default settings
index <- build_analog_index(climate_data)

# Build with finer geographic resolution than default
index <- build_analog_index(climate_data, geog_res_adj = 2)

# Build with downsampling for large datasets
index <- build_analog_index(
  large_climate_data,
  geog_res_adj = 1, env_res_adj = 1,
  downsample = 0.1,  # Reduce max bin size to 10%
  seed = 123         # Reproducible sampling
)

# Query the index multiple times
v1 <- analog_velocity(sites1, pool = index, env = kernel(max = 0.5))
v2 <- analog_velocity(sites2, pool = index, env = kernel(max = 0.3))
a1 <- analog_availability(sites3, pool = index, env = kernel(max = 0.5), geog = kernel(max = 100))
} # }
```
