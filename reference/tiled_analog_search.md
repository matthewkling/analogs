# Tiled analog search for memory-constrained queries

Performs analog searches on large datasets by dividing the focal region
into tiles and processing each tile separately. This reduces memory
usage, and with `method = "fft"` it also keeps each tile's FFT grid
small. Works with any `analog_*()` function.

## Usage

``` r
tiled_analog_search(
  x,
  pool,
  n_tiles,
  fun,
  geog,
  env = NULL,
  y = NULL,
  x_cov = NULL,
  weight = NULL,
  cell_area_weight = "auto",
  coord_type = c("auto", "lonlat", "projected"),
  ...,
  output_file = NULL,
  progress = TRUE
)
```

## Arguments

- x:

  Focal locations: a SpatRaster with environmental variable layers, or a
  matrix/data.frame with columns x, y, and environmental variables.

- pool:

  Reference locations (potential analog pool): a SpatRaster, or a
  matrix/data.frame with columns x, y, and environmental variables. An
  `analog_index` is not accepted, since each tile needs its own pool
  subset. Point tables carry no CRS, so `x` and `pool` must share a
  coordinate system whenever either is a point table (not checked), and
  `coord_type` should usually be given explicitly.

- n_tiles:

  Approximate number of tiles. The function will find a grid close to
  this number that creates square-ish tiles. Choosing larger values for
  n_tiles will reduce memory usage, but will also reduce computational
  efficiency. Choose the smallest n_tiles that fits your memory
  constraints.

- fun:

  An analog\_\* function to apply to each tile (e.g., analog_velocity,
  analog_impact).

- geog:

  A
  [`kernel()`](https://matthewkling.github.io/analogs/reference/kernel.md)
  object giving the geographic search treatment. Its `max` (the
  geographic search radius) is required and finite: it both sizes the
  tile halos (so edge focals still see all in-range analogs) and is
  forwarded to `fun`. Tiling is only beneficial when the geographic
  search is bounded, hence `geog` with a finite `max` is required.

- env:

  Optional
  [`kernel()`](https://matthewkling.github.io/analogs/reference/kernel.md)
  object giving the environmental search treatment, forwarded to `fun`.
  `NULL` (default) means no environmental kernel is passed.

- y:

  Optional values to aggregate across analogs. A SpatRaster matching
  `pool`'s CRS and extent when `pool` is a raster; a vector, factor, or
  matrix/data.frame with one row per row of `pool` when `pool` is a
  point table.

- x_cov:

  Optional focal covariances. A SpatRaster matching `x`'s CRS and extent
  when `x` is a raster; a matrix/data.frame with one row per row of `x`
  when `x` is a point table.

- weight:

  Optional per-pool-site weights. A single-layer SpatRaster matching
  `pool`'s CRS and extent when `pool` is a raster; a numeric vector with
  one value per row of `pool` when `pool` is a point table. See
  [`analog_search()`](https://matthewkling.github.io/analogs/reference/analog_search.md)
  for details.

- cell_area_weight:

  Controls cell-area weighting. One of `"auto"` (default; on for raster
  pools, off for point pools), `TRUE`, or `FALSE`. `TRUE` requires a
  raster pool. Unlike a non-tiled query, where weights are normalized to
  mean 1 over whatever pool is passed, here weights are normalized to
  mean 1 over the *full pool raster* once at the outer level, and the
  resulting raster is sliced per tile. This keeps absolute magnitudes of
  weighted stats (e.g. `sum_weights`) consistent across tiles.

- coord_type:

  One of `"auto"` (default), `"lonlat"`, or `"projected"`; see
  [`analog_search()`](https://matthewkling.github.io/analogs/reference/analog_search.md).
  `"auto"` reads the CRS of whichever of `x` and `pool` is a raster, and
  falls back to a coordinate-range check only when both are point
  tables. That check reads small projected coordinates (any within
  +/-180, +/-90) as lon/lat, so set this explicitly for point tables.
  Required to be `"projected"` for geographic kernels with
  `method = "fft"`.

- ...:

  Additional arguments passed to fun (e.g. `select`, `stat`, `k`). May
  include `covariates`, aligned with `pool` in the same way as `y`
  (cropped or subset per tile alongside `y`). May also include
  `normalize` for helpers that support it (e.g.
  [`analog_density()`](https://matthewkling.github.io/analogs/reference/analog_density.md));
  when normalization is requested, the global mean cell area is computed
  once over the full pool and propagated to each tile so per-tile
  `D_max` values are consistent.

- output_file:

  Optional filename for disk-based output. If specified and fun returns
  a SpatRaster, tiles are written to temporary files during processing
  and merged to output_file at the end. This is useful when results are
  too large to fit in memory. Ignored for data.frame results.

- progress:

  Logical indicating whether to show progress bar.

## Value

Same type as fun returns (SpatRaster or data.frame). If output_file is
specified, returns a disk-backed SpatRaster. When `x` is a point table,
rows come back in the order of `x`, and the `index` column refers to
rows of `x`, exactly as in a non-tiled call.

## Details

Tiled analog searches work by splitting x into a number of smaller tiles
and calling the requested analog function on each tile, using an analog
pool that is the size of the tile buffered by the geographic radius
(`geog`'s `max`). This buffer is necessary for correctness but increases
compute time, particularly if the radius is large. Raster results for
each tile are merged into a single results raster once all tiles have
processed (via temporary files when `output_file` is set).

Raster inputs are cropped to each tile and halo; point-table inputs are
subset by row. Each focal point belongs to exactly one tile, and focal
rows with missing coordinates are passed with the first tile, so they
come back as `NA` rows just as they would without tiling.

The function requires `geog` to have a finite `max`, as tiling is only
beneficial when geographic distance constraints limit the reference pool
size for each focal point. The function will warn if the radius is so
large that tiling provides minimal memory benefit.

If index_res is specified in ..., all tiles will use the same lattice
resolution. If index_res is not specified, each tile will independently
auto-tune its lattice resolution based on local data characteristics.
This adaptive behavior is generally fine and can even be beneficial when
environmental distributions vary substantially across the landscape
(e.g., mountains vs plains).
