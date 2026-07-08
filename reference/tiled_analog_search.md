# Tiled analog search for memory-constrained queries

Performs analog searches on large raster datasets by dividing the focal
region into tiles and processing each tile separately. This reduces
memory usage at the cost of increasing compute time. Works with any
`analog_*()` function.

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
  ...,
  output_file = NULL,
  progress = TRUE
)
```

## Arguments

- x:

  SpatRaster with focal locations (points to find analogs for).

- pool:

  SpatRaster with reference locations (potential analog pool).

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

  Optional SpatRaster with values to aggregate across analogs. Must have
  spatial properties matching pool.

- x_cov:

  Optional SpatRaster with covariates for focal points. Must have
  spatial properties matching x.

- weight:

  Optional single-layer SpatRaster of per-pool-cell weights. Must have
  the same CRS and extent as `pool`. See
  [`analog_search()`](https://matthewkling.github.io/analogs/reference/analog_search.md)
  for details.

- cell_area_weight:

  Controls cell-area weighting. One of `"auto"` (default; on for raster
  pools, which is always the case here), `TRUE`, or `FALSE`. Unlike a
  non-tiled query, where weights are normalized to mean 1 over whatever
  pool is passed, here weights are normalized to mean 1 over the *full
  pool raster* once at the outer level, and the resulting raster is
  sliced per tile. This keeps absolute magnitudes of weighted stats
  (e.g. `sum_weights`) consistent across tiles.

- ...:

  Additional arguments passed to fun (e.g. `select`, `stat`, `k`). May
  include `covariates` as a SpatRaster matching `pool`'s CRS and extent
  (cropped per tile alongside `y`). May also include `normalize` for
  helpers that support it (e.g.
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
specified, returns a disk-backed SpatRaster.

## Details

Tiled analog searches work by splitting x into a number of smaller tiles
and calling the requested analog function on each tile, using an analog
pool that is the size of the tile buffered by the geographic radius
(`geog`'s `max`). This buffer is necessary for correctness but increases
compute time, particularly if the radius is large. The results for each
tile are temporarily written to disk, and are merged into a single
results raster once all tiles have processed.

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
