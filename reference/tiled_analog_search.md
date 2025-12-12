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
  max_geog,
  values = NULL,
  x_cov = NULL,
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

- max_geog:

  Maximum geographic distance constraint (default: NULL = no geographic
  constraint). When specified, only reference locations within this
  distance are considered. Radius units should be specified in
  kilometers if `coord_type = "lonlat"`, or in projected coordinate
  units if `coord_type = "projected"`.

- values:

  Optional SpatRaster with values to aggregate across analogs. Must have
  spatial properties matching pool.

- x_cov:

  Optional SpatRaster with covariates for focal points. Must have
  spatial properties matching x.

- ...:

  Additional arguments passed to fun. Must include max_geog.

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
pool that is the size of the tile buffered by max_geog. This buffer is
necessary for correctness but increases compute time, particularly if
max_geog is large. The results for each tile are temporarily written to
disk, and are merged into a single results raster once all tiles have
processed.

The function requires max_geog to be specified, as tiling is only
beneficial when geographic distance constraints limit the reference pool
size for each focal point. The function will warn if max_geog is so
large that tiling provides minimal memory benefit.

If index_res is specified in ..., all tiles will use the same lattice
resolution. If index_res is not specified, each tile will independently
auto-tune its lattice resolution based on local data characteristics.
This adaptive behavior is generally fine and can even be beneficial when
climate distributions vary substantially across the landscape (e.g.,
mountains vs plains).
