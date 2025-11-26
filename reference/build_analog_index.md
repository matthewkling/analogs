# Build Analog Index

Pre-builds a reusable lattice index from reference climate data. The
index can be queried multiple times with different focal points and
parameters, avoiding the need to rebuild the lattice for each query.

## Usage

``` r
build_analog_index(
  pool,
  coord_type = c("auto", "lonlat", "projected"),
  index_res = 16
)
```

## Arguments

- pool:

  The reference dataset to search for analogs. Either:

  - Matrix/data.frame with columns x, y, and climate variables, or
    SpatRaster with climate variable layers, OR

  - An `analog_index` object created by `build_analog_index` (for
    repeated queries).

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

## Value

An S3 object of class `"analog_index"` containing:

- The compiled lattice index (internal C++ structure)

- Reference data

- Metadata: coordinate type, dimensions, ranges, resolution

- Diagnostics: bin counts and occupancy statistics

## Details

The lattice index is built over both geographic and climate dimensions,
allowing efficient spatial queries regardless of the constraint values
used at query time. For lon/lat coordinates, the index uses ECEF
(Earth-Centered Earth-Fixed) space internally for optimal performance.

Index resolution (`index_res`) controls the granularity of spatial
binning. The optimal value depends on your data size and query patterns.
Use
[`tune_index_res`](https://matthewkling.github.io/analogs/reference/tune_index_res.md)
to find the best resolution for your use case, or accept the default of
16 which works well for most applications.

## Examples

``` r
if (FALSE) { # \dontrun{
# Build index with default settings
index <- build_analog_index(climate_data)

# Build with explicit resolution
index <- build_analog_index(climate_data, index_res = 20)

# Query the index multiple times
v1 <- analog_velocity(sites1, pool = index, max_clim = 0.5)
v2 <- analog_velocity(sites2, pool = index, max_clim = 0.3)
a1 <- analog_availability(sites3, pool = index, max_clim = 0.5, max_geog = 100)
} # }
```
