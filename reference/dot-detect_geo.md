# Detect coordinate system from data ranges

Resolves `coord_type = "auto"` to either `"lonlat"` or `"projected"`.
When the input has CRS metadata (SpatRaster, SpatVector), that's
authoritative. For raw matrix / data.frame inputs we fall back to a
coordinate-magnitude heuristic: if all xy values fit within
`[-180, 180] x [-90, 90]`, treat as lonlat; otherwise projected.

## Usage

``` r
.detect_geo(xy, input = NULL)
```

## Details

The magnitude heuristic is unavoidably ambiguous (a small projected
region in meters can fit in lonlat bounds). Users with matrix-form
projected data in that range should pass `coord_type` explicitly.
