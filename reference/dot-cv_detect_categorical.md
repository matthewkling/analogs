# Detect categorical CV output by column structure

Categorical results from
[`analog_cv()`](https://matthewkling.github.io/analogs/reference/analog_cv.md)
with `stat = "tabulate"` carry `obs[_<name>]`, `primary[_<name>]`, and
`brier[_<name>]` columns. If any of these are present, this is a
categorical run; we extract the per-variable label vectors here so the
caller doesn't need to repeat the column-name parsing.

## Usage

``` r
.cv_detect_categorical(x)
```

## Details

Note: in raster output, character columns (obs, primary) are dropped by
`.cv_to_raster()`, so categorical CV is only fully usable with
data.frame output. We error informatively if the user passes a
SpatRaster that lacks the labels.
