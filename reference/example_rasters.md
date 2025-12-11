# Example spatiotemporal climate rasters

Load example climate data rasters covering a portion of North America.
Requires the `terra` package. The data are for two time periods (recent
decades and future: year 2041-2070, scenario SSP585, model gfdl-esm4).
Each period has data for two scaled climate variables (climatic water
deficit (CWD) and actual evapotranspiration (AET)), which were derived
using raw data from CHELSA v2 (Karger et al. 2017).

## Usage

``` r
example_rasters()
```

## Value

a list of two SpatRasters, each with two layers (CWD, AET).
