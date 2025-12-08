# Validate and normalize query parameters

Validates select/stat/k/weight/theta/x_cov/values combinations and
normalizes values. Returns a list with normalized parameters.

## Usage

``` r
.validate_query_params(
  focal = NULL,
  ref = NULL,
  x_cov = NULL,
  values = NULL,
  max_clim,
  max_geog,
  select,
  k,
  stat,
  weight,
  theta
)
```
