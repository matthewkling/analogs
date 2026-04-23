# Validate and normalize query parameters

Validates select/stat/k/kernel/theta/x_cov/y/covariates/lambda/se
combinations and normalizes values. Returns a list with normalized
parameters.

## Usage

``` r
.validate_query_params(
  focal = NULL,
  ref = NULL,
  x_cov = NULL,
  y = NULL,
  covariates = NULL,
  max_clim,
  max_geog,
  select,
  k,
  stat,
  kernel,
  theta,
  lambda = 0,
  se = "none"
)
```
