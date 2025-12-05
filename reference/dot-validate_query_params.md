# Validate and normalize query parameters

Validates select/stat/k/weight/theta/x_cov combinations and normalizes
values. Returns a list with normalized parameters.

## Usage

``` r
.validate_query_params(
  select,
  stat,
  k,
  weight,
  theta,
  x_cov = NULL,
  focal_mm = NULL
)
```
