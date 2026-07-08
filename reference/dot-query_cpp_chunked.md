# Internal helper: Query C++ with optional chunking and progress tracking

Wraps query_analog_index_cpp with optional chunking for progress bars.
Handles merging of chunked results.

## Usage

``` r
.query_cpp_chunked(
  index,
  focal,
  ref_mm,
  k,
  max_env,
  max_geog,
  select_code,
  aggregate_codes,
  env_kernel_code,
  geog_kernel_code,
  theta_env,
  theta_geog,
  x_cov,
  y,
  covariates,
  lambda,
  se_code,
  n_classes_per_var = integer(0),
  area_weight = NULL,
  user_weight = NULL,
  exclude_self = FALSE,
  show_progress = FALSE
)
```
