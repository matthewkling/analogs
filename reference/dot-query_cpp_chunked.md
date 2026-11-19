# Internal helper: Query C++ with optional chunking and progress tracking

Wraps query_analog_index_cpp. When `show_progress` is TRUE, focal rows
are processed in chunks from the main R thread (the C++ worker threads
cannot call R) and a progressr step is signalled after each chunk.
Results are identical to a single call: each focal row is computed
independently, and self-exclusion is carried per chunk through
`self_rows`.

## Usage

``` r
.query_cpp_chunked(
  index,
  focal,
  ref_mm,
  k,
  max_env,
  max_geog,
  min_geog = 0,
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
  self_rows = NULL,
  show_progress = FALSE
)
```
