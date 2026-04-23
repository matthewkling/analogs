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
  max_clim,
  max_geog,
  select_code,
  aggregate_codes,
  kernel_code,
  theta,
  x_cov,
  y,
  covariates,
  lambda,
  se_code,
  show_progress = FALSE
)
```
