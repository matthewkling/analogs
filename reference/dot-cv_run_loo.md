# Run the LOO path: single call with exclude_self = TRUE

For fun = analog_regression, passes x_covariates = cov_mat so
predictions come back in the result (LOO: focals = pool, so focal-side
covariates are just cov_mat).

## Usage

``` r
.cv_run_loo(fun, fun_name, pool_mat, y_mat, cov_mat, extra)
```
