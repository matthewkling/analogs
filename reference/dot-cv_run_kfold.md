# Run the k-fold path: loop over folds, rebuilding the index per fold

For fun = analog_regression, passes x_covariates = cov_matfocal_rows, so
per-fold predictions come back in the result, avoiding downstream
re-computation.

## Usage

``` r
.cv_run_kfold(fun, fun_name, pool_mat, y_mat, cov_mat, folds, extra)
```
