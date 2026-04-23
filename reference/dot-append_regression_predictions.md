# Append `pred` column(s) to an analog_regression result

Evaluates fitted values at focal locations by combining per-focal
regression coefficients in `result` with focal-side covariate values in
`x_covariates`. Handles both data.frame and SpatRaster `result`, both
single-y and multi-y regression, and both data.frame/matrix and
SpatRaster `x_covariates` inputs. Uses
[`.predict_from_coefs()`](https://matthewkling.github.io/analogs/reference/dot-predict_from_coefs.md)
for the arithmetic.

## Usage

``` r
.append_regression_predictions(result, x_covariates, y, covariates)
```

## Details

Name resolution uses the package's existing validator helpers to pull
`y_names` / `cov_names` from the original `y` / `covariates` arguments,
avoiding any need to parse coefficient column names.
