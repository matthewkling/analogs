# Determine the prediction target column/mechanism for a CV run

Returns one of "weighted_mean", "regression", or "none". Errors on
ambiguous analog_search stats (both weighted_mean and regression).

## Usage

``` r
.cv_determine_prediction_target(fun_name, extra, cov_mat)
```
