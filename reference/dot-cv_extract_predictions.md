# Extract predictions (one column per y variable) from a CV result data.frame

Dispatch:

- pred_target == "weighted_mean": read weighted_mean_yname columns.

- pred_target == "regression" AND fun_name == "analog_regression": read
  pred_yname columns that analog_regression() attached via x_covariates.

- pred_target == "regression" AND fun_name == "analog_search": compute
  via .predict_from_coefs(), since analog_search does not accept
  x_covariates.

## Usage

``` r
.cv_extract_predictions(
  pred_target,
  fun_name,
  result_df,
  y_names,
  cov_mat,
  cov_names
)
```
