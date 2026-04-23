# Predict from regression coefficients

Shared helper that evaluates fitted values from an
[`analog_regression()`](https://matthewkling.github.io/analogs/reference/analog_regression.md)
/ `analog_search(stat = "regression")` output by multiplying coefficient
columns with a matrix of covariate values (one row per focal).

## Usage

``` r
.predict_from_coefs(coefs_df, covariates_matrix, y_names, cov_names)
```

## Arguments

- coefs_df:

  A data.frame with `coef_intercept` and `coef_{covname}` columns
  (single-y case), or `coef_intercept_{yname}` and
  `coef_{covname}_{yname}` (multi-y case).

- covariates_matrix:

  Matrix with one row per focal and one column per covariate. Column
  order must match `cov_names`.

- y_names:

  Character vector of y variable names.

- cov_names:

  Character vector of covariate names, matching the order of columns in
  `covariates_matrix`.

## Value

A numeric matrix with `nrow(covariates_matrix)` rows and
`length(y_names)` columns, named by `y_names`.

## Details

Handles both single-y and multi-y coefficient layouts and returns a
matrix of predictions with n_focal rows and n_y columns.

This helper exists so that residual computation in
[`analog_cv()`](https://matthewkling.github.io/analogs/reference/analog_cv.md)
and any future user-facing prediction helper share the same arithmetic.
