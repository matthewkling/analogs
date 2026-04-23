# Extract obs and residual matrices from CV output

Handles both data.frame and SpatRaster input, and both single-y layout
(`obs`, `residual`) and multi-y layout (`obs_{name}`,
`residual_{name}`).

## Usage

``` r
.cv_extract_obs_residual(x)
```

## Value

a list with `obs` (matrix, n_cells x n_vars), `residual` (matrix, same
dims), and `var_names` (character vector).
