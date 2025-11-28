# Transform Climate Data for Global Mahalanobis Distance

Transforms climate data to decorrelated, unit-variance space, enabling
the use of Euclidean distance as a global Mahalanobis distance. Useful
for climate analog analysis when you want to account for spatial
correlation structure and/or standardize variables with different units.

## Usage

``` r
mahalanobis_transform(x, pool = NULL, center = TRUE, scale = TRUE)
```

## Arguments

- x:

  Climate data for focal/query points. Can be:

  - Matrix/data.frame with columns x, y, and climate variables

  - SpatRaster with climate variable layers

- pool:

  Optional reference climate data to transform jointly with `x`. When
  provided, both datasets are transformed using the same transformation
  (computed from their combined covariance structure). Must have the
  same number of climate variables as `x`. Same format options as `x`.

- center:

  Logical; if TRUE (default), center each variable to mean zero before
  transformation.

- scale:

  Logical; if TRUE (default), standardize each variable to unit variance
  before transformation (correlation-based). If FALSE, use
  covariance-based transformation which preserves relative variance
  magnitudes. Generally TRUE is recommended when climate variables have
  different units (e.g., temperature in °C vs precipitation in mm).

## Value

If `pool = NULL`: Returns transformed version of `x` in the same format
(matrix/data.frame/SpatRaster).

If `pool` is provided: Returns a list with two components:

- `$x`: Transformed version of `x`

- `$pool`: Transformed version of `pool`

Both use the same format as their respective inputs.

## Details

The transformation uses eigendecomposition of the covariance (or
correlation) matrix to decorrelate variables and scale to unit variance
in all directions. After transformation, Euclidean distance in the
transformed space is equivalent to Mahalanobis distance in the original
space.

The transformation is: X_transformed = (X - μ)

When `pool` is provided, the covariance structure is computed from the
combined dataset to ensure both are transformed to the same space.

If the covariance matrix is near-singular, a small regularization term
(1e-6 \* max eigenvalue) is added to stabilize the inversion.

## Examples

``` r
if (FALSE) { # \dontrun{
# Single dataset transformation
clim_transformed <- mahalanobis_transform(climate_data)

# Joint transformation for analog analysis
transformed <- mahalanobis_transform(
  x = baseline_climate,
  pool = future_climate
)

# Use transformed data with Euclidean distance (= global Mahalanobis)
analogs <- analog_search(
  x = transformed$x,
  pool = transformed$pool,
  mode = "knn_geog",
  max_clim = 2,  # Now in standardized units
  k = 1
)

# Covariance-based transformation (preserve relative variances)
transformed_cov <- mahalanobis_transform(
  x = climate_data,
  scale = FALSE
)
} # }
```
