# Get metadata from analog search results

Returns the metadata attributes attached to results from
[`analog_search()`](https://matthewkling.github.io/analogs/reference/analog_search.md),
or any `analog_*()` wrapper function that calls it, as a named list.
These same values can be accessed individually via `attr(x, "<name>")`;
`metadata()` is a convenience that bundles them into a consistent list
structure independent of the underlying attribute layout. Use it to
inspect or programmatically reconstruct the parameterization that
produced a result. Each attribute records one user-facing parameter or
input-data property.

## Usage

``` r
metadata(x)
```

## Arguments

- x:

  A result from any `analog_*()` function (data.frame or SpatRaster).

## Value

A named list. Elements depend on which attributes are present on `x`;
missing attributes are simply absent from the returned list (no `NA`
placeholders). Possible elements include:

### Selection and aggregation parameters

- `select`:

  Selection strategy: `"all"`, `"knn_env"`, or `"knn_geog"`.

- `k`:

  Number of neighbors for kNN selection modes; `NULL` for
  `select = "all"`.

- `stat`:

  Aggregation statistic(s) requested. Character vector; `"none"` for
  pair-mode results.

- `kernel_env`:

  Environmental weighting kernel shape: `"uniform"`, `"gaussian"`, or
  `"inverse"`. `"uniform"` when environment is not distance-weighted.

- `kernel_geog`:

  Geographic weighting kernel shape: `"uniform"`, `"gaussian"`, or
  `"inverse"`. `"uniform"` when geography is not distance-weighted.

- `theta_env`:

  Scale parameter for the environmental kernel (Gaussian bandwidth or
  inverse half-weight distance); `NULL` if the environmental kernel is
  uniform.

- `theta_geog`:

  Scale parameter for the geographic kernel; `NULL` if the geographic
  kernel is uniform.

- `lambda`:

  Ridge penalty for `stat = "regression"`; `0` for ordinary weighted
  least squares.

- `se`:

  Standard-error framing applied to SE-supporting stats: `"none"`,
  `"ess"`, or `"design"`.

- `max_env`:

  Environmental-distance threshold used for analog selection. Scalar for
  Euclidean / Mahalanobis radius, or length-`n_env` for per-variable
  thresholds.

- `max_geog`:

  Geographic-distance threshold used for analog selection. Units are km
  when `coord_type = "lonlat"`, projection units otherwise.

- `exclude_self`:

  Logical: was each focal's own pool row excluded from its analog
  neighborhood? `TRUE` for results from
  [`analog_cv()`](https://matthewkling.github.io/analogs/reference/analog_cv.md)
  (LOO) or from `analog_search(exclude_self = TRUE)`.

### Input-data properties

- `n_x`:

  Number of focal locations in the input `x`.

- `n_pool`:

  Number of reference locations in the input `pool`.

- `n_env`:

  Number of environmental variables.

- `coord_type`:

  Coordinate system: `"lonlat"` (great-circle distances in km) or
  `"projected"` (planar distances in projection units).

- `x_cov_provided`:

  Logical: was a per-focal covariance matrix (`x_cov`) supplied? If
  `TRUE`, environmental distances were computed in Mahalanobis units
  using location-specific covariance; otherwise Euclidean.

- `downsample_actual`:

  Actual downsampling rate applied to the reference pool (1.0 if no
  downsampling). May exceed the requested rate when maintaining at least
  one point per occupied bin requires it.

### Index diagnostics (lattice-related)

These attributes describe the internal lattice index used to accelerate
the search; they're useful for performance debugging but not needed to
reproduce the result. See
[`build_analog_index()`](https://matthewkling.github.io/analogs/reference/build_analog_index.md)
for detail.

- `binning_method`:

  Internal indexing method; one of `"lattice"` or `"lattice_ecef"`.

- `total_bins`, `n_bins_nonempty`, `min_bin_occupancy`,
  `max_bin_occupancy`, `avg_bin_occupancy`,
  `avg_nonempty_bin_occupancy`:

  Lattice occupancy summaries.

- `env_res_adj`, `geog_res_adj`:

  Per-family resolution adjustments used to build the index: each scales
  its family's bin count relative to a data-dependent default (`1` =
  default, `0` = deactivated). See
  [`build_analog_index()`](https://matthewkling.github.io/analogs/reference/build_analog_index.md).

- `env_target`, `geo_target`:

  Realized per-family bin-count targets passed to the lattice builder
  (the absolute values the `*_res_adj` resolve to).

- `bins_per_axis`:

  Integer vector of realized bins per axis (geographic axes first, then
  environment), showing how the budget was distributed.

### Cross-validation metadata

Present only when `x` is a result from
[`analog_cv()`](https://matthewkling.github.io/analogs/reference/analog_cv.md).

- `cv_method`:

  `"loo"` or `"kfold"`.

- `cv_fun`:

  Name of the analog function being cross-validated (e.g.
  `"analog_impact"`).

- `cv_n_folds`:

  Number of folds (k-fold only).

- `cv_pred_target`:

  Mechanism used to extract predictions for residuals:
  `"weighted_mean"`, `"regression"`, or `"none"`.

## Details

`metadata()` is the documented entry point for reading these attributes;
the attribute names themselves are also part of the public contract, so
`attr(x, "select")` etc. are equally valid for individual values. The
function's main advantages over raw
[`attributes()`](https://rdrr.io/r/base/attributes.html) are: it filters
out R-internal attributes (`names`, `class`, `dim`, ...), and the
`?metadata` help page gives a single canonical reference for every
attribute the package attaches.

## Examples

``` r
if (FALSE) { # \dontrun{
result <- analog_velocity(
  x = future_climate,
  pool = current_climate,
  env = kernel(max = 0.5)
)

# Get the full parameterization as a list
meta <- metadata(result)
meta$select
meta$max_env
meta$n_pool

# Or read individual attributes directly
attr(result, "select")
} # }
```
