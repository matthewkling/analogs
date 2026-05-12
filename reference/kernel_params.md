# Kernel parameter recommendations

Returns recommended values for the bandwidth (`theta`) and hard
truncation distance (`max`) of a kernel operating on Euclidean distances
in `d`-dimensional space. Gives theoretical answers to the questions,
"How big should `theta` be in order for my kernel to capture a given
fraction of pool sites for the typical focal site?" and "How big should
`max_clim` or `max_geog` be to truncate only a given percentage of
kernel weight for the typical focal site?" Accounts for the effects of
analog space multidimensionality on pairwise distance distributions,
which can result in one-dimensional intuitions being incorrect. Supports
the kernel types and data distributions used by the `analogs` package.

## Usage

``` r
kernel_params(
  fraction = NULL,
  theta = NULL,
  d,
  loss = NULL,
  kernel = c("gaussian", "uniform", "inverse_distance"),
  data_dist = c("mvn", "uniform")
)
```

## Arguments

- fraction:

  Target fraction of pool sites captured (in weight-proportional terms)
  by the kernel for the typical focal. Use this OR `theta`. Requires
  `data_dist = "mvn"`.

- theta:

  Bandwidth value. Use this OR `fraction`. For Gaussian kernels, this is
  the standard bandwidth parameter; for uniform kernels, this is the
  cutoff radius (also returned as `max`); for inverse-distance kernels,
  this is the epsilon regularization that determines the half-weight
  scale.

- d:

  Dimensionality of the space (e.g., number of climate variables after
  Mahalanobis transformation, or 2 for geographic).

- loss:

  Fraction of aggregate kernel weight to discard at the truncation
  distance `max`. If `NULL` (default), `max` is not computed.

- kernel:

  One of `"gaussian"` (default), `"uniform"`, or `"inverse_distance"`.

- data_dist:

  Distribution of cells in space. Either `"mvn"` (multivariate standard
  normal; default; appropriate for Mahalanobis-transformed climate data)
  or `"uniform"` (appropriate for geographic space).

## Value

A named list. For Gaussian and inverse-distance kernels: element
`theta`, and `max` if `loss` is specified. For uniform kernels: element
`max` (the single cutoff radius, which serves as both bandwidth and
truncation distance; supplied in the `analogs` package as `max_clim` or
`max_geog`).

## Details

Either `fraction` or `theta` should be provided. When `fraction` is
given, the function returns the `theta` that calibrates the kernel to
capture, on average, that fraction of pool sites (where partial capture
is in proportion to kernel weight). Switching kernel shapes at fixed
`fraction` holds the expected total kernel weight constant, so weighted
aggregate statistics (e.g. `sum_weights`) remain comparable across
kernels.

For uniform data, `fraction` is not meaningful because "fraction of
space" depends on landscape extent; `theta` must be supplied directly
(e.g. a dispersal-derived bandwidth for a geographic kernel).

When `loss` is specified, the function additionally returns `max`: the
truncation distance beyond which less than `loss` of aggregate kernel
weight is discarded. Useful for computational efficiency.

Recommendations are averages over the distribution of focal cells;
specific focal cells experience effective neighborhoods that vary around
these averages, with cells in dense climate regions seeing more
neighbors than cells in sparse regions.

## Examples

``` r
# Climate kernel: niche fraction of 5% in 4 climate variables
kernel_params(fraction = 0.05, d = 4, loss = 0.01)
#> $theta
#> [1] 0.536663
#> 
#> $max
#> [1] 1.723009
#> 

# Geographic kernel: 500 km dispersal-based bandwidth
kernel_params(theta = 500, d = 2, data_dist = "uniform", loss = 0.01)
#> $theta
#> [1] 500
#> 
#> $max
#> [1] 1517.427
#> 

# Switching kernels at fixed niche fraction (matched expected weight)
kernel_params(fraction = 0.05, d = 4, kernel = "gaussian")
#> $theta
#> [1] 0.536663
#> 
kernel_params(fraction = 0.05, d = 4, kernel = "uniform")
#> $max
#> [1] 0.8430439
#> 
kernel_params(fraction = 0.05, d = 4, kernel = "inverse_distance")
#> $theta
#> [1] 18.14302
#> 
```
