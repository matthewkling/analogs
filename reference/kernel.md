# Define the distance treatment for one analog dimension

Constructs a specification of how a single distance family — environment
or geography — is filtered and weighted in an analog search. A `kernel`
bundles the per-family choices: the hard distance thresholds (`max` and
`min`), the weighting kernel shape (`weight`), and that kernel's scale
parameter (`theta`). Pass one to the `env` argument of
[`analog_search()`](https://matthewkling.github.io/analogs/reference/analog_search.md)
and/or one to the `geog` argument.

## Usage

``` r
kernel(weight = NULL, theta = NULL, max = NULL, min = NULL)
```

## Arguments

- weight:

  Kernel shape for this family. One of `"uniform"` (constant weight 1;
  the default when `NULL`), `"gaussian"` (`exp(-d^2 / (2 theta^2))`), or
  `"inverse"` (`1 / (1 + d / theta)`, a heavy-tailed inverse-distance
  kernel). The family (environmental vs geographic) is determined by
  whether the kernel is passed as `env` or `geog`, so the shape name is
  unqualified.

- theta:

  Scale parameter for the `weight` kernel. For `"gaussian"` it is the
  bandwidth (sigma); for `"inverse"` it is the half-weight distance (the
  weight is 1/2 at `d = theta`). Ignored for `"uniform"`. Defaults to
  `NULL`, which lets downstream code apply a default of 1. See
  [`kernel_params()`](https://matthewkling.github.io/analogs/reference/kernel_params.md)
  for help choosing a value calibrated to a target coverage fraction.

- max:

  Hard upper distance threshold for this family: candidates beyond `max`
  (in this family's distance) are excluded. `NULL` (default) means no
  upper threshold. Usually a single radius. For the environmental
  family, `max` may also be a vector of per-variable absolute-difference
  thresholds (length equal to the number of environmental variables);
  the geographic family uses a single radius. Supplied to the search as
  `max_env` / `max_geog`.

- min:

  Hard lower distance threshold for this family: candidates closer than
  `min` (in this family's distance) are excluded, so the retained
  candidates form an annulus `min <= d <= max`. A single positive
  scalar, or `NULL` (default) for no lower threshold. Currently
  supported only for the **geographic** family; setting `min` on an
  environmental kernel is an error. The primary use case is buffered
  spatial cross-validation (e.g. via
  [`analog_cv()`](https://matthewkling.github.io/analogs/reference/analog_cv.md)),
  where excluding geographically near-duplicate candidates around each
  focal gives a less optimistic estimate of predictive skill than plain
  leave-one-out. Supplied to the search as `min_geog`.

## Value

An object of class `"analog_kernel"`: a list with elements `weight`,
`theta`, `max`, and `min` (each possibly `NULL`).

## Details

The overall kernel weight for a candidate is the product of the two
families' weights, so the families are specified independently and may
use different shapes (e.g. an inverse-distance environmental kernel
together with a Gaussian geographic kernel). A family with
`weight = "uniform"` (or `NULL`) contributes a constant weight of 1,
i.e. it filters (if `max`/`min` are set) but does not down-weight by
distance.

All four components are optional. Which combinations are valid depends
on the operation and is checked by
[`analog_search()`](https://matthewkling.github.io/analogs/reference/analog_search.md)
downstream (for example, climate velocity requires an environmental
`max`; a weighted statistic requires a non-uniform `weight` on at least
one family). A bare `NULL` passed as the `env` or `geog` argument is
equivalent to `kernel()` with all components unset: no threshold and no
weighting for that family.

## See also

[`analog_search()`](https://matthewkling.github.io/analogs/reference/analog_search.md),
[`kernel_params()`](https://matthewkling.github.io/analogs/reference/kernel_params.md)

## Examples

``` r
# Environment: keep analogs within 2 environmental-distance units, Gaussian-weighted
kernel(weight = "gaussian", theta = 0.5, max = 2)
#> <analog_kernel>
#>   weight: gaussian
#>   theta:  0.5
#>   max:    2
#>   min:    none

# Geography: hard 100 km cutoff, no distance weighting (uniform)
kernel(max = 100)
#> <analog_kernel>
#>   weight: uniform
#>   max:    100
#>   min:    none

# Geography: annulus keeping analogs between 5 km and 100 km of each focal
# (e.g. a 5 km buffer for spatial cross-validation)
kernel(max = 100, min = 5)
#> <analog_kernel>
#>   weight: uniform
#>   max:    100
#>   min:    5

# Inverse-distance environmental weighting, no hard cutoff
kernel("inverse", theta = 1)
#> <analog_kernel>
#>   weight: inverse
#>   theta:  1
#>   max:    none
#>   min:    none
```
