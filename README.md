
<!-- README.md is generated from README.Rmd. Please edit that file -->

# analogs

<!-- badges: start -->
<!-- badges: end -->

The **analogs** package provides fast, flexible tools for climate analog
analyses. It efficiently searches large reference datasets to find
locations that best match query sites based on climate similarity and
geographic proximity.

## Core Functions

The package supports a range of analyses. The core workhorse function is
`analog_search()`. Common query types are supported via dedicated helper
functions:

- `analog_impact()` predicts ecological state variables using analog
  impact models (AIMs)
- `analog_velocity()` finds the k nearest geographic analogs under a
  hard climate constraint
- `analog_similarity()` finds the k nearest climate analogs within a
  geographic constraint  
- `analog_availability()` counts analogs that meet both climate and
  distance thresholds
- `analog_intensity()` computes distance-weighted aggregations of analog
  properties

Geographic distance computations support projected coordinates as well
as lon-lat coordinates via great circle distances. Multivariate climate
distances are Euclidean by default, but you can use
`mahalanobis_transform()` to implement Mahalanobis distance based on
global spatial covariance, or supply `x_cov` to any `analog_*()`
function to account for site-specific climate covariance patterns.

## Computation

**analogs** uses optimized C++ algorithms and parallel processing to
handle computationally intensive pairwise distance calculations across
millions of locations. The package is built around a lattice search
index. While the above functions are sufficient for many analyses,
further optimization is possible by using `build_analog_index()` to
pre-build a reusable search index, and by using `tune_index_res()` to
find optimal index parameters for your specific data.

## Installation

You can install the development version of analogs from
[GitHub](https://github.com/) with:

``` r
pak::pak("matthewkling/analogs")
```

## Documentation

See the function documentation linked above for detailed usage examples,
or browse the [full reference
documentation](https://matthewkling.github.io/analogs/reference/index.html).
