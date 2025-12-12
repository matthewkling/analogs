
<!-- README.md is generated from README.Rmd. Please edit that file -->

# analogs

<!-- badges: start -->
<!-- badges: end -->

The **analogs** package provides fast, flexible tools for climate analog
analyses. It efficiently searches large reference datasets to find
locations that best match query sites based on climate similarity and
geographic proximity.

## Core Functions

The package supports a range of analysis types. The core workhorse
function is `analog_search()`, which supports a range of methods for
selecting and summarizing analogs for each query location. Common query
types are supported via simplified wrapper functions:

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

## Distance metrics

Geographic distance computations support projected coordinates as well
as lon-lat coordinates via great circle distances.

Climate similarity is measured using Euclidean or Mahalanobis distance
in climate space. Multivariate climate distances are Euclidean by
default, but you can use `mahalanobis_transform()` to implement
Mahalanobis distance based on global spatial covariance, or supply
`x_cov` to any analog function to account for site-specific covariance
patterns based on historic temporal climate variability.

## Computational optimization

Various aspects of the search architecture are designed to optimize
performance:

- The core search architecture is built in optimized C++.
- Parallel processing is available via the `n_threads` parameter.
- The package uses a lattice search index to efficiently query large
  candidate pools.
- For repeated queries against the same candidate pool,
  `build_analog_index()` can be used to pre-build a reusable search
  index.
- Use the `tune_index_res()` function to identify optimal index
  parameters for your specific data.
- Memory-safe queries for large raster datasets are available via
  `tiled_analog_search()`.
- To increase speed at the cost of some precision, you can optionally
  `downsample` your reference data.

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
