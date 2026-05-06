# analogs

The **analogs** package implements a general framework for
distance-based neighborhood models, where analogs are selected using
climate and/or geographic distance constraints, then summarized via
weighted statistics or local regression. Common methods like climate
analog impact models, climate velocity, and geographically weighted
regression are all specific configurations of this framework.

Every analysis follows the same two-stage pattern: (1) select a
neighborhood of analog locations from a reference pool based on climate
similarity, geographic proximity, or both; and (2) summarize the
neighborhood using counts, weighted means, regression coefficients, or
other statistics.

## Core Functions

The core engine is
[`analog_search()`](https://matthewkling.github.io/analogs/reference/analog_search.md),
which exposes the full flexibility of the framework. Simplified wrapper
functions configure it for common analysis types:

- [`analog_impact()`](https://matthewkling.github.io/analogs/reference/analog_impact.md)
  predicts ecological state variables using analog impact models (AIMs)
- [`analog_regression()`](https://matthewkling.github.io/analogs/reference/analog_regression.md)
  fits local weighted regressions across analog neighborhoods
- [`analog_velocity()`](https://matthewkling.github.io/analogs/reference/analog_velocity.md)
  finds the nearest geographic analogs under a climate constraint
- [`analog_similarity()`](https://matthewkling.github.io/analogs/reference/analog_similarity.md)
  finds the nearest climate analogs within a geographic constraint
- [`analog_availability()`](https://matthewkling.github.io/analogs/reference/analog_availability.md)
  counts analogs that meet climate and distance thresholds
- [`analog_intensity()`](https://matthewkling.github.io/analogs/reference/analog_intensity.md)
  computes distance-weighted aggregations of analog properties

## Distance metrics

Geographic distance computations support projected coordinates as well
as lon-lat coordinates via great circle distances.

Climate similarity is measured using Euclidean or Mahalanobis distance
in climate space. Multivariate climate distances are Euclidean by
default, but you can use
[`mahalanobis_transform()`](https://matthewkling.github.io/analogs/reference/mahalanobis_transform.md)
to implement Mahalanobis distance based on global spatial covariance, or
supply `x_cov` to any analog function to account for site-specific
covariance patterns based on historic temporal climate variability.

## Computational optimization

The package is designed to maximize performance in various ways:

- The core search architecture is built in optimized C++.
- Parallel processing is available via the `n_threads` parameter.
- The package uses a lattice search index to efficiently query large
  candidate pools.
- For repeated queries against the same candidate pool,
  [`build_analog_index()`](https://matthewkling.github.io/analogs/reference/build_analog_index.md)
  can be used to pre-build a reusable search index.
- Optimal index parameters for your specific data can be identified with
  [`tune_index_res()`](https://matthewkling.github.io/analogs/reference/tune_index_res.md),
  which is used internally by default.
- Memory-safe queries for large raster datasets are available via
  [`tiled_analog_search()`](https://matthewkling.github.io/analogs/reference/tiled_analog_search.md).
- To increase speed at the cost of some precision, you can adaptively
  `downsample` your reference data.

## Installation

You can install the development version of analogs from
[GitHub](https://github.com/) with:

``` r

pak::pak("matthewkling/analogs")
```

## Documentation

See the [package
vignette](https://matthewkling.github.io/analogs/articles/analogs.html)
for an overview of functionality with usage examples, browse the
function documentation linked above for details, or browse the
[reference
documentation](https://matthewkling.github.io/analogs/reference/index.html)
for a complete index of functions.
