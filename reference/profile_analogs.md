# Profile find_analogs Performance

Runs find_analogs with detailed C++ profiling instrumentation and
returns timing breakdowns for optimization analysis.

## Usage

``` r
profile_analogs(..., report = TRUE, plot = TRUE)
```

## Arguments

- ...:

  Arguments passed to find_analogs

- report:

  Logical; if TRUE, print formatted report

- plot:

  Logical; if TRUE, create visualization

## Value

List containing: - result: the normal find_analogs output - profile:
data.frame with timing details - counters: data.frame with operation
counts
