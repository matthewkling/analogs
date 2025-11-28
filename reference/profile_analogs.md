# Profile analog_search Performance

Runs analog_search with detailed C++ profiling instrumentation and
returns timing breakdowns for optimization analysis.

## Usage

``` r
profile_analogs(..., report = TRUE, plot = TRUE)
```

## Arguments

- ...:

  Arguments passed to analog_search

- report:

  Logical; if TRUE, print formatted report

- plot:

  Logical; if TRUE, create visualization

## Value

List containing: - result: the normal analog_search output - profile:
data.frame with timing details - counters: data.frame with operation
counts
