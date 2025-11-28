# Compare Performance Across Configurations

Runs multiple profiling tests with different settings to compare
performance.

## Usage

``` r
profile_compare(focal, ref, configs, labels = NULL)
```

## Arguments

- focal:

  Focal dataset

- ref:

  Reference dataset

- configs:

  List of configuration lists, each containing analog_search parameters

- labels:

  Character vector of labels for each config

## Value

data.frame comparing performance metrics
