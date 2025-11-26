# Validate Analog Index

Internal function to validate that an analog index is well-formed and
compatible with query data.

## Usage

``` r
.validate_analog_index(index, query_data = NULL, validate_ranges = FALSE)
```

## Arguments

- index:

  An `analog_index` object

- query_data:

  Optional matrix of query points to validate against

- validate_ranges:

  Logical; if TRUE and query_data provided, check that query points fall
  within index bounds

## Value

Invisible TRUE if valid, otherwise throws an error
