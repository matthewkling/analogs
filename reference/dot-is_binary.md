# Is a numeric vector binary-valued?

TRUE iff all finite values are in 0, 1 AND both classes are present. A
single-class vector is not classified as binary because key metrics
(AUC, TSS) are undefined.

## Usage

``` r
.is_binary(obs)
```
