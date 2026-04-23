# Optimized TSS: scan candidate thresholds and return maximizing threshold

TSS = sensitivity + specificity - 1, with positive classification when
pred \>= threshold. Returns list(tss, threshold). NA when only one class
is present.

## Usage

``` r
.tss_optim_weighted(obs, pred, w)
```
