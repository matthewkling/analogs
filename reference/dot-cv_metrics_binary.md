# Binary-outcome metrics: n, auc, tss, tss_threshold, brier

Assumes obs is already restricted to finite values. Degrades gracefully
(NA metrics) if only one class is present after filtering.

## Usage

``` r
.cv_metrics_binary(obs, pred, w)
```
