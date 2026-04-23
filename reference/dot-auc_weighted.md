# Weighted AUC via Mann-Whitney U, handling ties by averaging ranks

With unit weights this is the standard AUC. With nonuniform weights,
this is the weighted rank-sum version: each pair (pos, neg) contributes
w_pos \* w_neg, with ties contributing 0.5 \* w_pos \* w_neg.

## Usage

``` r
.auc_weighted(obs, pred, w)
```

## Details

Returns NA when only one class is present.
