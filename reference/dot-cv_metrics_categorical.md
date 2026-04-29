# Categorical-outcome metrics: n, accuracy, brier, n_classes, plus K^2 confusion-matrix rows

For each y variable in the CV result, returns a data.frame with the
headline scalars followed by `K * K` confusion rows, one per (observed,
predicted) class pair, with zero counts included. Confusion-row metric
names follow the format `confusion[<obs>|<pred>]`.

## Usage

``` r
.cv_metrics_categorical(obs_labels, primary_labels, brier_vec, levels, w)
```

## Details

Inputs already validated/extracted by the categorical path:

- obs_labels: character vector of observed class labels (NA for missing
  observation).

- primary_labels: character vector of predicted class labels (NA when no
  analog had a non-NA class).

- brier_vec: numeric vector of per-focal Brier scores, computed upstream
  on row-normalized votes (NA when no usable analogs or no observed
  class).

- levels: character vector of all class labels (the global level set,
  which determines K).

- w: per-location weights (or NULL for unweighted).
