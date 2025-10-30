# Helper: brute-force climate matches (scalar Euclidean radius)
bf_clim_radius <- function(focal, ref, radius) {
      stopifnot(ncol(focal) == ncol(ref))
      out <- vector("list", nrow(focal))
      for (i in seq_len(nrow(focal))) {
            d <- sqrt(rowSums((t(t(ref) - focal[i, ]))^2))
            out[[i]] <- which(d <= radius)
      }
      out
}

# Helper: brute-force climate kNN (Euclidean), return indices only
bf_clim_knn <- function(focal, ref, k) {
      stopifnot(ncol(focal) == ncol(ref))
      out <- vector("list", nrow(focal))
      for (i in seq_len(nrow(focal))) {
            d <- sqrt(rowSums((t(t(ref) - focal[i, ]))^2))
            ord <- order(d, seq_along(d))  # stable tie-break by index
            out[[i]] <- head(ord, k)
      }
      out
}

# Helper: brute-force per-var band pass
bf_clim_band <- function(focal, ref, band) {
      stopifnot(ncol(focal) == ncol(ref), length(band) == ncol(ref))
      out <- vector("list", nrow(focal))
      for (i in seq_len(nrow(focal))) {
            pass <- rep(TRUE, nrow(ref))
            for (j in seq_len(ncol(ref))) {
                  pass <- pass & (abs(ref[, j] - focal[i, j]) <= band[j])
            }
            out[[i]] <- which(pass)
      }
      out
}
