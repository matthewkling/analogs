#' Performance metrics for cross-validation output
#'
#' Computes per-variable prediction-performance metrics from the output of
#' [analog_cv()]. Handles both tabular (data.frame) and raster (SpatRaster)
#' CV output, both single-y and multi-y configurations, and both continuous
#' and binary outcomes.
#'
#' @section Output format:
#' Results are returned in long format with one row per
#' (variable, metric) pair and columns:
#'
#' - `variable`: name of the response variable.
#' - `type`: outcome type (`"continuous"` or `"binary"`).
#' - `metric`: metric name (see below).
#' - `value`: numeric metric value.
#'
#' This format is designed to accommodate additional outcome types and
#' metrics in the future without schema changes.
#'
#' @section Continuous metrics:
#' For variables classified as continuous:
#'
#' - `n`: number of locations with finite observed and predicted values.
#' - `rmse`: root mean squared error of held-out predictions.
#' - `mae`: mean absolute error.
#' - `bias`: mean signed residual (`mean(obs - pred)`). Positive values
#'   indicate systematic under-prediction.
#' - `r2`: out-of-sample R² (a.k.a. "predicted R²"), computed as
#'   `1 - SS_res / SS_tot` using the held-out residuals. Can be negative
#'   when predictions are worse than simply predicting the overall mean.
#'
#' @section Binary metrics:
#' For variables classified as binary (observed values all in `[0, 1]` with
#' both classes present):
#'
#' - `n`: number of locations with finite observed and predicted values.
#' - `auc`: area under the ROC curve (via Mann-Whitney U; handles ties).
#'   A threshold-independent measure of rank-based discrimination.
#' - `tss`: true skill statistic (`sensitivity + specificity - 1`) at the
#'   threshold that maximizes it over unique predicted values.
#' - `tss_threshold`: the threshold used for `tss`.
#' - `brier`: Brier score (mean squared error between prediction and 0/1
#'   outcome). A proper scoring rule; most interpretable when predictions
#'   are in `[0, 1]`.
#'
#' @section Outcome-type detection:
#' When `outcome_type = "auto"` (the default), each variable is classified as:
#'
#' - `"binary"` if observed values are all in `[0, 1]` after removing NAs,
#'   with both classes present.
#' - `"continuous"` otherwise.
#'
#' Users can override classification by passing a scalar type name (applies
#' to all variables) or a named character vector (one entry per variable).
#'
#' @param x Output from [analog_cv()]. Either a data.frame or a SpatRaster
#'   produced with `include_residuals = TRUE`.
#' @param outcome_type Controls outcome-type classification:
#'
#'   - `"auto"` (default): auto-detect as described above.
#'   - `"continuous"` or `"binary"`: force all variables to this type.
#'   - A named character vector with one entry per variable, giving the
#'     type for each (e.g., `c(biomass = "continuous", presence = "binary")`).
#' @param weights Optional numeric vector of per-location weights (one value
#'   per row/cell of `x`). Used for weighted versions of all metrics. Default
#'   `NULL` gives unweighted metrics.
#'
#' @return A data.frame in long format with columns `variable`, `type`,
#'   `metric`, and `value`. One row per (variable, metric) pair.
#'
#' @examples
#' \dontrun{
#' # Continuous outcome
#' cv <- analog_cv(
#'   fun      = analog_impact,
#'   pool     = sites,
#'   y        = sites$biomass,
#'   max_clim = 0.5,
#'   max_geog = 100,
#'   kernel   = "gaussian_clim",
#'   theta    = 0.2
#' )
#' cv_performance(cv)
#'
#' # Binary outcome (presence/absence)
#' cv_bin <- analog_cv(
#'   fun      = analog_impact,
#'   pool     = sites,
#'   y        = sites$presence,   # 0/1 values
#'   max_clim = 0.5,
#'   max_geog = 100,
#'   kernel   = "gaussian_clim",
#'   theta    = 0.2
#' )
#' cv_performance(cv_bin)
#'
#' # Parameter tuning via AUC
#' thetas <- c(0.1, 0.2, 0.3, 0.5)
#' auc <- sapply(thetas, function(th) {
#'   cv <- analog_cv(
#'     fun = analog_impact, pool = sites, y = sites$presence,
#'     max_clim = 0.5, max_geog = 100,
#'     kernel = "gaussian_clim", theta = th
#'   )
#'   perf <- cv_performance(cv)
#'   perf$value[perf$metric == "auc"]
#' })
#' thetas[which.max(auc)]
#' }
#'
#' @seealso [analog_cv()].
#'
#' @export
cv_performance <- function(x, outcome_type = "auto", weights = NULL) {

      # Extract obs/residual columns ------------------------------------------
      cols <- .cv_extract_obs_residual(x)
      obs_mat <- cols$obs
      res_mat <- cols$residual
      var_names <- cols$var_names
      n_vars <- length(var_names)

      # Validate weights ------------------------------------------------------
      if (!is.null(weights)) {
            if (!is.numeric(weights) || length(weights) != nrow(obs_mat)) {
                  stop("`weights` must be a numeric vector with one value per ",
                       "row/cell of `x` (expected length ", nrow(obs_mat), ").",
                       call. = FALSE)
            }
            if (any(weights < 0, na.rm = TRUE)) {
                  stop("`weights` must be non-negative.", call. = FALSE)
            }
      }

      # Resolve outcome type per variable -------------------------------------
      types <- .cv_resolve_outcome_types(outcome_type, obs_mat, var_names)

      # Dispatch to per-type metric computation -------------------------------
      parts <- vector("list", n_vars)
      for (j in seq_len(n_vars)) {
            obs <- obs_mat[, j]
            res <- res_mat[, j]
            pred <- obs - res  # recover held-out prediction from residual

            # Finite mask
            valid <- is.finite(obs) & is.finite(pred)
            if (!is.null(weights)) {
                  valid <- valid & is.finite(weights) & weights > 0
            }
            w <- if (is.null(weights)) NULL else weights[valid]

            metrics_df <- switch(
                  types[j],
                  continuous = .cv_metrics_continuous(obs[valid], pred[valid], w),
                  binary     = .cv_metrics_binary(obs[valid], pred[valid], w),
                  stop("Internal error: unhandled outcome type '", types[j], "'.",
                       call. = FALSE)
            )

            parts[[j]] <- data.frame(
                  variable = var_names[j],
                  type     = types[j],
                  metric   = metrics_df$metric,
                  value    = metrics_df$value,
                  stringsAsFactors = FALSE
            )
      }

      out <- do.call(rbind, parts)
      rownames(out) <- NULL
      out
}


# Internal helpers -----------------------------------------------------------

#' Continuous-outcome metrics: n, rmse, mae, bias, r2
#'
#' Returns a two-column data.frame (metric, value) in canonical order.
#'
#' @keywords internal
.cv_metrics_continuous <- function(obs, pred, w) {
      if (length(obs) == 0L) {
            return(data.frame(
                  metric = c("n", "rmse", "mae", "bias", "r2"),
                  value  = c(0, NA_real_, NA_real_, NA_real_, NA_real_),
                  stringsAsFactors = FALSE
            ))
      }

      if (is.null(w)) w <- rep(1.0, length(obs))
      sum_w <- sum(w)
      r <- obs - pred

      rmse <- sqrt(sum(w * r^2) / sum_w)
      mae  <- sum(w * abs(r)) / sum_w
      bias <- sum(w * r) / sum_w

      obs_wmean <- sum(w * obs) / sum_w
      ss_tot <- sum(w * (obs - obs_wmean)^2)
      ss_res <- sum(w * r^2)
      r2 <- if (ss_tot > 0) 1 - ss_res / ss_tot else NA_real_

      data.frame(
            metric = c("n", "rmse", "mae", "bias", "r2"),
            value  = c(length(obs), rmse, mae, bias, r2),
            stringsAsFactors = FALSE
      )
}


#' Binary-outcome metrics: n, auc, tss, tss_threshold, brier
#'
#' Assumes obs is already restricted to finite values. Degrades gracefully
#' (NA metrics) if only one class is present after filtering.
#'
#' @keywords internal
.cv_metrics_binary <- function(obs, pred, w) {
      n <- length(obs)
      if (n == 0L) {
            return(data.frame(
                  metric = c("n", "auc", "tss", "tss_threshold", "brier"),
                  value  = c(0, NA_real_, NA_real_, NA_real_, NA_real_),
                  stringsAsFactors = FALSE
            ))
      }

      if (is.null(w)) w <- rep(1.0, n)

      # Validate that forced-binary obs actually is binary
      if (!all(obs %in% c(0, 1))) {
            stop("Binary outcome type requested, but obs contains values ",
                 "outside {0, 1}.", call. = FALSE)
      }

      auc   <- .auc_weighted(obs, pred, w)
      tss_t <- .tss_optim_weighted(obs, pred, w)
      brier <- sum(w * (pred - obs)^2) / sum(w)

      data.frame(
            metric = c("n", "auc", "tss", "tss_threshold", "brier"),
            value  = c(n, auc, tss_t$tss, tss_t$threshold, brier),
            stringsAsFactors = FALSE
      )
}


#' Resolve outcome_type arg to a per-variable type vector
#'
#' @keywords internal
.cv_resolve_outcome_types <- function(outcome_type, obs_mat, var_names) {
      n_vars <- length(var_names)
      valid_types <- c("continuous", "binary")

      # Named vector: one entry per variable ---------------------------------
      # Trigger on any named input, including length-1 like c(y = "binary"),
      # so user intent is honored and name-based validation runs.
      if (!is.null(names(outcome_type))) {
            missing_names <- setdiff(var_names, names(outcome_type))
            if (length(missing_names) > 0L) {
                  stop("`outcome_type` is missing entries for: ",
                       paste(missing_names, collapse = ", "), ".",
                       call. = FALSE)
            }
            types <- unname(outcome_type[var_names])
            bad <- setdiff(types, valid_types)
            if (length(bad) > 0L) {
                  stop("`outcome_type` values must be one of: ",
                       paste(valid_types, collapse = ", "),
                       ". Got: ", paste(unique(bad), collapse = ", "), ".",
                       call. = FALSE)
            }
            return(types)
      }

      # Length > 1 but unnamed: error (ambiguous positional meaning)
      if (length(outcome_type) > 1L) {
            stop("`outcome_type` of length > 1 must be a named character ",
                 "vector with one entry per variable.", call. = FALSE)
      }

      # Scalar: "auto", "continuous", or "binary" ----------------------------
      if (!is.character(outcome_type) || length(outcome_type) != 1L) {
            stop("`outcome_type` must be 'auto', 'continuous', 'binary', or a ",
                 "named character vector with one entry per variable.",
                 call. = FALSE)
      }

      if (outcome_type %in% valid_types) {
            return(rep(unname(outcome_type), n_vars))
      }

      if (outcome_type == "auto") {
            return(vapply(seq_len(n_vars), function(j) {
                  if (.is_binary(obs_mat[, j])) "binary" else "continuous"
            }, character(1)))
      }

      stop("`outcome_type` must be 'auto', 'continuous', 'binary', or a ",
           "named character vector. Got: '", outcome_type, "'.",
           call. = FALSE)
}


#' Is a numeric vector binary-valued?
#'
#' TRUE iff all finite values are in `[0, 1]` AND both classes are present.
#' A single-class vector is not classified as binary because key metrics
#' (AUC, TSS) are undefined.
#'
#' @keywords internal
.is_binary <- function(obs) {
      obs_f <- obs[is.finite(obs)]
      if (length(obs_f) == 0L) return(FALSE)
      if (!all(obs_f %in% c(0, 1))) return(FALSE)
      length(unique(obs_f)) == 2L
}


#' Weighted AUC via Mann-Whitney U, handling ties by averaging ranks
#'
#' With unit weights this is the standard AUC. With nonuniform weights,
#' this is the weighted rank-sum version: each pair (pos, neg) contributes
#' w_pos * w_neg, with ties contributing 0.5 * w_pos * w_neg.
#'
#' Returns NA when only one class is present.
#'
#' @keywords internal
.auc_weighted <- function(obs, pred, w) {
      is_pos <- obs == 1
      is_neg <- obs == 0
      if (!any(is_pos) || !any(is_neg)) return(NA_real_)

      all_pred <- c(pred[is_pos], pred[is_neg])
      all_w    <- c(w[is_pos],    w[is_neg])
      all_lbl  <- c(rep(1L, sum(is_pos)), rep(0L, sum(is_neg)))

      ord <- order(all_pred)
      all_pred <- all_pred[ord]
      all_w    <- all_w[ord]
      all_lbl  <- all_lbl[ord]

      n <- length(all_pred)
      numerator <- 0.0
      cum_neg_before <- 0.0

      i <- 1L
      while (i <= n) {
            j <- i
            while (j < n && all_pred[j + 1L] == all_pred[i]) j <- j + 1L
            # Block [i..j] shares a prediction value
            block_pos_w <- sum(all_w[i:j][all_lbl[i:j] == 1L])
            block_neg_w <- sum(all_w[i:j][all_lbl[i:j] == 0L])

            if (block_pos_w > 0) {
                  numerator <- numerator +
                        block_pos_w * cum_neg_before +
                        0.5 * block_pos_w * block_neg_w
            }

            cum_neg_before <- cum_neg_before + block_neg_w
            i <- j + 1L
      }

      denom <- sum(w[is_pos]) * sum(w[is_neg])
      if (denom <= 0) return(NA_real_)
      numerator / denom
}


#' Optimized TSS: scan candidate thresholds and return maximizing threshold
#'
#' TSS = sensitivity + specificity - 1, with positive classification when
#' pred >= threshold. Returns list(tss, threshold). NA when only one class
#' is present.
#'
#' @keywords internal
.tss_optim_weighted <- function(obs, pred, w) {
      is_pos <- obs == 1
      is_neg <- obs == 0
      sum_w_pos <- sum(w[is_pos])
      sum_w_neg <- sum(w[is_neg])
      if (sum_w_pos == 0 || sum_w_neg == 0) {
            return(list(tss = NA_real_, threshold = NA_real_))
      }

      ord <- order(pred)
      p_sorted <- pred[ord]
      w_sorted <- w[ord]
      pos_sorted <- is_pos[ord]

      # Start at threshold = -Inf: everything classified positive.
      #   TP = sum_w_pos, FN = 0, FP = sum_w_neg, TN = 0
      TP <- sum_w_pos; FN <- 0
      FP <- sum_w_neg; TN <- 0
      best_tss <- 0
      best_thr <- -Inf

      n <- length(p_sorted)
      i <- 1L
      while (i <= n) {
            j <- i
            while (j < n && p_sorted[j + 1L] == p_sorted[i]) j <- j + 1L

            block_w_pos <- sum(w_sorted[i:j][pos_sorted[i:j]])
            block_w_neg <- sum(w_sorted[i:j][!pos_sorted[i:j]])

            # Flip whole tied block from predicted-positive to predicted-negative
            TP <- TP - block_w_pos; FN <- FN + block_w_pos
            FP <- FP - block_w_neg; TN <- TN + block_w_neg

            sens <- TP / sum_w_pos
            spec <- TN / sum_w_neg
            tss  <- sens + spec - 1

            if (is.finite(tss) && tss > best_tss) {
                  best_tss <- tss
                  next_val <- if (j < n) p_sorted[j + 1L] else p_sorted[j] + 1e-12
                  best_thr <- 0.5 * (p_sorted[j] + next_val)
            }

            i <- j + 1L
      }

      list(tss = best_tss, threshold = best_thr)
}


#' Extract obs and residual matrices from CV output
#'
#' Handles both data.frame and SpatRaster input, and both single-y layout
#' (`obs`, `residual`) and multi-y layout (`obs_{name}`, `residual_{name}`).
#'
#' @return a list with `obs` (matrix, n_cells x n_vars), `residual`
#'   (matrix, same dims), and `var_names` (character vector).
#'
#' @keywords internal
.cv_extract_obs_residual <- function(x) {

      is_raster <- inherits(x, "SpatRaster")

      if (is_raster) {
            if (!requireNamespace("terra", quietly = TRUE)) {
                  stop("Package 'terra' is required for SpatRaster input.",
                       call. = FALSE)
            }
            nms <- terra::names(x)
            vals <- as.matrix(terra::as.data.frame(x, xy = FALSE, na.rm = FALSE))
            colnames(vals) <- nms
            get_col <- function(cn) vals[, cn, drop = TRUE]
            col_names <- nms
      } else if (is.data.frame(x)) {
            col_names <- names(x)
            get_col <- function(cn) x[[cn]]
      } else {
            stop("`x` must be a data.frame or SpatRaster (output of analog_cv()).",
                 call. = FALSE)
      }

      has_single <- "obs" %in% col_names && "residual" %in% col_names
      obs_multi <- grep("^obs_", col_names, value = TRUE)
      res_multi <- grep("^residual_", col_names, value = TRUE)

      if (has_single) {
            obs <- matrix(get_col("obs"), ncol = 1L)
            res <- matrix(get_col("residual"), ncol = 1L)
            colnames(obs) <- colnames(res) <- "y"
            return(list(obs = obs, residual = res, var_names = "y"))
      }

      if (length(obs_multi) > 0L) {
            names_from_obs <- sub("^obs_", "", obs_multi)
            names_from_res <- sub("^residual_", "", res_multi)
            paired <- intersect(names_from_obs, names_from_res)
            if (length(paired) == 0L) {
                  stop("Could not pair any obs_* columns with residual_* columns. ",
                       "Did this come from analog_cv(include_residuals = TRUE)?",
                       call. = FALSE)
            }
            paired <- names_from_obs[names_from_obs %in% paired]

            obs <- do.call(cbind, lapply(paired,
                                         function(nm) get_col(paste0("obs_", nm))))
            res <- do.call(cbind, lapply(paired,
                                         function(nm) get_col(paste0("residual_", nm))))
            colnames(obs) <- colnames(res) <- paired

            return(list(obs = obs, residual = res, var_names = paired))
      }

      stop(
            "No obs/residual columns found in `x`. This function expects the ",
            "output of analog_cv() with `include_residuals = TRUE`.",
            call. = FALSE
      )
}
