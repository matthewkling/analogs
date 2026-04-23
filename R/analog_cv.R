#' Cross-validate an analog function
#'
#' Runs an analog impact or regression analysis in cross-validation mode,
#' generating held-out predictions and residuals for each location in `pool`.
#' Each location is predicted using only neighbors that exclude itself,
#' providing an honest assessment of how well the specified configuration
#' predicts observed `y` values. Supports leave-one-out (LOO) and k-fold
#' cross-validation methods.
#'
#' `analog_cv()` supports two CV methods:
#'
#' - `"loo"` (leave-one-out): Each focal location excludes its own
#'   pool row from its neighborhood. Implemented as a single call with
#'   self-exclusion. Fast and the most granular form of CV.
#' - `"kfold"`: Pool is partitioned into `n_folds` folds (or user-
#'   supplied `fold_id`). Each fold's locations are predicted using the
#'   remaining folds as the pool. Implemented as `k` separate calls with
#'   the index rebuilt per fold. Reduces optimism from spatial
#'   autocorrelation by holding out larger contiguous sets of locations
#'   (if folds are spatially blocked).
#'
#' Supported functions: [analog_impact()], [analog_regression()], and
#' [analog_search()]. Other `analog_*()` functions have no `y` input and
#' thus no prediction to validate.
#'
#' When `fun = analog_search`, residuals are computed against:
#'
#' - the `weighted_mean` column if `stat` includes `"weighted_mean"` but not
#'   `"regression"`;
#' - fitted values from regression coefficients if `stat` includes
#'   `"regression"` but not `"weighted_mean"`;
#'
#' If `stat` includes both, the prediction target is ambiguous and
#' `analog_cv()` will error. If it includes neither, residuals are skipped
#' and only the underlying search columns are returned.
#'
#' @param fun An analog function to cross-validate. Must be one of
#'   [analog_impact()], [analog_regression()], or [analog_search()] (passed
#'   as a function object, not a string).
#' @param pool The reference dataset. Matrix/data.frame with columns x, y,
#'   and climate variables, or a SpatRaster with climate variable layers.
#'   Pre-built `analog_index` objects are not supported; `analog_cv()`
#'   builds indices internally per fold (for k-fold) or once (for LOO).
#' @param y Response variable(s). Numeric vector, matrix, data.frame, or
#'   SpatRaster. Must have exactly the same number of rows/cells as `pool`.
#' @param covariates Predictor variables (required for regression; must be
#'   supplied whenever `fun` will fit local regressions).
#'   Matrix, data.frame, or SpatRaster. Must have exactly the same number of
#'   rows/cells as `pool`.
#' @param cv_method One of `"loo"` (default) or `"kfold"`.
#' @param n_folds Integer number of folds for k-fold CV. Pool rows are
#'   randomly assigned to folds. Ignored when `cv_method = "loo"` or when
#'   `fold_id` is supplied.
#' @param fold_id Optional integer vector of length `nrow(pool)` giving a
#'   fold assignment for each pool row. Overrides `n_folds`. Can be used to
#'   manually specify nonrandom folds, such as for spatial block cross-validation.
#' @param include_residuals Logical; if `TRUE` (default), the output includes
#'   `obs[_{yname}]` and `residual[_{yname}]` columns for each `y` variable
#'   (when a prediction target can be identified).
#' @param ... Additional arguments passed to `fun` (e.g., `max_clim`,
#'   `max_geog`, `kernel`, `theta`, `k`, `lambda`, `select`, `se`). Note:
#'   `fun` must accept `exclude_self` (directly or via `...`); [analog_search()]
#'   accepts it as a named parameter, and the wrapper helpers forward it via
#'   their own `...`.
#'
#' @return A data.frame or SpatRaster (matching the format of `pool`) with
#'   one row per pool location, containing all variables that `fun` would
#'   return, plus:
#'
#'   - `obs` / `obs_{yname}`: observed y value at this location (when
#'     `include_residuals = TRUE` and a prediction target is identified).
#'   - `residual` / `residual_{yname}`: observed minus held-out prediction.
#'   - `fold`: fold assignment (k-fold only).
#'
#'   Rows are ordered to match `pool`'s input row order.
#'
#' @examples
#' \dontrun{
#' # LOO for AIM
#' cv <- analog_cv(
#'   fun      = analog_impact,
#'   pool     = sites,
#'   y        = sites$biomass,
#'   max_clim = 0.5,
#'   max_geog = 100,
#'   kernel   = "gaussian_clim",
#'   theta    = 0.2
#' )
#' rmse <- sqrt(mean(cv$residual^2, na.rm = TRUE))
#'
#' # 10-fold CV for local regression
#' cv_reg <- analog_cv(
#'   fun         = analog_regression,
#'   pool        = sites,
#'   y           = sites$income,
#'   covariates  = data.frame(education = sites$edu),
#'   select      = "knn_geog",
#'   k           = 50,
#'   kernel      = "gaussian_geog",
#'   theta       = 20,
#'   cv_method   = "kfold",
#'   n_folds     = 10
#' )
#'
#' # Power-user CV via analog_search with a custom stat
#' cv_custom <- analog_cv(
#'   fun      = analog_search,
#'   pool     = sites,
#'   y        = sites$biomass,
#'   select   = "knn_clim",
#'   k        = 10,
#'   stat     = c("count", "weighted_mean"),
#'   kernel   = "gaussian_clim",
#'   theta    = 0.3
#' )
#' }
#'
#' @seealso [analog_search()], [analog_impact()], [analog_regression()].
#'
#' @export
analog_cv <- function(fun,
                      pool,
                      y,
                      covariates = NULL,
                      cv_method = c("loo", "kfold"),
                      n_folds = NULL,
                      fold_id = NULL,
                      include_residuals = TRUE,
                      ...) {

      cv_method <- match.arg(cv_method)

      # Validate fun ----------------------------------------------------------
      if (!is.function(fun)) {
            stop("`fun` must be a function (e.g., analog_impact, ",
                 "analog_regression, or analog_search).",
                 call. = FALSE)
      }
      supported <- list(
            analog_impact     = analog_impact,
            analog_regression = analog_regression,
            analog_search     = analog_search
      )
      fun_name <- NULL
      for (nm in names(supported)) {
            if (identical(fun, supported[[nm]])) { fun_name <- nm; break }
      }
      if (is.null(fun_name)) {
            stop(
                  "`fun` must be analog_impact, analog_regression, or analog_search. ",
                  "Other analog_* functions do not have a response variable to ",
                  "validate against.",
                  call. = FALSE
            )
      }

      # Validate pool / y / covariates ----------------------------------------
      if (is_analog_index(pool)) {
            stop("`pool` must not be a pre-built analog_index; pass raw pool data.",
                 call. = FALSE)
      }

      if (missing(y) || is.null(y)) {
            stop("`y` is required for cross-validation.", call. = FALSE)
      }

      if (fun_name == "analog_regression" && is.null(covariates)) {
            stop("`covariates` is required when `fun = analog_regression`.",
                 call. = FALSE)
      }

      # Coerce to matrix form -------------------------------------------------
      # Consistent across LOO and k-fold paths; also capture raster template
      # for output reassembly.
      pool_mat <- .format_data(pool)
      raster_template <- attr(pool_mat, "template")
      n_pool <- nrow(pool_mat)

      y_info <- .cv_coerce_values(y, "y", n_pool)
      y_mat <- y_info$mat
      y_names <- y_info$names

      cov_mat <- NULL
      cov_names <- NULL
      if (!is.null(covariates)) {
            cov_info <- .cv_coerce_values(covariates, "covariates", n_pool)
            cov_mat <- cov_info$mat
            cov_names <- cov_info$names
      }

      # Validate fold args / build fold assignments ---------------------------
      if (cv_method == "loo") {
            if (!is.null(n_folds) || !is.null(fold_id)) {
                  stop("`n_folds` and `fold_id` must be NULL when `cv_method = 'loo'`.",
                       call. = FALSE)
            }
            folds <- NULL
      } else {
            if (!is.null(fold_id)) {
                  if (length(fold_id) != n_pool) {
                        stop("`fold_id` must have length equal to nrow(pool).",
                             call. = FALSE)
                  }
                  folds <- as.integer(fold_id)
            } else {
                  if (is.null(n_folds)) {
                        stop("Either `n_folds` or `fold_id` must be supplied for kfold CV.",
                             call. = FALSE)
                  }
                  if (!is.numeric(n_folds) || length(n_folds) != 1L ||
                      n_folds < 2L || n_folds > n_pool) {
                        stop("`n_folds` must be an integer between 2 and nrow(pool).",
                             call. = FALSE)
                  }
                  n_folds <- as.integer(n_folds)
                  folds <- rep(seq_len(n_folds), length.out = n_pool)
                  folds <- sample(folds)
            }
      }

      # Disallow progress -----------------------------------------------------
      # exclude_self is incompatible with chunking, so progress bars are not
      # available under CV.
      extra <- list(...)
      if (isTRUE(extra$progress)) {
            stop("`progress = TRUE` is not supported inside analog_cv().",
                 call. = FALSE)
      }
      extra$progress <- FALSE

      # Determine prediction target -------------------------------------------
      # For analog_search as fun, this runs before dispatch so ambiguous
      # configurations error early rather than after running the search.
      pred_target <- .cv_determine_prediction_target(fun_name, extra, cov_mat)

      # Dispatch to LOO or k-fold path ----------------------------------------
      # For fun = analog_regression, we pass x_covariates through so that the
      # returned data.frame already carries `pred[_{yname}]` columns, avoiding
      # the need to redo the arithmetic on the R side afterward.
      if (cv_method == "loo") {
            res_df <- .cv_run_loo(
                  fun = fun,
                  fun_name = fun_name,
                  pool_mat = pool_mat,
                  y_mat = y_mat,
                  cov_mat = cov_mat,
                  extra = extra
            )
      } else {
            res_df <- .cv_run_kfold(
                  fun = fun,
                  fun_name = fun_name,
                  pool_mat = pool_mat,
                  y_mat = y_mat,
                  cov_mat = cov_mat,
                  folds = folds,
                  extra = extra
            )
      }

      # Append obs and residual columns ---------------------------------------
      # Only when a prediction target is known; otherwise skip with a message.
      if (isTRUE(include_residuals)) {
            if (pred_target == "none") {
                  message("analog_cv: no prediction-producing stat found in ",
                          "`stat`; skipping obs/residual columns. Pass stat ",
                          "including 'weighted_mean' or 'regression' to enable ",
                          "residuals.")
            } else {
                  preds <- .cv_extract_predictions(
                        pred_target = pred_target,
                        fun_name    = fun_name,
                        result_df   = res_df,
                        y_names     = y_names,
                        cov_mat     = cov_mat,
                        cov_names   = cov_names
                  )

                  n_y <- length(y_names)
                  for (v in seq_len(n_y)) {
                        if (n_y == 1L) {
                              obs_col <- "obs"
                              res_col <- "residual"
                        } else {
                              obs_col <- paste0("obs_", y_names[v])
                              res_col <- paste0("residual_", y_names[v])
                        }
                        res_df[[obs_col]] <- y_mat[, v]
                        res_df[[res_col]] <- y_mat[, v] - preds[, v]
                  }
            }
      }

      # Attach CV metadata as attributes ---------------------------------------
      # analog_summary() reads these to describe the CV run alongside the
      # normal search parameter summary. Attached to data.frame; the raster
      # branch below drops attrs by design.
      attr(res_df, "cv_method") <- cv_method
      if (cv_method == "kfold") {
            attr(res_df, "n_folds") <- length(unique(folds))
      }
      attr(res_df, "cv_fun") <- fun_name
      attr(res_df, "cv_pred_target") <- pred_target

      # Raster-in / raster-out ------------------------------------------------
      # Reassemble if pool was a SpatRaster.
      if (!is.null(raster_template)) {
            return(.cv_to_raster(res_df, raster_template))
      }

      res_df
}




# Internal helpers -----------------------------------------------------------

#' Determine the prediction target column/mechanism for a CV run
#'
#' Returns one of "weighted_mean", "regression", or "none".
#' Errors on ambiguous analog_search stats (both weighted_mean and regression).
#'
#' @keywords internal
.cv_determine_prediction_target <- function(fun_name, extra, cov_mat) {
      if (fun_name == "analog_impact") {
            # analog_impact default stat includes weighted_mean.
            stat_val <- extra$stat %||% c("count", "sum_weights", "weighted_mean")
            if (!"weighted_mean" %in% stat_val) {
                  stop("analog_cv with fun = analog_impact requires 'weighted_mean' ",
                       "to be in `stat` (it is by default).", call. = FALSE)
            }
            return("weighted_mean")
      }

      if (fun_name == "analog_regression") {
            # analog_regression always forces 'regression' into stat.
            if (is.null(cov_mat)) {
                  stop("Internal error: analog_regression requires covariates.",
                       call. = FALSE)
            }
            return("regression")
      }

      # analog_search: infer from user-provided stat
      stat_val <- extra$stat
      if (is.null(stat_val)) stat_val <- character(0)
      has_wm  <- "weighted_mean" %in% stat_val
      has_reg <- "regression"    %in% stat_val

      if (has_wm && has_reg) {
            stop(
                  "analog_cv with fun = analog_search: stat includes both ",
                  "'weighted_mean' and 'regression'. The prediction target for ",
                  "residuals is ambiguous; pass stat with only one of them, ",
                  "or use fun = analog_impact / analog_regression.",
                  call. = FALSE
            )
      }

      if (has_reg) {
            if (is.null(cov_mat)) {
                  stop("analog_cv with stat = 'regression' requires covariates.",
                       call. = FALSE)
            }
            return("regression")
      }

      if (has_wm) {
            return("weighted_mean")
      }

      "none"
}


#' Coerce a y or covariates argument to a numeric matrix
#' @keywords internal
.cv_coerce_values <- function(obj, arg_name, n_pool) {
      if (inherits(obj, "SpatRaster")) {
            if (!requireNamespace("terra", quietly = TRUE)) {
                  stop("Package 'terra' is required for SpatRaster inputs",
                       call. = FALSE)
            }
            nms <- terra::names(obj)
            mat <- as.matrix(terra::as.data.frame(obj, xy = FALSE, na.rm = FALSE))
      } else if (is.data.frame(obj)) {
            nms <- names(obj)
            mat <- as.matrix(obj)
      } else if (is.matrix(obj)) {
            nms <- colnames(obj)
            mat <- obj
      } else if (is.numeric(obj) && is.null(dim(obj))) {
            nms <- arg_name
            mat <- matrix(obj, ncol = 1L)
            colnames(mat) <- nms
      } else {
            stop("`", arg_name, "` must be a numeric vector, matrix, data.frame, ",
                 "or SpatRaster.", call. = FALSE)
      }

      storage.mode(mat) <- "double"
      if (nrow(mat) != n_pool) {
            stop("`", arg_name, "` must have exactly ", n_pool,
                 " rows/cells to match pool.", call. = FALSE)
      }
      if (is.null(nms) || any(!nzchar(nms))) {
            nms <- paste0(arg_name, seq_len(ncol(mat)))
            colnames(mat) <- nms
      }
      list(mat = mat, names = nms)
}


#' Run the LOO path: single call with exclude_self = TRUE
#'
#' For fun = analog_regression, passes x_covariates = cov_mat so predictions
#' come back in the result (LOO: focals = pool, so focal-side covariates
#' are just cov_mat).
#'
#' @keywords internal
.cv_run_loo <- function(fun, fun_name, pool_mat, y_mat, cov_mat, extra) {
      args <- c(
            list(
                  x    = pool_mat,
                  pool = pool_mat,
                  y    = y_mat,
                  exclude_self = TRUE
            ),
            extra
      )
      if (!is.null(cov_mat)) args$covariates <- cov_mat
      if (fun_name == "analog_regression" && !is.null(cov_mat)) {
            args$x_covariates <- cov_mat
      }
      do.call(fun, args)
}


#' Run the k-fold path: loop over folds, rebuilding the index per fold
#'
#' For fun = analog_regression, passes x_covariates = cov_mat[focal_rows, ]
#' so per-fold predictions come back in the result, avoiding downstream
#' re-computation.
#'
#' @keywords internal
.cv_run_kfold <- function(fun, fun_name, pool_mat, y_mat, cov_mat, folds, extra) {
      unique_folds <- sort(unique(folds))

      per_fold <- vector("list", length(unique_folds))
      for (i in seq_along(unique_folds)) {
            fid <- unique_folds[i]
            focal_rows <- which(folds == fid)
            train_rows <- which(folds != fid)

            args <- c(
                  list(
                        x    = pool_mat[focal_rows, , drop = FALSE],
                        pool = pool_mat[train_rows, , drop = FALSE],
                        y    = y_mat[train_rows, , drop = FALSE]
                  ),
                  extra
            )
            if (!is.null(cov_mat)) {
                  args$covariates <- cov_mat[train_rows, , drop = FALSE]
            }
            if (fun_name == "analog_regression" && !is.null(cov_mat)) {
                  args$x_covariates <- cov_mat[focal_rows, , drop = FALSE]
            }

            res_i <- do.call(fun, args)
            if (!is.data.frame(res_i)) {
                  stop("Internal error: expected data.frame from ", fun_name,
                       " in k-fold path; got ", class(res_i)[1],
                       ". Did a SpatRaster slip through?", call. = FALSE)
            }

            res_i$`.orig_row` <- focal_rows
            res_i$fold <- as.integer(fid)
            per_fold[[i]] <- res_i
      }

      # rbind and reorder. All per-fold data.frames should have the same columns.
      all_cols <- unique(unlist(lapply(per_fold, names)))
      per_fold <- lapply(per_fold, function(df) {
            missing_cols <- setdiff(all_cols, names(df))
            for (mc in missing_cols) df[[mc]] <- NA_real_
            df[, all_cols, drop = FALSE]
      })
      combined <- do.call(rbind, per_fold)

      # Reorder to match pool input
      combined <- combined[order(combined$`.orig_row`), , drop = FALSE]

      # Replace fold-local `index` with global row numbers and drop helper column
      combined$index <- combined$`.orig_row`
      combined$`.orig_row` <- NULL
      rownames(combined) <- NULL

      combined
}


#' Extract predictions (one column per y variable) from a CV result data.frame
#'
#' Dispatch:
#'   - pred_target == "weighted_mean": read weighted_mean[_{yname}] columns.
#'   - pred_target == "regression" AND fun_name == "analog_regression":
#'       read pred[_{yname}] columns that analog_regression() attached
#'       via x_covariates.
#'   - pred_target == "regression" AND fun_name == "analog_search":
#'       compute via .predict_from_coefs(), since analog_search does not
#'       accept x_covariates.
#'
#' @keywords internal
.cv_extract_predictions <- function(pred_target, fun_name, result_df,
                                    y_names, cov_mat, cov_names) {
      n_focal <- nrow(result_df)
      n_y <- length(y_names)
      preds <- matrix(NA_real_, nrow = n_focal, ncol = n_y)
      colnames(preds) <- y_names

      if (pred_target == "weighted_mean") {
            for (v in seq_len(n_y)) {
                  col <- if (n_y == 1L) "weighted_mean" else
                        paste0("weighted_mean_", y_names[v])
                  if (!col %in% names(result_df)) {
                        stop("Expected prediction column '", col,
                             "' not found in CV output. ",
                             "Was `weighted_mean` excluded from `stat`?",
                             call. = FALSE)
                  }
                  preds[, v] <- result_df[[col]]
            }
            return(preds)
      }

      if (pred_target == "regression" && fun_name == "analog_regression") {
            # analog_regression supplies pred columns via x_covariates.
            for (v in seq_len(n_y)) {
                  col <- if (n_y == 1L) "pred" else paste0("pred_", y_names[v])
                  if (!col %in% names(result_df)) {
                        stop("Expected prediction column '", col,
                             "' not found in CV output. Internal error: ",
                             "x_covariates was not forwarded to analog_regression.",
                             call. = FALSE)
                  }
                  preds[, v] <- result_df[[col]]
            }
            return(preds)
      }

      if (pred_target == "regression") {
            # analog_search path: compute predictions from coefficient columns.
            if (is.null(cov_mat) || is.null(cov_names)) {
                  stop("Internal error: regression prediction requires covariates.",
                       call. = FALSE)
            }
            cov_for_focals <- cov_mat[result_df$index, , drop = FALSE]
            preds_mat <- .predict_from_coefs(
                  coefs_df          = result_df,
                  covariates_matrix = cov_for_focals,
                  y_names           = y_names,
                  cov_names         = cov_names
            )
            preds[, ] <- preds_mat[, , drop = FALSE]
            return(preds)
      }

      stop("Internal error: unsupported pred_target '", pred_target, "'.",
           call. = FALSE)
}


#' Reassemble CV result as a SpatRaster when pool was a SpatRaster
#' @keywords internal
.cv_to_raster <- function(res_df, template) {
      if (!requireNamespace("terra", quietly = TRUE)) {
            stop("Package 'terra' required for SpatRaster output.", call. = FALSE)
      }
      drop_cols <- intersect(c("index", "x", "y"), names(res_df))
      num_cols <- setdiff(names(res_df), drop_cols)
      num_cols <- num_cols[vapply(res_df[num_cols], is.numeric, logical(1))]

      layers <- lapply(num_cols, function(cn) {
            r <- terra::setValues(template, res_df[[cn]])
            names(r) <- cn
            r
      })
      do.call(c, layers)
}
