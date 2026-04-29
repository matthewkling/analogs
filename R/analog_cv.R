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
#'   `"regression"` or `"tabulate"`;
#' - fitted values from regression coefficients if `stat` includes
#'   `"regression"` but not `"weighted_mean"` or `"tabulate"`;
#' - per-class weighted votes if `stat` includes `"tabulate"` (categorical
#'   y; produces Brier score and primary-class label rather than a numeric
#'   residual).
#'
#' If `stat` includes more than one of these, the prediction target is
#' ambiguous and `analog_cv()` will error. If it includes none, residuals
#' are skipped and only the underlying search columns are returned.
#'
#' For categorical CV (`stat = "tabulate"`), `y` must be a factor or
#' coercible-to-factor input (character, integer codes, single-layer
#' categorical SpatRaster). The output uses different residual-equivalent
#' columns: see `@return` below.
#'
#' @param fun An analog function to cross-validate. Must be one of
#'   [analog_impact()], [analog_regression()], or [analog_search()] (passed
#'   as a function object, not a string).
#' @param pool The reference dataset. Matrix/data.frame with columns x, y,
#'   and climate variables, or a SpatRaster with climate variable layers.
#'   Pre-built `analog_index` objects are not supported; `analog_cv()`
#'   builds indices internally per fold (for k-fold) or once (for LOO).
#' @param y Response variable(s). For continuous prediction targets
#'   (`weighted_mean`, `regression`): numeric vector, matrix, data.frame,
#'   or SpatRaster. For categorical (`tabulate`): factor or coercible-to-
#'   factor vector / character / integer / matrix / data.frame /
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
#'   per-focal residual-equivalent columns (see `@return`).
#' @param ... Additional arguments passed to `fun` (e.g., `max_clim`,
#'   `max_geog`, `kernel`, `theta`, `k`, `lambda`, `select`, `se`). Note:
#'   `fun` must accept `exclude_self` (directly or via `...`); [analog_search()]
#'   accepts it as a named parameter, and the wrapper helpers forward it via
#'   their own `...`.
#'
#' @return A data.frame or SpatRaster (matching the format of `pool`) with
#'   one row per pool location, containing all variables that `fun` would
#'   return, plus the following residual-equivalent columns when
#'   `include_residuals = TRUE` and a prediction target is identified.
#'
#'   For continuous prediction targets (`weighted_mean`, `regression`):
#'
#'   - `obs` / `obs_{yname}`: observed y value at this location.
#'   - `residual` / `residual_{yname}`: observed minus held-out prediction.
#'
#'   For categorical prediction target (`tabulate`):
#'
#'   - `obs` / `obs_{yname}`: observed class label (character).
#'   - `primary` / `primary_{yname}`: predicted (modal) class label
#'     (character; argmax across the per-class vote columns).
#'   - `brier` / `brier_{yname}`: per-focal Brier score, computed on
#'     row-normalized vote shares (range 0-2).
#'
#'   The underlying `n_<level>` vote columns from the analog search are
#'   also retained, so users can compute additional metrics (entropy,
#'   top-k accuracy, custom losses) in postprocessing.
#'
#'   Always present:
#'
#'   - `fold`: fold assignment (k-fold only).
#'
#'   Rows are ordered to match `pool`'s input row order. For SpatRaster
#'   output, character columns (`obs`, `primary`) are dropped with a
#'   message; pass `pool` as a data.frame to retain them.
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
#' # LOO categorical CV (vegetation projection)
#' cv_veg <- analog_cv(
#'   fun      = analog_impact,
#'   pool     = sites,
#'   y        = factor(sites$vegetation_type),
#'   stat     = "tabulate",
#'   max_clim = 0.5,
#'   max_geog = 100,
#'   kernel   = "gaussian_clim",
#'   theta    = 0.2
#' )
#' # Per-focal Brier score and primary class are in cv_veg$brier and
#' # cv_veg$primary; the n_<level> vote columns are also retained.
#' accuracy <- mean(cv_veg$primary == cv_veg$obs, na.rm = TRUE)
#' mean_brier <- mean(cv_veg$brier, na.rm = TRUE)
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

      # Pre-parse extra args (before y coercion) so we can detect tabulate
      # mode and route y through categorical coercion if needed.
      extra <- list(...)
      if (isTRUE(extra$progress)) {
            stop("`progress = TRUE` is not supported inside analog_cv().",
                 call. = FALSE)
      }
      extra$progress <- FALSE

      # Detect whether this is a tabulate-mode CV run. Different y semantics
      # (categorical vs continuous), different residual machinery downstream.
      is_tabulate <- .cv_is_tabulate_run(fun_name, extra)

      # y coercion: branched on tabulate mode.
      # Continuous path produces a numeric matrix; tabulate path produces
      # a list of factored vectors plus a captured global level set.
      y_levels_list <- NULL
      if (is_tabulate) {
            y_info <- .cv_coerce_categorical_y(y, n_pool)
            # For inner calls, package factored y into a data.frame (multi-y)
            # or use the bare factor (single-y). The analog_search-side
            # categorical helper accepts either.
            y_for_inner <- if (y_info$n_y == 1L) {
                  y_info$factored[[1L]]
            } else {
                  stats::setNames(
                        as.data.frame(y_info$factored, stringsAsFactors = FALSE),
                        y_info$names
                  )
            }
            y_names <- y_info$names
            y_levels_list <- y_info$levels_list
            # y_mat is also built as integer codes for downstream obs lookup
            # (one row per pool location, integer factor codes, NA preserved).
            y_mat <- do.call(cbind, lapply(y_info$factored,
                                           function(f) as.integer(f)))
            colnames(y_mat) <- y_names
            storage.mode(y_mat) <- "double"
      } else {
            y_info <- .cv_coerce_values(y, "y", n_pool)
            y_mat <- y_info$mat
            y_names <- y_info$names
            y_for_inner <- y_mat
      }

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
                  y_for_inner = y_for_inner,
                  cov_mat = cov_mat,
                  extra = extra
            )
      } else {
            res_df <- .cv_run_kfold(
                  fun = fun,
                  fun_name = fun_name,
                  pool_mat = pool_mat,
                  y_for_inner = y_for_inner,
                  y_levels_list = y_levels_list,
                  y_names = y_names,
                  cov_mat = cov_mat,
                  folds = folds,
                  extra = extra
            )
      }

      # Append obs and residual columns ---------------------------------------
      # Only when a prediction target is known; otherwise skip with a message.
      # Tabulate uses a different shape than continuous (Brier + primary
      # rather than residual), so it gets its own branch.
      if (isTRUE(include_residuals)) {
            if (pred_target == "none") {
                  message("analog_cv: no prediction-producing stat found in ",
                          "`stat`; skipping obs/residual columns. Pass stat ",
                          "including 'weighted_mean', 'regression', or ",
                          "'tabulate' to enable residuals.")

            } else if (pred_target == "tabulate") {
                  res_df <- .cv_append_tabulate_columns(
                        res_df       = res_df,
                        y_mat        = y_mat,
                        y_names      = y_names,
                        y_levels_list = y_levels_list
                  )

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
      # analog_metadata() reads these to describe the CV run alongside the
      # normal search parameter summary. Attached to data.frame; the raster
      # branch below drops attrs by design.
      attr(res_df, "cv_method") <- cv_method
      if (cv_method == "kfold") {
            attr(res_df, "cv_n_folds") <- length(unique(folds))
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

# Lightweight check: does this CV run involve stat = "tabulate"?
# Used early in analog_cv() to route y through categorical coercion before
# the inner calls. Mirrors the tabulate-detection inside
# .cv_determine_prediction_target() but doesn't run the full prediction-target
# resolution (which depends on cov_mat etc., not yet built at this point).
.cv_is_tabulate_run <- function(fun_name, extra) {
      if (fun_name == "analog_impact") {
            stat_val <- extra$stat %||% c("count", "sum_weights", "weighted_mean")
      } else if (fun_name == "analog_regression") {
            return(FALSE)
      } else {
            stat_val <- extra$stat %||% character(0)
      }
      "tabulate" %in% stat_val
}


# Determine the prediction target column/mechanism for a CV run
#
# Returns one of "weighted_mean", "regression", "tabulate", or "none".
# Errors on ambiguous analog_search stats (more than one of weighted_mean,
# regression, tabulate).
.cv_determine_prediction_target <- function(fun_name, extra, cov_mat) {
      if (fun_name == "analog_impact") {
            stat_val <- extra$stat %||% c("count", "sum_weights", "weighted_mean")
            has_tab <- "tabulate" %in% stat_val
            has_wm  <- "weighted_mean" %in% stat_val
            if (has_tab && has_wm) {
                  # Should have been caught by .validate_query_params, but re-check
                  stop("analog_cv: stat cannot include both 'tabulate' and ",
                       "'weighted_mean'. Use one or the other.",
                       call. = FALSE)
            }
            if (has_tab) return("tabulate")
            if (!has_wm) {
                  stop("analog_cv with fun = analog_impact requires 'weighted_mean' ",
                       "or 'tabulate' to be in `stat` (weighted_mean is default).",
                       call. = FALSE)
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
      has_tab <- "tabulate" %in% stat_val
      has_wm  <- "weighted_mean" %in% stat_val
      has_reg <- "regression"    %in% stat_val

      n_targets <- sum(has_tab, has_wm, has_reg)
      if (n_targets > 1L) {
            stop(
                  "analog_cv with fun = analog_search: stat includes more than ",
                  "one prediction-producing aggregator (among 'tabulate', ",
                  "'weighted_mean', 'regression'). The prediction target for ",
                  "residuals is ambiguous; pass stat with only one of them, ",
                  "or use the corresponding wrapper (analog_impact / ",
                  "analog_regression).",
                  call. = FALSE
            )
      }

      if (has_tab) return("tabulate")

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


# Coerce a y or covariates argument to a numeric matrix
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


# Coerce a categorical y argument for tabulate-mode CV.
#
# Accepts factor / character / integer / matrix / data.frame / SpatRaster
# inputs (anything analog_search's tabulate path accepts). Returns a list
# with:
#   factored      list of factor vectors (one per y variable), with global
#                 levels intact - this is what gets passed to inner calls
#                 so per-fold droplevels mostly preserves the level set.
#   names         character vector of y variable names.
#   levels_list   list of character vectors giving each variable's levels.
#                 Used downstream to compute primary class labels and to
#                 align tabulate-vote columns across folds.
#   n_y           number of y variables.
.cv_coerce_categorical_y <- function(obj, n_pool) {
      cols <- NULL
      var_names <- NULL

      if (inherits(obj, "SpatRaster")) {
            if (!requireNamespace("terra", quietly = TRUE)) {
                  stop("Package 'terra' is required for SpatRaster `y`.",
                       call. = FALSE)
            }
            df <- terra::as.data.frame(obj, xy = FALSE, na.rm = FALSE)
            if (nrow(df) != n_pool) {
                  stop("`y` SpatRaster has ", nrow(df), " cells but pool has ",
                       n_pool, " rows. They must match.", call. = FALSE)
            }
            var_names <- names(df)
            cols <- as.list(df)

      } else if (is.factor(obj)) {
            if (length(obj) != n_pool) {
                  stop("`y` length is ", length(obj), " but pool has ", n_pool,
                       " rows. They must match.", call. = FALSE)
            }
            cols <- list(obj)
            var_names <- "y"

      } else if (is.data.frame(obj)) {
            if (nrow(obj) != n_pool) {
                  stop("`y` has ", nrow(obj), " rows but pool has ", n_pool,
                       " rows. They must match.", call. = FALSE)
            }
            var_names <- names(obj)
            cols <- as.list(obj)

      } else if (is.matrix(obj)) {
            if (nrow(obj) != n_pool) {
                  stop("`y` has ", nrow(obj), " rows but pool has ", n_pool,
                       " rows. They must match.", call. = FALSE)
            }
            var_names <- colnames(obj)
            cols <- lapply(seq_len(ncol(obj)), function(j) obj[, j])

      } else if (is.atomic(obj) && is.null(dim(obj))) {
            if (length(obj) != n_pool) {
                  stop("`y` length is ", length(obj), " but pool has ", n_pool,
                       " rows. They must match.", call. = FALSE)
            }
            cols <- list(obj)
            var_names <- "y"

      } else {
            stop("`y` must be a factor, character/integer vector, matrix, ",
                 "data.frame, or SpatRaster for stat = 'tabulate' CV.",
                 call. = FALSE)
      }

      n_y <- length(cols)

      if (is.null(var_names) || any(!nzchar(var_names))) {
            var_names <- if (n_y == 1L) "y" else paste0("y", seq_len(n_y))
      }

      factored <- vector("list", n_y)
      levels_list <- vector("list", n_y)
      for (v in seq_len(n_y)) {
            f <- if (is.factor(cols[[v]])) {
                  droplevels(cols[[v]])  # match the analog_search-side helper
            } else {
                  factor(cols[[v]])
            }
            if (length(levels(f)) == 0L) {
                  stop("`y` column ", var_names[v], " has no non-NA values. ",
                       "Cannot run categorical CV on an all-NA variable.",
                       call. = FALSE)
            }
            factored[[v]] <- f
            levels_list[[v]] <- levels(f)
      }

      list(factored = factored, names = var_names,
           levels_list = levels_list, n_y = n_y)
}


# Run the LOO path: single call with exclude_self = TRUE
#
# For fun = analog_regression, passes x_covariates = cov_mat so predictions
# come back in the result (LOO: focals = pool, so focal-side covariates
# are just cov_mat).
.cv_run_loo <- function(fun, fun_name, pool_mat, y_for_inner, cov_mat, extra) {
      args <- c(
            list(
                  x    = pool_mat,
                  pool = pool_mat,
                  y    = y_for_inner,
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


# Run the k-fold path: loop over folds, rebuilding the index per fold
#
# For fun = analog_regression, passes `x_covariates = cov_mat[focal_rows, ]`
# so per-fold predictions come back in the result, avoiding downstream
# re-computation.
.cv_run_kfold <- function(fun, fun_name, pool_mat, y_for_inner,
                          y_levels_list, y_names, cov_mat, folds, extra) {
      unique_folds <- sort(unique(folds))

      # Helper: subset y_for_inner by row, preserving type. Three shapes
      # are possible:
      #   - numeric matrix (continuous y)
      #   - factor (categorical, single y)
      #   - data.frame of factors (categorical, multi-y)
      subset_y <- function(idx) {
            if (is.factor(y_for_inner)) {
                  y_for_inner[idx]
            } else if (is.data.frame(y_for_inner)) {
                  y_for_inner[idx, , drop = FALSE]
            } else {
                  y_for_inner[idx, , drop = FALSE]   # numeric matrix
            }
      }

      # Pre-compute the expected tabulate vote-column names for column
      # alignment after rbind. Only relevant for tabulate runs (signaled
      # by y_levels_list != NULL); empty for continuous CV.
      tabulate_vote_cols <- character(0)
      if (!is.null(y_levels_list)) {
            n_y <- length(y_names)
            for (v in seq_len(n_y)) {
                  lev <- y_levels_list[[v]]
                  prefix <- if (n_y == 1L && identical(y_names[v], "y")) {
                        "n_"   # bare unnamed single y matches the
                        # analog_search/tabulate column-naming rule
                  } else {
                        paste0(y_names[v], "_n_")
                  }
                  tabulate_vote_cols <- c(tabulate_vote_cols,
                                          paste0(prefix, lev))
            }
      }

      per_fold <- vector("list", length(unique_folds))
      for (i in seq_along(unique_folds)) {
            fid <- unique_folds[i]
            focal_rows <- which(folds == fid)
            train_rows <- which(folds != fid)

            args <- c(
                  list(
                        x    = pool_mat[focal_rows, , drop = FALSE],
                        pool = pool_mat[train_rows, , drop = FALSE],
                        y    = subset_y(train_rows)
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

      # rbind and reorder. Before rbind we align columns: per-fold dropped
      # tabulate levels (rare but possible when a fold's training subset is
      # missing a level entirely) get filled with 0 (no analogs in that
      # class for those focals); other missing columns get NA.
      all_cols <- unique(unlist(lapply(per_fold, names)))
      # Union with the expected tabulate columns from the global level set.
      # This catches the case where a level is missing from EVERY fold.
      all_cols <- unique(c(all_cols, tabulate_vote_cols))
      per_fold <- lapply(per_fold, function(df) {
            missing_cols <- setdiff(all_cols, names(df))
            for (mc in missing_cols) {
                  df[[mc]] <- if (mc %in% tabulate_vote_cols) 0 else NA_real_
            }
            df[, all_cols, drop = FALSE]
      })

      # Capture inner attributes from the first per-fold result before rbind,
      # which strips them. All folds share the same parameterization, so the
      # first fold's attrs are representative for the combined result.
      inner_attrs <- attributes(per_fold[[1]])
      inner_attrs[c("names", "row.names", "class")] <- NULL

      combined <- do.call(rbind, per_fold)

      # Reorder to match pool input
      combined <- combined[order(combined$`.orig_row`), , drop = FALSE]

      # Replace fold-local `index` with global row numbers and drop helper column
      combined$index <- combined$`.orig_row`
      combined$`.orig_row` <- NULL
      rownames(combined) <- NULL

      # Restore inner attributes (select, stat, kernel, theta, n_x, etc.).
      # CV-specific attributes (cv_method, cv_fun, ...) are written by the
      # outer analog_cv() wrapper after this returns.
      for (nm in names(inner_attrs)) {
            attr(combined, nm) <- inner_attrs[[nm]]
      }

      # n_x and n_pool from the first fold reflect that fold's sizes, not the
      # totals across the combined CV result. Override with totals.
      attr(combined, "n_x") <- nrow(pool_mat)
      attr(combined, "n_pool") <- nrow(pool_mat)

      combined
}


# Extract predictions (one column per y variable) from a CV result data.frame
#
# Dispatch:
#   - pred_target == "weighted_mean": read `weighted_mean[_{yname}]` columns.
#   - pred_target == "regression" AND fun_name == "analog_regression":
#       read `pred[_{yname}]` columns that analog_regression() attached
#       via x_covariates.
#   - pred_target == "regression" AND fun_name == "analog_search":
#       compute via .predict_from_coefs(), since analog_search does not
#       accept x_covariates.
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


# Append per-focal categorical-CV columns to the result data.frame.
#
# Adds:
#   obs_<y>     character: the observed class label per focal (NA for NA y).
#   primary_<y> character: the predicted class label (argmax of vote
#               columns; NA when all vote columns are zero / NA).
#   brier_<y>   numeric: per-focal Brier score, computed on row-normalized
#               vote shares. NA when the row sum of votes is 0 (no analogs
#               in any class) or the observed value is NA.
#
# Single-y uses bare names (obs, primary, brier); multi-y uses _<varname>
# suffixes, matching the convention used by the continuous pred_target paths.
#
# y_mat carries integer factor codes (1-based, NA preserved). y_levels_list
# carries the level names for each y variable so codes can be back-mapped
# to labels.
.cv_append_tabulate_columns <- function(res_df, y_mat, y_names, y_levels_list) {
      n_y <- length(y_names)

      for (v in seq_len(n_y)) {
            lev <- y_levels_list[[v]]
            K_v <- length(lev)

            # Identify this variable's vote columns. Naming follows the
            # query_analog_index.R convention: bare unnamed single y uses
            # n_<level>, otherwise <varname>_n_<level>. We probe for both
            # in case the analog_search-side helper used a different name.
            prefix_unnamed <- "n_"
            prefix_named   <- paste0(y_names[v], "_n_")

            cand_named <- paste0(prefix_named, lev)
            cand_unnamed <- paste0(prefix_unnamed, lev)

            if (all(cand_named %in% names(res_df))) {
                  vote_cols <- cand_named
            } else if (all(cand_unnamed %in% names(res_df))) {
                  vote_cols <- cand_unnamed
            } else {
                  # Diagnostic: neither full set is present. Fail loudly so
                  # the schema mismatch is obvious rather than producing
                  # silently wrong results.
                  stop("Internal error: tabulate vote columns for variable '",
                       y_names[v], "' not found in CV result. Expected ",
                       "either ", paste(sQuote(cand_named), collapse = ", "),
                       " or ", paste(sQuote(cand_unnamed), collapse = ", "), ".",
                       call. = FALSE)
            }

            votes <- as.matrix(res_df[, vote_cols, drop = FALSE])
            # Coerce in case rbind alignment introduced character-y artifacts
            storage.mode(votes) <- "double"

            # Row-normalize to probabilities. Rows summing to 0 (no analogs
            # in any class) are non-finite and become NA per-class, which
            # propagates through Brier as NA - documented as the right
            # behavior for empty neighborhoods.
            row_sums <- rowSums(votes, na.rm = FALSE)
            probs <- votes / row_sums  # NaN where row_sums == 0

            # Predicted (primary) class: argmax across vote columns. Rows
            # with all-NA or all-zero vote sums get NA primary.
            obs_codes <- y_mat[, v]   # integer codes, 1-based, possibly NA
            primary_idx <- max.col(votes, ties.method = "first")
            all_zero_or_na <- !is.finite(row_sums) | row_sums == 0
            primary_label <- ifelse(all_zero_or_na, NA_character_,
                                    lev[primary_idx])

            # Per-focal Brier score: sum_k (p_k - I[k == obs])^2.
            # Standard textbook form on row-normalized probabilities;
            # ranges 0 to 2.
            brier <- rep(NA_real_, nrow(res_df))
            valid <- !all_zero_or_na & !is.na(obs_codes)
            if (any(valid)) {
                  obs_code_valid <- obs_codes[valid]
                  probs_valid <- probs[valid, , drop = FALSE]
                  # Build indicator: rows where column k == obs_code -> 1
                  ind <- matrix(0, nrow = sum(valid), ncol = K_v)
                  ind[cbind(seq_len(sum(valid)), obs_code_valid)] <- 1
                  brier[valid] <- rowSums((probs_valid - ind)^2)
            }

            # Observed class label as character; NA for NA codes.
            obs_label <- ifelse(is.na(obs_codes), NA_character_,
                                lev[obs_codes])

            if (n_y == 1L) {
                  obs_col <- "obs"
                  prim_col <- "primary"
                  brier_col <- "brier"
            } else {
                  obs_col <- paste0("obs_", y_names[v])
                  prim_col <- paste0("primary_", y_names[v])
                  brier_col <- paste0("brier_", y_names[v])
            }
            res_df[[obs_col]] <- obs_label
            res_df[[prim_col]] <- primary_label
            res_df[[brier_col]] <- brier
      }

      res_df
}


# Reassemble CV result as a SpatRaster when pool was a SpatRaster
.cv_to_raster <- function(res_df, template) {
      if (!requireNamespace("terra", quietly = TRUE)) {
            stop("Package 'terra' required for SpatRaster output.", call. = FALSE)
      }
      drop_cols <- intersect(c("index", "x", "y"), names(res_df))
      data_cols <- setdiff(names(res_df), drop_cols)

      # Raster layers are numeric; character / factor columns get dropped
      # silently by `terra::setValues` so flag them up front.
      char_cols <- data_cols[vapply(res_df[data_cols],
                                    function(z) is.character(z) || is.factor(z),
                                    logical(1))]
      if (length(char_cols) > 0L) {
            message("analog_cv: dropping ", length(char_cols), " non-numeric ",
                    "column(s) from raster output: ",
                    paste(char_cols, collapse = ", "),
                    ". Pass `pool` as a data.frame to retain these columns ",
                    "(e.g., observed and primary class labels for tabulate CV).")
      }

      num_cols <- data_cols[vapply(res_df[data_cols], is.numeric, logical(1))]

      layers <- lapply(num_cols, function(cn) {
            r <- terra::setValues(template, res_df[[cn]])
            names(r) <- cn
            r
      })
      do.call(c, layers)
}
