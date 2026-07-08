#' Local weighted regression across analog neighborhoods
#'
#' Fits a weighted local regression of `y` on `covariates` within
#' each focal location's analog neighborhood. Analog neighborhoods are
#' defined by environmental similarity, geographic proximity, or both, while
#' covariates capture additional predictors that influence outcomes within
#' each neighborhood but are not captured by to the search dimensions.
#' This generalizes the weighted mean — which averages over all within-neighborhood
#' variation — by resolving variation driven by these auxiliary predictors.
#' Supports ordinary and ridge-penalized weighted least squares. With purely geographic
#' neighborhoods, this is equivalent to geographically weighted regression
#' (GWR); with environment-based neighborhoods, it extends the analog impact
#' model (AIM) framework to incorporate local covariate effects. This function is
#' a wrapper that calls [analog_search()] with `"regression"` included in `stat`.
#'
#' @param x Focal locations for which regressions will be fit. Should be a
#'   matrix/data.frame with columns x, y, and environmental variables, or a
#'   SpatRaster with environmental variable layers.
#' @param y Response variable(s) to model via local regression.
#'   Can be a numeric vector, matrix, data.frame, or SpatRaster.
#'   Must have exactly the same number of rows/cells as `pool`.
#'   A separate regression is fit for each variable.
#' @param covariates Predictor variables for local regression, supplied at
#'   pool locations. Can be a numeric vector (single covariate), matrix,
#'   data.frame, or SpatRaster. Must have exactly the same number of
#'   rows/cells as `pool`. Column/layer names carry through to output.
#'   These variables are NOT used for the analog search itself -- only for
#'   regression within each neighborhood.
#' @param x_covariates Optional predictor values at focal (`x`) locations, at
#'   which to evaluate the local regressions. When supplied, the output
#'   includes one additional column (`pred`) or one column per response
#'   variable (`pred_{yname}` for multi-y), giving the fitted value at the
#'   focal. Must have exactly the same number of rows/cells as `x`, and the
#'   same column/layer names as `covariates`. Default `NULL` returns only
#'   coefficients.
#' @param lambda Ridge penalty parameter (default: 0, giving ordinary
#'   weighted least squares). Higher values shrink high-variance coefficients
#'   toward zero, causing the intercept to approach the weighted mean as
#'   `lambda -> Inf`. Useful when some neighborhoods have few analogs
#'   relative to the number of covariates, or when covariates are strongly
#'   inter-correlated.
#' @param stat Statistic(s) to compute. `"regression"` is always included.
#'   Additional stats like `"count"`, `"ess"`, and `"weighted_mean"` can
#'   be requested alongside regression coefficients. Default includes
#'   `"count"` and `"ess"` for diagnostics.
#' @inheritParams analog_search
#' @param ... additional arguments passed to `analog_search()` (usually unneeded)
#'
#' @return A data.frame (or SpatRaster if `x` is a SpatRaster) with one
#'   row per focal location containing:
#'
#'   \itemize{
#'     \item `index`, `x`, `y`: Focal location identifiers
#'     \item Columns for any additional stats requested (e.g., `count`, `ess`)
#'     \item `coef_intercept`: Regression intercept
#'     \item `coef_{name}`: Regression slope for each covariate
#'     \item When `se != "none"`: `se_intercept` and `se_{name}` giving
#'       standard errors of each coefficient
#'     \item With multiple `y` variables: columns are named
#'       `coef_{coeff}_{varname}` (and `se_{coeff}_{varname}` when SEs
#'       are returned)
#'     \item When `x_covariates` is supplied: `pred` (single-y) or
#'       `pred_{varname}` (multi-y) giving the fitted value at each focal.
#'   }

#'   See [analog_search()] for column-naming conventions across stats
#'   and [metadata()] for attached metadata attributes.
#'
#' @details
#'
#' ## Method
#'
#' For each focal location, the function:
#'
#' \enumerate{
#'   \item Selects analog pool locations based on `select`, `env`, `geog`, and `k`
#'   \item Computes distance-based kernel weights for each analog (via the
#'     `env` / `geog` kernels)
#'   \item Fits a weighted least squares regression of `y` on `covariates`
#'     using these weights, with optional ridge penalty `lambda`
#'   \item Returns the regression coefficients (intercept + slopes), and
#'     optionally their standard errors (see `se` in [analog_search()]).
#' }
#'
#' The math: `beta = (X'WX + lambda * I_p)^{-1} X'Wy`, where `W` is diagonal
#' weights, `X` is the design matrix (intercept + covariates), and `I_p` penalizes
#' covariate coefficients only (not the intercept).
#'
#' ## Relationship to Weighted Mean
#'
#' With `lambda -> Inf`, covariate coefficients shrink to zero and the
#' intercept converges to the weighted mean. With `lambda = 0` (the default),
#' the full local regression is used. If covariates are centered (weighted
#' mean zero within each neighborhood), the intercept equals the weighted
#' mean at any lambda. The `lambda` parameter thus provides smooth
#' interpolation between a simple weighted average and a full local regression.
#'
#' ## Common Configurations
#'
#' \itemize{
#'   \item **Geographically weighted regression** (`select = "all"` or
#'     `"knn_geog"`, with `max_geog` and geographic weights): Local spatial
#'     regression using geographic proximity to define and weight neighborhoods,
#'     equivalent to GWR.
#'   \item **Environmental-nearest regression** (`select = "knn_env"`, with
#'     `max_geog` and `k`): Fixed-size neighborhoods of the k most similar
#'     environments within geographic range.
#' }
#'
#' ## Prediction
#'
#' To evaluate fitted values at the focal locations, pass covariate values at
#' those locations via `x_covariates`. The output will then include a `pred`
#' column (single-y case) or `pred_{yname}` columns (multi-y case) alongside
#' the usual coefficient columns.
#'
#' If you only need coefficients, leave `x_covariates = NULL`.
#'
#' @examples
#' \dontrun{
#' # GWR-style spatial regression
#' gwr_result <- analog_regression(
#'   x = sites,
#'   pool = sites,
#'   y = sites$income,
#'   covariates = data.frame(education = sites$edu, access = sites$access),
#'   select = "knn_geog",
#'   k = 50,
#'   env = NULL,
#'   geog = kernel("gaussian", theta = 20),
#'   se = "ess"
#' )
#'
#' # With focal-side covariates, to obtain fitted predictions
#' pred_result <- analog_regression(
#'   x = future_sites,
#'   pool = sites,
#'   y = sites$income,
#'   covariates   = data.frame(education = sites$edu),       # at pool
#'   x_covariates = data.frame(education = future_sites$edu) # at focals
#' )
#' head(pred_result[, c("coef_intercept", "coef_education", "pred")])
#' }
#'
#' @seealso [analog_search()] for the underlying flexible analog search function;
#'   [analog_impact()] for the standard AIM workflow;
#'   [analog_cv()] for cross-validation of regression fits.
#'
#' @export
analog_regression <- function(
            x,
            pool,
            y,
            weight = NULL,
            covariates,
            x_covariates = NULL,
            geog = NULL,
            env = kernel("gaussian"),
            select = "all",
            k = NULL,
            lambda = 0,
            stat = c("count", "ess", "regression"),
            se = c("none", "ess", "design"),
            normalize = "auto",
            x_cov = NULL,
            coord_type = "auto",
            env_res_adj = "auto",
            geog_res_adj = "auto",
            cell_area_weight = "auto",
            n_threads = NULL,
            progress = FALSE,
            ...
) {
      se <- match.arg(se)

      # Ensure "regression" is always in stat
      if (!"regression" %in% stat) {
            stat <- unique(c(stat, "regression"))
      }

      result <- analog_search(
            x           = x,
            pool        = pool,
            select      = select,
            stat        = stat,
            y           = y,
            weight      = weight,
            covariates  = covariates,
            env         = env,
            geog        = geog,
            k           = k,
            lambda      = lambda,
            se          = se,
            normalize   = normalize,
            x_cov       = x_cov,
            coord_type  = coord_type,
            env_res_adj= env_res_adj,
            geog_res_adj = geog_res_adj,
            cell_area_weight = cell_area_weight,
            n_threads   = n_threads,
            progress    = progress,
            ...
      )

      # Append predictions at focal locations if x_covariates was supplied
      if (!is.null(x_covariates)) {
            result <- .append_regression_predictions(
                  result       = result,
                  x_covariates = x_covariates,
                  y            = y,
                  covariates   = covariates
            )
      }

      result
}


# Internal helpers -----------------------------------------------------------

#' Append `pred` column(s) to an analog_regression result
#'
#' Evaluates fitted values at focal locations by combining per-focal
#' regression coefficients in `result` with focal-side covariate values
#' in `x_covariates`. Handles both data.frame and SpatRaster `result`,
#' both single-y and multi-y regression, and both data.frame/matrix and
#' SpatRaster `x_covariates` inputs. Uses [.predict_from_coefs()] for the
#' arithmetic.
#'
#' Name resolution uses the package's existing validator helpers to pull
#' `y_names` / `cov_names` from the original `y` / `covariates` arguments,
#' avoiding any need to parse coefficient column names.
#'
#' @keywords internal
.append_regression_predictions <- function(result, x_covariates, y, covariates) {

      # Resolve y_names and cov_names from the inputs (reuses existing validators
      # to accept the same input shapes that analog_search() accepts).
      # A dummy single-row "ref" satisfies row-count validation, since we only
      # need the names here.
      dummy_ref <- matrix(0, nrow = nrow_of(y), ncol = 1L)
      y_info    <- .validate_and_format_values(y, dummy_ref)
      y_names   <- y_info$names

      dummy_ref_c <- matrix(0, nrow = nrow_of(covariates), ncol = 1L)
      cov_info    <- .validate_and_format_covariates(covariates, dummy_ref_c)
      cov_names   <- cov_info$names

      # Coerce x_covariates and harvest raster template (if any) --------------
      is_raster_result <- inherits(result, "SpatRaster")

      xc <- .coerce_xcov_to_matrix(x_covariates, cov_names)
      xc_mat <- xc$mat
      xc_names <- xc$names

      # Validate row/cell count
      n_focals <- if (is_raster_result) {
            if (!requireNamespace("terra", quietly = TRUE)) {
                  stop("Package 'terra' is required for SpatRaster result.",
                       call. = FALSE)
            }
            terra::ncell(result)
      } else {
            nrow(result)
      }
      if (nrow(xc_mat) != n_focals) {
            stop("`x_covariates` has ", nrow(xc_mat), " rows/cells but `x` has ",
                 n_focals, ". These must match.", call. = FALSE)
      }

      # Validate column name alignment
      missing_covs <- setdiff(cov_names, xc_names)
      if (length(missing_covs) > 0L) {
            stop("`x_covariates` is missing columns for covariates: ",
                 paste(missing_covs, collapse = ", "),
                 ". Column/layer names in `x_covariates` must match those in ",
                 "`covariates`.", call. = FALSE)
      }
      # Reorder to cov_names
      xc_mat <- xc_mat[, cov_names, drop = FALSE]
      storage.mode(xc_mat) <- "double"

      # Get coefficients as data.frame for the arithmetic ---------------------
      coefs_df <- if (is_raster_result) {
            terra::as.data.frame(result, xy = FALSE, na.rm = FALSE)
      } else {
            as.data.frame(result)
      }

      preds <- .predict_from_coefs(
            coefs_df          = coefs_df,
            covariates_matrix = xc_mat,
            y_names           = y_names,
            cov_names         = cov_names
      )

      # Attach predictions to result ------------------------------------------
      n_y <- length(y_names)
      pred_colnames <- if (n_y == 1L) "pred" else paste0("pred_", y_names)

      if (is_raster_result) {
            pred_layers <- lapply(seq_len(n_y), function(j) {
                  stats::setNames(terra::setValues(result[[1]], preds[, j]),
                                  pred_colnames[j])
            })
            # Preserve existing attributes by appending the new layer(s)
            att <- attributes(result)
            out <- c(result, do.call(c, pred_layers))
            # c(...) on SpatRaster drops most custom attributes; re-attach.
            attributes(out) <- append(attributes(out),
                                      att[setdiff(names(att), names(attributes(out)))])
            return(out)
      }

      # data.frame path
      for (j in seq_len(n_y)) {
            result[[pred_colnames[j]]] <- preds[, j]
      }
      result
}


#' Coerce x_covariates input to numeric matrix with column names
#'
#' Accepts numeric vector (single covariate), matrix, data.frame, or
#' SpatRaster. Returns list(mat, names).
#'
#' @keywords internal
.coerce_xcov_to_matrix <- function(x_covariates, cov_names) {
      if (inherits(x_covariates, "SpatRaster")) {
            if (!requireNamespace("terra", quietly = TRUE)) {
                  stop("Package 'terra' is required for SpatRaster x_covariates.",
                       call. = FALSE)
            }
            nms <- terra::names(x_covariates)
            mat <- as.matrix(terra::as.data.frame(x_covariates,
                                                  xy = FALSE, na.rm = FALSE))
            return(list(mat = mat, names = nms))
      }
      if (is.data.frame(x_covariates)) {
            return(list(mat = as.matrix(x_covariates),
                        names = names(x_covariates)))
      }
      if (is.matrix(x_covariates)) {
            return(list(mat = x_covariates, names = colnames(x_covariates)))
      }
      if (is.numeric(x_covariates) && is.null(dim(x_covariates))) {
            # Single-covariate shorthand: valid only when there is one covariate
            if (length(cov_names) != 1L) {
                  stop("`x_covariates` given as a plain numeric vector is only ",
                       "allowed when `covariates` has a single variable.",
                       call. = FALSE)
            }
            mat <- matrix(x_covariates, ncol = 1L)
            colnames(mat) <- cov_names
            return(list(mat = mat, names = cov_names))
      }
      stop("`x_covariates` must be a numeric vector, matrix, data.frame, or ",
           "SpatRaster.", call. = FALSE)
}


#' Get row count of an input that may be a vector, matrix, data.frame, or SpatRaster
#'
#' @keywords internal
nrow_of <- function(obj) {
      if (inherits(obj, "SpatRaster")) {
            if (!requireNamespace("terra", quietly = TRUE)) {
                  stop("Package 'terra' is required for SpatRaster inputs.",
                       call. = FALSE)
            }
            return(terra::ncell(obj))
      }
      if (is.null(dim(obj))) return(length(obj))
      nrow(obj)
}
