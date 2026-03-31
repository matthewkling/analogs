#' Local weighted regression across analog neighborhoods
#'
#' Fits a weighted local regression of `y` on `covariates` within
#' each focal location's analog neighborhood. Analog neighborhoods are
#' defined by climatic similarity, geographic proximity, or both, while
#' covariates capture additional predictors that influence outcomes within
#' each neighborhood but are not captured by to the search dimensions.
#' This generalizes the weighted mean — which averages over all within-neighborhood
#' variation — by resolving variation driven by these auxiliary predictors.
#' Supports ordinary and ridge-penalized weighted least squares. With purely geographic
#' neighborhoods, this is equivalent to geographically weighted regression
#' (GWR); with climate-based neighborhoods, it extends the analog impact
#' model (AIM) framework to incorporate local covariate effects.
#'
#' This function is a wrapper that calls [analog_search()] with
#' `"regression"` included in `stat`.
#'
#' @param x Focal locations for which regressions will be fit. Should be a
#'   matrix/data.frame with columns x, y, and climate variables, or a
#'   SpatRaster with climate variable layers.
#' @param pool The reference dataset to search for analogs. Either:
#'   \itemize{
#'     \item Matrix/data.frame with columns x, y, and climate variables,
#'       or SpatRaster with climate variable layers, OR
#'     \item An `analog_index` object created by
#'       [build_analog_index()] (for repeated queries).
#'   }
#' @param y Response variable(s) to model via local regression.
#'   Can be a numeric vector, matrix, data.frame, or SpatRaster.
#'   Must have exactly the same number of rows/cells as `pool`.
#'   A separate regression is fit for each variable.
#' @param covariates Predictor variables for local regression. Can be a
#'   numeric vector (single covariate), matrix, data.frame, or SpatRaster.
#'   Must have exactly the same number of rows/cells as `pool`. Column/layer
#'   names carry through to output. These variables are NOT used for the
#'   analog search itself -- only for regression within each neighborhood.
#' @param lambda Ridge penalty parameter (default: 0, giving ordinary
#'   weighted least squares). Higher values shrink covariate coefficients
#'   toward zero, with the intercept approaching the weighted mean as
#'   `lambda -> Inf`. Useful when some neighborhoods have few analogs
#'   relative to the number of covariates.
#' @param stat Statistic(s) to compute. `"regression"` is always included.
#'   Additional stats like `"count"`, `"ess"`, and `"weighted_mean"` can
#'   be requested alongside regression coefficients. Default includes
#'   `"count"` and `"ess"` for diagnostics.
#' @inheritParams analog_search
#'
#' @return A data.frame (or SpatRaster if `x` is a SpatRaster) with one
#'   row per focal location containing:
#'   \itemize{
#'     \item `index`, `x`, `y`: Focal location identifiers
#'     \item Columns for any additional stats requested (e.g., `count`, `ess`)
#'     \item `intercept`: Regression intercept (predicted value when all
#'       covariates equal zero)
#'     \item One column per covariate: regression slope coefficients
#'     \item With multiple `y` variables: columns are named
#'       `{coeff}_{varname}` (e.g., `intercept_biomass`, `slope_biomass`)
#'   }
#'
#' @details
#' ## Method
#'
#' For each focal location, the function:
#' 1. Selects analog pool locations based on `select`, `max_clim`, `max_geog`, and `k`
#' 2. Computes distance-based weights for each analog (via `weight` and `theta`)
#' 3. Fits a weighted least squares regression of `y` on `covariates`
#'    using these weights, with optional ridge penalty `lambda`
#' 4. Returns the regression coefficients (intercept + slopes)
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
#'   \item **Climate-nearest regression** (`select = "knn_clim"`, with
#'     `max_geog` and `k`): Fixed-size neighborhoods of the k most similar
#'     climates within geographic range.
#' }
#'
#' ## Prediction
#'
#' The function returns coefficients only. Prediction at new covariate values
#' (e.g., a fine-resolution topography grid) can be done via regression algebra,
#' e.g.:
#'
#' \preformatted{
#' prediction <- intercept + beta_1 * covariate_1 + beta_2 * covariate_2
#' }
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
#'   max_clim = NULL,
#'   weight = "gaussian_geog",
#'   theta = 20
#' )
#' }
#'
#' @seealso [analog_search()] for the underlying flexible analog search function;
#'   [analog_impact()] for the standard AIM workflow.
#'
#' @export
analog_regression <- function(
            x,
            pool,
            y,
            covariates,
            max_geog = NULL,
            max_clim = NULL,
            select = "all",
            k = NULL,
            weight = c("gaussian_clim", "inverse_clim",
                       "gaussian_geog", "inverse_geog",
                       "gaussian_joint", "inverse_joint",
                       "uniform"),
            theta = NULL,
            lambda = 0,
            stat = c("count", "ess", "regression"),
            x_cov = NULL,
            coord_type = "auto",
            index_res = "auto",
            n_threads = NULL,
            progress = FALSE
) {
      weight <- match.arg(weight)

      # Ensure "regression" is always in stat
      if (!"regression" %in% stat) {
            stat <- unique(c(stat, "regression"))
      }

      analog_search(
            x           = x,
            pool        = pool,
            select      = select,
            stat        = stat,
            y           = y,
            covariates  = covariates,
            max_clim    = max_clim,
            max_geog    = max_geog,
            k           = k,
            weight      = weight,
            theta       = theta,
            lambda      = lambda,
            x_cov       = x_cov,
            coord_type  = coord_type,
            index_res   = index_res,
            n_threads   = n_threads,
            progress    = progress
      )
}
