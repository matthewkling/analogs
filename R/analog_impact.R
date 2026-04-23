#' Climate impact assessment via analog impact model
#'
#' Assesses potential climate change impacts using the analog impact model (AIM)
#' methodology. For each focal location's future climate, identifies locations
#' with similar baseline climates within a specified geographic range, then
#' aggregates their ecological characteristics weighted by climate similarity.
#' This quantifies what ecosystem conditions are likely accessible via dispersal
#' as climate changes.
#'
#' @param x Focal locations (generally with future climate conditions). Should be a
#'   matrix/data.frame with columns x, y, and climate variables, or a
#'   SpatRaster with climate variable layers.
#' @param pool The reference dataset (generally representing baseline climate conditions).
#' Either:
#'   \itemize{
#'     \item Matrix/data.frame with columns x, y, and climate variables,
#'       or SpatRaster with climate variable layers, OR
#'     \item An [analog_index()] object created by
#'       [build_analog_index()] (for repeated queries).
#'   }
#' @param y Ecological or environmental variable(s) for the same era as `pool`, to
#'   aggregate across climate analogs. Examples include occupancy of focal species, species
#'   richness, biomass, or any other ecological state variable. Can be a numeric vector
#'   (single variable), matrix or data.frame with numeric columns (multiple
#'   variables), or a SpatRaster with one or more numeric layers. Must have
#'   exactly the same number of reference locations as `pool`.
#' @param stat Statistic(s) to compute across analogs (default: c("count",
#'   "sum_weights", "weighted_mean")). See [analog_search()] for options.
#'   The default statistics provide a complete picture:
#'   \itemize{
#'     \item `"count"`: Analog availability (number of analogs)
#'     \item `"sum_weights"`: Analog intensity
#'     \item `"weighted_mean"`: Expected ecological state
#'   }
#' @param kernel Kernel decay function for weighting analogs during aggregation. Only
#'   weighting options that are based on *climate* are allowed:
#'   `"inverse_clim"` (default), `"gaussian_clim"`, `"inverse_joint"`,
#'   `"gaussian_joint"`. See [analog_search()] for details.
#' @inheritParams analog_search
#'
#' @return A data.frame with one row per focal location containing:
#'   \itemize{
#'     \item `index`: Row number from input `x` data
#'     \item `x, y`: Coordinates of focal location
#'     \item One column per requested statistic
#'     \item For value statistics with multiple variables: `{stat}_{varname}`
#'       (e.g., `weighted_mean_habitat_quality`)
#'     \item For `"regression"`: `coef_intercept` and `coef_{covariate}`
#'       coefficient columns
#'     \item When `se != "none"`: SE columns for SE-supporting stats
#'       (e.g., `se_weighted_mean`)
#'   }
#'
#' @details
#' ## The Analog Impact Model (AIM) Framework
#'
#' This function implements the "reverse analog" approach from the climate
#' change ecology literature. It addresses the question, "For a location's
#' future climate, what ecological conditions exist in current locations with
#' similar climates that are within dispersal range?"
#'
#' The methodology:
#' 1. For each focal location's future climate conditions
#' 2. Find all current locations with similar climates (within `max_clim`)
#' 3. Constrain to dispersal-reachable distance (within `max_geog`)
#' 4. Weight each analog by climate similarity (via `kernel` function)
#' 5. Aggregate ecosystem characteristics across these weighted analogs
#'
#' Unlike traditional AIM implementations that select k nearest climate neighbors,
#' this function uses all analogs within thresholds combined with climate-distance-based
#' kernel weighting. This approach eliminates arbitrary choice of k, provides smoother,
#' more continuous results, and lets the kernel function (via `theta`) naturally
#' control influence. (Note that the traditional version can be implemented via
#' `analog_search(select = "knn_clim", stat = "mean", ...))`.)
#'
#' ## Choosing Parameters
#'
#' \itemize{
#'   \item `max_geog`: Set based on species dispersal ability (e.g., 5-500 km)
#'   \item `max_clim`: Defines what counts as an "analog"
#'   \item `theta`: Controls kernel decay. The weight should decay to a small value
#'     at the `max_clim`/`max_geog` boundary. If `theta` is too large relative to
#'     thresholds, the hard cutoffs do most of the filtering and weighting becomes
#'     nearly uniform. For Gaussian kernels with three or fewer climate variables,
#'     a reasonable rule of thumb is to set `theta` to `max_* / 3`.
#' }
#'
#' ## Interpreting Results
#'
#' \itemize{
#'   \item `count`: How many analogs exist within max_clim and max_geog? Low counts
#'     indicate limited analog availability, while zero counts indicate climates
#'     that are novel within the geographic search radius.
#'   \item `sum_weights`: Total analog intensity. Low values indicate sparse
#'     or distant climate matches. This metric captures both the number and quality
#'     of analogs. Interpretation details vary based on the `kernel` parameter.
#'   \item `weighted_mean`: Expected ecosystem state if colonized by species
#'     from analog locations.
#' }
#'
#' @seealso [analog_search()] for the underlying flexible analog search function.
#'
#' @examples
#' \dontrun{
#' # Basic climate impact assessment
#' impact <- analog_impact(
#'   x = future_climate,
#'   pool = current_climate,
#'   y = current$habitat,
#'   max_geog = 100,    # 100 km dispersal range
#'   max_clim = 0.5     # Climate analog threshold
#' )
#'
#' # With uncertainty quantification on weighted_mean
#' impact_se <- analog_impact(
#'   x = future_climate,
#'   pool = current_climate,
#'   y = current$habitat,
#'   max_geog = 100,
#'   max_clim = 0.5,
#'   se = "ess"
#' )
#' }
#'
#' @export
analog_impact <- function(
            x,
            pool,
            y,
            covariates = NULL,
            max_geog = NULL,
            max_clim = 1.0,
            kernel = c("gaussian_clim", "inverse_clim",
                       "gaussian_joint", "inverse_joint"),
            theta = .25,
            stat = c("count", "sum_weights", "weighted_mean"),
            lambda = 0,
            se = c("none", "ess", "design"),
            x_cov = NULL,
            coord_type = "auto",
            index_res = "auto",
            n_threads = NULL,
            progress = FALSE
) {
      kernel <- match.arg(kernel)
      se <- match.arg(se)

      analog_search(
            x           = x,
            pool        = pool,
            select      = "all",
            stat        = stat,
            y           = y,
            covariates  = covariates,
            max_clim    = max_clim,
            max_geog    = max_geog,
            k           = NULL,
            kernel      = kernel,
            theta       = theta,
            lambda      = lambda,
            se          = se,
            x_cov       = x_cov,
            coord_type  = coord_type,
            index_res   = index_res,
            n_threads   = n_threads,
            progress    = progress
      )
}
