#' Climate velocity: geographically nearest climate analogs
#'
#' Computes, for each focal location, the geographic nearest neighbor(s) in a
#' reference dataset that satisfy a specified maximum climate distance threshold.
#' This helper wraps \code{\link{analog_search}} using \code{select = "knn_geog"}
#' and is used for estimating analog-based climate velocity.
#'
#' @inheritParams analog_search
#' @inherit analog_search return
#'
#' @examples
#' \dontrun{
#' # One-shot query
#' v <- analog_velocity(
#'   x = clim$clim1,
#'   pool = clim$clim2,
#'   max_clim = 0.5,
#'   k = 1
#' )
#'
#' # With pre-built index (for repeated queries)
#' index <- build_analog_index(clim$clim2)
#' v1 <- analog_velocity(x = sites1, pool = index, max_clim = 0.5, k = 1)
#' v2 <- analog_velocity(x = sites2, pool = index, max_clim = 0.3, k = 1)
#'
#' # With focal-specific covariance matrices
#' v_mahal <- analog_velocity(
#'   x = clim$clim1,
#'   pool = clim$clim2,
#'   x_cov = baseline_covariances,
#'   max_clim = 2,  # In Mahalanobis distance units
#'   k = 1
#' )
#' }
#'
#' @export
analog_velocity <- function(
            x,
            pool,
            x_cov = NULL,
            values = NULL,
            coord_type = "auto",

            max_clim,
            max_geog = NULL,
            k = 1,

            index_res = "auto",
            n_threads = NULL
) {
      analog_search(
            x           = x,
            pool        = pool,
            select      = "knn_geog",
            stat        = "none",  # Returns pairs
            max_clim    = max_clim,
            max_geog    = max_geog,
            x_cov       = x_cov,
            values = values,
            k           = k,
            weight      = NULL,
            theta       = NULL,
            coord_type  = coord_type,
            index_res   = index_res,
            n_threads   = n_threads
      )
}
