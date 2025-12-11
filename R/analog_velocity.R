#' Climate velocity: geographically nearest climate analogs
#'
#' Computes, for each focal location, the geographic nearest neighbor(s) in a
#' reference dataset that satisfy a specified maximum climate distance threshold.
#' This implements analog-based climate velocity (Hamann et al. 2015;
#' Dobrowski and Parks 2016).
#'
#' This function is a wrapper that calls [analog_search()] using `select = "knn_geog"`.
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
#' @references
#' Hamann A, Roberts DR, Barber QE, Carroll C, Nielsen SE (2015). "Velocity of
#' climate change algorithms for guiding conservation and management."
#' \emph{Global Change Biology}, \strong{21}(2), 997-1004.
#' \doi{10.1111/gcb.12736}
#'
#' Dobrowski SZ, Parks SA (2016). "Climate change velocity underestimates
#' climate change exposure in mountainous regions." \emph{Nature Communications},
#' \strong{7}, 12349. \doi{10.1038/ncomms12349}
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
            n_threads = NULL,
            downsample = 1.0,
            seed = NULL,
            progress = FALSE
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
            n_threads   = n_threads,
            downsample  = downsample,
            seed        = seed,
            progress    = progress
      )
}
