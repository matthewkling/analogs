#' Climate velocity: nearest geographic analogs within a climate envelope
#'
#' Computes, for each focal location, the geographic nearest neighbor(s) in a
#' reference dataset that satisfy a specified climate distance threshold. This
#' helper wraps \code{\link{analog_search}} using \code{select = "knn_geog"} and is
#' most commonly used for estimating climate velocity (the rate and direction
#' at which organisms would have to move to track constant climate conditions).
#'
#' @inheritParams analog_search
#'
#' @details
#' For each focal point, this function:
#' \enumerate{
#'   \item Identifies all reference points satisfying the climate (and optional
#'         geographic) threshold(s).
#'   \item Among those, selects the \code{k} nearest in \emph{geographic}
#'         distance.
#' }
#'
#' This is the classical operation needed for estimating \emph{climate
#' velocity}: the minimum relocation distance needed to maintain similar
#' climatic conditions under temporal change.
#'
#' @return A data set with a record for every focal site in `x` for each of
#' its `k` analogs from the `pool`, with the following variables:
#'   \itemize{
#'     \item `index`, `x`, `y`: reference variables identifying the focal site's
#'       row number (index) and geographic coordinates from input dataset `x`
#'     \item `analog_index`, `analog_x`, `analog_y`: results identifying an analog
#'       for each focal site, including its index in `pool` and its coordinates
#'     \item `clim_dist` Euclidean climate distance between the focal and analog sites
#'     \item `geog_dist` Geographic distance between the focal and analog sites,
#'       either in km or projected units.
#'   }
#'
#' If `x` is a SpatRaster AND `k = 1`, the result is a SpatRaster and the reference
#' variables are omitted; otherwise it is a data.frame.
#'
#' Diagnostic attributes (e.g., binning statistics) from the underlying spatial
#' index are preserved.
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
