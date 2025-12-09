#' Climate impact: nearest climate analogs within a geographic envelope
#'
#' Computes, for each focal location, the climate–nearest neighbor(s) in a
#' reference dataset that satisfy a specified geographic distance threshold.
#' This helper wraps \code{\link{analog_search}} using \code{select = "knn_clim"}.
#'
#' It is useful for estimating the potential ecological impact of local climate
#' change: e.g., how climate conditions at a site compare to those available
#' within a species' dispersal range.
#'
#' @inheritParams analog_search
#' @inherit analog_search return
#'
#' @details
#' For each focal location, \code{analog_impact()}:
#' \enumerate{
#'   \item Identifies all reference points within \code{max_geog} km (and
#'         optional climate filter).
#'   \item Selects the \code{k} closest in \emph{climate} distance.
#' }
#'
#' This is the natural "inverse" of \code{\link{analog_velocity}}: instead of finding
#' where the focal climate moves geographically, it finds the closest climatically
#' similar conditions that are geographically reachable.
#'
#' @examples
#' \dontrun{
#' # One-shot query
#' im <- analog_impact(
#'   x = clim$clim1,
#'   pool = clim$clim2,
#'   max_geog = 100,
#'   k = 20
#' )
#'
#' # With pre-built index (for repeated queries)
#' index <- build_analog_index(clim$clim2)
#' i1 <- analog_impact(x = sites1, pool = index, max_geog = 100, k = 20)
#' i2 <- analog_impact(x = sites2, pool = index, max_geog = 50, k = 10)
#' }
#'
#' @export
analog_impact <- function(
            x,
            pool,
            x_cov = NULL,
            values = NULL,
            coord_type = "auto",

            max_geog,
            max_clim = NULL,
            k = 20,

            index_res = "auto",
            n_threads = NULL
) {
      analog_search(
            x           = x,
            pool        = pool,
            select      = "knn_clim",
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
