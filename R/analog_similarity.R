#' Analog similarity: best climate analogs within a geographic envelope
#'
#' Finds, for each focal location, the climate–nearest neighbor(s) in a
#' reference dataset that satisfy a specified geographic distance threshold.
#' Among other uses, this operation is often the first step in a traditional
#' analog impact modeling (AIM) analysis.
#'
#' This function is a wrapper that calls [analog_search()] using `select = "knn_clim"`.
#'
#' @inheritParams analog_search
#' @inherit analog_search return
#'
#' @details
#' For each focal location, \code{analog_similarity()}:
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
#' im <- analog_similarity(
#'   x = clim$clim1,
#'   pool = clim$clim2,
#'   max_geog = 100,
#'   k = 20
#' )
#'
#' # With pre-built index (for repeated queries)
#' index <- build_analog_index(clim$clim2)
#' i1 <- analog_similarity(x = sites1, pool = index, max_geog = 100, k = 20)
#' i2 <- analog_similarity(x = sites2, pool = index, max_geog = 50, k = 10)
#' }
#'
#' @export
analog_similarity <- function(
            x,
            pool,
            x_cov = NULL,
            values = NULL,
            coord_type = "auto",

            max_geog,
            max_clim = NULL,
            k = 20,

            index_res = "auto",
            n_threads = NULL,
            downsample = 1.0,
            seed = NULL,
            progress = FALSE
) {
      analog_search(
            x           = x,
            pool        = pool,
            select      = "knn_clim",
            stat        = "none",  # Returns pairs
            max_clim    = max_clim,
            max_geog    = max_geog,
            x_cov       = x_cov,
            values      = values,
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
