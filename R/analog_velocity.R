#' Climate velocity: nearest geographic analogs within a climate envelope
#'
#' Computes, for each focal location, the geographic nearest neighbor(s) in a
#' reference dataset that satisfy a specified climate distance threshold. This
#' helper wraps \code{\link{find_analogs}} using \code{mode = "knn_geog"} and is
#' most commonly used for estimating climate velocity (the rate and direction
#' at which organisms would have to move to track constant climate conditions).
#'
#' @inheritParams find_analogs
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
#' @return A data.frame with one row per focal–analog pair, including:
#'   \itemize{
#'     \item \code{focal_index}, \code{analog_index}
#'     \item \code{focal_x}, \code{focal_y}, \code{analog_x}, \code{analog_y}
#'     \item \code{clim_dist}, \code{geog_dist} (if \code{report_dist = TRUE})
#'   }
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
#' }
#'
#' @export
analog_velocity <- function(
            x,
            pool,
            max_clim,
            k = 1,
            max_geog = NULL,
            coord_type = "auto",
            report_dist = TRUE,
            index_res = "auto",
            n_threads = NULL
) {
      find_analogs(
            x          = x,
            pool       = pool,
            mode       = "knn_geog",
            max_clim   = max_clim,
            max_geog   = max_geog,
            k          = k,
            weight     = NULL,
            theta      = NULL,
            coord_type = coord_type,
            report_dist = report_dist,
            index_res  = index_res,
            n_threads  = n_threads
      )
}
