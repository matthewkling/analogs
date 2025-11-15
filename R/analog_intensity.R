#' Analog intensity: weighted sum of analogs within climate/geographic limits
#'
#' Computes, for each focal location, the weighted sum of all reference locations
#' that satisfy the supplied climate and geographic constraints. The weights are
#' controlled by the \code{weight} and \code{theta} arguments and are applied
#' after filtering.
#'
#' @param focal Data.frame/matrix/SpatRaster of focal points.
#' @param ref   Data.frame/matrix/SpatRaster of reference points.
#' @param max_clim Maximum climate distance (scalar or vector), or NULL for no
#'   climate constraint.
#' @param max_geog Maximum geographic distance in km, or NULL for no geographic
#'   constraint.
#' @param weight Weighting kernel, one of:
#'   \itemize{
#'     \item \code{"uniform"}: weight = 1 for all matches.
#'     \item \code{"inverse_clim"}: weight = 1 / (climate_distance + theta),
#'       with \code{theta} as epsilon (if NULL, a small default is used).
#'     \item \code{"inverse_geog"}: weight = 1 / (geographic_distance + theta),
#'       with \code{theta} as epsilon (if NULL, a small default is used).
#'   }
#' @param theta Optional numeric hyperparameter for the weighting kernel
#'   (epsilon term for inverse kernels). See \code{weight} description.
#' @param coord_type "auto", "lonlat", or "projected".
#' @param n_threads Number of parallel compute threads to use.
#'
#' @return A data.frame with one row per focal location:
#'   \itemize{
#'     \item \code{focal_index}
#'     \item \code{focal_x}, \code{focal_y}
#'     \item \code{value}: weighted sum over all analogs
#'   }
#'
#' @export
analog_intensity <- function(
            focal,
            ref,
            max_clim   = NULL,
            max_geog   = NULL,
            weight     = c("uniform", "inverse_clim", "inverse_geog"),
            theta      = NULL,
            coord_type = "auto",
            n_threads = NULL
) {
      weight <- match.arg(weight)

      find_analogs(
            focal      = focal,
            ref        = ref,
            mode       = "sum",
            max_clim   = max_clim,
            max_geog   = max_geog,
            k          = NULL,        # required to be NULL for sum/mean
            weight     = weight,
            theta      = theta,
            coord_type = coord_type,
            report_dist = FALSE,       # no pairwise distances needed for sums
            n_threads = n_threads
      )
}
