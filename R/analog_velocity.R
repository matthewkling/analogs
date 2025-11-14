#' Climate velocity: nearest geographic analogs within a climate envelope
#'
#' Computes, for each focal location, the geographic nearest neighbor(s) in a
#' reference dataset that satisfy a specified climate distance threshold. This
#' helper wraps \code{\link{find_analogs}} using \code{mode = "knn_geog"} and is
#' most commonly used for estimating climate velocity (the rate and direction
#' at which organisms would have to move to track constant climate conditions).
#'
#' @param focal Data.frame, matrix, or SpatRaster containing focal locations.
#'   Must include columns \code{x}, \code{y}, and one or more climate variables.
#'
#' @param ref Data.frame, matrix, or SpatRaster containing reference locations.
#'   Must be structured like \code{focal}.
#'
#' @param max_clim Maximum allowable climate distance between a focal and an
#'   analog. May be:
#'   \itemize{
#'     \item a scalar: an Euclidean radius in climate space, or
#'     \item a vector: per-variable absolute differences.
#'   }
#'   Only reference locations within this climate envelope are considered.
#'
#' @param k Number of nearest geographic analogs to return for each focal
#'   (default = 1).
#'
#' @param max_geog Optional additional geographic constraint in kilometers.
#'   If provided, analogs must be within \code{max_geog} km of the focal.
#'
#' @param coord_type One of \code{"auto"}, \code{"lonlat"}, or
#'   \code{"projected"}, indicating how geographic distance should be computed.
#'   \code{"auto"} attempts to infer this from coordinate ranges.
#'
#' @param report_dist Logical; if TRUE (default), include climate and geographic
#'   distance columns in the output.
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
#' v <- velocity(
#'   focal = clim$clim1,
#'   ref   = clim$clim2,
#'   max_clim = 0.5,
#'   k = 1
#' )
#' }
#'
#' @export
analog_velocity <- function(
            focal,
            ref,
            max_clim,
            k = 1,
            max_geog = NULL,
            coord_type = "auto",
            report_dist = TRUE
) {
      find_analogs(
            focal      = focal,
            ref        = ref,
            mode       = "knn_geog",
            max_clim   = max_clim,
            max_geog   = max_geog,
            k          = k,
            weight     = NULL,
            theta      = NULL,
            coord_type = coord_type,
            report_dist = report_dist
      )
}
