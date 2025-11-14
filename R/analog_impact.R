#' Climate impact: nearest climate analogs within a geographic envelope
#'
#' Computes, for each focal location, the climate–nearest neighbor(s) in a
#' reference dataset that satisfy a specified geographic distance threshold.
#' This helper wraps \code{\link{find_analogs}} using \code{mode = "knn_clim"}.
#'
#' It is useful for estimating the potential ecological impact of local climate
#' change: e.g., how climate conditions at a site compare to those available
#' within a species' dispersal range.
#'
#' @param focal Data.frame, matrix, or SpatRaster containing focal locations
#'   with columns \code{x}, \code{y}, and climate variables.
#'
#' @param ref Data.frame, matrix, or SpatRaster containing reference locations.
#'
#' @param max_geog Maximum allowable geographic distance in kilometers. Only
#'   reference locations within \code{max_geog} km are considered.
#'
#' @param k Number of nearest climate analogs to return for each focal
#'   (default = 20).
#'
#' @param max_clim Optional additional climate constraint (scalar or vector).
#'   Useful for restricting comparisons to very similar climate neighborhoods.
#'
#' @param coord_type One of \code{"auto"}, \code{"lonlat"}, or
#'   \code{"projected"}. \code{"auto"} attempts to detect the appropriate mode.
#'
#' @param report_dist Logical; if TRUE (default), include climate and geographic
#'   distance columns in the output.
#'
#' @details
#' For each focal location, \code{impact()}:
#' \enumerate{
#'   \item Identifies all reference points within \code{max_geog} km (and
#'         optional climate filter).
#'   \item Selects the \code{k} closest in \emph{climate} distance.
#' }
#'
#' This is the natural “inverse” of \code{\link{velocity}}: instead of finding
#' where the focal climate moves geographically, it finds the closest climatically
#' similar conditions that are geographically reachable.
#'
#' @return A data.frame with one row per focal–analog pair, including:
#'   \itemize{
#'     \item \code{focal_index}, \code{analog_index}
#'     \item \code{focal_x}, \code{focal_y}, \code{analog_x}, \code{analog_y}
#'     \item \code{clim_dist}, \code{geog_dist} (if \code{report_dist = TRUE})
#'   }
#' Diagnostic attributes from the underlying spatial index are preserved.
#'
#' @examples
#' \dontrun{
#' im <- impact(
#'   focal = clim$clim1,
#'   ref   = clim$clim2,
#'   max_geog = 100,
#'   k = 20
#' )
#' }
#'
#' @export
analog_impact <- function(
            focal,
            ref,
            max_geog,
            k = 20,
            max_clim = NULL,
            coord_type = "auto",
            report_dist = TRUE
) {
      find_analogs(
            focal      = focal,
            ref        = ref,
            mode       = "knn_clim",
            max_clim   = max_clim,
            max_geog   = max_geog,
            k          = k,
            weight     = NULL,
            theta      = NULL,
            coord_type = coord_type,
            report_dist = report_dist
      )
}
