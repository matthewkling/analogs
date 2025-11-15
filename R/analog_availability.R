#' Analog availability: count of all analogs within climate/geographic limits
#'
#' Computes, for each focal location, how many reference locations satisfy
#' the supplied climate and geographic constraints. This is useful for
#' mapping analog "availability" or environmental similarity density.
#'
#' @param focal Data.frame/matrix/SpatRaster of focal points.
#' @param ref   Data.frame/matrix/SpatRaster of reference points.
#' @param max_clim Maximum climate distance (scalar or vector), or NULL.
#' @param max_geog Maximum geographic distance in km, or NULL.
#' @param coord_type "auto", "lonlat", or "projected".
#' @param n_threads Number of parallel compute threads to use.
#'
#' @return A data.frame with columns:
#'   - focal_index
#'   - focal_x, focal_y
#'   - value (the count of analogs)
#'
#' @export
analog_availability <- function(
            focal,
            ref,
            max_clim = NULL,
            max_geog = NULL,
            coord_type = "auto",
            n_threads = NULL
) {
      find_analogs(
            focal      = focal,
            ref        = ref,
            mode       = "count",
            max_clim   = max_clim,
            max_geog   = max_geog,
            k          = NULL,   # required to be NULL
            weight     = NULL,   # required to be NULL
            theta      = NULL,   # required to be NULL
            coord_type = coord_type,
            report_dist = FALSE,  # no pair distances in aggregate mode
            n_threads = n_threads
      )
}
