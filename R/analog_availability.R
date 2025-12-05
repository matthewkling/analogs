#' Analog availability: count of all analogs within climate/geographic limits
#'
#' Computes, for each focal location, how many reference locations satisfy
#' the supplied climate and geographic constraints. This is useful for
#' mapping analog "availability" or environmental similarity density.
#'
#' @inheritParams analog_search
#'
#' @return A data.frame with columns:
#'   - focal_index
#'   - focal_x, focal_y
#'   - value (the count of analogs)
#'
#' @examples
#' \dontrun{
#' # One-shot query
#' avail <- analog_availability(
#'   x = sites,
#'   pool = climate_data,
#'   max_clim = 0.5,
#'   max_geog = 100
#' )
#'
#' # With pre-built index (for repeated queries)
#' index <- build_analog_index(climate_data)
#' a1 <- analog_availability(x = sites1, pool = index, max_clim = 0.5, max_geog = 100)
#' a2 <- analog_availability(x = sites2, pool = index, max_clim = 0.3, max_geog = 50)
#' }
#'
#' @export
analog_availability <- function(
            x,
            pool,
            max_clim = NULL,
            max_geog = NULL,
            x_cov = NULL,
            coord_type = "auto",
            index_res = "auto",
            n_threads = NULL
) {
      analog_search(
            x           = x,
            pool        = pool,
            select      = "all",
            aggregate   = "count",
            max_clim    = max_clim,
            max_geog    = max_geog,
            k           = NULL,
            weight      = NULL,
            theta       = NULL,
            x_cov       = x_cov,
            coord_type  = coord_type,
            report_dist = FALSE,  # no pair distances in aggregate mode
            index_res   = index_res,
            n_threads   = n_threads
      )
}
