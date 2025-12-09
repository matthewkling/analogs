#' Analog intensity: weighted sum of analogs within climate/geographic limits
#'
#' Computes, for each focal location, the weighted sum of all reference locations
#' that satisfy the supplied climate and geographic constraints. The weights are
#' controlled by the \code{weight} and \code{theta} arguments and are applied
#' after filtering.
#'
#' @inheritParams analog_search
#' @inherit analog_search return
#'
#' @examples
#' \dontrun{
#' # One-shot query with inverse weighting
#' intens <- analog_intensity(
#'   x = sites,
#'   pool = climate_data,
#'   max_clim = 0.5,
#'   max_geog = 100,
#'   weight = "inverse_clim"
#' )
#'
#' # Gaussian weighting by climate distance
#' intens_gauss <- analog_intensity(
#'   x = sites,
#'   pool = climate_data,
#'   max_clim = 0.5,
#'   max_geog = 100,
#'   weight = "gaussian_clim",
#'   theta = 0.2  # bandwidth parameter
#' )
#'
#' # Joint Gaussian weighting (both climate and geography)
#' intens_joint <- analog_intensity(
#'   x = sites,
#'   pool = climate_data,
#'   max_clim = 0.5,
#'   max_geog = 100,
#'   weight = "gaussian_joint",
#'   theta = c(0.2, 50)  # c(clim_bandwidth, geog_bandwidth)
#' )
#'
#' # With pre-built index (for repeated queries)
#' index <- build_analog_index(climate_data)
#' i1 <- analog_intensity(x = sites1, pool = index, max_clim = 0.5,
#'                        weight = "inverse_clim")
#' i2 <- analog_intensity(x = sites2, pool = index, max_geog = 100,
#'                        weight = "inverse_geog")
#' }
#'
#' @export
analog_intensity <- function(
            x,
            pool,
            x_cov      = NULL,
            values = NULL,
            coord_type = "auto",

            max_clim   = NULL,
            max_geog   = NULL,

            weight     = c("uniform", "inverse_clim", "inverse_geog",
                           "gaussian_clim", "gaussian_geog",
                           "gaussian_joint", "inverse_joint"),
            theta      = NULL,

            index_res = "auto",
            n_threads = NULL
) {
      weight <- match.arg(weight)

      analog_search(
            x           = x,
            pool        = pool,
            select      = "all",
            stat        = "sum_weights",
            max_clim    = max_clim,
            max_geog    = max_geog,
            k           = NULL,
            weight      = weight,
            theta       = theta,
            x_cov       = x_cov,
            values = values,
            coord_type  = coord_type,
            index_res   = index_res,
            n_threads   = n_threads
      )
}
