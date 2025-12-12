#' Analog availability: count of all analogs within climate/geographic limits
#'
#' Computes, for each focal location, how many reference locations satisfy
#' the supplied climate and geographic constraints. This is useful for
#' mapping analog "availability" or environmental similarity density.
#'
#' This function is a wrapper that calls [analog_search()] using
#' `select = "all"` and `stat = "count"`.
#'
#' @inheritParams analog_search
#' @inherit analog_search return
#'
#' @references
#' Stralberg D, Carroll C, Pedlar JH, Wilsey CB, McKenney DW, Nielsen SE (2018).
#' "Macrorefugia for North American trees and songbirds: Climatic limiting
#' factors and multi-scale topographic influences." \emph{Global Ecology and
#' Biogeography}, \strong{27}(6), 690-703. \doi{10.1111/geb.12731}
#'
#' Carroll C, Lawler JJ, Roberts DR, Hamann A (2015). "Biotic and climatic
#' velocity identify contrasting areas of vulnerability to climate change."
#' \emph{PLOS ONE}, \strong{10}(10), e0140486.
#' \doi{10.1371/journal.pone.0140486}
#'
#' Parks SA, Holsinger LM, Abatzoglou JT, Littlefield CE, Zeller KA (2023).
#' "Protected areas not likely to serve as steppingstones for species
#' undergoing climate-induced range shifts." \emph{Global Change Biology},
#' \strong{29}(10), 2681-2696. \doi{10.1111/gcb.16629}
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
#' @seealso [analog_search()] for the underlying flexible analog search function;
#'   [tiled_analog_search()] for memory-safe searches on large raster datasets.
#'
#' @export
analog_availability <- function(
            x,
            pool,
            x_cov = NULL,
            coord_type = "auto",

            max_clim = NULL,
            max_geog = NULL,

            index_res = "auto",
            n_threads = NULL,
            downsample = 1.0,
            seed = NULL,
            progress = FALSE
) {
      analog_search(
            x           = x,
            pool        = pool,
            select      = "all",
            stat        = "count",
            max_clim    = max_clim,
            max_geog    = max_geog,
            k           = NULL,
            weight      = NULL,
            theta       = NULL,
            x_cov       = x_cov,
            values      = NULL,
            coord_type  = coord_type,
            index_res   = index_res,
            n_threads   = n_threads,
            downsample  = downsample,
            seed        = seed,
            progress    = progress
      )
}
