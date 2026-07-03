#' Analog similarity: best climate analogs within a geographic envelope
#'
#' Finds, for each focal location, the climate–nearest neighbor(s) in a
#' reference dataset that satisfy a specified geographic distance threshold.
#' This function is a wrapper that calls [analog_search()] using `select = "knn_clim"`.
#'
#' @inheritParams analog_search
#'
#' @return A data.frame, or a SpatRaster when `x` is one and `k = 1`.
#'   Contains one row per focal-analog pair with `index`, `x`, `y`,
#'   `analog_index`, `analog_x`, `analog_y`, `clim_dist`, and
#'   `geog_dist`. See [analog_search()] for full column conventions
#'   and [metadata()] for attached metadata attributes.
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
#' Among other uses, this operation is often the first step in a traditional
#' analog impact modeling (AIM) analysis -- though see [analog_impact()] for a
#' more complete AIM implementation.
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
#' @seealso [analog_search()] for the underlying flexible analog search function;
#'   [tiled_analog_search()] for memory-safe searches on large raster datasets.
#'
#' @export
analog_similarity <- function(
            x,
            pool,
            x_cov = NULL,
            y = NULL,
            weight = NULL,
            coord_type = "auto",

            max_geog,
            max_clim = NULL,
            k = 20,

            clim_res_adj = "auto",
            geog_res_adj = "auto",
            cell_area_weight = "auto",
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
            y           = y,
            weight      = weight,
            k           = k,
            kernel      = NULL,
            theta       = NULL,
            coord_type  = coord_type,
            clim_res_adj= clim_res_adj,
            geog_res_adj = geog_res_adj,
            cell_area_weight = cell_area_weight,
            n_threads   = n_threads,
            downsample  = downsample,
            seed        = seed,
            progress    = progress
      )
}
