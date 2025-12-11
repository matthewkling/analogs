
#' Climate impact assessment via analog impact model
#'
#' Assesses potential climate change impacts using the analog impact model (AIM)
#' methodology. For each focal location's future climate, identifies locations
#' with similar baseline climates within a specified geographic range, then
#' aggregates their ecological characteristics weighted by climate similarity.
#' This quantifies what ecosystem conditions are likely accessible via dispersal
#' as climate changes.
#'
#' @param x Focal locations (generally with future climate conditions). Should be a
#'   matrix/data.frame with columns x, y, and climate variables, or a
#'   SpatRaster with climate variable layers.
#' @param pool The reference dataset (generally representing baseline climate conditions).
#' Either:
#'   \itemize{
#'     \item Matrix/data.frame with columns x, y, and climate variables,
#'       or SpatRaster with climate variable layers, OR
#'     \item An [analog_index()] object created by
#'       [build_analog_index()] (for repeated queries).
#'   }
#' @param values Ecological or environmental variable(s) for the same era as `pool`, to
#'   aggregate across climate analogs. Must have exactly `nrow(pool)` rows.
#'   Examples include occupancy of focal species, species richness, biomass,
#'   or any other ecological state variable.
#' @param stat Statistic(s) to compute across analogs (default: c("count",
#'   "sum_weights", "weighted_mean")). See [analog_search()] for options.
#'   The default statistics provide a complete picture:
#'   \itemize{
#'     \item `"count"`: Analog availability (number of analogs)
#'     \item `"sum_weights"`: Analog intensity
#'     \item `"weighted_mean"`: Expected ecological state
#'   }
#' @param weight Function for weighting analogs during aggregation. Only
#'   weight options that are based on *climate* are allowed:
#'   `"inverse_clim"` (default), `"gaussian_clim"`, `"inverse_joint"`,
#'   `"gaussian_joint"`. See [analog_search()] for details.
#' @inheritParams analog_search
#'
#' @return A data.frame with one row per focal location containing:
#'   \itemize{
#'     \item `index`: Row number from input `x` data
#'     \item `x, y`: Coordinates of focal location
#'     \item One column per requested statistic
#'     \item For value statistics with multiple variables: `{stat}_{varname}`
#'       (e.g., `weighted_mean_habitat_quality`)
#'   }
#'
#' @details
#' ## The Analog Impact Model (AIM) Framework
#'
#' This function implements the "reverse analog" approach from the climate
#' change ecology literature. It addresses the question, "For a location's
#' future climate, what ecological conditions exist in current locations with
#' similar climates that are within dispersal range?"
#'
#' The methodology:
#' 1. For each focal location's future climate conditions
#' 2. Find all current locations with similar climates (within `max_clim`)
#' 3. Constrain to dispersal-reachable distance (within `max_geog`)
#' 4. Weight each analog by climate similarity (via `weight` function)
#' 5. Aggregate ecosystem characteristics across these weighted analogs
#'
#' Unlike traditional AIM implementations that select k nearest climate neighbors,
#' this function uses all analogs within thresholds combined with climate-distance-based
#' weighting. This approach eliminates arbitrary choice of k, provides smoother,
#' more continuous results, and lets the weight function (via `theta`) naturally
#' control influence. (Note that the traditional version can be implemented via
#' `analog_search(select = "knn_clim", stat = "mean", ...))`.)
#'
#' ## Choosing Parameters
#'
#' \itemize{
#'   \item `max_geog`: Set based on species dispersal ability (e.g., 5-500 km)
#'   \item `max_clim`: Defines what counts as an "analog"
#'   \item `theta`: Controls weight decay. The weight should decay to a small value
#'     at the `max_clim`/`max_geog` boundary. If `theta` is too large relative to
#'     thresholds, the hard cutoffs do most of the filtering and weighting becomes
#'     nearly uniform. For Gaussian weights with three or fewer climate variables,
#'     a reasonable rule of thumb is to set `theta` to `max_* / 3`.
#' }
#'
#' ## Interpreting Results
#'
#' \itemize{
#'   \item `count`: How many analogs exist within max_clim and max_geog? Low counts
#'     indicate limited analog availability, while zero counts indicate climates
#'     that are novel within the geographic search radius.
#'   \item `sum_weights`: Total analog intensity. Low values indicate sparse
#'     or distant climate matches. This metric captures both the number and quality
#'     of analogs. Interpretation details vary based on the `weight` parameter.
#'   \item `weighted_mean`: Expected ecosystem state if colonized by species
#'     from analog locations.
#' }
#'
#' @seealso [analog_search()] for the underlying flexible analog search function.
#'
#' @references
#' Parks SA, Holsinger LM, Miller C, Parisien MA (2018). "Analog-based fire
#' regime and vegetation shifts in mountainous regions of the western US."
#' \emph{Ecography}, \strong{41}(6), 910-921. \doi{10.1111/ecog.03378}
#'
#' Parks SA, Dobrowski SZ, Shaw JD, Miller C (2019). "Living on the edge:
#' Trailing edge forests at risk of fire-facilitated conversion to non-forest."
#' \emph{Ecosphere}, \strong{10}(3), e02651. \doi{10.1002/ecs2.2651}
#'
#' Holsinger L, Parks SA, Parisien MA, Miller C, Batllori E, Moritz MA (2019).
#' "Climate change likely to reshape vegetation in North America's largest
#' protected areas." \emph{Conservation Science and Practice}, \strong{1}(7),
#' e50. \doi{10.1111/csp2.50}
#'
#' Dobrowski SZ, Littlefield CE, Lyons DS, Hollenberg CH, Carroll C, Parks SA,
#' Abatzoglou JT, Hegewisch K, Gage J (2021). "Protected-area targets could be
#' undermined by climate change-driven shifts in ecoregions and biomes."
#' \emph{Communications Earth & Environment}, \strong{2}(1), 198.
#' \doi{10.1038/s43247-021-00270-z}
#'
#' @examples
#' \dontrun{
#' # Basic climate impact assessment
#' impact <- analog_impact(
#'   x = future_climate,
#'   pool = current_climate,
#'   values = current$habitat,
#'   max_geog = 100,    # 100 km dispersal range
#'   max_clim = 0.5     # Climate analog threshold
#' )
#'
#' # Multiple ecosystem variables
#' values_df <- data.frame(
#'   habitat_quality = current$habitat,
#'   species_richness = current$richness,
#'   forest_cover = current$forest
#' )
#' impact_multi <- analog_impact(
#'   x = future_climate,
#'   pool = current_climate,
#'   values = values_df,
#'   max_geog = 150,
#'   max_clim = 0.4,
#'   theta = 0.25
#' )
#'
#' # Custom statistics and weighting
#' impact_custom <- analog_impact(
#'   x = future_climate,
#'   pool = current_climate,
#'   values = current$biomass,
#'   stat = c("count", "weighted_mean", "weighted_sum"),
#'   max_clim = 0.6,
#'   max_geog = 200,
#'   weight = "gaussian_joint",    # Weight by both climate and geography
#'   theta = c(0.2, 50)           # Climate and geographic decay
#' )
#'
#' # With pre-built index for multiple scenarios
#' current_index <- build_analog_index(current_climate)
#'
#' impact_current <- analog_impact(current_climate, current_index,
#'                                  values = current$quality,
#'                                  max_geog = 100)
#' impact_ssp126 <- analog_impact(future_ssp126, current_index,
#'                                  values = current$quality,
#'                                  max_geog = 100)
#' impact_ssp585 <- analog_impact(future_ssp585, current_index,
#'                                  values = current$quality,
#'                                  max_geog = 100)
#' }
#'
#' @export
analog_impact <- function(
            x,
            pool,
            values,
            max_geog = NULL,
            max_clim = 1.0,
            weight = c("gaussian_clim", "inverse_clim",
                       "gaussian_joint", "inverse_joint"),
            theta = .25,
            stat = c("count", "sum_weights", "weighted_mean"),
            x_cov = NULL,
            coord_type = "auto",
            index_res = "auto",
            n_threads = NULL,
            downsample = 1.0,
            seed = NULL,
            progress = FALSE
) {
      weight <- match.arg(weight)

      analog_search(
            x           = x,
            pool        = pool,
            select      = "all",
            stat        = stat,
            values      = values,
            max_clim    = max_clim,
            max_geog    = max_geog,
            k           = NULL,
            weight      = weight,
            theta       = theta,
            x_cov       = x_cov,
            coord_type  = coord_type,
            index_res   = index_res,
            n_threads   = n_threads,
            downsample  = downsample,
            seed        = seed,
            progress    = progress
      )
}
