#' Analog density: kernel-weighted analog count within climate/geographic limits
#'
#' Computes, for each focal location, the kernel-weighted sum of all reference
#' locations that satisfy the supplied climate and geographic constraints.
#'
#' This function is a wrapper that calls [analog_search()] using
#' `select = "all"` and `stat = "sum_weights"`. It also supports `"count"`,
#' `"mean_weights"`, and `"ess"` as additional stats.
#'
#' By default (`normalize = "auto"`), the returned `sum_weights` column is
#' divided by a global scalar `D_max` whenever `pool` meets the prerequisites
#' (raster `pool` with cell-area weighting active). This gives values on roughly
#' `[0, 1]`, interpretable as the fraction of theoretical maximum analog density
#' achieved at each focal. See [analog_search()] for more details.
#'
#' @param stat Character vector of one or more density-style summary statistics.
#'   Allowed options: `"sum_weights"` (default), `"mean_weights"`, `"count"`,
#'   `"ess"`. See [analog_search()] for detailed stat descriptions.
#'
#' @param ... Additional arguments forwarded to [analog_search()] (e.g.
#'   `mean_cell_area` for tiled workflows).
#'
#' @inheritParams analog_search
#'
#' @return A data.frame, or a SpatRaster when `x` is one. Contains
#'   `index`, `x`, `y`, and `sum_weights`. When `normalize = TRUE`, the
#'   value of `D_max` is attached as a result attribute. See [metadata()]
#'   for attached metadata attributes.
#'
#' @references
#' Mahony CR, Cannon AJ, Wang T, Aitken SN (2017). "A closer look at novel
#' climates: New methods and insights at continental to landscape scales."
#' \emph{Global Change Biology}, \strong{23}(9), 3934-3955.
#' \doi{10.1111/gcb.13645}
#'
#' Abatzoglou JT, Dobrowski SZ, Parks SA (2020). "Multivariate climate
#' departures have outpaced univariate changes across global lands."
#' \emph{Scientific Reports}, \strong{10}(1), 3891.
#' \doi{10.1038/s41598-020-60270-5}
#'
#' Williams JW, Jackson ST, Kutzbach JE (2007). "Projected distributions of
#' novel and disappearing climates by 2100 AD." \emph{Proceedings of the
#' National Academy of Sciences}, \strong{104}(14), 5738-5742.
#' \doi{10.1073/pnas.0606292104}
#'
#' @examples
#' \dontrun{
#' # Normalized density (default): returns sum_weights on roughly [0, 1]
#' dens <- analog_density(
#'   x = sites,
#'   pool = climate_data,
#'   clim = kernel("gaussian", theta = 0.2, max = 0.5),
#'   geog = kernel(max = 100)
#' )
#' attr(dens, "D_max")  # the global denominator used
#'
#' # Raw (unnormalized) kernel-weighted sums
#' dens_raw <- analog_density(
#'   x = sites,
#'   pool = climate_data,
#'   clim = kernel("gaussian", theta = 0.2, max = 0.5),
#'   geog = kernel(max = 100),
#'   normalize = FALSE
#' )
#'
#' # Joint Gaussian weighting (both climate and geography)
#' intens_joint <- analog_density(
#'   x = sites,
#'   pool = climate_data,
#'   clim = kernel("gaussian", theta = 0.2, max = 0.5),
#'   geog = kernel("gaussian", theta = 50, max = 100)
#' )
#'
#' # With pre-built index (for repeated queries)
#' index <- build_analog_index(climate_data)
#' i1 <- analog_density(x = sites1, pool = index,
#'                      clim = kernel("inverse", max = 0.5),
#'                      geog = kernel(max = 100))
#' i2 <- analog_density(x = sites2, pool = index,
#'                      clim = kernel("gaussian", theta = 0.3),
#'                      geog = kernel(max = 100))
#' }
#'
#' @seealso [analog_search()] for the underlying flexible analog search function;
#'   [tiled_analog_search()] for memory-safe searches on large raster datasets.
#'
#' @export
analog_density <- function(
            x,
            pool,
            x_cov      = NULL,
            weight     = NULL,
            coord_type = "auto",

            stat        = c("sum_weights"),
            clim        = NULL,
            geog        = NULL,

            normalize  = "auto",

            clim_res_adj = "auto",
            geog_res_adj = "auto",
            cell_area_weight = "auto",
            n_threads = NULL,
            downsample = 1.0,
            seed = NULL,
            progress = FALSE,
            ...
) {
      # Validate
      allowed_density_stats <- c("count", "sum_weights", "mean_weights", "ess")
      bad <- setdiff(stat, allowed_density_stats)
      if (length(bad) > 0L) {
            stop(
                  "`stat` value(s) not supported by analog_density(): ",
                  paste(sprintf("'%s'", bad), collapse = ", "),
                  ". Allowed: ",
                  paste(sprintf("'%s'", allowed_density_stats), collapse = ", "),
                  ". For y-based stats (mean, weighted_mean, regression), ",
                  "use analog_impact() or analog_regression() instead.",
                  call. = FALSE
            )
      }

      analog_search(
            x           = x,
            pool        = pool,
            select      = "all",
            stat        = stat,
            k           = NULL,
            clim        = clim,
            geog        = geog,
            x_cov       = x_cov,
            y           = NULL, # not relevant for density
            weight      = weight,
            normalize   = normalize,
            coord_type  = coord_type,
            clim_res_adj= clim_res_adj,
            geog_res_adj = geog_res_adj,
            cell_area_weight = cell_area_weight,
            n_threads   = n_threads,
            downsample  = downsample,
            seed        = seed,
            progress    = progress,
            ...
      )
}
