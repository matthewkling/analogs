#' Find Climate Analogs
#'
#' Identifies locations in a reference dataset that are climatically similar to
#' focal locations, with optional constraints on climate distance and geographic
#' distance. This function uses a two-stage approach: first selecting analogs
#' based on specified criteria, then optionally aggregating the results.
#'
#' The function uses a spatial indexing structure (lattice-based) to quickly
#' search through large reference datasets. Climate similarity is measured
#' using Euclidean distance in climate space (ideally pre-whitened; see Details).
#' Geographic distance can be computed for lon/lat coordinates (great-circle
#' distance) or projected coordinates (planar distance).
#'
#' @param x Focal locations for which analogs will be found. Should be a
#'   matrix/data.frame with columns x, y, and climate variables, or a
#'   SpatRaster with climate variable layers.
#'
#' @param pool The reference dataset to search for analogs. Either:
#'   \itemize{
#'     \item Matrix/data.frame with columns x, y, and climate variables,
#'       or SpatRaster with climate variable layers, OR
#'     \item An \code{\link{analog_index}} object created by
#'       \code{\link{build_analog_index}} (for repeated queries).
#'   }
#'
#' @param select Character string specifying the analog selection strategy.
#'   One of:
#'   \itemize{
#'     \item \code{"all"} (default): Select all analogs that satisfy the
#'       \code{max_clim} and \code{max_geog} constraints.
#'     \item \code{"knn_clim"}: For each focal, select up to \code{k} analogs
#'       with smallest climate distance, subject to filters.
#'     \item \code{"knn_geog"}: For each focal, select up to \code{k} analogs
#'       with smallest geographic distance, subject to filters.
#'   }
#'
#' @param stat Statistic used to aggregate selected analogs. Either:
#'   \itemize{
#'     \item \code{NULL} or \code{"none"} (default): Return all selected
#'       analog pairs as a data.frame.
#'     \item \code{"count"}: For each focal, count the number of selected analogs.
#'     \item \code{"sum_weights"}: For each focal, sum the weights of selected
#'       analogs (see \code{weight} and \code{theta}).
#'     \item \code{"mean_weights"}: For each focal, mean of weights of selected
#'       analogs.
#'   }
#'
#' @param max_clim Maximum climate distance constraint (default: NULL = no
#'   climate constraint). Can be either:
#'   \itemize{
#'     \item A scalar: Euclidean radius in climate space (e.g., 0.5)
#'     \item A vector: Per-variable absolute differences (length must equal
#'       number of climate variables)
#'   }
#'   Only reference locations within this climate distance are considered.
#'   When \code{x_cov} is provided, scalar thresholds are interpreted in
#'   Mahalanobis distance units.
#'
#' @param max_geog Maximum geographic distance constraint (default:
#'   NULL = no geographic constraint). When specified, only reference locations
#'   within this distance are considered. Radius units should be specified in
#'   kilometers if \code{coord_type = "lonlat"}, or in projected coordinate units
#'   if \code{coord_type = "projected"}.
#'
#' @param x_cov Optional focal-specific covariance matrices for Mahalanobis
#'   distance calculations. Should be a matrix or data.frame with one row per
#'   focal location and one column per unique covariance component. For n climate
#'   variables, there are n*(n+1)/2 unique components, ordered as: variances
#'   first (diagonals), then covariances (upper triangle by row).
#'
#' @param k Number of nearest analogs to return per focal location for kNN
#'   selection modes. Required when \code{select} is \code{"knn_geog"} or
#'   \code{"knn_clim"}; must be \code{NULL} for \code{select = "all"}.
#'
#' @param weight Weighting function for matches, used only when
#'   \code{stat} is \code{"sum_weights"} or \code{"mean_weights"}. One of:
#'   \itemize{
#'     \item \code{"uniform"}: All matches weighted equally (weight = 1.0).
#'     \item \code{"inverse_clim"}: Inverse climate distance,
#'       weight = 1 / (climate_distance + eps), with epsilon given by \code{theta}.
#'     \item \code{"inverse_geog"}: Inverse geographic distance,
#'       weight = 1 / (geographic_distance + eps), with epsilon given by \code{theta}.
#'     \item \code{"gaussian_clim"}: Gaussian kernel on climate distance,
#'       weight = exp(-climate_distance^2 / (2*sigma^2)), with sigma given by \code{theta}.
#'     \item \code{"gaussian_geog"}: Gaussian kernel on geographic distance,
#'       weight = exp(-geographic_distance^2 / (2*sigma^2)), with sigma given by \code{theta}.
#'     \item \code{"gaussian_joint"}: Joint Gaussian kernel,
#'       weight = exp(-climate_distance^2/(2*sigma_c^2) - geographic_distance^2/(2*sigma_g^2)),
#'       with sigma values given by \code{theta} as a 2-element vector c(sigma_clim, sigma_geog).
#'     \item \code{"inverse_joint"}: Joint inverse distance,
#'       weight = 1 / sqrt((climate_distance + eps_clim)^2 + (geographic_distance + eps_geog)^2),
#'       with epsilon values given by \code{theta} as a 2-element vector c(eps_clim, eps_geog).
#'   }
#'   For \code{stat} in \code{NULL}, \code{"pairs"}, or \code{"count"},
#'   \code{weight} must be \code{NULL}.
#'
#' @param theta Optional numeric parameter used by weighting functions
#'   when \code{stat} is \code{"sum_weights"} or \code{"mean_weights"} and
#'   \code{weight} is not \code{"uniform"}. Interpretation depends on \code{weight}:
#'   \itemize{
#'     \item For \code{"inverse_clim"} or \code{"inverse_geog"}: epsilon value
#'       added to distances (scalar; default: 1e-12 for climate, 1e-6 for geography).
#'     \item For \code{"gaussian_clim"} or \code{"gaussian_geog"}: sigma bandwidth
#'       parameter (scalar; larger values = slower decay with distance).
#'     \item For \code{"gaussian_joint"}: 2-element vector c(sigma_clim, sigma_geog)
#'       giving bandwidths for climate and geographic dimensions.
#'     \item For \code{"inverse_joint"}: 2-element vector c(eps_clim, eps_geog)
#'       giving epsilon values for climate and geographic dimensions.
#'   }
#'   If \code{theta} is \code{NULL}, sensible defaults are used for single-parameter
#'   weights. For \code{weight = "uniform"} or for non-weighted aggregations, \code{theta}
#'   must be \code{NULL}.
#'
#' @param report_dist Logical; if TRUE (default), include distance columns in
#'   output when \code{stat} is \code{NULL} or \code{"pairs"}. Set to FALSE
#'   for more compact output.
#'
#' @param coord_type Coordinate system type (default: "auto"):
#'   \itemize{
#'     \item \code{"auto"}: Automatically detect from coordinate ranges.
#'     \item \code{"lonlat"}: Unprojected lon/lat coordinates (uses great-circle distance;
#'       assumes \code{max_geog} is in km).
#'     \item \code{"projected"}: Projected XY coordinates (uses planar distance;
#'       assumes \code{max_geog} is in projection units).
#'   }
#'
#' @param index_res Tuning parameter giving the number of bins per dimension
#'   of the internally-used lattice search index. Either:
#'   \itemize{
#'     \item A positive integer.
#'     \item \code{"auto"} (the default): Automatically tune the index resolution
#'       by optimizing compute time on a subsample of focal points. If focal has
#'       relatively few rows, auto-tuning is skipped and a default resolution of
#'       16 is used.
#'   }
#'   Ignored if \code{pool} is an \code{analog_index} (uses index's resolution).
#'
#' @param n_threads Optional integer number of threads to use for the
#'   computation. If \code{NULL} (default), the global RcppParallel setting
#'   is used (see \code{RcppParallel::setThreadOptions}).
#'
#' @return
#' The return value depends on the \code{stat} parameter:
#'
#' **For stat = NULL or "pairs"**:
#' A data.frame with one row per focal-analog pair:
#' \itemize{
#'   \item \code{index}: Index of focal location (1-based).
#'   \item \code{x, y}: Coordinates of focal location.
#'   \item \code{analog_index}: Index of analog location in reference dataset (1-based).
#'   \item \code{analog_x, analog_y}: Coordinates of analog location.
#'   \item \code{clim_dist}: Climate distance (if \code{report_dist = TRUE}).
#'   \item \code{geog_dist}: Geographic distance in km (if \code{report_dist = TRUE}).
#' }
#'
#' **For stat = "sum_weights", "mean_weights", or "count"**:
#' A data.frame with one row per focal location:
#' \itemize{
#'   \item \code{index}: Index of focal location (1-based).
#'   \item \code{x, y}: Coordinates of focal location.
#'   \item \code{value}: Aggregated value (count, sum of weights, or mean of weights).
#' }
#'
#' All outputs include diagnostic attributes propagated from the C++ core.
#'
#' @examples
#' \dontrun{
#' # Basic pair queries
#' analog_search(x = focal, pool = ref, select = "all", max_clim = 0.5)
#' analog_search(x = focal, pool = ref, select = "knn_geog", max_clim = 0.5, k = 1)
#'
#' # Aggregated queries
#' analog_search(x = focal, pool = ref, select = "all", stat = "count", max_clim = 0.5)
#' analog_search(x = focal, pool = ref, select = "all", stat = "mean_weights",
#'               max_clim = 0.5, weight = "gaussian_clim", theta = 0.1)
#'
#' # With pre-built index (for repeated queries)
#' index <- build_analog_index(ref)
#' analog_search(x = focal1, pool = index, select = "knn_geog", max_clim = 0.5, k = 1)
#' analog_search(x = focal2, pool = index, select = "all", stat = "count", max_clim = 0.3)
#' }
#'
#' @export
analog_search <- function(
            x,
            pool,
            select = "all",
            stat = NULL,
            max_clim = NULL,
            max_geog = NULL,
            x_cov = NULL,
            k = NULL,
            weight = NULL,
            theta = NULL,
            report_dist = TRUE,
            coord_type = c("auto", "lonlat", "projected"),
            index_res = "auto",
            n_threads = NULL
) {

      # Validate and normalize query parameters
      params <- .validate_query_params(select, stat, k, weight, theta)
      select <- params$select
      stat <- params$stat
      k <- params$k
      weight <- params$weight
      theta <- params$theta

      # Validate coord_type
      coord_type <- match.arg(coord_type)

      # Check if pool is already an index
      if (is_analog_index(pool)) {
            # Pool is pre-built index - use it directly
            index <- pool

      } else {
            # Pool is raw data - need to build index

            # Tune resolution if needed
            if (identical(index_res, "auto")) {
                  # Use tune_index_res for auto-tuning
                  index_res_int <- tune_index_res(
                        x = x,
                        pool = pool,
                        x_cov = x_cov,
                        select = select,
                        stat = stat,
                        max_clim = max_clim,
                        max_geog = max_geog,
                        k = k,
                        weight = weight,
                        theta = theta,
                        coord_type = coord_type,
                        n_threads = n_threads,
                        verbose = FALSE
                  )
            } else if (is.numeric(index_res)) {
                  index_res_int <- as.integer(index_res)
            } else {
                  index_res_int <- 16L  # Default
            }

            # Build index from raw pool data
            index <- build_analog_index(
                  pool = pool,
                  coord_type = coord_type,
                  index_res = index_res_int
            )
      }

      # Query the index
      return(query_analog_index(
            x = x,
            index = index,
            select = select,
            stat = stat,
            max_clim = max_clim,
            max_geog = max_geog,
            x_cov = x_cov,
            k = k,
            weight = weight,
            theta = theta,
            report_dist = report_dist,
            n_threads = n_threads
      ))
}
