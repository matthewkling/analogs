#' Find Climate Analogs
#'
#' Identifies locations in a reference dataset that are climatically similar to
#' focal locations, with optional constraints on climate distance and geographic
#' distance. This function supports multiple use cases including climate velocity
#' analysis, analog availability mapping, and climate impact assessment.
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
#' @param mode Character string specifying the analog search mode. One of:
#'   \itemize{
#'     \item \code{"knn_clim"}: For each focal, return up to \code{k} analogs
#'       with smallest climate distance, subject to \code{max_clim} and
#'       \code{max_geog} filters.
#'     \item \code{"knn_geog"}: For each focal, return up to \code{k} analogs
#'       with smallest geographic distance, subject to \code{max_clim} and
#'       \code{max_geog} filters.
#'     \item \code{"all"}: Return all analogs that satisfy the filters.
#'     \item \code{"count"}: For each focal, count how many analogs satisfy
#'       the filters.
#'     \item \code{"sum"}: For each focal, sum weights of all analogs that
#'       satisfy the filters (see \code{weight} and \code{theta}).
#'     \item \code{"mean"}: For each focal, mean of weights of all analogs that
#'       satisfy the filters.
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
#'   first (diagonals), then covariances (upper triangle by row). For example:
#'   \itemize{
#'     \item 2 variables: c(var1, var2, cov12)
#'     \item 3 variables: c(var1, var2, var3, cov12, cov13, cov23)
#'   }
#'   When provided, all climate distances are computed as Mahalanobis distances
#'   using each focal's covariance structure. For focal-specific variances only
#'   (no covariance), set off-diagonal covariances to zero. Default is NULL
#'   (Euclidean climate distance).
#'
#' @param k Number of nearest analogs to return per focal location for kNN
#'   modes. Required when \code{mode} is \code{"knn_geog"} or \code{"knn_clim"};
#'   must be \code{NULL} for other modes.
#'
#' @param weight Weighting function for matches, used only when
#'   \code{mode} is \code{"sum"} or \code{"mean"}. One of:
#'   \itemize{
#'     \item \code{"uniform"}: All matches weighted equally (weight = 1.0).
#'     \item \code{"inverse_clim"}: Weight = 1 / (climate_distance + epsilon),
#'       with epsilon given by \code{theta} (or a small default if \code{theta}
#'       is \code{NULL}).
#'     \item \code{"inverse_geog"}: Weight = 1 / (geographic_distance + epsilon),
#'       with epsilon given by \code{theta} (or a small default if \code{theta}
#'       is \code{NULL}).
#'     \item \code{"gaussian_clim"}: Gaussian kernel on climate distance,
#'       weight = exp(-climate_distance^2 / (2*sigma^2)), with sigma (bandwidth)
#'       given by \code{theta}.
#'     \item \code{"gaussian_geog"}: Gaussian kernel on geographic distance,
#'       weight = exp(-geographic_distance^2 / (2*sigma^2)), with sigma (bandwidth)
#'       given by \code{theta}.
#'     \item \code{"gaussian_joint"}: Bivariate Gaussian kernel (product of
#'       independent Gaussians over climate and geographic distances), with
#'       bandwidths given by \code{theta} as a 2-element vector c(sigma_clim, sigma_geog).
#'     \item \code{"inverse_joint"}: Joint inverse distance,
#'       weight = 1 / sqrt((climate_distance + eps_clim)^2 + (geographic_distance + eps_geog)^2),
#'       with epsilon values given by \code{theta} as a 2-element vector c(eps_clim, eps_geog).
#'   }
#'   For \code{mode} in \code{"knn_geog"}, \code{"knn_clim"}, \code{"count"},
#'   or \code{"all"}, \code{weight} must be \code{NULL}.
#'
#' @param theta Optional numeric parameter used by weighting functions
#'   when \code{mode} is \code{"sum"} or \code{"mean"} and \code{weight} is
#'   not \code{"uniform"}. Interpretation depends on \code{weight}:
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
#'   weights. For \code{weight = "uniform"} or for non-aggregating modes, \code{theta}
#'   must be \code{NULL}.
#'
#' @param report_dist Logical; if TRUE (default), include distance columns in
#'   output when \code{mode} is \code{"knn_geog"}, \code{"knn_clim"} or
#'   \code{"all"}. Set to FALSE for more compact output.
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
#'     \item \code{"auto"} (the default): Automatically tune the intex resolution
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
#' The return value depends on the \code{mode} parameter:
#'
#' **For mode = "knn_geog", "knn_clim" or "all"**:
#' A data.frame with one row per focal-analog pair:
#' \itemize{
#'   \item \code{focal_index}: Index of focal location (1-based).
#'   \item \code{focal_x, focal_y}: Coordinates of focal location.
#'   \item \code{analog_index}: Index of analog location in reference dataset (1-based).
#'   \item \code{analog_x, analog_y}: Coordinates of analog location.
#'   \item \code{clim_dist}: Climate distance (if \code{report_dist = TRUE}).
#'   \item \code{geog_dist}: Geographic distance in km (if \code{report_dist = TRUE}).
#' }
#'
#' **For mode = "sum", "mean", or "count"**:
#' A data.frame with one row per focal location:
#' \itemize{
#'   \item \code{focal_index}: Index of focal location (1-based).
#'   \item \code{focal_x, focal_y}: Coordinates of focal location.
#'   \item \code{value}: Aggregated value (count, sum of weights, or mean of weights).
#' }
#'
#' All outputs include diagnostic attributes propagated from the C++ core,
#' including:
#' \itemize{
#'   \item \code{total_bins}: Number of spatial bins in the lattice index.
#'   \item \code{avg_bin_occupancy}: Average points per bin.
#'   \item \code{min_bin_occupancy, max_bin_occupancy}: Range of bin occupancy.
#'   \item \code{binning_method}: Method used ("multi_dim_lattice" or "none").
#'   \item \code{n_ref, n_clim}: Size of reference dataset and number of climate variables.
#' }
#'
#' @details
#' **Common Use Cases:**
#'
#' \strong{Climate Velocity} (nearest geographic neighbor with similar climate):
#' \preformatted{
#' analog_search(
#'   x        = clim$clim1,
#'   pool     = clim$clim2,
#'   mode     = "knn_geog",
#'   max_clim = 0.5,
#'   max_geog = NULL,
#'   k        = 1
#' )
#' }
#'
#' \strong{Climate Impact} (climatically similar locations within dispersal range):
#' \preformatted{
#' analog_search(
#'   x        = clim$clim1,
#'   pool     = clim$clim2,
#'   mode     = "knn_clim",
#'   max_clim = 0.5,
#'   max_geog = 100,
#'   k        = 20
#' )
#' }
#'
#' \strong{Analog Availability} (count of suitable locations):
#' \preformatted{
#' analog_search(
#'   x        = clim$clim1,
#'   pool     = clim$clim1,
#'   mode     = "count",
#'   max_clim = 0.5,
#'   max_geog = 100
#' )
#' }
#'
#' \strong{Weighted Analog Intensity} (e.g., distance-weighted availability):
#' \preformatted{
#' analog_search(
#'   x        = clim$clim1,
#'   pool     = clim$clim1,
#'   mode     = "sum",
#'   max_clim = 0.5,
#'   max_geog = 100,
#'   weight   = "inverse_geog",
#'   theta    = 1e-6
#' )
#' }
#'
#' \strong{Using a Pre-built Index} (for repeated queries):
#' \preformatted{
#' # Build index once
#' index <- build_analog_index(clim$clim2, index_res = 16)
#'
#' # Query multiple times
#' v1 <- analog_search(x = sites1, pool = index, mode = "knn_geog", max_clim = 0.5, k = 1)
#' v2 <- analog_search(x = sites2, pool = index, mode = "knn_geog", max_clim = 0.3, k = 1)
#' }
#'
#' \strong{Focal-specific Mahalanobis Distance}:
#' \preformatted{
#' # With focal-specific covariance matrices
#' analog_search(
#'   x        = clim$clim1,
#'   pool     = clim$clim2,
#'   x_cov    = focal_covariances,  # n_focal x 3 matrix for 2 climate vars
#'   mode     = "knn_geog",
#'   max_clim = 2,  # In Mahalanobis distance units
#'   k        = 1
#' )
#' }
#'
#' @export
analog_search <- function(
            x,
            pool,
            mode,
            max_clim = NULL,
            max_geog = NULL,
            x_cov = NULL,
            k = NULL,
            weight = NULL,
            theta = NULL,
            report_dist = TRUE,
            coord_type = c("auto", "lonlat", "projected"),
            index_res = "auto",
            resolutions = NULL, n_sample = NULL, n_reps = NULL,
            n_threads = NULL
) {

      # Validate and normalize query parameters
      params <- .validate_query_params(mode, k, weight, theta)
      mode <- params$mode
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
                        mode = mode,
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

      # Query the index (validation of x_cov happens here)
      return(query_analog_index(
            x = x,
            index = index,
            mode = mode,
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
