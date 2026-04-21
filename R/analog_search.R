#' Find Climate Analogs
#'
#' Identifies locations in a reference dataset that are climatically similar
#' and/or geographically proximal to focal locations. Analog searches use a
#' two-stage approach: first selecting analogs based on specified criteria,
#' then optionally aggregating the results.
#'
#' ## Parameter categories
#'
#' \itemize{
#'     \item *Data parameters*
#'       (`x`, `pool`, `x_cov`, `y`, `covariates`, `coord_type`)
#'       give attributes of the data on which to operate.
#'     \item *Selection parameters*
#'       (`select`, `max_clim`, `max_geog`, `k`)
#'       define which analogs to `select` from the `pool` for each `x`.
#'     \item *Aggregation parameters*
#'       (`stat`, `weight`, `theta`, `lambda`)
#'       control how selected analogs are summarized.
#'     \item *Computation parameters*
#'       (`n_threads`, `index_res`, `downsample`, `seed`, `progress`)
#'       specify behavior for controlling compute performance.
#' }
#'
#' ## Distance metrics
#'
#' Geographic distance can be computed for lon/lat coordinates (great-circle
#' distance) or projected coordinates (planar distance).
#'
#' Climate similarity is measured using Euclidean or Mahalanobis distance in
#' climate space. In general, when multiple climate variables are used, it is
#' recommended to use pre-whitened (scaled) climate data, to avoid major artifacts
#' from climate variables with different units. Pre-whitening can be done using
#' `scale()` for dataset-wide Euclidean distances, or `mahalanobis_transform()`
#' for dataset-wide Mahalanobis distances.
#'
#' The function also supports climate distance calculations based on
#' *local temporal* covariance structure at focal locations, via the `x_cov`
#' parameter. These local covariance values need to be pre-calculated.
#'
#'
#' ## Computational optimization
#'
#' The analog search architecture is designed with compute performance in mind:
#' \itemize{
#'     \item All internal computations are done in C++.
#'     \item Searches use a lattice-based indexing structure to efficiently
#'       search through large reference datasets. By default, the lattice
#'       resolution is tuned for optimal performance.
#'     \item Parallel processing is available via the `threads` parameter.
#'     \item You can `downsample` prohibitively large reference pool datasets
#'       to improve speed and memory, using a stratified sampling
#'       scheme that reduces loss of precision relative to random sampling.
#'     \item For large datasets, enable `progress = TRUE` to display
#'       a progress bar during computation.
#'     \item For raster datasets that are too large to fit in memory,
#'       `tiled_analog_search()` offers a memory-safe option.
#' }
#'
#' @param x Focal locations for which analogs will be found. Should be a
#'   matrix/data.frame with columns x, y, and climate variables, or a
#'   SpatRaster with climate variable layers.
#'
#' @param pool The reference dataset to search for analogs. Either:
#'   \itemize{
#'     \item Matrix/data.frame with columns x, y, and climate variables,
#'       or SpatRaster with climate variable layers, OR
#'     \item An `analog_index` object created by
#'       [build_analog_index()] (for repeated queries).
#'   }
#'
#' @param x_cov Optional focal-specific covariance matrices for Mahalanobis
#'   distance calculations. Should be a matrix or data.frame with one row per
#'   focal location and one column per unique covariance component, or a
#'   SpatRaster with a layer for each component. For n climate variables,
#'   there are n*(n+1)/2 unique components, ordered as: variances first
#'   (diagonals), then covariances (upper triangle by row).
#'
#' @param y Optional user-defined variables for each reference location
#'   in `pool` to aggregate across selected analogs. Can be a numeric vector
#'   (single variable), matrix or data.frame with numeric columns (multiple
#'   variables), or a SpatRaster with one or more numeric layers. Must have
#'   exactly the same number of reference locations as `pool`.
#'
#'   When provided, enables value-based aggregation stats `"sum"`, `"mean"`,
#'   `"weighted_sum"`, `"weighted_mean"`, and `"regression"`. For stat =
#'   NULL/"none" (pairs mode), y value columns are included in output for each
#'   analog pair.
#'
#' @param covariates Optional auxiliary predictor variables for each reference
#'   location in `pool`, used with `stat = "regression"`. Can be a numeric
#'   vector (single covariate), matrix or data.frame with numeric columns
#'   (multiple covariates), or a SpatRaster with one or more numeric layers.
#'   Must have exactly the same number of rows/cells as `pool`. Column names
#'   are used for output layer naming. These variables are NOT used for the
#'   analog search itself -- only for local regression within each analog
#'   neighborhood.
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
#' @param select Character string specifying the analog selection strategy.
#'   One of:
#'   \itemize{
#'     \item `"all"` (default): Select all analogs that satisfy the
#'       `max_clim` and `max_geog` constraints.
#'     \item `"knn_clim"`: For each focal, select up to `k` analogs
#'       with smallest climate distance, subject to filters.
#'     \item `"knn_geog"`: For each focal, select up to `k` analogs
#'       with smallest geographic distance, subject to filters.
#'   }
#'
#' @param k Number of nearest analogs to return per focal location for kNN
#'   selection modes. Required when \code{select} is \code{"knn_geog"} or
#'   \code{"knn_clim"}; must be \code{NULL} for \code{select = "all"}.
#'
#' @param stat Statistic(s) used to aggregate selected analogs. Either:
#'   \itemize{
#'     \item \code{NULL} or \code{"none"}: Return all selected analog pairs as a data.frame.
#'     \item \code{"count"}: For each focal, count the number of selected analogs.
#'     \item \code{"sum_weights"}: For each focal, sum the weights of selected
#'       analogs (see \code{weight} and \code{theta}).
#'     \item \code{"mean_weights"}: For each focal, mean of weights of selected
#'       analogs.
#'     \item \code{"sum"}: Sum of \code{y} values across analogs (requires \code{y}).
#'     \item \code{"mean"}: Mean of \code{y} values across analogs (requires \code{y}).
#'     \item \code{"weighted_sum"}: Sum of (\code{y} × weight) across analogs
#'       (requires \code{y} and \code{weight}).
#'     \item \code{"weighted_mean"}: Weighted mean of \code{y} values across analogs
#'       (requires \code{y} and \code{weight}).
#'     \item \code{"ess"}: Kish's effective sample size (ESS), computed as the
#'       squared sum of weights divided by the sum of squared weights
#'       (requires \code{weight}).
#'     \item \code{"regression"}: Weighted least squares (or ridge) regression
#'       of \code{y} on \code{covariates} within each analog neighborhood.
#'       Returns intercept and slope coefficients. Requires \code{y},
#'       \code{covariates}, and \code{weight}. See \code{lambda} for
#'       regularization.
#'     \item A character vector combining multiple stats (e.g.,
#'       \code{c("count", "weighted_mean", "regression")}).
#'       Note: \code{"none"} cannot be combined with other stats.
#'   }
#'
#' @param weight Weighting function for matches, used only when
#'   \code{stat} includes \code{"sum_weights"} or \code{"mean_weights"}. One of:
#'   \itemize{
#'     \item \code{"uniform"}: All matches weighted equally (weight = 1.0).
#'     \item \code{"inverse_clim"}: Inverse climate distance,
#'       weight = 1 / (climate_distance + eps), with epsilon given by \code{theta}.
#'     \item \code{"inverse_geog"}: Inverse geographic distance,
#'       weight = 1 / (geographic_distance + eps), with epsilon given by \code{theta}.
#'     \item \code{"gaussian_clim"}: Gaussian kernel on climate distance,
#'       weight = exp(-climate_distance^2 / (2 sigma^2)), with sigma given by \code{theta}.
#'     \item \code{"gaussian_geog"}: Gaussian kernel on geographic distance,
#'       weight = exp(-geographic_distance^2 / (2 sigma^2)), with sigma given by \code{theta}.
#'     \item \code{"gaussian_joint"}: Gaussian kernel on combined distance,
#'       weight = exp(-(clim_dist^2 / (2 sigma_clim^2) + geog_dist^2 / (2 sigma_geog^2))),
#'       with sigmas given by \code{theta}.
#'     \item \code{"inverse_joint"}: Inverse joint distance,
#'       weight = 1 / (sqrt(clim_dist^2 + geog_dist^2) + eps), with epsilon given by \code{theta}.
#'   }
#'
#' @param theta Optional numeric parameter used by weighting functions
#'   when \code{stat} includes \code{"sum_weights"} or \code{"mean_weights"} and
#'   \code{weight} is not \code{"uniform"}. Interpretation depends on \code{weight}:
#'   \itemize{
#'     \item For \code{"inverse_clim"} or \code{"inverse_geog"}: epsilon value
#'       added to distances (scalar; default: 1e-12 for climate, 1e-6 for geography).
#'     \item For \code{"gaussian_clim"} or \code{"gaussian_geog"}: sigma bandwidth
#'       parameter (scalar; larger values = slower decay with distance).
#'     \item For \code{"gaussian_joint"} or \code{"inverse_joint"}: 2-element vector
#'       \code{c(theta_clim, theta_geog)} (defaults: 1 for climate, 1 for geography).
#'   }
#'
#' @param lambda Ridge penalty parameter for \code{stat = "regression"}
#'   (default: 0, giving ordinary weighted least squares). Higher values
#'   shrink covariate coefficients toward zero, with the intercept
#'   approaching the weighted mean as \code{lambda -> Inf}. Ignored when
#'   \code{"regression"} is not in \code{stat}.
#'
#' @param coord_type Coordinate system type:
#'   \itemize{
#'     \item \code{"auto"} (default): Automatically detect from coordinate ranges.
#'     \item \code{"lonlat"}: Unprojected lon/lat coordinates (uses great-circle distance;
#'       assumes \code{max_geog} is in km).
#'     \item \code{"projected"}: Projected XY coordinates (uses planar distance;
#'       assumes \code{max_geog} is in projection units).
#'   }
#'
#' @param downsample Optional downsampling rate (0-1) for the reference pool,
#'   indicating the proportion of points to retain. Values < 1 reduce memory
#'   and improve speed at some cost to precision. Default is 1.0 (no downsampling).
#'   Ignored if \code{pool} is a pre-built index.
#'
#' @param seed Optional random seed for reproducible downsampling. If \code{NULL}
#'   (default), uses current R random state. Ignored if \code{pool} is a pre-built
#'   index or \code{downsample = 1}.
#'
#' @param index_res Tuning parameter giving the number of bins per dimension
#'   of the internally-used lattice search index. Either:
#'   \itemize{
#'     \item A positive integer.
#'     \item `"auto"` (the default): Automatically tune the index resolution
#'       by optimizing compute time on a subsample of focal points. If focal has
#'       relatively few rows, auto-tuning is skipped and a default resolution of
#'       16 is used.
#'   }
#'   Ignored if `pool` is an `analog_index` (uses index's resolution).
#'
#' @param n_threads Optional integer number of threads to use for the
#'   computation. If `NULL` (default), the global RcppParallel setting
#'   is used (see `RcppParallel::setThreadOptions`).
#'
#' @param progress Logical; if `TRUE`, display a progress bar during
#'   computation. Progress tracking works by splitting the focal dataset into
#'   chunks and processing them sequentially. Useful for large datasets. Default
#'   is `FALSE`.
#'
#' @return Return type depends on input format and query mode.
#'
#' Returns a data.frame, unless `x` is a SpatRaster and results have exactly one record per
#' input cell (aggregation mode, or pairwise with `k = 1`), in which case returns a
#' SpatRaster with one layer per output variable.
#'
#' Pairwise mode (`stat = NULL` or `"none"`) returns one row per focal-analog pair,
#' with the following variables:
#' \itemize{
#'     \item `index`, `x`, `y`: Focal location (1-based index and coordinates) corresponding to input `x`
#'     \item `analog_index`, `analog_x`, `analog_y`: Analog location corresponding to input `pool`
#'     \item `clim_dist`: Climate distance (Euclidean or Mahalanobis)
#'     \item `geog_dist`: Geographic distance (km for lonlat, projection units otherwise)
#'     \item Value columns (if `y` provided): one per variable
#' }
#'
#' Aggregation mode (one or more `stat` values) returns one row per focal location,
#' with the following variables:
#' \itemize{
#'     \item `index`, `x`, `y`: Focal location
#'     \item One column per requested statistic. For `stat` with single `y` variable:
#'       column named by stat (e.g., `sum`, `mean`). For `stat` with multiple `y`
#'       variables: columns named `{stat}_{varname}` (e.g., `sum_biomass`, `mean_richness`)
#'     \item For `stat = "regression"`: columns for `coef_intercept` and `coef_{covariate}`,
#'       or `coef_intercept_{varname}` and  `coef_{covariate}_{varname}` with multiple `y`
#'       variables.
#' }
#'
#' All results include metadata attributes (`select`, `stat`, `weight`, etc.).
#' Use [analog_summary()] to view a formatted summary.
#'
#' @references
#' Mahony CR, Cannon AJ, Wang T, Aitken SN (2017). "A closer look at novel
#' climates: new methods and insights at continental to landscape scales."
#' \emph{Global Change Biology}, \strong{23}(9), 3934-3955.
#' \doi{10.1111/gcb.13645}
#'
#' Hamann A, Roberts DR, Barber QE, Carroll C, Nielsen SE (2015). "Velocity of
#' climate change algorithms for guiding conservation and management."
#' \emph{Global Change Biology}, \strong{21}(2), 997-1004.
#' \doi{10.1111/gcb.12736}
#'
#' Grenier P, Parent A-C, Huard D, Anctil F, Chaumont D (2013). "An assessment
#' of six dissimilarity metrics for climate analogs." \emph{Journal of Applied
#' Meteorology and Climatology}, \strong{52}(4), 733-752.
#' \doi{10.1175/JAMC-D-12-0170.1}
#'
#' @seealso [tiled_analog_search()] offers memory-safe searches on large raster
#'   datasets. Helper functions such as [analog_impact()], [analog_velocity()],
#'   and [analog_intensity()] offer simpler interfaces for common search types.
#'
#' @export
analog_search <- function(

      # data
      x,
      pool,
      x_cov = NULL,
      y = NULL,
      covariates = NULL,

      # candidate filtering
      max_clim = NULL,
      max_geog = NULL,
      select = "all",
      k = NULL,

      # analog aggregation
      stat = NULL,
      weight = NULL,
      theta = NULL,
      lambda = 0,

      # args passed to build_lattice_index
      coord_type = c("auto", "lonlat", "projected"),
      downsample = 1.0,
      seed = NULL,
      index_res = "auto",

      n_threads = NULL,
      progress = FALSE
) {

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
                        y = y,
                        covariates = covariates,
                        select = select,
                        stat = stat,
                        max_clim = max_clim,
                        max_geog = max_geog,
                        k = k,
                        weight = weight,
                        theta = theta,
                        lambda = lambda,
                        coord_type = coord_type,
                        n_threads = n_threads,
                        verbose = FALSE,
                        downsample = downsample,
                        seed = seed
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
                  index_res = index_res_int,
                  downsample = downsample,
                  seed = seed
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
            y = y,
            covariates = covariates,
            k = k,
            weight = weight,
            theta = theta,
            lambda = lambda,
            n_threads = n_threads,
            show_progress = progress
      ))
}
