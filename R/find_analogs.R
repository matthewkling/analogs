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
#'
#' @param max_geog Maximum geographic distance constraint (default:
#'   NULL = no geographic constraint). When specified, only reference locations
#'   within this distance are considered. Radius units should be specified in
#'   kilometers if \code{coord_type = "lonlat"}, or in projected coordinate units
#'   if \code{coord_type = "projected"}.
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
#'   }
#'   For \code{mode} in \code{"knn_geog"}, \code{"knn_clim"}, \code{"count"},
#'   or \code{"all"}, \code{weight} must be \code{NULL}.
#'
#' @param theta Optional numeric parameter used by some weighting kernels
#'   when \code{mode} is \code{"sum"} or \code{"mean"} and \code{weight} is
#'   not \code{"uniform"}. Currently interpreted as:
#'   \itemize{
#'     \item For \code{"gaussian_clim"}: sigma bandwidth for climate distance.
#'     \item For \code{"gaussian_geog"}: sigma bandwidth for geographic distance.
#'     \item For \code{"gaussian_joint"}: length-two vector of sigma bandwidths
#'       for climate and geographic distances, respectively.
#'     \item For \code{"inverse_clim"}: epsilon added to climate distance.
#'     \item For \code{"inverse_geog"}: epsilon added to geographic distance.
#'   }
#'   If \code{theta} is \code{NULL}, a small default epsilon is used. For
#'   \code{weight = "uniform"} or for non-aggregating modes, \code{theta}
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
#' @param resolution,n_sample,n_reps Parameters passed to \code{tune_index_res()}.
#'   Ignored if \code{pool} is an \code{analog_index} or if \code{index_res} isn't
#'   \code{"auto"}.
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
#' find_analogs(
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
#' find_analogs(
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
#' find_analogs(
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
#' find_analogs(
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
#' v1 <- find_analogs(x = sites1, pool = index, mode = "knn_geog", max_clim = 0.5, k = 1)
#' v2 <- find_analogs(x = sites2, pool = index, mode = "knn_geog", max_clim = 0.3, k = 1)
#' }
#'
#' @export
find_analogs <- function(
            x,
            pool,
            mode = c("knn_clim", "knn_geog", "count", "sum", "mean", "all"),
            max_clim = NULL,
            max_geog = NULL,
            k = NULL,
            weight = c("uniform",
                       "gaussian_clim", "gaussian_geog", "gaussian_joint",
                       "inverse_clim", "inverse_geog", "inverse_joint"),
            theta = NULL,
            report_dist = TRUE,
            coord_type = c("auto", "lonlat", "projected"),
            index_res = "auto",
            resolutions = NULL, n_sample = NULL, n_reps = NULL,
            n_threads = NULL
) {

      # [placeholder: throw warnings if build-specific params and prebuilt index are both supplied]


      # Build index (if needed) --------------------------------------

      if(!is_analog_index(pool)){

            # Tune resolution (if needed) ---------------------------

            if(identical(index_res, "auto")) {
                  # Use tune_index_res for auto-tuning
                  index_res_int <- tune_index_res(
                        x = focal_mm,
                        pool = ref_mm,
                        mode = mode,
                        max_clim = max_clim,
                        max_geog = max_geog,
                        k = k,
                        weight = weight,
                        theta = theta,
                        coord_type = coord_type,
                        resolutions = resolutions,
                        n_sample = n_sample,
                        n_reps = n_reps,
                        verbose = FALSE
                  )
            } else if (is.numeric(index_res)) {
                  index_res_int <- as.integer(index_res)
            } else {
                  index_res_int <- 16L  # Default
            }

            # Build index
            index <- build_analog_index(
                  pool = ref_mm,
                  coord_type = coord_type,
                  index_res = index_res_int
            )
      }


      # Query index --------------------------------------------------

      return(query_analog_index(
            x = x,
            index = pool,
            mode = mode,
            max_clim = max_clim,
            max_geog = max_geog,
            k = k,
            weight = weight,
            theta = theta,
            report_dist = report_dist,
            n_threads = n_threads
      ))

}



