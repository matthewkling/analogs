#' Build Analog Index
#'
#' Pre-builds a reusable lattice index from reference climate data. The index
#' can be queried multiple times with different focal points and parameters,
#' avoiding the need to rebuild the lattice for each query.
#'
#' @param pool The reference dataset to search for analogs. Should be a
#'  matrix/data.frame with columns x, y, and climate variables, or a SpatRaster
#'  with climate variable layers.
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
#' @param downsample Optional downsampling rate (0-1) indicating the proportion of
#'   points in `pool` to retain. Downsampling reduces memory use and improves
#'   query speed at the cost of some precision; adaptive stratified sampling is
#'   used to minimize loss of precision. The default is 1.0 (no downsampling).
#'   See Details for more info.
#'
#' @param seed Optional random seed for reproducible downsampling. If `NULL`
#'   (default), uses current R random state.
#'
#' @return An S3 object of class \code{"analog_index"} containing:
#'   \itemize{
#'     \item The compiled lattice index (internal C++ structure)
#'     \item Reference data
#'     \item Metadata: coordinate type, dimensions, ranges, resolution
#'     \item Diagnostics: bin counts, occupancy statistics, and downsampling info
#'   }
#'
#' @details
#' The lattice index is a multidimensional grid of bins, built over both geographic
#' and climate dimensions. This structure enables efficient analog searches by
#' first filtering and sorting bins of similar points before computing exact
#' results. For lon/lat coordinates, the index uses ECEF (Earth-Centered Earth-Fixed)
#' space internally for optimal performance.
#'
#' Index resolution (`index_res`) controls the granularity of spatial
#' binning. The optimal value depends on your data size and query patterns.
#' Use [tune_index_res()] to find the best resolution for your use case,
#' or accept the default of 16 which works well for many applications.
#'
#' ## Downsampling
#'
#' For very large datasets, downsampling can significantly improve memory usage
#' and query speed, at the cost of some precision. The `downsample` parameter controls
#' the target fraction of the data points in `pool` that are retained in the index.
#' Downsampling uses an adaptive stratified approach: densely-packed bins are thinned more
#' aggressively while sparse bins are preserved, which helps reduce imprecision
#' in sparse regions compared to fully random sampling. Note: The actual rate may be
#' higher than requested if maintaining at least one point per occupied bin requires
#' it (common with sparse data or fine-grained binning); check `index$downsample_actual`.
#'
#' Each remaining analog in the downsampled pool gets a `sample_weight` indicating
#' the number of points it represents in the original pool; this weight is the inverse
#' of the sampling rate in the analog's index bin. For pair queries (`stat = "none"`),
#' results include each analog's `sample_weight`. For aggregation stats (count, sum,
#' mean, etc.), sampling weights are used internally to automatically correct for
#' the downsampling bias.
#'
#' @examples
#' \dontrun{
#' # Build index with default settings
#' index <- build_analog_index(climate_data)
#'
#' # Build with explicit resolution
#' index <- build_analog_index(climate_data, index_res = 20)
#'
#' # Build with downsampling for large datasets
#' index <- build_analog_index(
#'   large_climate_data,
#'   index_res = 16,
#'   downsample = 0.1,  # Reduce max bin size to 10%
#'   seed = 123         # Reproducible sampling
#' )
#'
#' # Query the index multiple times
#' v1 <- analog_velocity(sites1, pool = index, max_clim = 0.5)
#' v2 <- analog_velocity(sites2, pool = index, max_clim = 0.3)
#' a1 <- analog_availability(sites3, pool = index, max_clim = 0.5, max_geog = 100)
#' }
#'
#' @export
build_analog_index <- function(pool,
                               coord_type = c("auto", "lonlat", "projected"),
                               index_res = 16,
                               downsample = 1.0,
                               seed = NULL) {

      # Validate inputs
      coord_type <- match.arg(coord_type)

      if (!is.numeric(index_res) || length(index_res) != 1L || index_res <= 0) {
            stop("index_res must be a positive integer")
      }
      index_res <- as.integer(index_res)

      if (!is.numeric(downsample) || length(downsample) != 1L ||
          downsample <= 0 || downsample > 1) {
            stop("downsample must be a number between 0 and 1 (exclusive of 0)")
      }

      # Handle seed
      if (!is.null(seed)) {
            if (!is.numeric(seed) || length(seed) != 1L) {
                  stop("seed must be a single integer or NULL")
            }
            seed <- as.integer(seed)
      } else {
            # Use R's current RNG state to generate a seed
            seed <- as.integer(sample.int(.Machine$integer.max, 1))
      }

      # Normalize pool to standard matrix format
      ref_mm <- .format_data(pool)

      # Detect coordinate system if auto
      if (coord_type == "auto") {
            coord_type <- .detect_geo(ref_mm[, 1:2, drop = FALSE])
      }

      # Call C++ to build index
      index_cpp <- build_analog_index_cpp(
            ref_mm = ref_mm,
            coord_type = coord_type,
            index_res = index_res,
            downsample = downsample,
            seed = seed
      )

      # Construct S3 object
      obj <- structure(
            list(
                  # C++ components
                  lattice_xptr = index_cpp$lattice_xptr,
                  ecef_xptr = index_cpp$ecef_xptr,

                  # Reference data (needed for distance calculations)
                  ref_data = ref_mm,

                  # Metadata
                  coord_type = coord_type,
                  n_ref = index_cpp$n_ref,
                  n_clim = index_cpp$n_clim,
                  index_res = index_res,

                  # Coordinate ranges
                  coord_mins = index_cpp$coord_mins,
                  coord_maxs = index_cpp$coord_maxs,

                  # Climate ranges
                  clim_mins = index_cpp$clim_mins,
                  clim_maxs = index_cpp$clim_maxs,

                  # Internal flags
                  use_ecef = index_cpp$use_ecef,

                  # Diagnostics
                  total_bins = index_cpp$total_bins,
                  n_cells_nonempty = index_cpp$n_cells_nonempty,
                  min_cell_occ = index_cpp$min_cell_occ,
                  max_cell_occ = index_cpp$max_cell_occ,

                  # Downsampling info
                  downsample_target = downsample,
                  downsample_actual = index_cpp$downsample_actual,
                  seed = seed
            ),
            class = "analog_index"
      )

      return(obj)
}


#' Print Method for Analog Index
#'
#' @param x An \code{analog_index} object
#' @param ... Additional arguments (ignored)
#'
#' @export
print.analog_index <- function(x, ...) {
      cat("Analog Index\n")
      cat("============\n\n")

      cat("Reference data:\n")
      cat(sprintf("  %d locations\n", x$n_ref))
      cat(sprintf("  %d climate variables\n", x$n_clim))
      cat(sprintf("  Coordinate type: %s\n", x$coord_type))

      cat("\nCoordinate ranges:\n")
      cat(sprintf("  X: [%.3f, %.3f]\n", x$coord_mins[1], x$coord_maxs[1]))
      cat(sprintf("  Y: [%.3f, %.3f]\n", x$coord_mins[2], x$coord_maxs[2]))

      cat("\nClimate ranges:\n")
      for (i in seq_along(x$clim_mins)) {
            cat(sprintf("  Variable %d: [%.3f, %.3f]\n",
                        i, x$clim_mins[i], x$clim_maxs[i]))
      }

      cat("\nIndex structure:\n")
      cat(sprintf("  Resolution: %d bins per dimension\n", x$index_res))
      cat(sprintf("  Total bins: %s\n", format(x$total_bins, big.mark = ",")))
      cat(sprintf("  Non-empty bins: %s (%.1f%%)\n",
                  format(x$n_cells_nonempty, big.mark = ","),
                  100 * x$n_cells_nonempty / x$total_bins))

      if (x$n_cells_nonempty > 0) {
            avg_occ <- x$n_ref / x$n_cells_nonempty
            cat(sprintf("  Bin occupancy: %d-%d (avg: %.1f) points per bin\n",
                        x$min_cell_occ, x$max_cell_occ, avg_occ))
      }

      # Show downsampling info if applied
      if (!is.null(x$downsample_actual) && x$downsample_actual < 1.0) {
            cat(sprintf("\nDownsampling:\n"))
            cat(sprintf("  Target rate: %.1f%%\n", x$downsample_target * 100))
            cat(sprintf("  Actual rate: %.1f%% of points retained\n",
                        x$downsample_actual * 100))
            cat(sprintf("  Seed: %d\n", x$seed))
      }

      if (x$use_ecef) {
            cat("\n  Note: Using ECEF coordinates internally for lon/lat data\n")
      }

      invisible(x)
}


#' Check if Object is an Analog Index
#'
#' @param x Object to test
#' @return Logical indicating whether \code{x} is an \code{analog_index}
#'
#' @export
is_analog_index <- function(x) {
      inherits(x, "analog_index")
}


#' Validate Analog Index
#'
#' Internal function to validate that an analog index is well-formed and
#' compatible with query data.
#'
#' @param index An \code{analog_index} object
#' @param query_data Optional matrix of query points to validate against
#' @param validate_ranges Logical; if TRUE and query_data provided, check
#'   that query points fall within index bounds
#'
#' @return Invisible TRUE if valid, otherwise throws an error
#' @keywords internal
.validate_analog_index <- function(index,
                                   query_data = NULL,
                                   validate_ranges = FALSE) {

      # Check class
      if (!is_analog_index(index)) {
            stop("Object is not an analog_index")
      }

      # Check required components
      required <- c("lattice_xptr", "ref_data", "coord_type",
                    "n_ref", "n_clim", "index_res")
      missing <- setdiff(required, names(index))
      if (length(missing) > 0) {
            stop("Invalid analog_index: missing components: ",
                 paste(missing, collapse = ", "))
      }

      # Check Xptr validity
      if (is.null(index$lattice_xptr)) {
            stop("Invalid analog_index: lattice pointer is NULL")
      }

      # Validate against query data if provided
      if (!is.null(query_data)) {
            query_ncol <- ncol(query_data)
            expected_ncol <- 2 + index$n_clim

            if (query_ncol != expected_ncol) {
                  stop(sprintf(
                        "Query data has %d columns but index expects %d (2 coords + %d climate)",
                        query_ncol, expected_ncol, index$n_clim
                  ))
            }

            # Optionally check if query points are within index bounds
            if (validate_ranges) {
                  query_x <- query_data[, 1]
                  query_y <- query_data[, 2]

                  x_out <- query_x < index$coord_mins[1] | query_x > index$coord_maxs[1]
                  y_out <- query_y < index$coord_mins[2] | query_y > index$coord_maxs[2]

                  n_out <- sum(x_out | y_out)
                  if (n_out > 0) {
                        warning(sprintf(
                              "%d query point(s) fall outside index coordinate bounds",
                              n_out
                        ))
                  }

                  # Check climate bounds (warnings only)
                  for (i in seq_len(index$n_clim)) {
                        clim_col <- query_data[, 2 + i]
                        clim_out <- clim_col < index$clim_mins[i] |
                              clim_col > index$clim_maxs[i]
                        n_clim_out <- sum(clim_out)
                        if (n_clim_out > 0) {
                              warning(sprintf(
                                    "%d query point(s) have climate variable %d outside index bounds",
                                    n_clim_out, i
                              ))
                        }
                  }
            }
      }

      invisible(TRUE)
}


#' Detect coordinate system from data ranges
#' @keywords internal
.detect_geo <- function(xy) {
      lon_rng <- range(xy[, 1], na.rm = TRUE)
      lat_rng <- range(xy[, 2], na.rm = TRUE)

      if (
            all(is.finite(c(lon_rng, lat_rng))) &&
            lon_rng[1] >= -180 &&
            lon_rng[2] <= 180 &&
            lat_rng[1] >= -90 &&
            lat_rng[2] <= 90
      ) {
            "lonlat"
      } else {
            "projected"
      }
}
