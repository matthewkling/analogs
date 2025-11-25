#' Build Analog Index
#'
#' Pre-builds a reusable lattice index from reference climate data. The index
#' can be queried multiple times with different focal points and parameters,
#' avoiding the need to rebuild the lattice for each query.
#'
#' @inheritParams find_analogs
#'
#' @return An S3 object of class \code{"analog_index"} containing:
#'   \itemize{
#'     \item The compiled lattice index (internal C++ structure)
#'     \item Reference data
#'     \item Metadata: coordinate type, dimensions, ranges, resolution
#'     \item Diagnostics: bin counts and occupancy statistics
#'   }
#'
#' @details
#' The lattice index is built over both geographic and climate dimensions,
#' allowing efficient spatial queries regardless of the constraint values used
#' at query time. For lon/lat coordinates, the index uses ECEF (Earth-Centered
#' Earth-Fixed) space internally for optimal performance.
#'
#' Index resolution (\code{index_res}) controls the granularity of spatial
#' binning. The optimal value depends on your data size and query patterns.
#' Use \code{\link{tune_index_res}} to find the best resolution for your use case,
#' or accept the default of 16 which works well for most applications.
#'
#' @examples
#' \dontrun{
#' # Build index with default settings
#' index <- build_analog_index(climate_data)
#'
#' # Build with explicit resolution
#' index <- build_analog_index(climate_data, index_res = 20)
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
                               index_res = 16) {

      # Validate inputs
      coord_type <- match.arg(coord_type)

      if (!is.numeric(index_res) || length(index_res) != 1L || index_res <= 0) {
            stop("index_res must be a positive integer")
      }
      index_res <- as.integer(index_res)

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
            index_res = index_res
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
                  max_cell_occ = index_cpp$max_cell_occ
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

                  # Check climate bounds
                  for (k in seq_len(index$n_clim)) {
                        clim_col <- query_data[, 2 + k]
                        clim_out <- clim_col < index$clim_mins[k] |
                              clim_col > index$clim_maxs[k]
                        n_clim_out <- sum(clim_out)
                        if (n_clim_out > 0) {
                              warning(sprintf(
                                    "%d query point(s) fall outside index bounds for climate variable %d",
                                    n_clim_out, k
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
