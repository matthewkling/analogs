# ---- Internal Helper Functions ---------------------------------------------

#' Null-coalescing operator
#' @keywords internal
`%||%` <- function(x, y) {
      if (is.null(x)) y else x
}

#' Validate and normalize query parameters
#'
#' Validates mode/k/weight/theta combinations and normalizes values.
#' Returns a list with normalized parameters.
#'
#' @keywords internal
.validate_query_params <- function(mode, k, weight, theta) {

      # Validate mode
      mode <- match.arg(mode, c("knn_clim", "knn_geog", "count", "sum", "mean", "all"))

      # Validate and normalize weight
      # Note: Only inverse_clim and inverse_geog are currently implemented
      # Other weights in the list are placeholders for future implementation
      if (!is.null(weight)) {
            weight <- match.arg(weight, c("uniform", "gaussian_clim", "gaussian_geog",
                                          "gaussian_joint", "inverse_clim", "inverse_geog",
                                          "inverse_joint"))
            if (!weight %in% c("uniform", "inverse_clim", "inverse_geog")) {
                  weight <- NULL
            }
      }

      # Validate mode/k/weight/theta combinations
      if (mode %in% c("knn_clim", "knn_geog")) {
            # kNN modes
            if (is.null(k)) k <- 1L
            k <- as.integer(k)
            if (length(k) != 1L || k <= 0L) {
                  stop("For mode '", mode, "', k must be a positive integer.")
            }
            if (!is.null(weight)) {
                  stop("For mode '", mode, "', weight must be NULL.")
            }
            if (!is.null(theta)) {
                  stop("For mode '", mode, "', theta must be NULL.")
            }

      } else if (mode %in% c("all", "count")) {
            # Filter modes
            if (!is.null(k)) {
                  stop("For mode '", mode, "', k must be NULL.")
            }
            if (!is.null(weight)) {
                  stop("For mode '", mode, "', weight must be NULL.")
            }
            if (!is.null(theta)) {
                  stop("For mode '", mode, "', theta must be NULL.")
            }
            # Leave k as NULL (converted to 0L later in query_analog_index if needed)

      } else if (mode %in% c("sum", "mean")) {
            # Aggregate modes
            if (!is.null(k)) {
                  stop("For mode '", mode, "', k must be NULL.")
            }
            # Leave k as NULL (converted to 0L later in query_analog_index if needed)

            valid_weights <- c("uniform", "inverse_clim", "inverse_geog")
            if (is.null(weight)) weight <- "uniform"
            if (!weight %in% valid_weights) {
                  stop("For mode '", mode, "', weight must be one of: ",
                       paste(valid_weights, collapse = ", "))
            }

            if (identical(weight, "uniform")) {
                  if (!is.null(theta)) {
                        stop("For weight = 'uniform', theta must be NULL.")
                  }
            } else {
                  # inverse_clim or inverse_geog
                  if (!is.null(theta)) {
                        if (!is.numeric(theta) || length(theta) != 1L || theta <= 0) {
                              stop("theta must be a single positive numeric value, or NULL.")
                        }
                  }
            }
      }

      # Return normalized parameters
      list(
            mode = mode,
            k = k,
            weight = weight,
            theta = theta
      )
}

#' Extract coordinates and climate data from input
#' @keywords internal
.select_xy_climate <- function(obj) {
      nm <- colnames(obj)

      # Try to find x,y columns by name
      if (!is.null(nm) && all(c("x", "y") %in% nm)) {
            xy_idx <- match(c("x", "y"), nm)
      } else {
            # Default to first two columns
            xy_idx <- 1:2
      }

      coords <- as.matrix(obj[, xy_idx, drop = FALSE])
      climate <- as.matrix(obj[,
                               setdiff(seq_len(ncol(obj)), xy_idx),
                               drop = FALSE
      ])

      storage.mode(coords) <- "double"
      storage.mode(climate) <- "double"

      if (ncol(coords) != 2L) {
            stop("Coordinate data must have exactly 2 columns (x, y)")
      }
      if (ncol(climate) < 1L) {
            stop(
                  "No climate variable columns found after extracting coordinates"
            )
      }

      cbind(coords, climate)
}

#' Normalize input to standard format
#' @keywords internal
.format_data <- function(r) {
      if (inherits(r, "SpatRaster")) {
            # Convert SpatRaster to data.frame
            if (!requireNamespace("terra", quietly = TRUE)) {
                  stop("Package 'terra' is required for SpatRaster inputs")
            }
            df <- terra::as.data.frame(r, xy = TRUE, na.rm = FALSE)
            .select_xy_climate(df)
      } else if (is.matrix(r) || is.data.frame(r)) {
            .select_xy_climate(r)
      } else {
            stop("Input must be a data.frame, matrix, or SpatRaster")
      }
}
