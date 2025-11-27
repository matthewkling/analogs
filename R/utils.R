# ---- Internal Helper Functions ---------------------------------------------

#' Null-coalescing operator
#' @keywords internal
`%||%` <- function(x, y) {
      if (is.null(x)) y else x
}

#' Validate and normalize query parameters
#'
#' Validates mode/k/weight/theta/x_cov combinations and normalizes values.
#' Returns a list with normalized parameters.
#'
#' @keywords internal
.validate_query_params <- function(mode, k, weight, theta, x_cov = NULL, focal_mm = NULL) {

      # Validate mode
      mode <- match.arg(mode, c("knn_clim", "knn_geog", "count", "sum", "mean", "all"))

      # Validate and normalize weight
      if (!is.null(weight)) {
            weight <- match.arg(weight, c("uniform", "gaussian_clim", "gaussian_geog",
                                          "gaussian_joint", "inverse_clim", "inverse_geog",
                                          "inverse_joint"))
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

            valid_weights <- c("uniform", "gaussian_clim", "gaussian_geog",
                               "gaussian_joint", "inverse_clim", "inverse_geog",
                               "inverse_joint")
            if (is.null(weight)) weight <- "uniform"
            if (!weight %in% valid_weights) {
                  stop("For mode '", mode, "', weight must be one of: ",
                       paste(valid_weights, collapse = ", "))
            }

            # Validate theta based on weight type
            if (identical(weight, "uniform")) {
                  if (!is.null(theta)) {
                        stop("For weight = 'uniform', theta must be NULL.")
                  }
            } else if (weight %in% c("gaussian_joint", "inverse_joint")) {
                  # Joint weights require 2-element theta
                  if (is.null(theta)) {
                        stop("For weight = '", weight, "', theta must be a numeric vector of length 2.")
                  }
                  if (!is.numeric(theta) || length(theta) != 2L) {
                        stop("For weight = '", weight, "', theta must be a numeric vector of length 2: ",
                             "c(theta_clim, theta_geog)")
                  }
                  if (any(theta <= 0)) {
                        stop("For weight = '", weight, "', both theta values must be positive.")
                  }
            } else {
                  # Single-dimension weights (gaussian_clim, gaussian_geog, inverse_clim, inverse_geog)
                  if (!is.null(theta)) {
                        if (!is.numeric(theta) || length(theta) != 1L || theta <= 0) {
                              stop("For weight = '", weight, "', theta must be a single positive numeric value, or NULL.")
                        }
                  }
            }
      }

      # Validate and format x_cov if provided
      x_cov_mat <- NULL
      if (!is.null(x_cov)) {
            if (is.null(focal_mm)) {
                  stop("Internal error: focal_mm required for x_cov validation")
            }
            x_cov_mat <- .validate_and_format_x_cov(x_cov, focal_mm)
      }

      # Return normalized parameters
      list(
            mode = mode,
            k = k,
            weight = weight,
            theta = theta,
            x_cov = x_cov_mat
      )
}

#' Validate and format x_cov parameter
#' @keywords internal
.validate_and_format_x_cov <- function(x_cov, focal_mm) {

      # focal_mm is already formatted matrix with coords + climate
      n_focal <- nrow(focal_mm)
      n_clim <- ncol(focal_mm) - 2

      # Expected number of covariance components
      n_cov_components <- n_clim * (n_clim + 1) / 2

      # Convert to matrix if needed
      if (inherits(x_cov, "SpatRaster")) {
            if (!requireNamespace("terra", quietly = TRUE)) {
                  stop("Package 'terra' is required for SpatRaster x_cov", call. = FALSE)
            }
            x_cov <- terra::as.matrix(x_cov, wide = TRUE)
      } else if (is.data.frame(x_cov)) {
            x_cov <- as.matrix(x_cov)
      } else if (!is.matrix(x_cov)) {
            stop("x_cov must be a matrix, data.frame, or SpatRaster")
      }

      # Validate dimensions
      if (nrow(x_cov) != n_focal) {
            stop(sprintf(
                  "x_cov must have same number of rows as focal data (%d), but has %d rows",
                  n_focal, nrow(x_cov)
            ))
      }

      if (ncol(x_cov) != n_cov_components) {
            stop(sprintf(
                  "For %d climate variables, x_cov must have %d columns (n*(n+1)/2), but has %d",
                  n_clim, n_cov_components, ncol(x_cov)
            ))
      }

      # Check for non-finite values
      if (any(!is.finite(x_cov))) {
            stop("x_cov contains non-finite values")
      }

      # Basic positive-definiteness check on first focal's covariance matrix
      # (full check is expensive, just do sanity check)
      test_cov <- .reconstruct_cov_matrix(x_cov[1, ], n_clim)
      test_eig <- eigen(test_cov, symmetric = TRUE, only.values = TRUE)$values

      if (any(test_eig <= 0)) {
            warning("First focal's covariance matrix is not positive definite. ",
                    "This may cause issues. Check your covariance matrices.")
      }

      # Ensure storage mode is double
      storage.mode(x_cov) <- "double"

      return(x_cov)
}


#' Reconstruct symmetric covariance matrix from lower triangle
#' @keywords internal
.reconstruct_cov_matrix <- function(cov_vec, n_clim) {
      cov_mat <- matrix(0, n_clim, n_clim)

      # Fill diagonal (variances)
      for (i in seq_len(n_clim)) {
            cov_mat[i, i] <- cov_vec[i]
      }

      # Fill off-diagonals (covariances)
      if (n_clim > 1) {
            idx <- n_clim + 1
            for (i in seq_len(n_clim - 1)) {
                  for (j in (i + 1):n_clim) {
                        cov_mat[i, j] <- cov_vec[idx]
                        cov_mat[j, i] <- cov_vec[idx]  # Symmetric
                        idx <- idx + 1
                  }
            }
      }

      return(cov_mat)
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
