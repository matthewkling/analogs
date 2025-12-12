# ---- Internal Helper Functions ---------------------------------------------

#' Null-coalescing operator
#' @keywords internal
`%||%` <- function(x, y) {
      if (is.null(x)) y else x
}

#' Validate and normalize query parameters
#'
#' Validates select/stat/k/weight/theta/x_cov/values combinations and normalizes values.
#' Returns a list with normalized parameters.
#'
#' @keywords internal
.validate_query_params <- function(focal = NULL, ref = NULL,
                                   x_cov = NULL, values = NULL,
                                   max_clim, max_geog,
                                   select, k,
                                   stat, weight, theta) {

      # Validate select
      select <- match.arg(select, c("all", "knn_clim", "knn_geog"))

      # Validate max_clim, max_geog
      if(!is.null(max_clim)){
            if(!is.numeric(max_clim) ||
               min(max_clim) <= 0 ||
               !length(max_clim) %in% c(1, ncol(focal)-2)) {
                  stop("max_clim must be a non-negative numeric, with length 1 or",
                       "length matching the number of climate variables.")
            }
      }
      if(!is.null(max_geog)){
            if(!is.numeric(max_geog) ||
               min(max_geog) <= 0 ||
               length(max_geog) != 1) {
                  stop("max_geog must be a non-negative numeric value of length 1.")
            }
      }

      # Normalize stat (NULL becomes "none")
      if (is.null(stat)) {
            stat <- "none"
      } else if (is.character(stat)) {
            # Validate each stat value
            valid_stats <- c("none", "count", "sum_weights", "mean_weights",
                             "sum", "mean", "weighted_sum", "weighted_mean")
            invalid <- setdiff(stat, valid_stats)
            if (length(invalid) > 0) {
                  stop("Invalid stat value(s): ", paste(invalid, collapse = ", "),
                       ". Must be one of: ", paste(valid_stats, collapse = ", "))
            }

            # Check that "none" is not combined with others
            if ("none" %in% stat && length(stat) > 1) {
                  stop('stat = "none" cannot be combined with other aggregations')
            }
      } else {
            stop("stat must be NULL or a character vector")
      }

      # Validate and normalize weight
      if (!is.null(weight)) {
            weight <- match.arg(weight, c("uniform", "gaussian_clim", "gaussian_geog",
                                          "gaussian_joint", "inverse_clim", "inverse_geog",
                                          "inverse_joint"))
      }

      # Validate select/k combination
      if (select %in% c("knn_clim", "knn_geog")) {
            # kNN selection modes require k
            if (is.null(k)) k <- 1L
            k <- as.integer(k)
            if (length(k) != 1L || k <= 0L) {
                  stop("For select = '", select, "', k must be a positive integer.")
            }
      } else {
            # select = "all" doesn't use k
            if (!is.null(k)) {
                  stop("For select = '", select, "', k must be NULL.")
            }
      }

      # Check for value-based stats
      value_stats <- c("sum", "mean", "weighted_sum", "weighted_mean")
      has_value_stat <- any(stat %in% value_stats)

      # If value stats requested, values must be provided
      if (has_value_stat && is.null(values)) {
            requested_value_stats <- intersect(stat, value_stats)
            stop("stat includes ", paste(requested_value_stats, collapse = ", "),
                 " but 'values' parameter is NULL. ",
                 "These stats require 'values' to be provided.")
      }

      # Validate stat/weight/theta combinations
      has_weighted_stat <- any(stat %in% c("sum_weights", "mean_weights",
                                           "weighted_sum", "weighted_mean"))

      if (has_weighted_stat) {
            # Weighted aggregation modes require weight
            valid_weights <- c("uniform", "gaussian_clim", "gaussian_geog",
                               "gaussian_joint", "inverse_clim", "inverse_geog",
                               "inverse_joint")
            if (is.null(weight)) {
                  stop("For stat including weighted aggregations, weight must be specified. ",
                       "Valid options: ", paste(valid_weights, collapse = ", "))
            }
            if (!weight %in% valid_weights) {
                  stop("For stat including weighted aggregations, weight must be one of: ",
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

      } else {
            # Non-weighted aggregations (none, count, sum, mean)
            if (!is.null(weight)) {
                  stop("For stat = ", paste(stat, collapse = ", "), ", weight must be NULL.")
            }
            if (!is.null(theta)) {
                  stop("For stat = ", paste(stat, collapse = ", "), ", theta must be NULL.")
            }
      }

      # Validate and format x_cov if provided
      x_cov_mat <- NULL
      if (!is.null(x_cov)) {
            if (is.null(focal)) {
                  stop("Internal error: focal required for x_cov validation")
            }
            x_cov_mat <- .validate_and_format_x_cov(x_cov, focal)
      }

      # Validate and format values if provided
      values_mat <- NULL
      values_names <- NULL
      if (!is.null(values)) {
            if (is.null(ref)) {
                  stop("Internal error: ref required for values validation")
            }

            result <- .validate_and_format_values(values, ref)
            values_mat <- result$matrix
            values_names <- result$names
      }

      # Return normalized parameters
      list(
            select = select,
            stat = stat,
            k = k,
            weight = weight,
            theta = theta,
            x_cov = x_cov_mat,
            values = values_mat,
            values_names = values_names
      )
}

#' Validate and format values parameter
#' @keywords internal
.validate_and_format_values <- function(values, ref) {

      n_ref <- nrow(ref)

      # Handle SpatRaster input
      if (inherits(values, "SpatRaster")) {
            if (!requireNamespace("terra", quietly = TRUE)) {
                  stop("Package 'terra' is required for SpatRaster values", call. = FALSE)
            }

            # Convert to data.frame (keeps all cells including NA)
            df <- terra::as.data.frame(values, xy = FALSE, na.rm = FALSE)

            # Check dimensions match ref
            if (nrow(df) != n_ref) {
                  stop(sprintf(
                        "values SpatRaster has %d cells but reference data has %d rows. They must match.",
                        nrow(df), n_ref
                  ))
            }

            values_names <- names(df)
            values <- as.matrix(df)

      } else if (is.vector(values)) {
            # Convert vector to single-column matrix
            values <- matrix(values, ncol = 1)
            values_names <- "value_1"

      } else if (is.data.frame(values)) {
            values_names <- colnames(values)
            values <- as.matrix(values)

      } else if (is.matrix(values)) {
            values_names <- colnames(values)

      } else {
            stop("values must be a vector, matrix, data.frame, or SpatRaster")
      }

      # Validate dimensions
      if (nrow(values) != n_ref) {
            stop(sprintf(
                  "values must have same number of rows as reference data (%d), but has %d rows",
                  n_ref, nrow(values)
            ))
      }

      # Check for numeric
      if (!is.numeric(values)) {
            stop("values must be numeric")
      }

      # Generate names if missing
      if (is.null(values_names)) {
            n_vars <- ncol(values)
            values_names <- if (n_vars == 1) {
                  "value_1"
            } else {
                  paste0("value_", seq_len(n_vars))
            }
      }

      # Ensure storage mode is double
      storage.mode(values) <- "double"

      list(
            matrix = values,
            names = values_names
      )
}

#' Validate and format x_cov parameter
#' @keywords internal
.validate_and_format_x_cov <- function(x_cov, focal) {

      # focal is already formatted matrix with coords + climate
      n_focal <- nrow(focal)
      n_clim <- ncol(focal) - 2

      # Expected number of covariance components
      n_cov_components <- n_clim * (n_clim + 1) / 2

      # Convert to matrix if needed
      if (inherits(x_cov, "SpatRaster")) {
            if (!requireNamespace("terra", quietly = TRUE)) {
                  stop("Package 'terra' is required for SpatRaster x_cov", call. = FALSE)
            }
            # Convert to data.frame (keeps all cells including NA) to match .format_data behavior
            df <- terra::as.data.frame(x_cov, xy = FALSE, na.rm = FALSE)
            x_cov <- as.matrix(df)
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
            df <- .select_xy_climate(df)
            attr(df, "template") <- setNames(terra::setValues(r[[1]], NA), "raster")
            return(df)
      } else if (is.matrix(r) || is.data.frame(r)) {
            df <- .select_xy_climate(r)
            attr(df, "template") <- NULL
            return(df)
      } else {
            stop("Input must be a data.frame, matrix, or SpatRaster")
      }
}



.format_output <- function(out, x, stat, select, k, weight, theta, x_cov_mat){

      if(! requireNamespace("terra", quietly = TRUE) || # terra not available
         is.null(attr(x, "template")) || # x wasn't a raster
         (any(stat == "none") && (k != 1 || select == "all")) || # query not compatible
         nrow(out) != terra::ncell(attr(x, "template")) # data not compatible (e.g. called from tune_index_res)
      ){
            # raster not relevant -- return data.frame as is

      } else {
            # rasterize with template
            att <- attributes(out) # capture cpp attributes
            vars <- setdiff(names(out), c("x", "y", "index"))
            out <- terra::rast(
                  lapply(vars, function(v){
                        setNames(terra::setValues(attr(x, "template"), out[[v]]), v)
                  })
            )
            terra::varnames(out) <- vars  # set varnames to match layer names

            # add attributes
            attributes(out) <- append(attributes(out), att[setdiff(names(att), names(attributes(out)))])
      }

      attr(out, "select")    <- select
      attr(out, "stat")      <- stat
      attr(out, "k")         <- k
      attr(out, "weight")    <- weight
      attr(out, "theta")     <- theta
      attr(out, "x_cov")     <- !is.null(x_cov_mat)

      return(out)
}
