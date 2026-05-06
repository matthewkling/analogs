# Internal helper functions

# Validation helpers ------------------------------------------------

# Soft/hard caps on the number of distinct classes per y variable when
# stat = "tabulate". These thresholds gate misuse (e.g., passing a continuous
# variable as a categorical), since tabulate's output width grows linearly
# in K.
.tabulate_K_warn  <- 100L
.tabulate_K_error <- 1000L

# Validate and normalize query parameters
#
# Validates select/stat/k/kernel/theta/x_cov/y/covariates/lambda/se
# combinations and normalizes values. Returns a list with normalized parameters.
.validate_query_params <- function(focal = NULL, ref = NULL,
                                   x_cov = NULL, y = NULL,
                                   covariates = NULL,
                                   max_clim, max_geog,
                                   select, k,
                                   stat, kernel, theta,
                                   lambda = 0,
                                   se = "none") {

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
                             "sum", "mean", "weighted_sum", "weighted_mean",
                             "ess", "regression", "tabulate")
            invalid <- setdiff(stat, valid_stats)
            if (length(invalid) > 0) {
                  stop("Invalid stat value(s): ", paste(invalid, collapse = ", "),
                       ". Must be one of: ", paste(valid_stats, collapse = ", "))
            }

            # Check that "none" is not combined with others
            if ("none" %in% stat && length(stat) > 1) {
                  stop('stat = "none" cannot be combined with other aggregations')
            }

            # tabulate is mutually exclusive with continuous-y stats and
            # regression (different y semantics). It can be combined with
            # regular non-y stats (count, sum_weights, mean_weights, ess).
            if ("tabulate" %in% stat) {
                  incompatible <- intersect(stat,
                                            c("sum", "mean", "weighted_sum",
                                              "weighted_mean", "regression"))
                  if (length(incompatible) > 0L) {
                        stop("stat = 'tabulate' cannot be combined with ",
                             paste(sprintf("'%s'", incompatible), collapse = ", "),
                             ". 'tabulate' treats `y` as categorical, while ",
                             "the other listed stats treat `y` as continuous. ",
                             "Run separate queries if you need both.",
                             call. = FALSE)
                  }
            }
      } else {
            stop("stat must be NULL or a character vector")
      }

      # Validate and normalize kernel
      if (!is.null(kernel)) {
            kernel <- match.arg(kernel, c("uniform", "gaussian_clim", "gaussian_geog",
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

      # If y value stats requested, y must be provided
      if (has_value_stat && is.null(y)) {
            requested_value_stats <- intersect(stat, value_stats)
            stop("stat includes ", paste(requested_value_stats, collapse = ", "),
                 " but 'y' parameter is NULL. ",
                 "These stats require 'y' to be provided.")
      }

      # tabulate also requires y
      if ("tabulate" %in% stat && is.null(y)) {
            stop("stat includes 'tabulate' but 'y' parameter is NULL. ",
                 "Tabulate requires 'y' to be provided.",
                 call. = FALSE)
      }

      # Validate regression stat requirements
      if ("regression" %in% stat) {
            if (is.null(y)) {
                  stop("stat includes 'regression' but 'y' parameter is NULL. ",
                       "Regression requires 'y' to be provided.")
            }
            if (is.null(covariates)) {
                  stop("stat includes 'regression' but 'covariates' parameter is NULL. ",
                       "Regression requires 'covariates' to be provided.")
            }
            if (is.null(kernel)) {
                  stop("stat includes 'regression' but 'kernel' is NULL. ",
                       "Regression requires a kernel weighting function. ",
                       "Valid options: uniform, gaussian_clim, gaussian_geog, ",
                       "gaussian_joint, inverse_clim, inverse_geog, inverse_joint")
            }
      }

      # Validate lambda
      if (!is.numeric(lambda) || length(lambda) != 1L || lambda < 0) {
            stop("lambda must be a single non-negative numeric value.")
      }

      # Validate se
      # Accept a single string; normalize to one of "none"/"ess"/"design".
      se <- match.arg(se, c("none", "ess", "design"))

      # If se != "none" but no requested stat supports SE, warn.
      # SE-supporting stats: weighted_mean, regression.
      if (!identical(se, "none")) {
            se_stats <- c("weighted_mean", "regression")
            if (!any(stat %in% se_stats)) {
                  warning("`se = \"", se, "\"` was requested but no requested stat supports SE ",
                          "(SE is currently defined for: ",
                          paste(se_stats, collapse = ", "), "). ",
                          "No SE columns will be produced.")
            }
      }

      # Validate stat/kernel/theta combinations.
      # tabulate sums per-class WEIGHTS, so it's a weighted aggregation.
      has_weighted_stat <- any(stat %in% c("sum_weights", "mean_weights",
                                           "weighted_sum", "weighted_mean",
                                           "ess", "regression", "tabulate"))

      if (has_weighted_stat) {
            # Weighted aggregation modes require kernel
            valid_kernels <- c("uniform", "gaussian_clim", "gaussian_geog",
                               "gaussian_joint", "inverse_clim", "inverse_geog",
                               "inverse_joint")
            if (is.null(kernel)) {
                  stop("For stat including weighted aggregations, kernel must be specified. ",
                       "Valid options: ", paste(valid_kernels, collapse = ", "))
            }
            if (!kernel %in% valid_kernels) {
                  stop("For stat including weighted aggregations, kernel must be one of: ",
                       paste(valid_kernels, collapse = ", "))
            }

            # Validate theta based on kernel type
            if (identical(kernel, "uniform")) {
                  if (!is.null(theta)) {
                        stop("For kernel = 'uniform', theta must be NULL.")
                  }
            } else if (kernel %in% c("gaussian_joint", "inverse_joint")) {
                  # Joint kernels require 2-element theta
                  if (is.null(theta)) {
                        stop("For kernel = '", kernel, "', theta must be a numeric vector of length 2.")
                  }
                  if (!is.numeric(theta) || length(theta) != 2L) {
                        stop("For kernel = '", kernel, "', theta must be a numeric vector of length 2: ",
                             "c(theta_clim, theta_geog)")
                  }
                  if (any(theta <= 0)) {
                        stop("For kernel = '", kernel, "', both theta values must be positive.")
                  }
            } else {
                  # Single-dimension kernels (gaussian_clim, gaussian_geog, inverse_clim, inverse_geog)
                  if (!is.null(theta)) {
                        if (!is.numeric(theta) || length(theta) != 1L || theta <= 0) {
                              stop("For kernel = '", kernel, "', theta must be a single positive numeric value, or NULL.")
                        }
                  }
            }

      } else {
            # Non-weighted aggregations (none, count, sum, mean)
            if (!is.null(kernel)) {
                  stop("For stat = ", paste(stat, collapse = ", "), ", kernel must be NULL.")
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

      # Validate and format y values if provided.
      # Branch on whether tabulate was requested: tabulate factors each y
      # column/layer and passes 1-based integer codes through the same
      # numeric-matrix path used by continuous stats. The factor levels are
      # captured per variable for downstream column naming.
      values_mat <- NULL
      values_names <- NULL
      values_levels <- NULL   # list of character vectors, one per y column;
      # NULL when y is continuous
      if (!is.null(y)) {
            if (is.null(ref)) {
                  stop("Internal error: ref required for y validation")
            }

            if ("tabulate" %in% stat) {
                  result <- .validate_and_format_y_categorical(y, ref)
                  values_mat    <- result$matrix
                  values_names  <- result$names
                  values_levels <- result$levels
            } else {
                  # Existing continuous path. If y is a factor here, we error
                  # with a hint pointing the user at stat = "tabulate".
                  result <- .validate_and_format_values(y, ref)
                  values_mat   <- result$matrix
                  values_names <- result$names
            }
      }

      # Validate and format covariates if provided
      covariates_mat <- NULL
      covariates_names <- NULL
      if (!is.null(covariates)) {
            if (is.null(ref)) {
                  stop("Internal error: ref required for covariates validation")
            }

            result <- .validate_and_format_covariates(covariates, ref)
            covariates_mat <- result$matrix
            covariates_names <- result$names
      }

      # Return normalized parameters
      list(
            select = select,
            stat = stat,
            k = k,
            kernel = kernel,
            theta = theta,
            lambda = lambda,
            se = se,
            x_cov = x_cov_mat,
            y = values_mat,
            values_names = values_names,
            values_levels = values_levels,
            covariates = covariates_mat,
            covariates_names = covariates_names
      )
}


# Validate and format covariates parameter
.validate_and_format_covariates <- function(covariates, ref) {

      n_ref <- nrow(ref)

      # Handle SpatRaster input
      if (inherits(covariates, "SpatRaster")) {
            if (!requireNamespace("terra", quietly = TRUE)) {
                  stop("Package 'terra' is required for SpatRaster covariates", call. = FALSE)
            }

            df <- terra::as.data.frame(covariates, xy = FALSE, na.rm = FALSE)

            if (nrow(df) != n_ref) {
                  stop(sprintf(
                        "covariates SpatRaster has %d cells but pool has %d rows. They must match.",
                        nrow(df), n_ref
                  ))
            }

            cov_names <- names(df)
            mat <- as.matrix(df)

      } else if (is.data.frame(covariates)) {
            if (nrow(covariates) != n_ref) {
                  stop(sprintf(
                        "covariates has %d rows but pool has %d rows. They must match.",
                        nrow(covariates), n_ref
                  ))
            }

            # Check all columns are numeric
            non_numeric <- names(covariates)[!vapply(covariates, is.numeric, logical(1))]
            if (length(non_numeric) > 0) {
                  stop("All covariates columns must be numeric. Non-numeric: ",
                       paste(non_numeric, collapse = ", "))
            }

            cov_names <- names(covariates)
            mat <- as.matrix(covariates)

      } else if (is.matrix(covariates)) {
            if (nrow(covariates) != n_ref) {
                  stop(sprintf(
                        "covariates has %d rows but pool has %d rows. They must match.",
                        nrow(covariates), n_ref
                  ))
            }
            if (!is.numeric(covariates)) {
                  stop("covariates matrix must be numeric.")
            }

            cov_names <- colnames(covariates)
            mat <- covariates

      } else if (is.numeric(covariates) && is.null(dim(covariates))) {
            # Numeric vector — single covariate
            if (length(covariates) != n_ref) {
                  stop(sprintf(
                        "covariates has length %d but pool has %d rows. They must match.",
                        length(covariates), n_ref
                  ))
            }

            cov_names <- "covariate"
            mat <- matrix(covariates, ncol = 1)
      } else {
            stop("covariates must be a numeric vector, matrix, data.frame, or SpatRaster.")
      }

      # Ensure names exist
      if (is.null(cov_names)) {
            cov_names <- paste0("cov", seq_len(ncol(mat)))
      }

      list(matrix = mat, names = cov_names)
}

# Validate and format y parameter (continuous path)
.validate_and_format_values <- function(y, ref) {

      n_ref <- nrow(ref)

      # Catch factor input early with a clear hint.
      if (is.factor(y) ||
          (is.data.frame(y) && any(vapply(y, is.factor, logical(1))))) {
            stop("`y` is a factor, but the requested stat is not 'tabulate'. ",
                 "Pass `stat = \"tabulate\"` for categorical aggregation, ",
                 "or convert `y` to numeric for continuous stats.",
                 call. = FALSE)
      }

      # Handle SpatRaster input
      if (inherits(y, "SpatRaster")) {
            if (!requireNamespace("terra", quietly = TRUE)) {
                  stop("Package 'terra' is required for SpatRaster 'y'", call. = FALSE)
            }

            # Convert to data.frame (keeps all cells including NA)
            df <- terra::as.data.frame(y, xy = FALSE, na.rm = FALSE)

            # Check dimensions match ref
            if (nrow(df) != n_ref) {
                  stop(sprintf(
                        "y SpatRaster has %d cells but pool has %d rows. They must match.",
                        nrow(df), n_ref
                  ))
            }

            values_names <- names(df)
            y <- as.matrix(df)

      } else if (is.vector(y)) {
            # Convert vector to single-column matrix
            y <- matrix(y, ncol = 1)
            values_names <- "value_1"

      } else if (is.data.frame(y)) {
            values_names <- colnames(y)
            y <- as.matrix(y)

      } else if (is.matrix(y)) {
            values_names <- colnames(y)

      } else {
            stop("'y' must be a vector, matrix, data.frame, or SpatRaster")
      }

      # Validate dimensions
      if (nrow(y) != n_ref) {
            stop(sprintf(
                  "'y' must have same number of rows as pool (%d), but has %d rows",
                  n_ref, nrow(y)
            ))
      }

      # Check for numeric
      if (!is.numeric(y)) {
            stop("'y' must be numeric")
      }

      # Generate names if missing
      if (is.null(values_names)) {
            n_vars <- ncol(y)
            values_names <- if (n_vars == 1) {
                  "y1"
            } else {
                  paste0("y", seq_len(n_vars))
            }
      }

      # Ensure storage mode is double
      storage.mode(y) <- "double"

      list(
            matrix = y,
            names = values_names
      )
}

# Validate and format y parameter (categorical path, for stat = "tabulate")
#
# Accepts vector / matrix / data.frame / SpatRaster shapes (same as the
# continuous path). Each column / layer is independently coerced to a factor;
# the resulting 1-based integer codes are written into a numeric matrix to be
# passed verbatim through the existing values_sexp pathway. NA stays NA.
#
# Returns:
#   matrix  - n_ref x n_vars numeric matrix of 1-based codes (NA = missing)
#   names   - character vector of variable names
#   levels  - list (length n_vars) of character vectors giving each variable's
#             factor levels (i.e., the column/layer names per output bin)
.validate_and_format_y_categorical <- function(y, ref) {

      n_ref <- nrow(ref)

      # First, normalize y into a list of per-variable vectors plus a names
      # vector, so we can factor each column independently regardless of the
      # original container type.
      cols <- NULL
      values_names <- NULL

      if (inherits(y, "SpatRaster")) {
            if (!requireNamespace("terra", quietly = TRUE)) {
                  stop("Package 'terra' is required for SpatRaster 'y'", call. = FALSE)
            }
            df <- terra::as.data.frame(y, xy = FALSE, na.rm = FALSE)
            if (nrow(df) != n_ref) {
                  stop(sprintf(
                        "y SpatRaster has %d cells but pool has %d rows. They must match.",
                        nrow(df), n_ref
                  ))
            }
            values_names <- names(df)
            cols <- as.list(df)

      } else if (is.factor(y)) {
            # Bare factor vector
            if (length(y) != n_ref) {
                  stop(sprintf(
                        "'y' length is %d but pool has %d rows. They must match.",
                        length(y), n_ref
                  ))
            }
            cols <- list(y)
            values_names <- "y1"

      } else if (is.data.frame(y)) {
            if (nrow(y) != n_ref) {
                  stop(sprintf(
                        "'y' has %d rows but pool has %d rows. They must match.",
                        nrow(y), n_ref
                  ))
            }
            values_names <- names(y)
            cols <- as.list(y)

      } else if (is.matrix(y)) {
            if (nrow(y) != n_ref) {
                  stop(sprintf(
                        "'y' has %d rows but pool has %d rows. They must match.",
                        nrow(y), n_ref
                  ))
            }
            values_names <- colnames(y)
            cols <- lapply(seq_len(ncol(y)), function(j) y[, j])

      } else if (is.atomic(y) && is.null(dim(y))) {
            # Plain vector (numeric, integer, character, logical)
            if (length(y) != n_ref) {
                  stop(sprintf(
                        "'y' length is %d but pool has %d rows. They must match.",
                        length(y), n_ref
                  ))
            }
            cols <- list(y)
            values_names <- "y1"

      } else {
            stop("'y' must be a vector, factor, matrix, data.frame, or SpatRaster",
                 call. = FALSE)
      }

      n_vars <- length(cols)

      # Generate variable names if missing
      if (is.null(values_names) || any(!nzchar(values_names))) {
            values_names <- if (n_vars == 1L) "y1" else paste0("y", seq_len(n_vars))
      }

      # Factor each column, capture levels, build code matrix.
      # Always factor in R (per design): the user's stat = "tabulate"
      # declaration is what determines categorical behavior; we don't
      # second-guess input type. Sparse integer inputs are densely re-coded
      # by factor() so the C++ accumulator allocates only what's needed.
      levels_list <- vector("list", n_vars)
      code_mat <- matrix(NA_real_, nrow = n_ref, ncol = n_vars)

      total_K_warn  <- FALSE

      for (v in seq_len(n_vars)) {
            col <- cols[[v]]

            # Reject types that can't sensibly be factors (lists, dates with
            # tz semantics, raw, complex, etc.). Allow logical (factored as
            # FALSE/TRUE).
            if (is.list(col) ||
                inherits(col, "POSIXt") ||
                inherits(col, "Date") ||
                is.complex(col) ||
                is.raw(col)) {
                  stop("`y` column ", values_names[v],
                       " has type that cannot be coerced to a factor for ",
                       "stat = 'tabulate' (got class ",
                       paste(class(col), collapse = "/"), "). Convert to a ",
                       "character or factor first.", call. = FALSE)
            }

            f <- if (is.factor(col)) {
                  # Drop unused levels: even if the user pre-factored with a
                  # broader levels set, we only allocate for levels that
                  # actually appear (matches the "factor in R, dense in C++"
                  # contract).
                  droplevels(col)
            } else {
                  factor(col)
            }

            lev <- levels(f)
            K_v <- length(lev)

            if (K_v == 0L) {
                  stop("`y` column ", values_names[v],
                       " has no non-NA values. Cannot tabulate an all-NA ",
                       "variable.", call. = FALSE)
            }

            if (K_v > .tabulate_K_error) {
                  stop("`y` column ", values_names[v], " has ", K_v,
                       " distinct levels, exceeding the limit of ",
                       .tabulate_K_error, ". Tabulate is intended for ",
                       "categorical variables; this looks more like a ",
                       "continuous variable. If you really need this many ",
                       "levels, raise `analogs:::.tabulate_K_error`.",
                       call. = FALSE)
            }
            if (K_v > .tabulate_K_warn) {
                  total_K_warn <- TRUE
            }

            levels_list[[v]] <- lev
            # Integer codes (1-based, NA for missing) -> double for the
            # existing values_sexp matrix pathway. Round-trip is exact.
            code_mat[, v] <- as.numeric(as.integer(f))
      }

      if (isTRUE(total_K_warn)) {
            warning("At least one `y` column has more than ", .tabulate_K_warn,
                    " distinct levels. Output width grows linearly in the ",
                    "number of classes; ensure `y` is meant to be categorical.",
                    call. = FALSE)
      }

      list(
            matrix = code_mat,
            names  = values_names,
            levels = levels_list
      )
}

# Validate and format x_cov parameter
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

# Validate exclude_self parameter and its compatibility with other args
#
# Called from analog_search (and anywhere else that surfaces exclude_self).
# Enforces identical(x, pool), disallows pre-built indices, and disallows
# downsampling. Also disallows progress (chunking is incompatible with the
# simple j==i self-exclusion check).
.validate_exclude_self <- function(exclude_self, x, pool, downsample, progress) {
      if (!is.logical(exclude_self) || length(exclude_self) != 1L || is.na(exclude_self)) {
            stop("`exclude_self` must be TRUE or FALSE.", call. = FALSE)
      }
      if (!exclude_self) return(invisible(TRUE))

      if (is_analog_index(pool)) {
            stop(
                  "`exclude_self = TRUE` is not supported when `pool` is a pre-built ",
                  "analog_index. Pass the raw pool data instead.",
                  call. = FALSE
            )
      }

      if (!identical(x, pool)) {
            stop(
                  "`exclude_self = TRUE` requires `x` and `pool` to be the same object ",
                  "(checked via identical()). See `analog_cv()` for standard ",
                  "cross-validation workflows.",
                  call. = FALSE
            )
      }

      if (!is.null(downsample) && !isTRUE(all.equal(downsample, 1.0))) {
            stop(
                  "`exclude_self = TRUE` is not compatible with `downsample != 1`. ",
                  "Self-exclusion semantics are ill-defined under downsampling.",
                  call. = FALSE
            )
      }

      if (isTRUE(progress)) {
            stop(
                  "`exclude_self = TRUE` is not compatible with `progress = TRUE`. ",
                  "Disable progress or run without self-exclusion.",
                  call. = FALSE
            )
      }

      invisible(TRUE)
}


# Other helpers ------------------------------------------------

# Null-coalescing operator
`%||%` <- function(x, y) {
      if (is.null(x)) y else x
}

# Reconstruct symmetric covariance matrix from lower triangle
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

# Extract coordinates and climate data from input
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

# Normalize input to standard format
.format_data <- function(r) {
      if (inherits(r, "SpatRaster")) {
            # Convert SpatRaster to data.frame
            if (!requireNamespace("terra", quietly = TRUE)) {
                  stop("Package 'terra' is required for SpatRaster inputs")
            }
            df <- terra::as.data.frame(r, xy = TRUE, na.rm = FALSE)
            df <- .select_xy_climate(df)
            attr(df, "template") <- stats::setNames(terra::setValues(r[[1]], NA), "raster")
            return(df)
      } else if (is.matrix(r) || is.data.frame(r)) {
            df <- .select_xy_climate(r)
            attr(df, "template") <- NULL
            return(df)
      } else {
            stop("Input must be a data.frame, matrix, or SpatRaster")
      }
}

# Compute cell-area weights for a SpatRaster pool.
#
# Returns a length-ncell numeric vector of weights, normalized so the mean
# over finite (non-NA cellSize) values is 1.0. The first layer's NA mask is
# *not* applied; cells with NA climate values still get a finite area weight
# (they get filtered out via the climate-NA path during querying anyway).
#
# Used by build_analog_index() when cell_area_weight is "auto" (and pool is
# a raster) or TRUE.
.compute_cell_area_weights <- function(pool_raster) {
      if (!requireNamespace("terra", quietly = TRUE)) {
            stop("Package 'terra' is required to compute cell-area weights",
                 call. = FALSE)
      }
      if (!inherits(pool_raster, "SpatRaster")) {
            stop("Internal error: .compute_cell_area_weights requires a SpatRaster",
                 call. = FALSE)
      }

      # cellSize returns geodesic km^2 for lonlat; planar projection-units^2
      # otherwise. The unit here doesn't matter because we normalize.
      area <- terra::values(terra::cellSize(pool_raster, mask = FALSE,
                                            unit = "km"))
      area <- as.numeric(area)

      # Normalize to mean 1 over finite values
      finite <- is.finite(area)
      if (!any(finite)) {
            stop("Internal error: cellSize returned no finite values.",
                 call. = FALSE)
      }
      mean_area <- mean(area[finite])
      if (mean_area <= 0 || !is.finite(mean_area)) {
            stop("Internal error: non-positive mean cell area.", call. = FALSE)
      }
      area / mean_area
}

# Validate and format the user-supplied per-pool-point `weight` argument.
#
# Accepts numeric vector, single-column matrix/data.frame, or single-layer
# SpatRaster (mirroring the `y` and `covariates` input patterns). Returns a
# length-n_ref numeric vector with NAs converted to 0 (NA pool points are
# excluded from any aggregation).
#
# Returns NULL if `weight` is NULL.
.validate_and_format_weight <- function(weight, ref) {
      if (is.null(weight)) return(NULL)

      n_ref <- nrow(ref)

      if (inherits(weight, "SpatRaster")) {
            if (!requireNamespace("terra", quietly = TRUE)) {
                  stop("Package 'terra' is required for SpatRaster weight",
                       call. = FALSE)
            }
            if (terra::nlyr(weight) != 1L) {
                  stop("`weight` SpatRaster must have exactly one layer.",
                       call. = FALSE)
            }
            df <- terra::as.data.frame(weight, xy = FALSE, na.rm = FALSE)
            if (nrow(df) != n_ref) {
                  stop(sprintf(
                        "weight SpatRaster has %d cells but pool has %d rows. They must match.",
                        nrow(df), n_ref
                  ), call. = FALSE)
            }
            w <- as.numeric(df[[1L]])

      } else if (is.data.frame(weight)) {
            if (ncol(weight) != 1L) {
                  stop("`weight` data.frame must have exactly one column.",
                       call. = FALSE)
            }
            if (nrow(weight) != n_ref) {
                  stop(sprintf(
                        "weight has %d rows but pool has %d rows. They must match.",
                        nrow(weight), n_ref
                  ), call. = FALSE)
            }
            w <- as.numeric(weight[[1L]])

      } else if (is.matrix(weight)) {
            if (ncol(weight) != 1L) {
                  stop("`weight` matrix must have exactly one column.",
                       call. = FALSE)
            }
            if (nrow(weight) != n_ref) {
                  stop(sprintf(
                        "weight has %d rows but pool has %d rows. They must match.",
                        nrow(weight), n_ref
                  ), call. = FALSE)
            }
            w <- as.numeric(weight[, 1L])

      } else if (is.numeric(weight) && is.null(dim(weight))) {
            if (length(weight) != n_ref) {
                  stop(sprintf(
                        "weight has length %d but pool has %d rows. They must match.",
                        length(weight), n_ref
                  ), call. = FALSE)
            }
            w <- as.numeric(weight)

      } else {
            stop("`weight` must be a numeric vector, single-column matrix or ",
                 "data.frame, or single-layer SpatRaster.", call. = FALSE)
      }

      # Validate values: NAs allowed (converted to 0), other non-finite or
      # negative values are errors.
      bad <- !is.na(w) & (!is.finite(w) | w < 0)
      if (any(bad)) {
            stop("`weight` must contain only non-negative finite values (NAs allowed).",
                 call. = FALSE)
      }
      w[is.na(w)] <- 0

      w
}

.format_output <- function(out, x, stat, select, k, kernel, theta, x_cov_mat,
                           lambda, se, exclude_self, downsample_actual,
                           cell_area_weight_applied = FALSE,
                           weight_provided = FALSE){

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
                        stats::setNames(terra::setValues(attr(x, "template"), out[[v]]), v)
                  })
            )
            terra::varnames(out) <- vars  # set varnames to match layer names

            # add attributes
            attributes(out) <- append(attributes(out), att[setdiff(names(att), names(attributes(out)))])
      }

      attr(out, "select")            <- select
      attr(out, "stat")               <- stat
      attr(out, "k")                  <- k
      attr(out, "kernel")             <- kernel
      attr(out, "theta")              <- theta
      attr(out, "lambda")             <- lambda
      attr(out, "se")                 <- se
      attr(out, "exclude_self")       <- exclude_self
      attr(out, "x_cov_provided")     <- !is.null(x_cov_mat)
      attr(out, "downsample_actual")  <- downsample_actual
      attr(out, "cell_area_weight")   <- cell_area_weight_applied
      attr(out, "weight_provided")    <- weight_provided

      return(out)
}

#' Predict from regression coefficients
#'
#' Shared helper that evaluates fitted values from an `analog_regression()` /
#' `analog_search(stat = "regression")` output by multiplying coefficient
#' columns with a matrix of covariate values (one row per focal).
#'
#' Handles both single-y and multi-y coefficient layouts and returns a
#' matrix of predictions with n_focal rows and n_y columns.
#'
#' This helper exists so that residual computation in `analog_cv()` and any
#' future user-facing prediction helper share the same arithmetic.
#'
#' @param coefs_df A data.frame with `coef_intercept` and `coef_{covname}`
#'   columns (single-y case), or `coef_intercept_{yname}` and
#'   `coef_{covname}_{yname}` (multi-y case).
#' @param covariates_matrix Matrix with one row per focal and one column per
#'   covariate. Column order must match `cov_names`.
#' @param y_names Character vector of y variable names.
#' @param cov_names Character vector of covariate names, matching the order
#'   of columns in `covariates_matrix`.
#'
#' @return A numeric matrix with `nrow(covariates_matrix)` rows and
#'   `length(y_names)` columns, named by `y_names`.
#'
#' @keywords internal
.predict_from_coefs <- function(coefs_df, covariates_matrix, y_names, cov_names) {
      n_focal <- nrow(covariates_matrix)
      n_y <- length(y_names)
      n_cov <- length(cov_names)

      # Add a leading 1 column for the intercept
      design <- cbind(1.0, covariates_matrix)
      colnames(design) <- c("intercept", cov_names)

      out <- matrix(NA_real_, nrow = n_focal, ncol = n_y)
      colnames(out) <- y_names

      # Detect single-y vs multi-y layout by looking for the unsuffixed
      # `coef_intercept` column.
      single_y <- ("coef_intercept" %in% names(coefs_df)) && (n_y == 1L)

      for (j in seq_len(n_y)) {
            yn <- y_names[j]

            if (single_y) {
                  col_intercept <- "coef_intercept"
                  col_slopes <- paste0("coef_", cov_names)
            } else {
                  col_intercept <- paste0("coef_intercept_", yn)
                  col_slopes <- paste0("coef_", cov_names, "_", yn)
            }

            missing_cols <- setdiff(c(col_intercept, col_slopes), names(coefs_df))
            if (length(missing_cols) > 0) {
                  stop(
                        "Missing expected coefficient columns: ",
                        paste(missing_cols, collapse = ", "),
                        call. = FALSE
                  )
            }

            coef_mat <- cbind(
                  coefs_df[[col_intercept]],
                  do.call(cbind, lapply(col_slopes, function(cn) coefs_df[[cn]]))
            )
            if (is.null(dim(coef_mat))) {
                  coef_mat <- matrix(coef_mat, ncol = 1L)
            }

            out[, j] <- rowSums(design * coef_mat)
      }

      out
}
