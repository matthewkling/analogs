
#' Search lattice index for analogs
#' @keywords internal
query_analog_index <- function(x,
                               index,
                               mode,
                               max_clim,
                               max_geog,
                               k,
                               weight,
                               theta,
                               report_dist,
                               n_threads) {

      # Validate index
      .validate_analog_index(index)

      # Input validation (same as main function)
      mode <- match.arg(mode, c("knn_clim", "knn_geog", "count", "sum", "mean", "all"))

      weight <- match.arg(weight, c("uniform", "gaussian_clim", "gaussian_geog",
                                    "gaussian_joint", "inverse_clim", "inverse_geog",
                                    "inverse_joint"))
      if(!weight %in% c("inverse_clim", "inverse_geog")) weight <- NULL

      # Validate mode/k/weight/theta combinations (same logic as main function)
      if (mode %in% c("knn_clim", "knn_geog")) {
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
      } else {
            if (!is.null(k)) {
                  stop("For mode '", mode, "', k must be NULL.")
            }
            k <- 0L
      }

      if (mode %in% c("all", "count")) {
            if (!is.null(weight)) {
                  stop("For mode '", mode, "', weight must be NULL.")
            }
            if (!is.null(theta)) {
                  stop("For mode '", mode, "', theta must be NULL.")
            }
      }

      if (mode %in% c("sum", "mean")) {
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
                  if (!is.null(theta)) {
                        if (!is.numeric(theta) || length(theta) != 1L || theta <= 0) {
                              stop("theta must be a single positive numeric value, or NULL.")
                        }
                  }
            }
      } else {
            if (!is.null(weight)) {
                  stop("weight must be NULL when mode is not 'sum' or 'mean'.")
            }
            if (!is.null(theta)) {
                  stop("theta must be NULL when mode is not 'sum' or 'mean'.")
            }
      }

      # Format focal data
      focal_mm <- .format_data(x)

      # Validate compatibility
      .validate_analog_index(index, focal_mm, validate_ranges = FALSE)

      # Parse constraints
      max_geog_num <- if (is.null(max_geog)) Inf else as.numeric(max_geog)[1L]
      max_clim_val <- if (is.null(max_clim)) Inf else max_clim

      # Map mode/weight/theta for C++
      mode_code <- switch(
            mode,
            "knn_clim" = 0L,
            "knn_geog" = 1L,
            "count"    = 2L,
            "sum"      = 3L,
            "mean"     = 4L,
            "all"      = 5L
      )

      weight_code <- if (mode %in% c("sum","mean")) {
            switch(
                  weight,
                  "uniform"      = 1L,
                  "inverse_clim" = 2L,
                  "inverse_geog" = 3L
            )
      } else {
            0L
      }

      theta_num <- if (is.null(theta)) NA_real_ else as.numeric(theta)[1L]
      k_core <- if (mode %in% c("knn_clim","knn_geog")) as.integer(k) else 0L

      # Thread control
      if (!is.null(n_threads)) {
            if (!requireNamespace("RcppParallel", quietly = TRUE)) {
                  stop("Package 'RcppParallel' is required to control 'n_threads'.",
                       call. = FALSE)
            }
            RcppParallel::setThreadOptions(numThreads = as.integer(n_threads)[1L])
      }

      # Call C++ query function
      res <- query_analog_index_cpp(
            index_list = index,
            focal_mm = focal_mm,
            ref_mm = index$ref_data,
            k = k_core,
            max_clim = max_clim_val,
            max_geog = max_geog_num,
            mode_code = mode_code,
            weight_code = weight_code,
            theta = theta_num
      )

      # Capture diagnostic attributes
      cpp_attrs <- attributes(res)
      cpp_attrs$names <- NULL
      cpp_attrs$class <- NULL

      # Post-process results (same as main function)
      if (mode %in% c("knn_clim", "knn_geog", "all")) {
            out <- .emit_pairs_cpp(
                  res,
                  focal_mm,
                  index$ref_data,
                  report_dist = report_dist,
                  geo_mode = index$coord_type
            )
            for (nm in names(cpp_attrs)) {
                  attr(out, nm) <- cpp_attrs[[nm]]
            }
            attr(out, "mode")   <- mode
            attr(out, "weight") <- weight
            attr(out, "theta")  <- theta
            return(out)
      }

      if (mode %in% c("sum", "mean", "count")) {
            values <- as.numeric(res)
            if (length(values) != nrow(focal_mm)) {
                  stop("Internal error: aggregate result length does not match number of focals.")
            }

            out <- data.frame(
                  focal_index = seq_len(nrow(focal_mm)),
                  focal_x     = focal_mm[, 1],
                  focal_y     = focal_mm[, 2],
                  value       = values,
                  stringsAsFactors = FALSE
            )

            for (nm in names(cpp_attrs)) {
                  attr(out, nm) <- cpp_attrs[[nm]]
            }
            attr(out, "mode")   <- mode
            attr(out, "weight") <- weight
            attr(out, "theta")  <- theta
            return(out)
      }

      stop("Unreachable code - please report this bug")
}
