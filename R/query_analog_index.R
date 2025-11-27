#' Search lattice index for analogs
#' @keywords internal
query_analog_index <- function(x,
                               index,
                               mode,
                               max_clim,
                               max_geog,
                               x_cov,
                               k,
                               weight,
                               theta,
                               report_dist,
                               n_threads) {

      # Validate index
      .validate_analog_index(index)

      # Format focal data (needed for validation)
      focal_mm <- .format_data(x)

      # Validate compatibility
      .validate_analog_index(index, focal_mm, validate_ranges = FALSE)

      # Validate and normalize query parameters (including x_cov)
      params <- .validate_query_params(mode, k, weight, theta, x_cov, focal_mm)
      mode <- params$mode
      k <- params$k
      weight <- params$weight
      theta <- params$theta
      x_cov_mat <- params$x_cov

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
                  "uniform"        = 1L,
                  "inverse_clim"   = 2L,
                  "inverse_geog"   = 3L,
                  "gaussian_clim"  = 4L,
                  "gaussian_geog"  = 5L,
                  "gaussian_joint" = 6L,
                  "inverse_joint"  = 7L
            )
      } else {
            0L
      }

      # Handle theta: convert to numeric vector (length 1 or 2)
      # For joint weights, theta should already be length 2 from validation
      # For single weights, theta is length 1 or NULL (becomes NA_real_)
      theta_vec <- if (is.null(theta)) {
            NA_real_
      } else {
            as.numeric(theta)
      }

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
      # Pass x_cov_mat directly - C++ handles NULL properly
      res <- query_analog_index_cpp(
            index_list = index,
            focal_mm = focal_mm,
            ref_mm = index$ref_data,
            k = k_core,
            max_clim = max_clim_val,
            max_geog = max_geog_num,
            mode_code = mode_code,
            weight_code = weight_code,
            theta = theta_vec,
            x_cov = x_cov_mat  # NULL or matrix; C++ handles both
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
