#' Search lattice index for analogs
#' @keywords internal
query_analog_index <- function(x,
                               index,
                               select,
                               stat,
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
      params <- .validate_query_params(select, stat, k, weight, theta, x_cov, focal_mm)
      select <- params$select
      stat <- params$stat
      k <- params$k
      weight <- params$weight
      theta <- params$theta
      x_cov_mat <- params$x_cov

      # Parse constraints
      max_geog_num <- if (is.null(max_geog)) Inf else as.numeric(max_geog)[1L]
      max_clim_val <- if (is.null(max_clim)) Inf else max_clim

      # Map select/stat/weight for C++
      # Select codes: 0=knn_clim, 1=knn_geog, 2=all
      select_code <- switch(
            select,
            "knn_clim" = 0L,
            "knn_geog" = 1L,
            "all"      = 2L
      )

      # Aggregate codes: 0=pairs, 1=count, 2=sum_weights, 3=mean_weights
      aggregate_code <- switch(
            stat,
            "pairs"        = 0L,
            "count"        = 1L,
            "sum_weights"  = 2L,
            "mean_weights" = 3L
      )

      # Weight codes (only used when stat is sum_weights or mean_weights)
      weight_code <- if (stat %in% c("sum_weights", "mean_weights")) {
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
      theta_vec <- if (is.null(theta)) {
            NA_real_
      } else {
            as.numeric(theta)
      }

      k_core <- if (select %in% c("knn_clim", "knn_geog")) as.integer(k) else 0L

      # Thread control
      if (!is.null(n_threads)) {
            if (!requireNamespace("RcppParallel", quietly = TRUE)) {
                  stop("Package 'RcppParallel' is required to control 'n_threads'.",
                       call. = FALSE)
            }
            RcppParallel::setThreadOptions(numThreads = as.integer(n_threads)[1L])
      }

      # Call C++ query function with new parameter names
      res <- query_analog_index_cpp(
            index_list = index,
            focal_mm = focal_mm,
            ref_mm = index$ref_data,
            k = k_core,
            max_clim = max_clim_val,
            max_geog = max_geog_num,
            select_code = select_code,
            aggregate_code = aggregate_code,
            weight_code = weight_code,
            theta = theta_vec,
            x_cov = x_cov_mat
      )

      # Capture diagnostic attributes
      cpp_attrs <- attributes(res)
      cpp_attrs$names <- NULL
      cpp_attrs$class <- NULL

      # Post-process results based on aggregate type
      if (stat == "pairs") {
            out <- .emit_pairs_cpp(
                  res,
                  focal_mm,
                  index$ref_data,
                  report_dist = report_dist,
                  geo_mode = index$coord_type,
                  x_cov = x_cov_mat
            )
            names(out) <- gsub("focal_", "", names(out))
            for (nm in names(cpp_attrs)) {
                  attr(out, nm) <- cpp_attrs[[nm]]
            }
            attr(out, "select")    <- select
            attr(out, "stat")      <- stat
            attr(out, "weight")    <- weight
            attr(out, "theta")     <- theta
            return(out)
      }

      if (stat %in% c("sum_weights", "mean_weights", "count")) {
            values <- as.numeric(res)
            if (length(values) != nrow(focal_mm)) {
                  stop("Internal error: stat result length does not match number of focals.")
            }

            out <- data.frame(
                  index = seq_len(nrow(focal_mm)),
                  x     = focal_mm[, 1],
                  y     = focal_mm[, 2],
                  value       = values,
                  stringsAsFactors = FALSE
            )

            for (nm in names(cpp_attrs)) {
                  attr(out, nm) <- cpp_attrs[[nm]]
            }
            attr(out, "select")    <- select
            attr(out, "stat") <- stat
            attr(out, "weight")    <- weight
            attr(out, "theta")     <- theta
            return(out)
      }

      stop("Unreachable code - please report this bug")
}
