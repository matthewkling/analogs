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

      # Map select for C++
      # Select codes: 0=knn_clim, 1=knn_geog, 2=all
      select_code <- switch(
            select,
            "knn_clim" = 0L,
            "knn_geog" = 1L,
            "all"      = 2L
      )

      # Map stat(s) for C++
      # Aggregate codes: 0=none, 1=count, 2=sum_weights, 3=mean_weights
      stat_name_to_code <- c(
            "none"         = 0L,
            "count"        = 1L,
            "sum_weights"  = 2L,
            "mean_weights" = 3L
      )
      aggregate_codes <- stat_name_to_code[stat]
      names(aggregate_codes) <- NULL  # Remove names for C++

      # Weight code (only used when stat includes sum_weights or mean_weights)
      has_weighted_stat <- any(stat %in% c("sum_weights", "mean_weights"))
      weight_code <- if (has_weighted_stat) {
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

      # Call C++ query function with vector of aggregate codes
      res <- query_analog_index_cpp(
            index_list = index,
            focal_mm = focal_mm,
            ref_mm = index$ref_data,
            k = k_core,
            max_clim = max_clim_val,
            max_geog = max_geog_num,
            select_code = select_code,
            aggregate_codes = aggregate_codes,
            weight_code = weight_code,
            theta = theta_vec,
            x_cov = x_cov_mat
      )

      # Capture diagnostic attributes
      cpp_attrs <- attributes(res)
      cpp_attrs$names <- NULL
      cpp_attrs$class <- NULL
      cpp_attrs$dim <- NULL
      cpp_attrs$dimnames <- NULL

      # Post-process results based on aggregate type
      if (identical(stat, "none") || (length(stat) == 1 && stat[1] == "none")) {
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

      # Aggregation mode(s)
      # res is a matrix with n_focal rows and length(stat) columns
      if (!is.matrix(res)) {
            stop("Internal error: expected matrix result from C++ for aggregation stats")
      }

      if (nrow(res) != nrow(focal_mm)) {
            stop("Internal error: stat result rows do not match number of focals.")
      }

      if (ncol(res) != length(stat)) {
            stop("Internal error: stat result columns do not match number of stats.")
      }

      # Build output data.frame with named columns for each stat
      out <- data.frame(
            index = seq_len(nrow(focal_mm)),
            x     = focal_mm[, 1],
            y     = focal_mm[, 2],
            stringsAsFactors = FALSE
      )

      # Add each stat as a named column
      for (i in seq_along(stat)) {
            out[[stat[i]]] <- res[, i]
      }

      # Add attributes
      for (nm in names(cpp_attrs)) {
            attr(out, nm) <- cpp_attrs[[nm]]
      }
      attr(out, "select")    <- select
      attr(out, "stat")      <- stat
      attr(out, "weight")    <- weight
      attr(out, "theta")     <- theta

      return(out)
}
