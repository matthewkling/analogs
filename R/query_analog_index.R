#' Search lattice index for analogs
#' @keywords internal
query_analog_index <- function(x,
                               index,
                               select,
                               stat,
                               max_clim,
                               max_geog,
                               x_cov,
                               values,
                               k,
                               weight,
                               theta,
                               n_threads) {

      # Validate index
      .validate_analog_index(index)

      # Format focal data (needed for validation)
      focal <- .format_data(x)

      # Validate compatibility
      .validate_analog_index(index, focal, validate_ranges = FALSE)

      # Validate and normalize query parameters
      params <- .validate_query_params(focal, index$ref_data,
                                       x_cov, values,
                                       max_clim, max_geog,
                                       select, k,
                                       stat, weight, theta)
      select <- params$select
      stat <- params$stat
      k <- params$k
      weight <- params$weight
      theta <- params$theta
      x_cov_mat <- params$x_cov
      values_mat <- params$values
      values_names <- params$values_names

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
      # Aggregate codes: 0=none, 1=count, 2=sum_weights, 3=mean_weights,
      #                  4=sum, 5=mean, 6=weighted_sum, 7=weighted_mean
      stat_name_to_code <- c(
            "none"           = 0L,
            "count"          = 1L,
            "sum_weights"    = 2L,
            "mean_weights"   = 3L,
            "sum"            = 4L,
            "mean"           = 5L,
            "weighted_sum"   = 6L,
            "weighted_mean"  = 7L
      )
      aggregate_codes <- stat_name_to_code[stat]
      names(aggregate_codes) <- NULL  # Remove names for C++

      # Weight code (only used when stat includes weighted stats)
      has_weighted_stat <- any(stat %in% c("sum_weights", "mean_weights",
                                           "weighted_sum", "weighted_mean"))
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

      # Call C++ query function with vector of aggregate codes and values
      res <- query_analog_index_cpp(
            index_list = index,
            focal = focal,
            ref_mm = index$ref_data,
            k = k_core,
            max_clim = max_clim_val,
            max_geog = max_geog_num,
            select_code = select_code,
            aggregate_codes = aggregate_codes,
            weight_code = weight_code,
            theta = theta_vec,
            x_cov = x_cov_mat,
            values = values_mat
      )

      # Capture diagnostic attributes (before they get lost)
      cpp_attrs <- attributes(res)
      # Remove standard R attributes that we'll replace
      cpp_attrs$names <- NULL
      cpp_attrs$class <- NULL
      cpp_attrs$dim <- NULL
      cpp_attrs$dimnames <- NULL

      # Post-process results based on aggregate type
      if (identical(stat, "none") || (length(stat) == 1 && stat[1] == "none")) {
            # Pairs mode - res is a list
            out <- .emit_pairs_cpp(
                  res,
                  focal,
                  index$ref_data,
                  report_dist = TRUE,
                  geo_mode = index$coord_type,
                  x_cov = x_cov_mat,
                  values = values_mat,
                  values_names = values_names
            )
            names(out) <- gsub("focal_", "", names(out))

            # Add attributes to data.frame
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
      # res is a matrix with n_focal rows and (n_regular_stats + n_value_stats * n_vars) columns
      if (!is.matrix(res)) {
            stop("Internal error: expected matrix result from C++ for aggregation stats")
      }

      if (nrow(res) != nrow(focal)) {
            stop("Internal error: stat result rows do not match number of focals.")
      }

      # Build output data.frame with named columns for each stat
      out <- data.frame(
            index = seq_len(nrow(focal)),
            x     = focal[, 1],
            y     = focal[, 2],
            stringsAsFactors = FALSE
      )

      # Determine which stats are value-based
      value_stats <- c("sum", "mean", "weighted_sum", "weighted_mean")
      regular_stats <- c("count", "sum_weights", "mean_weights")

      regular_stat_list <- intersect(stat, regular_stats)
      value_stat_list <- intersect(stat, value_stats)

      n_regular <- length(regular_stat_list)
      n_value <- length(value_stat_list)
      n_vars <- if (!is.null(values_mat)) ncol(values_mat) else 0

      # Expected number of columns
      expected_cols <- n_regular + n_value * n_vars
      if (ncol(res) != expected_cols) {
            stop("Internal error: stat result columns (", ncol(res),
                 ") do not match expected (", expected_cols, ")")
      }

      col_idx <- 1

      # Add regular stats (apply to all focals, not variable-specific)
      for (s in regular_stat_list) {
            out[[s]] <- res[, col_idx]
            col_idx <- col_idx + 1
      }

      # Add value stats (one column per stat per variable)
      if (n_value > 0 && n_vars > 0) {
            for (var_idx in seq_len(n_vars)) {
                  var_name <- values_names[var_idx]
                  for (s in value_stat_list) {
                        col_name <- paste(s, var_name, sep = "_")
                        out[[col_name]] <- res[, col_idx]
                        col_idx <- col_idx + 1
                  }
            }
      }

      # Add attributes to data.frame
      for (nm in names(cpp_attrs)) {
            attr(out, nm) <- cpp_attrs[[nm]]
      }
      attr(out, "select")    <- select
      attr(out, "stat")      <- stat
      attr(out, "weight")    <- weight
      attr(out, "theta")     <- theta

      return(out)
}
