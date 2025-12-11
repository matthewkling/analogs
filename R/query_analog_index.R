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
                               n_threads,
                               show_progress = FALSE) {

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
      names(aggregate_codes) <- NULL

      # Map weight function for C++
      # Weight codes: 0=none, 1=uniform, 2=inverse_clim, 3=inverse_geog,
      #               4=gaussian_clim, 5=gaussian_geog,
      #               6=gaussian_joint, 7=inverse_joint
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

      # Convert theta to numeric vector (use NA_real_ for NULL so C++ can apply defaults)
      theta_vec <- if (is.null(theta)) {
            NA_real_
      } else {
            as.numeric(theta)
      }

      # Transform k for C++ (0L for "all", integer otherwise)
      k_core <- if (select %in% c("knn_clim", "knn_geog")) as.integer(k) else 0L

      # If n_threads explicitly provided, set RcppParallel threads
      # This needs to happen before calling C++
      if (!is.null(n_threads)) {
            if (!is.numeric(n_threads) || n_threads < 1) {
                  stop("`n_threads` must be a positive integer",
                       call. = FALSE)
            }
            RcppParallel::setThreadOptions(numThreads = as.integer(n_threads)[1L])
      }

      # Call C++ query function (with optional chunking/progress)
      res <- .query_cpp_chunked(
            index = index,
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
            values = values_mat,
            show_progress = show_progress
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

            return(.format_output(out, focal, stat, select,
                                  k, weight, theta, x_cov_mat))
      }

      # Aggregation mode(s)
      # res is a matrix with n_focal rows and (n_regular_stats + n_value_stats * n_vars) columns
      if (!is.matrix(res)) {
            stop("Internal error: expected matrix result from C++ for aggregation stats")
      }

      if (nrow(res) != nrow(focal)) {
            stop("Internal error: stat result rows do not match number of focals.")
      }

      # Build output data.frame
      out <- data.frame(
            index = seq_len(nrow(focal)),
            x = focal[, 1],
            y = focal[, 2]
      )

      # Determine column structure based on stats requested and values provided
      # Regular stats (count, sum_weights, mean_weights): 1 column each
      # Value-based stats (sum, mean, weighted_sum, weighted_mean): n_values columns each
      regular_stats <- c("count", "sum_weights", "mean_weights")
      value_stats <- c("sum", "mean", "weighted_sum", "weighted_mean")

      col_idx <- 1
      for (s in stat) {
            if (s %in% regular_stats) {
                  # Single column for this stat
                  out[[s]] <- res[, col_idx]
                  col_idx <- col_idx + 1
            } else if (s %in% value_stats) {
                  # One column per value variable
                  n_vals <- if (is.null(values_mat)) 0 else ncol(values_mat)
                  if (n_vals == 0) {
                        stop("Internal error: value stat requested but no values provided")
                  }

                  if (n_vals == 1) {
                        # Single value variable: column named by stat
                        out[[s]] <- res[, col_idx]
                        col_idx <- col_idx + 1
                  } else {
                        # Multiple value variables: columns named {stat}_{varname}
                        for (j in seq_len(n_vals)) {
                              var_name <- if (!is.null(values_names)) {
                                    values_names[j]
                              } else {
                                    paste0("var", j)
                              }
                              col_name <- paste0(s, "_", var_name)
                              out[[col_name]] <- res[, col_idx]
                              col_idx <- col_idx + 1
                        }
                  }
            }
      }

      # Add C++ diagnostic attributes
      for (nm in names(cpp_attrs)) {
            attr(out, nm) <- cpp_attrs[[nm]]
      }

      # Format and return
      return(.format_output(out, focal, stat, select, k, weight, theta, x_cov_mat))
}


#' Internal helper: Query C++ with optional chunking and progress tracking
#'
#' Wraps query_analog_index_cpp with optional chunking for progress bars.
#' Handles merging of chunked results.
#'
#' @keywords internal
.query_cpp_chunked <- function(index, focal, ref_mm, k, max_clim, max_geog,
                               select_code, aggregate_codes, weight_code, theta,
                               x_cov, values, show_progress = FALSE) {

      # If progress not requested or dataset too small, just call C++ directly
      if (!show_progress || nrow(focal) < 100) {
            return(query_analog_index_cpp(
                  index_list = index,
                  focal_mm = focal,
                  ref_mm = ref_mm,
                  k = k,
                  max_clim = max_clim,
                  max_geog = max_geog,
                  select_code = select_code,
                  aggregate_codes = aggregate_codes,
                  weight_code = weight_code,
                  theta = theta,
                  x_cov_sexp = x_cov,
                  values_sexp = values
            ))
      }

      # Split focal into chunks
      n_chunks <- min(100, nrow(focal))
      chunk_size <- ceiling(nrow(focal) / n_chunks)

      # Progress bar
      pb <- txtProgressBar(min = 0, max = n_chunks, style = 3)
      on.exit(close(pb), add = TRUE)

      # Collect results per chunk
      results <- vector("list", n_chunks)

      for (i in seq_len(n_chunks)) {
            start_idx <- (i - 1) * chunk_size + 1
            end_idx <- min(i * chunk_size, nrow(focal))
            focal_chunk <- focal[start_idx:end_idx, , drop = FALSE]

            # Get x_cov chunk if provided
            x_cov_chunk <- if (!is.null(x_cov)) {
                  x_cov[start_idx:end_idx, , drop = FALSE]
            } else {
                  NULL
            }

            # Query this chunk
            results[[i]] <- query_analog_index_cpp(
                  index_list = index,
                  focal_mm = focal_chunk,
                  ref_mm = ref_mm,
                  k = k,
                  max_clim = max_clim,
                  max_geog = max_geog,
                  select_code = select_code,
                  aggregate_codes = aggregate_codes,
                  weight_code = weight_code,
                  theta = theta,
                  x_cov_sexp = x_cov_chunk,
                  values_sexp = values  # Note: values are NOT chunked (they're for ref pool)
            )

            setTxtProgressBar(pb, i)
      }

      # Merge results based on type
      first <- results[[1]]

      if (is.list(first) && !is.matrix(first)) {
            # Pairs mode: concatenate vectors in list
            # Pre-allocate based on total lengths
            total_len <- sum(vapply(results, function(x) length(x[[1]]), integer(1)))
            merged <- list()

            for (nm in names(first)) {
                  merged[[nm]] <- vector(typeof(first[[nm]]), total_len)
                  pos <- 1
                  for (chunk_res in results) {
                        chunk_vec <- chunk_res[[nm]]
                        len <- length(chunk_vec)
                        merged[[nm]][pos:(pos + len - 1)] <- chunk_vec
                        pos <- pos + len
                  }
            }
      } else {
            # Aggregation mode: rbind matrices
            # Pre-allocate full result matrix
            total_rows <- sum(vapply(results, nrow, integer(1)))
            merged <- matrix(0, nrow = total_rows, ncol = ncol(first))
            colnames(merged) <- colnames(first)

            pos <- 1
            for (chunk_res in results) {
                  nrows <- nrow(chunk_res)
                  merged[pos:(pos + nrows - 1), ] <- chunk_res
                  pos <- pos + nrows
            }
      }

      # Preserve attributes from first chunk (diagnostics are constant across chunks)
      attrs <- attributes(first)
      attrs$dim <- NULL
      attrs$dimnames <- NULL
      attrs$names <- NULL
      for (nm in names(attrs)) {
            attr(merged, nm) <- attrs[[nm]]
      }
      attr(merged, "n_focal") <- nrow(focal)

      return(merged)
}
