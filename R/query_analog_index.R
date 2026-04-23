#' Search lattice index for analogs
#' @keywords internal
query_analog_index <- function(x,
                               index,
                               select,
                               stat,
                               max_clim,
                               max_geog,
                               x_cov,
                               y,
                               covariates,
                               k,
                               kernel,
                               theta,
                               lambda,
                               se,
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
                                       x_cov, y, covariates,
                                       max_clim, max_geog,
                                       select, k,
                                       stat, kernel, theta, lambda,
                                       se)
      select <- params$select
      stat <- params$stat
      k <- params$k
      kernel <- params$kernel
      theta <- params$theta
      lambda <- params$lambda
      se <- params$se
      x_cov_mat <- params$x_cov
      values_mat <- params$y
      values_names <- params$values_names
      covariates_mat <- params$covariates
      covariates_names <- params$covariates_names

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
      #                  4=sum, 5=mean, 6=weighted_sum, 7=weighted_mean,
      #                  8=ess, 9=regression
      stat_name_to_code <- c(
            "none"           = 0L,
            "count"          = 1L,
            "sum_weights"    = 2L,
            "mean_weights"   = 3L,
            "sum"            = 4L,
            "mean"           = 5L,
            "weighted_sum"   = 6L,
            "weighted_mean"  = 7L,
            "ess"            = 8L,
            "regression"     = 9L
      )
      aggregate_codes <- stat_name_to_code[stat]
      names(aggregate_codes) <- NULL

      # Map kernel function for C++
      # kernel codes: 0=none, 1=uniform, 2=inverse_clim, 3=inverse_geog,
      #               4=gaussian_clim, 5=gaussian_geog,
      #               6=gaussian_joint, 7=inverse_joint
      has_weighted_stat <- any(stat %in% c("sum_weights", "mean_weights",
                                           "weighted_sum", "weighted_mean",
                                           "ess", "regression"))
      kernel_code <- if (has_weighted_stat) {
            switch(
                  kernel,
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

      # Map se for C++
      # SE codes: 0=none, 1=ess, 2=design
      se_code <- switch(
            se,
            "none"   = 0L,
            "ess"    = 1L,
            "design" = 2L
      )

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
            kernel_code = kernel_code,
            theta = theta_vec,
            x_cov = x_cov_mat,
            y = values_mat,
            covariates = covariates_mat,
            lambda = lambda,
            se_code = se_code,
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
                                  k, kernel, theta, x_cov_mat))
      }

      # Aggregation mode(s)
      # res is a matrix with n_focal rows and variable number of columns
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
      # Regular stats (count, sum_weights, mean_weights, ess): 1 column each
      # Value-based stats (sum, mean, weighted_sum, weighted_mean): n_values columns each
      # weighted_mean SE (when se != "none"): one column per value var, emitted
      #   immediately after its weighted_mean column within each var's block.
      # Regression: (n_covariates + 1) columns per value variable for coefs;
      #   same again for SEs when se != "none".
      regular_stats <- c("count", "sum_weights", "mean_weights", "ess")
      value_stats <- c("sum", "mean", "weighted_sum", "weighted_mean")

      want_se <- !identical(se, "none")

      col_idx <- 1

      # IMPORTANT: C++ writes columns grouped by type (regular first, then value
      # stats grouped by variable, then regression). We must read in that order,
      # NOT in request order.

      # Read all regular stats first
      for (s in stat) {
            if (s %in% regular_stats) {
                  out[[s]] <- res[, col_idx]
                  col_idx <- col_idx + 1
            }
      }

      # Then read all value stats, grouped by variable (v0_stat0, v0_stat1, ...,
      # v1_stat0, ...). The order of stats within each variable's block matches
      # the order they appear in `stat`, filtered to value stats. When se != "none"
      # and weighted_mean is requested, an SE column is emitted right after the
      # weighted_mean column within each variable's block.
      n_vals <- if (is.null(values_mat)) 0 else ncol(values_mat)
      value_stats_in_order <- intersect(stat, value_stats)  # preserves user order
      has_weighted_mean <- "weighted_mean" %in% value_stats_in_order

      if (length(value_stats_in_order) > 0) {
            if (n_vals == 0) {
                  stop("Internal error: value stat requested but no y values provided")
            }

            for (j in seq_len(n_vals)) {
                  var_name <- if (!is.null(values_names)) {
                        values_names[j]
                  } else {
                        paste0("var", j)
                  }
                  for (s in value_stats_in_order) {
                        # Name for the stat column
                        col_name_stat <- if (n_vals == 1) s else paste0(s, "_", var_name)
                        out[[col_name_stat]] <- res[, col_idx]
                        col_idx <- col_idx + 1

                        # SE column follows immediately if applicable
                        if (want_se && s == "weighted_mean") {
                              col_name_se <- if (n_vals == 1) {
                                    "se_weighted_mean"
                              } else {
                                    paste0("se_weighted_mean_", var_name)
                              }
                              out[[col_name_se]] <- res[, col_idx]
                              col_idx <- col_idx + 1
                        }
                  }
            }
      }

      # Then read regression coefficients and (optionally) standard errors.
      # Layout in C++: for each variable v, write reg_dim coefficients, then
      # (if se != "none") a matching block of reg_dim SEs in a second pass
      # that runs after all variable coefficient blocks.
      if ("regression" %in% stat) {
            n_covs <- if (is.null(covariates_mat)) 0 else ncol(covariates_mat)
            reg_dim <- n_covs + 1  # intercept + covariates

            # Base names: intercept, then covariate names
            base_names <- c("intercept", covariates_names %||% paste0("cov", seq_len(n_covs)))

            if (n_vals == 1) {
                  # Single value variable
                  for (d in seq_len(reg_dim)) {
                        out[[paste0("coef_", base_names[d])]] <- res[, col_idx]
                        col_idx <- col_idx + 1
                  }
                  if (want_se) {
                        for (d in seq_len(reg_dim)) {
                              out[[paste0("se_", base_names[d])]] <- res[, col_idx]
                              col_idx <- col_idx + 1
                        }
                  }
            } else {
                  # Multiple value variables: {prefix}_{coeff}_{varname}
                  # First all coefficients (per var), then all SEs (per var) if requested
                  for (j in seq_len(n_vals)) {
                        var_name <- if (!is.null(values_names)) {
                              values_names[j]
                        } else {
                              paste0("var", j)
                        }
                        for (d in seq_len(reg_dim)) {
                              col_name <- paste0("coef_", base_names[d], "_", var_name)
                              out[[col_name]] <- res[, col_idx]
                              col_idx <- col_idx + 1
                        }
                  }
                  if (want_se) {
                        for (j in seq_len(n_vals)) {
                              var_name <- if (!is.null(values_names)) {
                                    values_names[j]
                              } else {
                                    paste0("var", j)
                              }
                              for (d in seq_len(reg_dim)) {
                                    col_name <- paste0("se_", base_names[d], "_", var_name)
                                    out[[col_name]] <- res[, col_idx]
                                    col_idx <- col_idx + 1
                              }
                        }
                  }
            }
      }

      # Add C++ diagnostic attributes
      for (nm in names(cpp_attrs)) {
            attr(out, nm) <- cpp_attrs[[nm]]
      }

      # Format and return
      return(.format_output(out, focal, stat, select, k, kernel, theta, x_cov_mat))
}


#' Internal helper: Query C++ with optional chunking and progress tracking
#'
#' Wraps query_analog_index_cpp with optional chunking for progress bars.
#' Handles merging of chunked results.
#'
#' @keywords internal
.query_cpp_chunked <- function(index, focal, ref_mm, k, max_clim, max_geog,
                               select_code, aggregate_codes, kernel_code, theta,
                               x_cov, y, covariates, lambda, se_code,
                               show_progress = FALSE) {

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
                  weight_code = kernel_code,
                  theta = theta,
                  x_cov_sexp = x_cov,
                  values_sexp = y,
                  covariates_sexp = covariates,
                  lambda = lambda,
                  se_code = se_code
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
                  weight_code = kernel_code,
                  theta = theta,
                  x_cov_sexp = x_cov_chunk,
                  values_sexp = y,     # Not chunked (ref pool)
                  covariates_sexp = covariates,  # Not chunked (ref pool)
                  lambda = lambda,
                  se_code = se_code
            )

            setTxtProgressBar(pb, i)
      }

      # Merge results based on type
      first <- results[[1]]

      if (is.list(first) && !is.matrix(first)) {
            # Pairs mode: concatenate vectors in list
            total_len <- sum(vapply(results, function(x) length(x[[1]]), integer(1)))

            merged <- lapply(seq_along(first), function(j) {
                  if (is.integer(first[[j]])) {
                        out <- integer(total_len)
                  } else {
                        out <- numeric(total_len)
                  }
                  pos <- 1L
                  for (r in results) {
                        n <- length(r[[j]])
                        out[pos:(pos + n - 1)] <- r[[j]]
                        pos <- pos + n
                  }
                  out
            })
            names(merged) <- names(first)

            # Adjust focal indices for chunks
            # (indices need offset for each chunk)
            offset <- 0L
            pos <- 1L
            for (ch in seq_len(n_chunks)) {
                  n <- length(results[[ch]][[1]])
                  chunk_start <- (ch - 1) * chunk_size
                  # Focal indices in results are 1-based within chunk; add offset
                  # (Already handled by emit_pairs_cpp returning per-chunk indices)
                  pos <- pos + n
            }

            # Copy attributes from first result
            for (nm in names(attributes(first))) {
                  if (!nm %in% c("names")) {
                        attr(merged, nm) <- attr(first, nm)
                  }
            }
            attr(merged, "n_focal") <- nrow(focal)

            return(merged)
      } else {
            # Aggregation mode: rbind matrices
            merged <- do.call(rbind, results)

            # Copy attributes from first result
            for (nm in names(attributes(first))) {
                  if (!nm %in% c("dim", "dimnames")) {
                        attr(merged, nm) <- attr(first, nm)
                  }
            }
            attr(merged, "n_focal") <- nrow(focal)

            return(merged)
      }
}
