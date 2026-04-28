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
                               exclude_self = FALSE,
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
      values_levels <- params$values_levels   # NULL unless tabulate
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
      #                  8=ess, 9=regression, 10=tabulate
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
            "regression"     = 9L,
            "tabulate"       = 10L
      )
      aggregate_codes <- stat_name_to_code[stat]
      names(aggregate_codes) <- NULL

      # Map kernel function for C++ (passed to C++ as weight_code).
      # Tabulate is a weighted aggregation (sums combined_weight per class),
      # so it triggers weighting alongside the existing weighted stats.
      has_weighted_stat <- any(stat %in% c("sum_weights", "mean_weights",
                                           "weighted_sum", "weighted_mean",
                                           "ess", "regression", "tabulate"))
      weight_code <- if (has_weighted_stat) {
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

      # Build n_classes_per_var to pass to C++. Length 0 (=integer(0)) when
      # tabulate not requested; length n_vars otherwise.
      n_classes_per_var <- if ("tabulate" %in% stat && !is.null(values_levels)) {
            as.integer(vapply(values_levels, length, integer(1L)))
      } else {
            integer(0)
      }

      # If n_threads explicitly provided, set RcppParallel threads
      if (!is.null(n_threads)) {
            if (!is.numeric(n_threads) || n_threads < 1) {
                  stop("`n_threads` must be a positive integer",
                       call. = FALSE)
            }
            RcppParallel::setThreadOptions(numThreads = as.integer(n_threads)[1L])
      }

      # exclude_self is incompatible with chunking/progress; upstream validation
      # should have caught this, but enforce again defensively.
      if (isTRUE(exclude_self) && isTRUE(show_progress)) {
            stop("`exclude_self = TRUE` is incompatible with `progress = TRUE`.",
                 call. = FALSE)
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
            y = values_mat,
            covariates = covariates_mat,
            lambda = lambda,
            se_code = se_code,
            n_classes_per_var = n_classes_per_var,
            exclude_self = exclude_self,
            show_progress = show_progress
      )

      # Capture diagnostic attributes (before they get lost)
      cpp_attrs <- attributes(res)
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

            for (nm in names(cpp_attrs)) {
                  attr(out, nm) <- cpp_attrs[[nm]]
            }

            return(.format_output(out, focal, stat, select,
                                  k, kernel, theta, x_cov_mat))
      }

      # Aggregation mode(s)
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

      col_idx <- 1L
      n_vals <- if (is.null(values_mat)) 0L else ncol(values_mat)
      want_se <- (se != "none")

      # Regular stats (count, sum_weights, mean_weights, ess): 1 col each
      regular_stats <- intersect(stat,
                                 c("count", "sum_weights", "mean_weights", "ess"))
      for (s in regular_stats) {
            out[[s]] <- res[, col_idx]
            col_idx <- col_idx + 1L
      }

      # Value-based stats: one column per value var
      value_stats <- intersect(stat,
                               c("sum", "mean", "weighted_sum", "weighted_mean"))

      if (length(value_stats) > 0L) {
            for (v in seq_len(n_vals)) {
                  var_name <- if (!is.null(values_names)) {
                        values_names[v]
                  } else {
                        paste0("var", v)
                  }

                  for (s in value_stats) {
                        col_name <- if (n_vals == 1L) s else paste0(s, "_", var_name)
                        out[[col_name]] <- res[, col_idx]
                        col_idx <- col_idx + 1L

                        # weighted_mean SE is emitted immediately after
                        if (s == "weighted_mean" && want_se) {
                              se_name <- if (n_vals == 1L) "se_weighted_mean" else
                                    paste0("se_weighted_mean_", var_name)
                              out[[se_name]] <- res[, col_idx]
                              col_idx <- col_idx + 1L
                        }
                  }
            }
      }

      # Tabulate: K_v columns per y variable v, contiguous in v-order.
      # Column naming follows the package convention:
      #   single y, no name        -> "n_<level>"
      #   y has a name (or multi-y)-> "<varname>_n_<level>"
      if ("tabulate" %in% stat) {
            for (v in seq_len(n_vals)) {
                  var_name <- if (!is.null(values_names)) {
                        values_names[v]
                  } else {
                        paste0("var", v)
                  }
                  lev <- values_levels[[v]]
                  K_v <- length(lev)

                  # Decide whether to prefix the variable name. We always
                  # prefix when y has a real name attached (e.g. raster layer
                  # or data.frame column) or when there are multiple y vars,
                  # to disambiguate. For a single bare vector with our
                  # auto-generated "y1", we drop the prefix.
                  user_named <- (n_vals > 1L) ||
                        (!identical(var_name, "y1") &&
                               !identical(var_name, "value_1"))

                  for (kk in seq_len(K_v)) {
                        col_name <- if (user_named) {
                              paste0(var_name, "_n_", lev[kk])
                        } else {
                              paste0("n_", lev[kk])
                        }
                        out[[col_name]] <- res[, col_idx]
                        col_idx <- col_idx + 1L
                  }
            }
      }

      # Regression: (n_covariates + 1) columns per value variable for coefs;
      # same again for SEs when se != "none".
      if ("regression" %in% stat) {
            n_covs <- if (is.null(covariates_mat)) 0 else ncol(covariates_mat)
            reg_dim <- n_covs + 1

            base_names <- c("intercept", covariates_names %||% paste0("cov", seq_len(n_covs)))

            if (n_vals == 1) {
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

      for (nm in names(cpp_attrs)) {
            attr(out, nm) <- cpp_attrs[[nm]]
      }

      return(.format_output(out, focal, stat, select, k, kernel, theta, x_cov_mat))
}


#' Internal helper: Query C++ with optional chunking and progress tracking
#'
#' Wraps query_analog_index_cpp with optional chunking for progress bars.
#' Handles merging of chunked results.
#'
#' @keywords internal
.query_cpp_chunked <- function(index, focal, ref_mm, k, max_clim, max_geog,
                               select_code, aggregate_codes, weight_code, theta,
                               x_cov, y, covariates, lambda, se_code,
                               n_classes_per_var = integer(0),
                               exclude_self = FALSE,
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
                  weight_code = weight_code,
                  theta = theta,
                  x_cov_sexp = x_cov,
                  values_sexp = y,
                  covariates_sexp = covariates,
                  lambda = lambda,
                  se_code = se_code,
                  n_classes_per_var = n_classes_per_var,
                  exclude_self = exclude_self
            ))
      }

      # Chunked path (progress bar). exclude_self is incompatible with chunking
      # because the worker's self-check uses j == i, which only holds for the
      # full focal set. Upstream validation should prevent this combination.
      if (isTRUE(exclude_self)) {
            stop("Internal error: exclude_self = TRUE reached chunked path. ",
                 "This combination should have been rejected upstream.",
                 call. = FALSE)
      }

      # Split focal into chunks
      n_chunks <- min(100, nrow(focal))
      chunk_size <- ceiling(nrow(focal) / n_chunks)

      # Progress bar
      pb <- utils::txtProgressBar(min = 0, max = n_chunks, style = 3)
      on.exit(close(pb), add = TRUE)

      # Collect results per chunk
      results <- vector("list", n_chunks)

      for (i in seq_len(n_chunks)) {
            start_idx <- (i - 1) * chunk_size + 1
            end_idx <- min(i * chunk_size, nrow(focal))
            focal_chunk <- focal[start_idx:end_idx, , drop = FALSE]

            x_cov_chunk <- if (!is.null(x_cov)) {
                  x_cov[start_idx:end_idx, , drop = FALSE]
            } else {
                  NULL
            }

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
                  values_sexp = y,
                  covariates_sexp = covariates,
                  lambda = lambda,
                  se_code = se_code,
                  n_classes_per_var = n_classes_per_var,
                  exclude_self = FALSE   # always FALSE in chunked path
            )

            utils::setTxtProgressBar(pb, i)
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

            for (nm in names(attributes(first))) {
                  if (!nm %in% c("dim", "dimnames")) {
                        attr(merged, nm) <- attr(first, nm)
                  }
            }
            attr(merged, "n_focal") <- nrow(focal)

            return(merged)
      }
}
