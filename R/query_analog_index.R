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
                               weight = NULL,
                               k,
                               kernel_clim,
                               kernel_geog,
                               theta_clim,
                               theta_geog,
                               lambda,
                               se,
                               normalize = "auto",
                               exclude_self = FALSE,
                               n_threads,
                               show_progress = FALSE) {

      # Validate index
      .validate_analog_index(index)

      # Capture the original focal size (before NA stripping below) so we
      # can reconstruct full-length output and validate user-supplied
      # focal-keyed inputs (x_cov) against the user's actual input shape.
      n_focal_original <- if (inherits(x, "SpatRaster")) {
            if (!requireNamespace("terra", quietly = TRUE)) {
                  stop("Package 'terra' is required for SpatRaster inputs",
                       call. = FALSE)
            }
            terra::ncell(x)
      } else {
            nrow(x)
      }

      # Format focal data. NA-containing rows are stripped here; row_map
      # records the original-focal row index of each kept row. Downstream
      # post-processing scatters C++ output back to original-focal indexing
      # (NA at stripped positions) and remaps the `index` column in pairs
      # mode. Identical mechanism to pool-side stripping in build_analog_index().
      focal <- .format_data(x)
      focal_row_map <- attr(focal, "row_map")  # NULL if no rows stripped

      # Validate compatibility
      .validate_analog_index(index, focal, validate_ranges = FALSE)

      # Validate and normalize the user-supplied per-pool weight, if any.
      # Returns a length-n_pool_used numeric vector with NAs zeroed out, or
      # NULL. User input is keyed to the original pool size; pool_row_map
      # translates to the post-NA-strip ordering used by ref_data and C++.
      user_weight_vec <- .validate_and_format_weight(
            weight,
            index$ref_data,
            pool_row_map = index$pool_row_map,
            n_pool_original = index$n_pool
      )

      # Permissive warning: a downsampled index has discarded some pool points,
      # so weights for those points are silently dropped. The user gets
      # weights only for the surviving sample.
      if (!is.null(user_weight_vec) &&
          !is.null(index$downsample_actual) &&
          index$downsample_actual < 1.0) {
            warning(
                  sprintf(
                        "`weight` was provided for %d points, but the index was built with downsample = %.3f, retaining ~%.1f%% of points. Weights for discarded points have no effect; weights for surviving points are used as supplied.",
                        length(user_weight_vec),
                        index$downsample_target %||% NA_real_,
                        100 * index$downsample_actual
                  ),
                  call. = FALSE
            )
      }

      # Validate and normalize query parameters. The four row_map / size
      # arguments let the per-pool validators (y, covariates) accept user
      # input keyed to the original pool size, and the focal-side validator
      # (x_cov) accept input keyed to the original focal size, then
      # translate each to the post-NA-strip ordering used by ref_data /
      # focal and the C++ side.
      params <- .validate_query_params(focal, index$ref_data,
                                       x_cov, y, covariates,
                                       max_clim, max_geog,
                                       select, k,
                                       stat,
                                       kernel_clim, kernel_geog,
                                       theta_clim, theta_geog,
                                       lambda,
                                       se,
                                       pool_row_map = index$pool_row_map,
                                       n_pool_original = index$n_pool,
                                       focal_row_map = focal_row_map,
                                       n_focal_original = n_focal_original)
      select <- params$select
      stat <- params$stat
      k <- params$k
      kernel_clim <- params$kernel_clim
      kernel_geog <- params$kernel_geog
      theta_clim <- params$theta_clim
      theta_geog <- params$theta_geog
      lambda <- params$lambda
      se <- params$se
      x_cov_mat <- params$x_cov
      values_mat <- params$y
      values_names <- params$values_names
      values_levels <- params$values_levels   # NULL unless tabulate
      covariates_mat <- params$covariates
      covariates_names <- params$covariates_names

      # Resolve `normalize = "auto"` to a concrete logical. `"auto"` does
      # the right thing whenever it can: TRUE iff every precondition is
      # met (raster-derived index with active cell-area weighting AND a
      # climate-aware kernel AND a finite max_geog), FALSE otherwise. This
      # silent fall-back lets the same call work across raster vs.
      # non-raster pools, with vs. without max_geog, and with any kernel
      # choice -- without the user having to manually toggle the argument.
      #
      # Explicit `normalize = TRUE` is stricter: any precondition failure
      # raises an error in .validate_normalize_compat() below. That gives
      # users who explicitly opted in clear feedback about why
      # normalization can't apply.
      if (identical(normalize, "auto")) {
            normalize <- .normalize_is_feasible(stat, max_geog, index)
      }

      # Validate normalize compatibility now that kernel/max_geog/etc. are
      # resolved and we have the index in hand. Throws on any
      # incompatibility (uniform/geo-only kernel, missing max_geog,
      # non-raster index, missing cell-area weighting, etc.).
      .validate_normalize_compat(normalize, stat, max_geog, index)

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

      # Map the per-family kernel shapes to C++ integer codes. The C++ side uses
      # a product model: combined weight = clim-family weight x geog-family
      # weight, each shape UNIFORM(0)/GAUSSIAN(1)/INVERSE(2). theta is resolved
      # to a per-family scalar (NA_real_ lets C++ apply its default of 1).
      .kernel_code <- function(k) {
            switch(k, "uniform" = 0L, "gaussian" = 1L, "inverse" = 2L,
                   stop("Internal error: unknown family kernel '", k, "'.",
                        call. = FALSE))
      }
      clim_kernel_code <- .kernel_code(kernel_clim)
      geog_kernel_code <- .kernel_code(kernel_geog)
      theta_clim_cpp <- if (is.null(theta_clim)) NA_real_ else as.numeric(theta_clim)
      theta_geog_cpp <- if (is.null(theta_geog)) NA_real_ else as.numeric(theta_geog)

      # Map se for C++
      # SE codes: 0=none, 1=ess, 2=design
      se_code <- switch(
            se,
            "none"   = 0L,
            "ess"    = 1L,
            "design" = 2L
      )

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

      # Pool-weight streams. Cell-area weights live on the index (build-time
      # property); user weights are query-time. Either may be NULL.
      area_weight_vec <- index$cell_area_weight  # NULL if not applied

      # Compute the global D_max scalar once per query, if normalization is
      # both requested AND would actually apply to a stat in this query.
      # Validation above ensures preconditions are met when we reach
      # .compute_global_dmax(). The stat-relevance gate mirrors the
      # validator's short-circuit so an explicit `normalize = TRUE` with
      # stat = "count" (no normalizable column) is a silent no-op and
      # doesn't fail trying to integrate against a missing kernel.
      do_dmax <- isTRUE(normalize) && any(stat %in% c("sum_weights", "tabulate"))
      D_max <- if (do_dmax) {
            .compute_global_dmax(
                  kernel_clim    = kernel_clim,
                  kernel_geog    = kernel_geog,
                  theta_clim     = theta_clim,
                  theta_geog     = theta_geog,
                  max_geog       = max_geog_num,
                  mean_cell_area = index$mean_cell_area
            )
      } else {
            NA_real_
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
            clim_kernel_code = clim_kernel_code,
            geog_kernel_code = geog_kernel_code,
            theta_clim = theta_clim_cpp,
            theta_geog = theta_geog_cpp,
            x_cov = x_cov_mat,
            y = values_mat,
            covariates = covariates_mat,
            lambda = lambda,
            se_code = se_code,
            n_classes_per_var = n_classes_per_var,
            area_weight = area_weight_vec,
            user_weight = user_weight_vec,
            exclude_self = exclude_self,
            show_progress = show_progress
      )

      # Capture diagnostic attributes (before they get lost)
      cpp_attrs <- attributes(res)
      cpp_attrs$names <- NULL
      cpp_attrs$class <- NULL
      cpp_attrs$dim <- NULL
      cpp_attrs$dimnames <- NULL

      # Emission flags for the three weight streams in pair mode. Each column
      # is surfaced only when its mechanism is actually active for this query.
      emit_sample_weight <- !is.null(index$downsample_actual) &&
            index$downsample_actual < 1.0
      emit_area_weight   <- !is.null(area_weight_vec)
      emit_user_weight   <- !is.null(user_weight_vec)

      # Post-process results based on aggregate type
      if (identical(stat, "none") || (length(stat) == 1 && stat[1] == "none")) {
            # Pairs mode - res is a list. Normalization does not apply to
            # pairs mode (per-pair kernel weights are reported as-is); if
            # normalize = TRUE was requested with stat = "none", the
            # validation above would have warned about no-op stat anyway.
            out <- .emit_pairs_cpp(
                  res,
                  focal,
                  index$ref_data,
                  report_dist = TRUE,
                  geo_mode = index$coord_type,
                  x_cov = x_cov_mat,
                  values = values_mat,
                  values_names = values_names,
                  emit_sample_weight = emit_sample_weight,
                  emit_area_weight   = emit_area_weight,
                  emit_user_weight   = emit_user_weight
            )
            names(out) <- gsub("focal_", "", names(out))

            # Remap analog_index from stripped-pool indexing (what the C++
            # workers see) back to original-pool indexing (what the user
            # passed in). NULL pool_row_map means identity, no remap needed.
            # Preserve NA values, which mark focals with no analog.
            if (!is.null(index$pool_row_map) && !is.null(out$analog_index)) {
                  ai <- out$analog_index
                  na_mask <- is.na(ai)
                  if (any(!na_mask)) {
                        ai[!na_mask] <- index$pool_row_map[ai[!na_mask]]
                  }
                  out$analog_index <- ai
            }

            # Focal-side remapping. Two cases when focal NA stripping
            # happened:
            # (1) k == 1 && select != "all": exactly one row per (stripped)
            #     focal. Reconstruct to n_focal_original rows, with NA in
            #     all columns at stripped focal positions, so the result
            #     can rasterize back onto the user's input grid.
            # (2) Otherwise: variable rows per focal. We don't reconstruct
            #     (it's not rasterizable), but we still remap the `index`
            #     column from stripped-focal indexing back to original-
            #     focal indexing so user-side joins work.
            if (!is.null(focal_row_map)) {
                  # When focal NA stripping happened, two cases for the
                  # output:
                  # (1) Single row per focal (k == 1 && select != "all"):
                  #     reconstruct to n_focal_original rows with NA at
                  #     stripped focal positions. This matches the shape
                  #     users would get without NA stripping (one row per
                  #     input cell, with NA pair data for NA focals), and
                  #     enables rasterization for SpatRaster inputs.
                  # (2) Variable rows per focal (k > 1 or select == "all"):
                  #     output has variable rows per focal. We don't
                  #     reconstruct (it's not rasterizable and the row
                  #     count is naturally smaller without NA focals);
                  #     just remap the `index` column so the user's
                  #     existing-input row indexing still works.
                  single_row_per_focal <- isTRUE(k == 1) &&
                        !identical(select, "all")

                  if (single_row_per_focal) {
                        # Reconstruct to n_focal_original rows, NA at stripped
                        # positions. After this, `out$index` is 1:n_focal_original
                        # (matches user-input row order).
                        n_out <- n_focal_original
                        scat <- function(col, na_val) {
                              full <- rep(na_val, n_out)
                              full[focal_row_map] <- col
                              full
                        }
                        out_full <- list()
                        for (nm in names(out)) {
                              col <- out[[nm]]
                              if (identical(nm, "index")) {
                                    out_full[[nm]] <- seq_len(n_out)
                              } else if (is.integer(col)) {
                                    out_full[[nm]] <- scat(col, NA_integer_)
                              } else {
                                    out_full[[nm]] <- scat(col, NA_real_)
                              }
                        }
                        out <- as.data.frame(out_full,
                                             stringsAsFactors = FALSE)
                  } else {
                        # Multi-row-per-focal; just remap the `index` column.
                        if (!is.null(out$index)) {
                              idx_col <- out$index
                              na_mask <- is.na(idx_col)
                              if (any(!na_mask)) {
                                    idx_col[!na_mask] <- focal_row_map[idx_col[!na_mask]]
                              }
                              out$index <- idx_col
                        }
                  }
            }

            for (nm in names(cpp_attrs)) {
                  attr(out, nm) <- cpp_attrs[[nm]]
            }
            out <- .attach_index_res_attrs(out, index)

            # Pair-mode normalization is not currently defined; attach
            # normalize attribute for transparency but no D_max.
            attr(out, "normalize") <- isTRUE(normalize)

            return(.format_output(out, focal, stat, select,
                                  k, kernel_clim, kernel_geog, theta_clim, theta_geog, x_cov_mat,
                                  lambda, se, exclude_self,
                                  index$downsample_actual,
                                  max_clim = max_clim, max_geog = max_geog,
                                  cell_area_weight_applied = emit_area_weight,
                                  weight_provided = emit_user_weight))
      }

      # Aggregation mode(s)
      if (!is.matrix(res)) {
            stop("Internal error: expected matrix result from C++ for aggregation stats")
      }

      if (nrow(res) != nrow(focal)) {
            stop("Internal error: stat result rows do not match number of focals.")
      }

      # Reconstruct output at original focal length when focal NA stripping
      # was applied. C++ produced n_focal_used rows; we scatter them into
      # n_focal_original rows with NA at stripped positions. This keeps the
      # output aligned to the user's input shape (essential for raster
      # output, where row order corresponds to cell index).
      if (is.null(focal_row_map)) {
            n_out <- nrow(focal)
            res_full <- res
            x_full <- focal[, 1]
            y_full <- focal[, 2]
      } else {
            n_out <- n_focal_original
            res_full <- matrix(NA_real_, nrow = n_out, ncol = ncol(res))
            res_full[focal_row_map, ] <- res
            x_full <- rep(NA_real_, n_out)
            y_full <- rep(NA_real_, n_out)
            x_full[focal_row_map] <- focal[, 1]
            y_full[focal_row_map] <- focal[, 2]
      }

      # Build output data.frame
      out <- data.frame(
            index = seq_len(n_out),
            x = x_full,
            y = y_full
      )

      col_idx <- 1L
      n_vals <- if (is.null(values_mat)) 0L else ncol(values_mat)
      want_se <- (se != "none")

      # Track tabulate column names as we assemble them, so the
      # post-aggregation normalization step can find them without having
      # to reconstruct the naming convention.
      tabulate_col_names <- character(0)

      # Regular stats (count, sum_weights, mean_weights, ess): 1 col each
      regular_stats <- intersect(stat,
                                 c("count", "sum_weights", "mean_weights", "ess"))
      for (s in regular_stats) {
            out[[s]] <- res_full[, col_idx]
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
                        out[[col_name]] <- res_full[, col_idx]
                        col_idx <- col_idx + 1L

                        # weighted_mean SE is emitted immediately after
                        if (s == "weighted_mean" && want_se) {
                              se_name <- if (n_vals == 1L) "se_weighted_mean" else
                                    paste0("se_weighted_mean_", var_name)
                              out[[se_name]] <- res_full[, col_idx]
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
                        out[[col_name]] <- res_full[, col_idx]
                        col_idx <- col_idx + 1L

                        # Track for downstream normalization.
                        tabulate_col_names <- c(tabulate_col_names, col_name)
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
                        out[[paste0("coef_", base_names[d])]] <- res_full[, col_idx]
                        col_idx <- col_idx + 1
                  }
                  if (want_se) {
                        for (d in seq_len(reg_dim)) {
                              out[[paste0("se_", base_names[d])]] <- res_full[, col_idx]
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
                              out[[col_name]] <- res_full[, col_idx]
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
                                    out[[col_name]] <- res_full[, col_idx]
                                    col_idx <- col_idx + 1
                              }
                        }
                  }
            }
      }

      for (nm in names(cpp_attrs)) {
            attr(out, nm) <- cpp_attrs[[nm]]
      }
      out <- .attach_index_res_attrs(out, index)

      # Apply D_max normalization to sum_weights and tabulate columns,
      # if requested AND if a normalizable stat is present. The
      # `do_dmax` predicate (computed once above) is the single source
      # of truth for "did we actually normalize anything?" -- it covers
      # the silent-no-op case where the user passed normalize = TRUE
      # explicitly but the requested stat doesn't include any
      # normalizable column.
      if (do_dmax) {
            out <- .apply_dmax_normalization(out, D_max, stat, tabulate_col_names)
      }

      # Record normalization metadata on the output. We report what the
      # user *requested* (resolved logical) for the `normalize` attribute,
      # and the actual D_max used (NA when no normalization was applied).
      attr(out, "normalize") <- isTRUE(normalize)
      attr(out, "D_max")     <- if (do_dmax) D_max else NA_real_

      return(.format_output(out, focal, stat, select, k, kernel_clim, kernel_geog, theta_clim, theta_geog, x_cov_mat,
                            lambda, se, exclude_self, index$downsample_actual,
                            max_clim = max_clim, max_geog = max_geog,
                            cell_area_weight_applied = emit_area_weight,
                            weight_provided = emit_user_weight))
}


#' Internal helper: attach index resolution parameters as result attributes
#'
#' Copies the lattice resolution knobs and realized bin layout from the
#' `analog_index` object onto a query result, so they surface via
#' [metadata()]. These are properties of the index (set at build time),
#' reused here rather than recomputed. Missing fields (e.g. on legacy index
#' objects) are simply skipped.
#'
#' @keywords internal
.attach_index_res_attrs <- function(out, index) {
      for (nm in c("geog_res_adj", "clim_res_adj", "geo_target", "clim_target", "bins_per_axis")) {
            val <- index[[nm]]
            if (!is.null(val)) attr(out, nm) <- val
      }
      out
}


#' Internal helper: Query C++ with optional chunking and progress tracking
#'
#' Wraps query_analog_index_cpp with optional chunking for progress bars.
#' Handles merging of chunked results.
#'
#' @keywords internal
.query_cpp_chunked <- function(index, focal, ref_mm, k, max_clim, max_geog,
                               select_code, aggregate_codes,
                               clim_kernel_code, geog_kernel_code,
                               theta_clim, theta_geog,
                               x_cov, y, covariates, lambda, se_code,
                               n_classes_per_var = integer(0),
                               area_weight = NULL,
                               user_weight = NULL,
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
                  clim_kernel_code = clim_kernel_code,
                  geog_kernel_code = geog_kernel_code,
                  theta_clim = theta_clim,
                  theta_geog = theta_geog,
                  x_cov_sexp = x_cov,
                  values_sexp = y,
                  covariates_sexp = covariates,
                  lambda = lambda,
                  se_code = se_code,
                  n_classes_per_var = n_classes_per_var,
                  area_weight_sexp = area_weight,
                  user_weight_sexp = user_weight,
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
                  clim_kernel_code = clim_kernel_code,
                  geog_kernel_code = geog_kernel_code,
                  theta_clim = theta_clim,
                  theta_geog = theta_geog,
                  x_cov_sexp = x_cov_chunk,
                  values_sexp = y,
                  covariates_sexp = covariates,
                  lambda = lambda,
                  se_code = se_code,
                  n_classes_per_var = n_classes_per_var,
                  area_weight_sexp = area_weight,
                  user_weight_sexp = user_weight,
                  exclude_self = FALSE   # always FALSE in chunked path
            )

            utils::setTxtProgressBar(pb, i)
      }

      # Merge results based on type
      first <- results[[1]]

      if (is.list(first) && !is.matrix(first)) {
            # Pairs mode. The C++ payload has four parallel sublists
            #   indices, sample_weights, area_weights, user_weights
            # each of length n_focal_chunk, concatenated end-to-end across
            # chunks. Plus two scalar boolean flags (has_area_weight,
            # has_user_weight) carried through unchanged from the first chunk.
            sublist_names <- c("indices", "sample_weights",
                               "area_weights", "user_weights")
            scalar_names  <- c("has_area_weight", "has_user_weight")

            merged <- list()
            for (nm in sublist_names) {
                  merged[[nm]] <- do.call(c, lapply(results, function(r) r[[nm]]))
            }
            for (nm in scalar_names) {
                  merged[[nm]] <- first[[nm]]
            }

            for (nm in names(attributes(first))) {
                  if (!nm %in% c("names")) {
                        attr(merged, nm) <- attr(first, nm)
                  }
            }
            attr(merged, "n_x") <- nrow(focal)

            return(merged)
      } else {
            # Aggregation mode: rbind matrices
            merged <- do.call(rbind, results)

            for (nm in names(attributes(first))) {
                  if (!nm %in% c("dim", "dimnames")) {
                        attr(merged, nm) <- attr(first, nm)
                  }
            }
            attr(merged, "n_x") <- nrow(focal)

            return(merged)
      }
}
