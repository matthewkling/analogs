#' Tune Index Resolution
#'
#' Automatically finds the optimal lattice index resolution for your data and
#' query pattern using adaptive bracketing search. Runs test queries with
#' different resolutions and recommends the one with the fastest compute speed.
#'
#' @inheritParams analog_search
#'
#' @param default_res Default resolution to use as starting point for search.
#'   Default is 16.
#' @param verbose Logical; if TRUE, print the selected resolution. Default is FALSE.
#'
#' @return An integer giving the recommended index resolution (bins per dimension).
#'
#' @details
#' The function uses an adaptive bracketing algorithm:
#' \enumerate{
#'   \item Starts with three resolutions: default/2, default, default*2
#'   \item Evaluates elapsed time for each
#'   \item If minimum is at an edge, expands search in that direction
#'   \item Returns resolution with lowest elapsed time
#' }
#'
#' This typically requires only 3-5 query evaluations total, making it much
#' faster than exhaustive grid search.
#'
#' The function only performs tuning for non-trivial problem sizes (>2000 focal
#' points). For smaller datasets, it returns the default resolution.
#'
#' A subsample of focal points is used for benchmarking to keep tuning fast
#' while still being representative of actual query performance.
#'
#' @examples
#' \dontrun{
#' # Find optimal resolution for velocity queries
#' optimal_res <- tune_index_res(
#'   x = sample_sites,
#'   pool = climate_data,
#'   mode = "knn_geog",
#'   max_clim = 0.5,
#'   k = 1
#' )
#'
#' # Use the optimized resolution
#' index <- build_analog_index(climate_data, index_res = optimal_res)
#' }
#'
#' @export
tune_index_res <- function(x,
                           pool,
                           mode = c("knn_clim", "knn_geog", "count", "sum", "mean", "all"),
                           max_clim = NULL,
                           max_geog = NULL,
                           k = NULL,
                           weight = NULL,
                           theta = NULL,
                           x_cov = NULL,
                           coord_type = c("auto", "lonlat", "projected"),
                           n_threads = NULL,
                           default_res = 16L,
                           verbose = FALSE) {

      # Helper: detect monotonic timings with a tolerance
      is_strict_monotonic <- function(x, tol = 0.15) {
            if (length(x) < 3L) return(FALSE)
            dx <- diff(x)
            # strictly decreasing or increasing, with relative tolerance
            all(dx <= -tol * abs(x[-length(x)])) ||
                  all(dx >=  tol * abs(x[-length(x)]))
      }

      # Validate and normalize query parameters
      # Note: We don't have focal_mm yet, but .validate_query_params handles x_cov=NULL gracefully
      # For x_cov validation, we'll need to format focal data first
      focal_mm <- .format_data(x)

      params <- .validate_query_params(mode, k, weight, theta, x_cov, focal_mm)
      mode <- params$mode
      k <- params$k
      weight <- params$weight
      theta <- params$theta
      x_cov_mat <- params$x_cov  # Already validated if provided

      # Validate coord_type
      coord_type <- match.arg(coord_type)

      n_focal <- nrow(focal_mm)

      # Only tune for non-trivial problem sizes
      if (n_focal <= 2000L) {
            if (verbose) {
                  message("Dataset too small for tuning (n=", n_focal,
                          "); using default resolution of ", default_res, ".")
            }
            return(as.integer(default_res))
      }

      # Subsample focal sites for faster benchmarking
      n_samp <- min(1000L, max(100L, as.integer(n_focal * 0.01)))
      idx    <- sample.int(n_focal, n_samp)
      focal_mm_samp <- focal_mm[idx, , drop = FALSE]

      # Subsample x_cov if provided
      x_cov_samp <- NULL
      if (!is.null(x_cov_mat)) {
            x_cov_samp <- x_cov_mat[idx, , drop = FALSE]
      }

      # Helper to evaluate timing for a given index_res
      eval_time <- function(r) {
            r <- as.integer(max(1L, r))  # ensure valid

            # Build index with this resolution
            index <- build_analog_index(
                  pool = pool,
                  coord_type = coord_type,
                  index_res = r
            )

            # Time the query
            st <- system.time({
                  result <- query_analog_index(
                        x = focal_mm_samp,
                        index = index,
                        mode = mode,
                        max_clim = max_clim,
                        max_geog = max_geog,
                        x_cov = x_cov_samp,
                        k = k,
                        weight = weight,
                        theta = theta,
                        report_dist = FALSE,  # Faster without distances
                        n_threads = n_threads
                  )
            })

            st[["elapsed"]]
      }

      # Start from heuristic center
      r0 <- as.integer(default_res)
      r_vals <- c(r0 %/% 2L, r0, r0 * 2L)
      r_vals <- unique(r_vals[r_vals > 0L])

      # Evaluate initial bracket
      if (verbose) {
            message("Evaluating initial bracket: {", paste(r_vals, collapse = ", "), "}")
      }
      times <- vapply(r_vals, eval_time, numeric(1))
      if (verbose) {
            message("  Times: {", paste(sprintf("%.3f", times), collapse = ", "), "} sec")
      }

      # One-step adaptive refinement
      best_idx <- which.min(times)
      best_r   <- r_vals[best_idx]

      if (best_idx == 1L && length(r_vals) > 1L) {
            # best is smallest → try halving again
            r_try <- max(1L, best_r %/% 2L)
            if (verbose) {
                  message("Minimum at lower edge; trying r=", r_try)
            }
            t_try <- eval_time(r_try)
            if (verbose) {
                  message("  Time: ", sprintf("%.3f", t_try), " sec")
            }
            r_vals <- c(r_try, r_vals)
            times  <- c(t_try, times)

      } else if (best_idx == length(r_vals) && length(r_vals) > 1L) {
            # best is largest → try doubling again
            r_try <- best_r * 2L
            if (verbose) {
                  message("Minimum at upper edge; trying r=", r_try)
            }
            t_try <- eval_time(r_try)
            if (verbose) {
                  message("  Time: ", sprintf("%.3f", t_try), " sec")
            }
            r_vals <- c(r_vals, r_try)
            times  <- c(times, t_try)
      }

      # Recompute best after any refinement
      best_idx <- which.min(times)
      best_r   <- r_vals[best_idx]

      # Sort evaluated points by r for monotonicity check
      o      <- order(r_vals)
      r_eval <- r_vals[o]
      t_eval <- times[o]

      # Check if we're in a k=1 velocity scenario (where tuning matters less)
      is_k1_velocity <- (mode == "knn_geog" &&
                               !is.null(k) &&
                               is.numeric(k) &&
                               k == 1L)

      # Optional: warn if timings are monotonic in a meaningful regime
      if (length(t_eval) >= 4L &&
          n_samp >= 300L &&
          !is_k1_velocity &&
          is_strict_monotonic(t_eval)) {

            warning(
                  "Auto-tuning of index_res did not detect an interior minimum; ",
                  "elapsed times were monotonic across tested values {",
                  paste(r_eval, collapse = ", "),
                  "}. The optimal index_res may lie outside this range. ",
                  "Consider manually specifying index_res."
            )
      }

      r <- as.integer(best_r)
      if (verbose) {
            message("\nSelected resolution: ", r, " (", sprintf("%.3f", times[best_idx]), " sec)")
      }

      return(r)
}
