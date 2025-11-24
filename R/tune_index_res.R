#' Tune Index Resolution
#'
#' Benchmarks different lattice index resolutions to find the optimal setting
#' for your data and query pattern. This function runs test queries with various
#' resolutions and recommends the one with the best performance.
#'
#' @inheritParams find_analogs
#'
#' @param resolutions Integer vector of resolutions to test. Default tests
#'   a range from 4 to 32 bins per dimension.
#' @param n_sample Number of focal points to use for benchmarking. If \code{x}
#'   has more rows than this, a random sample will be taken. Default is 500.
#' @param n_reps Number of repetitions for each resolution test. Default is 3.
#' @param verbose Logical; if TRUE, print progress and results. Default is FALSE.
#'
#' @return An integer giving the recommended index resolution (bins per dimension).
#'   Invisibly returns a data.frame with detailed benchmark results including
#'   columns: \code{resolution}, \code{mean_time_ms}, \code{sd_time_ms}.
#'
#' @details
#' The function works by:
#' \enumerate{
#'   \item Taking a sample of focal points (or using all if small dataset)
#'   \item Building an index at each test resolution
#'   \item Running the specified query multiple times
#'   \item Measuring elapsed time for each configuration
#'   \item Recommending the resolution with lowest mean time
#' }
#'
#' The optimal resolution depends on:
#' \itemize{
#'   \item Size of reference dataset (larger → higher resolution helpful)
#'   \item Dimensionality (more climate variables → consider lower resolution)
#'   \item Query selectivity (tight constraints → higher resolution helpful)
#'   \item Hardware (memory/cache considerations)
#' }
#'
#' For most applications, testing resolutions from 8 to 24 is sufficient.
#' Very large datasets (>100k points) may benefit from resolutions up to 32.
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
#'
#' # Test specific resolutions
#' res <- tune_index_res(
#'   x = sites,
#'   pool = climate_data,
#'   mode = "knn_clim",
#'   max_geog = 100,
#'   k = 20,
#'   resolutions = c(8, 12, 16, 20, 24)
#' )
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
                           coord_type = c("auto", "lonlat", "projected"),
                           resolutions = NULL,
                           n_sample = NULL,
                           n_reps = NULL,
                           verbose = FALSE) {

      # Defaults
      resolutions <- resolutions %||% c(4, 8, 12, 16, 20, 24, 32)
      n_sample <- n_sample %||% 500
      n_reps <- n_reps %||% 3

      # Validate inputs
      coord_type <- match.arg(coord_type)
      mode <- match.arg(mode)

      if (!is.numeric(resolutions) || any(resolutions <= 0)) {
            stop("resolutions must be a vector of positive integers")
      }
      resolutions <- as.integer(resolutions)
      resolutions <- sort(unique(resolutions))

      if (!is.numeric(n_sample) || length(n_sample) != 1 || n_sample <= 0) {
            stop("n_sample must be a positive integer")
      }
      n_sample <- as.integer(n_sample)

      if (!is.numeric(n_reps) || length(n_reps) != 1 || n_reps <= 0) {
            stop("n_reps must be a positive integer")
      }
      n_reps <- as.integer(n_reps)

      # Format data
      x_mm <- .format_data(x)
      pool_mm <- .format_data(pool)

      # Sample if needed
      n_x <- nrow(x_mm)
      if (n_x > n_sample) {
            idx <- sample.int(n_x, n_sample)
            x_mm <- x_mm[idx, , drop = FALSE]
            if (verbose) {
                  message(sprintf("Using %d of %d focal points for benchmarking",
                                  n_sample, n_x))
            }
      } else {
            if (verbose) {
                  message(sprintf("Using all %d focal points for benchmarking", n_x))
            }
      }

      if (verbose) {
            message(sprintf("\nTesting %d resolutions with %d repetitions each...\n",
                            length(resolutions), n_reps))
      }

      # Storage for results
      results <- data.frame(
            resolution = integer(),
            rep = integer(),
            time_ms = numeric(),
            stringsAsFactors = FALSE
      )

      # Benchmark each resolution
      for (res in resolutions) {

            if (verbose) {
                  cat(sprintf("  Resolution %2d: ", res))
            }

            # Build index once for this resolution
            build_start <- Sys.time()
            index <- build_analog_index(
                  pool = pool_mm,
                  coord_type = coord_type,
                  index_res = res
            )
            build_time <- as.numeric(difftime(Sys.time(), build_start, units = "secs"))

            # Run query n_reps times
            times <- numeric(n_reps)
            for (rep in seq_len(n_reps)) {
                  query_start <- Sys.time()

                  result <- find_analogs(
                        x = x_mm,
                        pool = index,
                        mode = mode,
                        max_clim = max_clim,
                        max_geog = max_geog,
                        k = k,
                        weight = weight,
                        theta = theta,
                        report_dist = FALSE  # Faster without distances
                  )

                  query_time <- as.numeric(difftime(Sys.time(), query_start, units = "secs"))
                  times[rep] <- query_time * 1000  # Convert to ms
            }

            # Store results
            for (rep in seq_len(n_reps)) {
                  results <- rbind(results, data.frame(
                        resolution = res,
                        rep = rep,
                        time_ms = times[rep],
                        stringsAsFactors = FALSE
                  ))
            }

            if (verbose) {
                  cat(sprintf("build=%.2fs  query=%.1f±%.1f ms\n",
                              build_time,
                              mean(times),
                              sd(times)))
            }
      }

      # Summarize by resolution
      summary_stats <- aggregate(
            time_ms ~ resolution,
            data = results,
            FUN = function(x) c(mean = mean(x), sd = sd(x))
      )

      summary_df <- data.frame(
            resolution = summary_stats$resolution,
            mean_time_ms = summary_stats$time_ms[, "mean"],
            sd_time_ms = summary_stats$time_ms[, "sd"],
            stringsAsFactors = FALSE
      )

      # Find best resolution
      best_idx <- which.min(summary_df$mean_time_ms)
      best_res <- summary_df$resolution[best_idx]
      best_time <- summary_df$mean_time_ms[best_idx]

      if (verbose) {
            message("\n--- Summary ---")
            message(sprintf("Optimal resolution: %d bins per dimension", best_res))
            message(sprintf("Query time: %.1f ms (mean over %d reps)",
                            best_time, n_reps))

            # Show relative performance
            message("\nRelative performance:")
            for (i in seq_len(nrow(summary_df))) {
                  rel_perf <- summary_df$mean_time_ms[i] / best_time
                  marker <- if (i == best_idx) " ← optimal" else ""
                  message(sprintf("  Resolution %2d: %.2fx slower%s",
                                  summary_df$resolution[i],
                                  rel_perf,
                                  marker))
            }
      }

      # Return recommended resolution
      invisible(summary_df)
      return(best_res)
}
