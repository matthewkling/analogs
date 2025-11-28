#' Profile analog_search Performance
#'
#' Runs analog_search with detailed C++ profiling instrumentation and returns
#' timing breakdowns for optimization analysis.
#'
#' @param ... Arguments passed to analog_search
#' @param report Logical; if TRUE, print formatted report
#' @param plot Logical; if TRUE, create visualization
#'
#' @return List containing:
#'   - result: the normal analog_search output
#'   - profile: data.frame with timing details
#'   - counters: data.frame with operation counts
#'
#' @export
profile_analogs <- function(..., report = TRUE, plot = TRUE) {

      # Capture arguments
      args <- list(...)

      # Call profiling version
      prof_result <- .Call(
            `_analogs_profile_find_analogs`,
            args$focal_mm,
            args$ref_mm,
            args$k,
            args$max_clim,
            args$max_geog,
            args$geo_mode,
            args$mode_code,
            args$weight_code,
            args$theta,
            args$lattice_res,
            TRUE  # enable_profiling
      )

      # Extract and format profiling data
      prof_data <- prof_result$profile

      profile_df <- data.frame(
            event = prof_data$event_names,
            time_ms = prof_data$event_times_ms,
            count = prof_data$event_counts,
            stringsAsFactors = FALSE
      )
      profile_df$time_per_call_ms <- profile_df$time_ms / pmax(profile_df$count, 1)
      profile_df$pct_total <- 100 * profile_df$time_ms / sum(profile_df$time_ms, na.rm = TRUE)
      profile_df <- profile_df[order(-profile_df$time_ms), ]

      counters_df <- data.frame(
            counter = prof_data$counter_names,
            value = prof_data$counter_values,
            stringsAsFactors = FALSE
      )
      counters_df <- counters_df[order(-counters_df$value), ]

      # Print report
      if (report) {
            cat("\n=== PERFORMANCE PROFILE ===\n\n")

            cat("Top time consumers:\n")
            print(head(profile_df[, c("event", "time_ms", "pct_total", "count")], 10))

            cat("\n\nOperation counts:\n")
            print(head(counters_df, 15))

            # Compute derived metrics
            n_focal <- counters_df$value[counters_df$counter == "focal_points_processed"]
            if (length(n_focal) > 0 && n_focal > 0) {
                  candidates_tested <- counters_df$value[counters_df$counter == "candidates_tested"]
                  candidates_passed <- counters_df$value[counters_df$counter == "candidates_passed"]

                  cat("\n\nDerived metrics:\n")
                  cat(sprintf("  Avg candidates tested per focal: %.1f\n",
                              candidates_tested / n_focal))
                  cat(sprintf("  Candidate pass rate: %.1f%%\n",
                              100 * candidates_passed / candidates_tested))

                  geog_tests <- counters_df$value[counters_df$counter == "geog_tests"]
                  clim_tests <- counters_df$value[counters_df$counter == "clim_tests"]
                  if (length(geog_tests) > 0) {
                        cat(sprintf("  Geog tests per focal: %.1f\n", geog_tests / n_focal))
                  }
                  if (length(clim_tests) > 0) {
                        cat(sprintf("  Clim tests per focal: %.1f\n", clim_tests / n_focal))
                  }
            }

            cat("\n")
      }

      # Create visualization
      if (plot && requireNamespace("ggplot2", quietly = TRUE)) {
            library(ggplot2)

            top_events <- head(profile_df, 12)

            p <- ggplot(top_events, aes(x = reorder(event, time_ms), y = time_ms)) +
                  geom_col(fill = "steelblue") +
                  geom_text(aes(label = sprintf("%.1f%%", pct_total)),
                            hjust = -0.1, size = 3) +
                  coord_flip() +
                  labs(
                        title = "Performance Profile: Time by Operation",
                        x = NULL,
                        y = "Time (ms)"
                  ) +
                  theme_minimal() +
                  theme(
                        panel.grid.major.y = element_blank(),
                        plot.title = element_text(face = "bold")
                  )

            print(p)
      }

      # Return structured results
      invisible(list(
            result = prof_result$result,
            profile = profile_df,
            counters = counters_df
      ))
}


#' Compare Performance Across Configurations
#'
#' Runs multiple profiling tests with different settings to compare performance.
#'
#' @param focal Focal dataset
#' @param ref Reference dataset
#' @param configs List of configuration lists, each containing analog_search parameters
#' @param labels Character vector of labels for each config
#'
#' @return data.frame comparing performance metrics
#'
#' @export
profile_compare <- function(focal, ref, configs, labels = NULL) {

      if (is.null(labels)) {
            labels <- paste0("config_", seq_along(configs))
      }

      results <- list()

      for (i in seq_along(configs)) {
            cat(sprintf("\nTesting %s...\n", labels[i]))

            config <- configs[[i]]
            config$focal <- focal
            config$ref <- ref
            config$report <- FALSE
            config$plot <- FALSE

            results[[i]] <- do.call(profile_analogs, config)
      }

      # Extract key metrics
      comparison <- data.frame(
            config = labels,
            total_time_ms = sapply(results, function(r) sum(r$profile$time_ms)),
            lattice_build_ms = sapply(results, function(r) {
                  idx <- which(r$profile$event == "lattice_build")
                  if (length(idx) > 0) r$profile$time_ms[idx] else NA
            }),
            worker_time_ms = sapply(results, function(r) {
                  idx <- which(r$profile$event == "parallel_worker_execution")
                  if (length(idx) > 0) r$profile$time_ms[idx] else NA
            }),
            candidates_per_focal = sapply(results, function(r) {
                  ct <- r$counters$value[r$counters$counter == "candidates_tested"]
                  nf <- r$counters$value[r$counters$counter == "focal_points_processed"]
                  if (length(ct) > 0 && length(nf) > 0) ct / nf else NA
            }),
            stringsAsFactors = FALSE
      )

      comparison$speedup <- comparison$total_time_ms[1] / comparison$total_time_ms

      print(comparison)

      invisible(list(
            comparison = comparison,
            detailed_results = results
      ))
}


#' Run Comprehensive Profiling Suite
#'
#' Tests different problem configurations to identify bottlenecks
#'
#' @export
profile_suite <- function() {

      cat("=== COMPREHENSIVE PROFILING SUITE ===\n\n")

      # Generate test data
      cat("Generating test data...\n")
      set.seed(123)

      # Small test
      focal_small <- matrix(rnorm(100 * 4), ncol = 4)
      ref_small <- matrix(rnorm(1000 * 4), ncol = 4)

      # Medium test
      focal_med <- matrix(rnorm(1000 * 4), ncol = 4)
      ref_med <- matrix(rnorm(10000 * 4), ncol = 4)

      # Large test (if memory allows)
      # focal_large <- matrix(rnorm(10000 * 4), ncol = 4)
      # ref_large <- matrix(rnorm(100000 * 4), ncol = 4)

      configs <- list(
            list(mode = "knn_geog", max_clim = 0.5, k = 1, lattice_res = 10),
            list(mode = "knn_geog", max_clim = 0.5, k = 1, lattice_res = 20),
            list(mode = "knn_clim", max_geog = 100, k = 5, lattice_res = 10),
            list(mode = "count", max_clim = 0.5, max_geog = 100, lattice_res = 10),
            list(mode = "sum", max_clim = 0.5, max_geog = 100,
                 weight = "inverse_clim", theta = 1e-6, lattice_res = 10)
      )

      labels <- c(
            "velocity_res10",
            "velocity_res20",
            "impact_k5",
            "availability",
            "intensity"
      )

      cat("\n--- Testing medium dataset (1k focal, 10k ref) ---\n")
      results_med <- profile_compare(focal_med, ref_med, configs, labels)

      return(invisible(results_med))
}
