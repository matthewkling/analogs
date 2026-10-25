# dev/fft_bench.R
#
# Developer-only benchmark of the FFT method against the exact engine. NOT part
# of the built package (see .Rbuildignore).
#
# For each scenario, times `method = "fft"` on all focals and the exact engine
# on a random subset of focals, extrapolating the exact time linearly (exact
# cost is proportional to the number of focals), and reports the FFT error on
# that subset. Scenarios use raster-like landscapes: a regular 1 km grid with
# smooth climate gradients plus noise, queried against a uniformly shifted
# "future" climate.
#
# Usage (from the package root, after installing or load_all()):
#
#   source("dev/fft_bench.R")
#   fft_bench()                      # default scenarios
#   fft_bench(sides = c(250, 500))   # smaller landscapes only
#   fft_bench(n_threads = 8)
#
# Note: the first call in a session also pays terra's load time; a warm-up
# query is run before timing.


# Simulated landscape: side x side cells at 1 km, two climate variables
.bench_landscape <- function(side, seed = 1) {
      set.seed(seed)
      g <- expand.grid(x = seq_len(side) - 0.5, y = seq_len(side) - 0.5)
      e1 <- 2 * sin(g$x / 150) + g$y / 200 + stats::rnorm(nrow(g), sd = 0.2)
      e2 <- cos(g$y / 120) + g$x / 300 + stats::rnorm(nrow(g), sd = 0.2)
      cbind(x = g$x, y = g$y, e1 = e1, e2 = e2)
}

# Default scenarios: name -> list(geog, env, stat)
.bench_scenarios <- function() {
      list(
            "climate only (gaussian)" = list(
                  geog = NULL,
                  env  = kernel("gaussian", theta = 0.3, max = 1),
                  stat = c("weighted_mean", "sum_weights")),
            "AIM: gaussian climate x 100 km disk" = list(
                  geog = kernel(max = 100),
                  env  = kernel("gaussian", theta = 0.3, max = 1),
                  stat = c("weighted_mean", "sum_weights")),
            "gaussian climate x gaussian geog" = list(
                  geog = kernel("gaussian", theta = 50, max = 150),
                  env  = kernel("gaussian", theta = 0.3, max = 1),
                  stat = c("weighted_mean", "sum_weights", "ess")),
            "regression, climate x geog" = list(
                  geog = kernel("gaussian", theta = 50, max = 150),
                  env  = kernel("gaussian", theta = 0.3, max = 1),
                  stat = c("regression", "count"))
      )
}

#' Benchmark the FFT method against the exact engine
#'
#' @param sides Landscape side lengths (km); each gives side^2 pool cells and
#'   as many focals.
#' @param scenarios Named list of scenarios (see `.bench_scenarios()`).
#' @param n_exact Number of focals the exact engine is timed on.
#' @param method The FFT method to benchmark (`"fft"` or an `fft_method()`).
#' @param n_threads Threads for both engines.
#' @return A data.frame of results (also printed as it goes).
fft_bench <- function(sides = c(250, 500, 1000),
                      scenarios = .bench_scenarios(),
                      n_exact = 5000, method = "fft", n_threads = NULL) {

      # Warm-up (loads terra, initializes the thread pool)
      invisible(analog_search(matrix(stats::rnorm(40), 10),
                              matrix(stats::rnorm(400), 100),
                              stat = "count", env = kernel(max = 1),
                              coord_type = "projected"))

      out <- list()
      for (side in sides) {
            pool <- .bench_landscape(side)
            fut <- pool
            fut[, "e1"] <- fut[, "e1"] + 0.7
            fut[, "e2"] <- fut[, "e2"] - 0.3
            y <- pool[, "e1"] + pool[, "e2"] + stats::rnorm(nrow(pool))
            covs <- cbind(c1 = stats::rnorm(nrow(pool)))
            rows <- sort(sample.int(nrow(fut), min(n_exact, nrow(fut))))

            for (nm in names(scenarios)) {
                  sc <- scenarios[[nm]]
                  args <- list(pool = pool, y = y, geog = sc$geog, env = sc$env,
                               stat = sc$stat, coord_type = "projected",
                               n_threads = n_threads)
                  if ("regression" %in% sc$stat) args$covariates <- covs

                  t_fft <- tryCatch(
                        system.time(ff <- do.call(analog_search,
                                                  c(list(x = fut), args,
                                                    list(method = method))))[["elapsed"]],
                        error = function(e) {
                              message(nm, " (", side, " km): ", conditionMessage(e))
                              NA_real_
                        })
                  if (is.na(t_fft)) next
                  t_exact <- system.time(
                        ex <- do.call(analog_search,
                                      c(list(x = fut[rows, , drop = FALSE]), args,
                                        list(method = "exact")))
                  )[["elapsed"]] * nrow(fut) / length(rows)

                  # Error of the main output column, relative to its spread
                  col <- if ("regression" %in% sc$stat) "coef_intercept" else "weighted_mean"
                  e <- ex[[col]]
                  f <- ff[[col]][rows]
                  ok <- is.finite(e) & is.finite(f)
                  err <- stats::median(abs(f[ok] - e[ok])) / stats::sd(e[ok])

                  grid <- attr(ff, "fft_grid")
                  res <- data.frame(
                        scenario = nm, cells = nrow(pool),
                        grid = paste(grid$n_nodes, collapse = "x"),
                        exact_s = round(t_exact, 1), fft_s = round(t_fft, 2),
                        speedup = round(t_exact / t_fft),
                        median_err_sd = signif(err, 2),
                        stringsAsFactors = FALSE)
                  print(res, row.names = FALSE)
                  out[[length(out) + 1L]] <- res
            }
      }
      invisible(do.call(rbind, out))
}
