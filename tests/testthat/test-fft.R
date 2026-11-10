# Tests for the binned-convolution (FFT) method. The FFT method approximates
# the exact engine, so comparisons use tolerances scaled to each column and
# restrict to focals with reasonably populated neighborhoods (sparse
# neighborhoods amplify binning error in ratio statistics like regression).

fft_test_data <- function(seed = 1, n = 3000, n_focal = 400, n_env = 2) {
      set.seed(seed)
      pool <- cbind(x = runif(n, 0, 100), y = runif(n, 0, 100),
                    matrix(rnorm(n * n_env), n,
                           dimnames = list(NULL, paste0("e", seq_len(n_env)))))
      focal <- pool[sample(n, n_focal), ]
      focal[, -(1:2)] <- focal[, -(1:2)] + rnorm(n_focal * n_env, sd = 0.3)
      covs <- cbind(c1 = rnorm(n), c2 = pool[, 3] + rnorm(n, sd = 0.5))
      y <- cbind(a = sin(pool[, 3]) + pool[, 4] + covs[, 1] + rnorm(n, sd = 0.3),
                 b = pool[, 3]^2 + rnorm(n))
      list(pool = pool, focal = focal, covs = covs, y = y)
}

# 90th-percentile error over well-populated focals, relative to the column's
# mean absolute value. (A high quantile rather than the max: a few focals whose
# weight is dominated by two or three analogs swing with any perturbation of
# their weights, binning included.)
fft_rel_err <- function(exact, approx, col, ok) {
      e <- exact[[col]][ok]
      f <- approx[[col]][ok]
      unname(stats::quantile(abs(f - e), 0.9)) / mean(abs(e))
}

run_both <- function(..., bin = 0.1) {
      list(exact = analog_search(..., method = "exact"),
           fft   = analog_search(..., method = fft_method(bin = bin)))
}


test_that("fft_method() validates its arguments and method resolves", {
      expect_s3_class(fft_method(), "analog_fft_method")
      expect_error(fft_method(bin = 0), "bin")
      expect_error(fft_method(bin = 2), "bin")
      expect_error(fft_method(check = -1), "check")
      expect_error(fft_method(max_cells = 0), "max_cells")
      expect_null(.resolve_method("exact"))
      expect_s3_class(.resolve_method("fft"), "analog_fft_method")
      expect_error(.resolve_method("fast"), "method")
})


test_that("fft matches exact for weighted aggregation stats", {
      d <- fft_test_data()
      r <- run_both(d$focal, d$pool, y = d$y, weight = runif(nrow(d$pool), 0.5, 2),
                    env = kernel("gaussian", theta = 0.4, max = 1.5),
                    stat = c("count", "sum_weights", "ess", "weighted_mean"),
                    se = "ess", coord_type = "projected")

      expect_identical(names(r$fft), names(r$exact))
      expect_identical(attr(r$fft, "binning_method"), "fft")

      ok <- r$exact$count >= 30
      expect_true(sum(ok) > 100)
      for (col in c("sum_weights", "ess", "weighted_mean_a", "weighted_mean_b",
                    "se_weighted_mean_a")) {
            expect_lt(fft_rel_err(r$exact, r$fft, col, ok), 0.02)
      }
      # count has hard kernel edges, so its error shrinks more slowly
      expect_lt(fft_rel_err(r$exact, r$fft, "count", ok), 0.05)
})


test_that("fft matches exact for regression, including missing responses", {
      d <- fft_test_data()
      y <- d$y
      y[sample(nrow(y), 200), "b"] <- NA
      r <- run_both(d$focal, d$pool, y = y, covariates = d$covs,
                    env = kernel("gaussian", theta = 0.4, max = 1.5),
                    stat = c("count", "regression"), se = "ess", lambda = 1,
                    coord_type = "projected")

      expect_identical(names(r$fft), names(r$exact))
      ok <- r$exact$count >= 50
      for (col in c("coef_intercept_a", "coef_c1_a", "coef_intercept_b",
                    "coef_c2_b", "se_intercept_a", "se_c1_b")) {
            expect_lt(fft_rel_err(r$exact, r$fft, col, ok), 0.03)
      }
})


test_that("fft matches exact for tabulate with a per-variable box limit", {
      d <- fft_test_data()
      cls <- sample(c("p", "q", "r"), nrow(d$pool), replace = TRUE)
      r <- run_both(d$focal, d$pool, y = cls,
                    env = kernel("inverse", theta = 0.4, max = c(1, 1.5)),
                    stat = c("sum_weights", "tabulate"), coord_type = "projected")

      expect_identical(names(r$fft), names(r$exact))
      ok <- r$exact$sum_weights > stats::quantile(r$exact$sum_weights, 0.1)
      for (col in c("n_p", "n_q", "n_r")) {
            expect_lt(fft_rel_err(r$exact, r$fft, col, ok), 0.05)
      }
})


test_that("fft leave-one-out removes each focal's own contribution", {
      d <- fft_test_data(n = 2000)
      y <- d$pool[, 3] + rnorm(nrow(d$pool), sd = 2)
      args <- list(d$pool, d$pool, y = y,
                   env = kernel("gaussian", theta = 0.2, max = 0.8),
                   stat = c("count", "weighted_mean"), coord_type = "projected")
      ex   <- do.call(analog_search, c(args, exclude_self = TRUE, method = "exact"))
      self <- do.call(analog_search, c(args, exclude_self = FALSE, method = "exact"))
      ff   <- do.call(analog_search, c(args, exclude_self = TRUE,
                                       method = list(fft_method(bin = 0.1))))

      ok <- ex$count >= 20
      err <- median(abs(ff$weighted_mean[ok] - ex$weighted_mean[ok]))
      loo_effect <- median(abs(self$weighted_mean[ok] - ex$weighted_mean[ok]))
      expect_lt(err, loo_effect / 10)
})


test_that("fft binning error shrinks with finer bins", {
      d <- fft_test_data(n = 2000)
      args <- list(d$focal, d$pool, y = d$y[, "a"],
                   env = kernel("gaussian", theta = 0.4, max = 1.5),
                   stat = c("count", "weighted_mean"), coord_type = "projected")
      ex <- do.call(analog_search, c(args, method = "exact"))
      ok <- ex$count >= 30
      err <- vapply(c(0.4, 0.1), function(b) {
            ff <- do.call(analog_search, c(args, method = list(fft_method(bin = b))))
            median(abs(ff$weighted_mean[ok] - ex$weighted_mean[ok]))
      }, numeric(1))
      # Smooth kernels: error ~ bin^2, so a 4x finer grid should cut it well
      # below half
      expect_lt(err[2], err[1] / 4)
})


test_that("fft handles three environmental variables", {
      d <- fft_test_data(n = 2000, n_focal = 200, n_env = 3)
      r <- run_both(d$focal, d$pool, y = d$pool[, 3] + d$pool[, 5],
                    env = kernel("gaussian", theta = 0.5),
                    stat = c("sum_weights", "weighted_mean"),
                    coord_type = "projected", bin = 0.2)
      ok <- r$exact$sum_weights > stats::quantile(r$exact$sum_weights, 0.1)
      expect_lt(fft_rel_err(r$exact, r$fft, "weighted_mean", ok), 0.05)
})


test_that("fft check attaches a named accuracy report", {
      d <- fft_test_data(n = 1000, n_focal = 100)
      res <- analog_search(d$focal, d$pool, y = d$y[, "a"],
                           env = kernel("gaussian", theta = 0.4, max = 1.5),
                           stat = c("sum_weights", "weighted_mean"),
                           coord_type = "projected",
                           method = fft_method(check = 25))
      chk <- attr(res, "fft_check")
      expect_s3_class(chk, "data.frame")
      expect_identical(chk$column, c("sum_weights", "weighted_mean"))
      expect_true(all(chk$n <= 25))
      expect_true(all(chk$rms_err < chk$rms_exact))
})


test_that("fft works through wrappers and analog_cv", {
      d <- fft_test_data(n = 1500)
      env <- kernel("gaussian", theta = 0.4, max = 1.5)

      reg <- analog_regression(d$focal, d$pool, y = d$y[, "a"], covariates = d$covs,
                               env = env, coord_type = "projected", method = "fft")
      expect_true(all(c("coef_intercept", "coef_c1", "coef_c2") %in% names(reg)))

      cv <- analog_cv(analog_impact, pool = d$pool, y = d$y[, "a"], env = env,
                      cv_method = "loo", coord_type = "projected", method = "fft")
      expect_true(all(c("obs", "residual") %in% names(cv)))
      expect_true(sum(is.finite(cv$residual)) > 1000)
})


test_that("fft rejects unsupported configurations", {
      d <- fft_test_data(n = 500, n_focal = 50)
      env <- kernel("gaussian", theta = 0.4, max = 1.5)
      base <- list(d$focal, d$pool, y = d$y[, "a"], env = env,
                   stat = "weighted_mean", coord_type = "projected", method = "fft")
      try_fft <- function(...) {
            a <- utils::modifyList(base, list(...))
            do.call(analog_search, a)
      }

      expect_error(try_fft(select = "knn_env", k = 5, stat = "weighted_mean"),
                   'select = "all"')
      expect_error(try_fft(stat = "none"), "stat")
      expect_error(try_fft(geog = kernel(max = 20), coord_type = "lonlat"),
                   "projected")
      expect_error(try_fft(downsample = 0.5, env_res_adj = 1, geog_res_adj = 1),
                   "downsample")
      expect_error(try_fft(stat = "regression", covariates = d$covs, se = "design"),
                   "design")

      four <- cbind(d$pool, e3 = rnorm(nrow(d$pool)), e4 = rnorm(nrow(d$pool)))
      expect_error(analog_search(four, four, y = d$y[, "a"], env = env,
                                 stat = "weighted_mean", coord_type = "projected",
                                 method = "fft"),
                   "at most 3")

      idx <- build_analog_index(d$pool, coord_type = "projected")
      expect_error(analog_search(d$focal, idx, y = d$y[, "a"], env = env,
                                 stat = "weighted_mean", method = "fft"),
                   "analog_index")
})


# Geographic (and combined geographic x environmental) kernels -----------------

fft_geo_data <- function(seed = 2, n = 6000, n_focal = 500) {
      set.seed(seed)
      pool <- cbind(x = runif(n, 0, 1000), y = runif(n, 0, 1000),
                    e1 = rnorm(n), e2 = rnorm(n))
      pool[, "e1"] <- pool[, "e1"] + pool[, "x"] / 300   # climate gradient in space
      focal <- pool[sample(n, n_focal), ]
      focal[, "e1"] <- focal[, "e1"] + 0.5                # shifted (e.g. future) climate
      list(pool = pool, focal = focal,
           y = pool[, "e1"] + pool[, "e2"] + rnorm(n))
}


test_that("fft matches exact for geographic-only kernels", {
      d <- fft_geo_data()

      r <- run_both(d$focal, d$pool, y = d$y,
                    geog = kernel("gaussian", theta = 50, max = 150),
                    stat = c("sum_weights", "weighted_mean"),
                    coord_type = "projected", bin = 0.25)
      expect_identical(names(r$fft), names(r$exact))
      ok <- r$exact$sum_weights > stats::quantile(r$exact$sum_weights, 0.05)
      for (col in c("sum_weights", "weighted_mean")) {
            expect_lt(fft_rel_err(r$exact, r$fft, col, ok), 0.02)
      }

      # Hard-edged disk: anti-aliased edges keep count error small
      r <- run_both(d$focal, d$pool, y = d$y, geog = kernel(max = 150),
                    stat = c("count", "mean"), coord_type = "projected", bin = 0.25)
      ok <- r$exact$count >= 30
      rel <- abs(r$fft$count[ok] - r$exact$count[ok]) / r$exact$count[ok]
      expect_lt(stats::median(rel), 0.02)
      expect_lt(fft_rel_err(r$exact, r$fft, "mean", ok), 0.1)
})


test_that("fft matches exact for combined geographic and environmental kernels", {
      d <- fft_geo_data()
      r <- run_both(d$focal, d$pool, y = d$y,
                    geog = kernel("inverse", theta = 100, max = 300),
                    env = kernel("gaussian", theta = 0.4, max = 1.5),
                    stat = c("count", "sum_weights", "ess", "weighted_mean"),
                    se = "ess", coord_type = "projected", bin = 0.25)
      expect_identical(names(r$fft), names(r$exact))
      ok <- r$exact$count >= 30
      expect_true(sum(ok) > 100)
      for (col in c("sum_weights", "ess", "weighted_mean", "se_weighted_mean")) {
            expect_lt(fft_rel_err(r$exact, r$fft, col, ok), 0.03)
      }
})


test_that("fft regression with geographic x environmental kernels matches exact", {
      d <- fft_geo_data()
      covs <- cbind(c1 = rnorm(nrow(d$pool)))
      r <- run_both(d$focal, d$pool, y = d$y, covariates = covs,
                    geog = kernel("gaussian", theta = 100, max = 300),
                    env = kernel("gaussian", theta = 0.5, max = 1.5),
                    stat = c("count", "regression"), se = "ess",
                    coord_type = "projected", bin = 0.2)
      ok <- r$exact$count >= 50
      for (col in c("coef_intercept", "coef_c1", "se_intercept")) {
            expect_lt(fft_rel_err(r$exact, r$fft, col, ok), 0.05)
      }
})


test_that("fft buffered leave-one-out (geographic annulus) matches exact", {
      d <- fft_geo_data(n = 3000)
      args <- list(d$pool, d$pool, y = d$y,
                   geog = kernel("gaussian", theta = 80, max = 250, min = 40),
                   env = kernel("gaussian", theta = 0.5, max = 1.5),
                   stat = c("count", "weighted_mean"), exclude_self = TRUE,
                   coord_type = "projected")
      ex <- do.call(analog_search, c(args, method = "exact"))
      ff <- do.call(analog_search, c(args, method = "fft"))
      ok <- ex$count >= 30
      expect_lt(fft_rel_err(ex, ff, "weighted_mean", ok), 0.05)
})


test_that("fft trims the grid to where focals can reach pool points", {
      d <- fft_geo_data()
      # Focals confined to one corner: the grid should cover that corner plus
      # the kernel's reach, not the whole pool extent.
      corner <- d$pool[d$pool[, "x"] < 200 & d$pool[, "y"] < 200, ][1:100, ]
      r <- run_both(corner, d$pool, y = d$y,
                    geog = kernel("gaussian", theta = 30, max = 90),
                    stat = c("sum_weights", "weighted_mean"),
                    coord_type = "projected", bin = 0.25)
      grid <- attr(r$fft, "fft_grid")
      span <- (grid$n_nodes[1:2] - 1) * grid$spacing[1:2]
      expect_true(all(span < 300))
      ok <- r$exact$sum_weights > 0
      expect_lt(fft_rel_err(r$exact, r$fft, "weighted_mean", ok), 0.02)
})


test_that("fft works with projected raster inputs and normalization", {
      skip_if_not_installed("terra")
      set.seed(4)
      r <- terra::rast(nrows = 40, ncols = 40, ext = terra::ext(0, 4e5, 0, 4e5),
                       crs = "EPSG:32633")
      xy <- terra::xyFromCell(r, seq_len(terra::ncell(r)))
      e1 <- r; terra::values(e1) <- xy[, 1] / 2e5 + rnorm(terra::ncell(r), sd = 0.2)
      e2 <- r; terra::values(e2) <- xy[, 2] / 2e5 + rnorm(terra::ncell(r), sd = 0.2)
      pool <- c(e1, e2); names(pool) <- c("e1", "e2")
      pool[[1]][1:20] <- NA                        # NA cells are stripped
      fut <- pool + 0.3

      args <- list(fut, pool, geog = kernel(max = 1e5),
                   env = kernel("gaussian", theta = 0.3, max = 1),
                   stat = c("count", "sum_weights"))
      ex <- do.call(analog_density, c(args, method = "exact"))
      ff <- do.call(analog_density, c(args, method = list(fft_method(bin = 0.2))))
      expect_s4_class(ff, "SpatRaster")
      expect_identical(names(ff), names(ex))
      expect_true(isTRUE(attr(ff, "normalize")))

      e <- terra::values(ex[["sum_weights"]])[, 1]
      f <- terra::values(ff[["sum_weights"]])[, 1]
      expect_identical(is.na(e), is.na(f))
      ok <- is.finite(e) & e > stats::quantile(e, 0.05, na.rm = TRUE)
      expect_lt(stats::median(abs(f[ok] - e[ok]) / e[ok]), 0.03)
})


test_that("fft matches exact for value stats with missing responses", {
      d <- fft_test_data()
      y <- d$y
      y[sample(nrow(y), 300), "b"] <- NA
      r <- run_both(d$focal, d$pool, y = y,
                    env = kernel("gaussian", theta = 0.4, max = 1.5),
                    stat = c("count", "sum", "mean", "weighted_mean"),
                    se = "ess", coord_type = "projected")
      ok <- r$exact$count >= 30
      for (col in c("weighted_mean_a", "weighted_mean_b", "se_weighted_mean_b")) {
            expect_lt(fft_rel_err(r$exact, r$fft, col, ok), 0.02)
      }
      for (col in c("sum_b", "mean_b")) {
            expect_lt(fft_rel_err(r$exact, r$fft, col, ok), 0.1)
      }
})


test_that("fft check works with leave-one-out", {
      # Narrow kernel and noisy y, so each focal's own row carries real weight:
      # if the exact rerun on the checked subset excluded the wrong pool rows,
      # its error would be on the order of the leave-one-out effect itself.
      set.seed(1)
      n <- 1500
      pool <- cbind(x = runif(n), y = runif(n), e1 = rnorm(n), e2 = rnorm(n))
      y <- pool[, 3] + rnorm(n, sd = 2)
      args <- list(pool, pool, y = y, env = kernel("gaussian", theta = 0.15, max = 0.6),
                   stat = c("count", "weighted_mean"), coord_type = "projected")
      loo  <- do.call(analog_search, c(args, exclude_self = TRUE))
      self <- do.call(analog_search, c(args, exclude_self = FALSE))
      ok <- is.finite(loo$weighted_mean) & is.finite(self$weighted_mean)
      loo_gap <- sqrt(mean((loo$weighted_mean[ok] - self$weighted_mean[ok])^2))

      res <- do.call(analog_search, c(args, exclude_self = TRUE,
                                      method = list(fft_method(check = 300))))
      chk <- attr(res, "fft_check")
      expect_identical(chk$column, c("count", "weighted_mean"))
      expect_lt(chk$rms_err[2], loo_gap / 10)
})



test_that("fft runs through tiled_analog_search and matches untiled fft", {
      skip_if_not_installed("terra")
      set.seed(4)
      r <- terra::rast(nrows = 80, ncols = 80, ext = terra::ext(0, 8e4, 0, 8e4),
                       crs = "EPSG:32633")
      xy <- terra::xyFromCell(r, seq_len(terra::ncell(r)))
      e1 <- r; terra::values(e1) <- 3 * xy[, 1] / 8e4 + rnorm(terra::ncell(r), sd = 0.2)
      e2 <- r; terra::values(e2) <- 2 * xy[, 2] / 8e4 + rnorm(terra::ncell(r), sd = 0.2)
      pool <- c(e1, e2); names(pool) <- c("e1", "e2")
      fut <- pool + 0.3
      y <- e1 + e2; names(y) <- "y"
      args <- list(geog = kernel(max = 1e4),
                   env = kernel("gaussian", theta = 0.15, max = 0.45),
                   stat = c("weighted_mean", "sum_weights"), method = "fft")

      untiled <- do.call(analog_impact, c(list(x = fut, pool = pool, y = y), args))
      tiled <- do.call(tiled_analog_search,
                       c(list(x = fut, pool = pool, y = y, n_tiles = 4,
                              fun = analog_impact, progress = FALSE), args))

      u <- terra::values(untiled[["weighted_mean"]])[, 1]
      t <- terra::values(tiled[["weighted_mean"]])[, 1]
      expect_identical(is.finite(u), is.finite(t))
      ok <- is.finite(u)
      # Tiles place grid nodes differently, so results differ by binning error only
      expect_lt(stats::median(abs(t[ok] - u[ok])) / stats::sd(u[ok]), 0.01)
})


test_that("fft runs through tiled_analog_search with point tables", {
      set.seed(5)
      n <- 4000
      pool <- data.frame(x = runif(n, 0, 8e4), y = runif(n, 0, 8e4))
      pool$e1 <- 3 * pool$x / 8e4 + rnorm(n, sd = 0.2)
      pool$e2 <- 2 * pool$y / 8e4 + rnorm(n, sd = 0.2)
      focal <- pool[sample(n, 1500), ]
      focal$e1 <- focal$e1 + 0.3
      focal[3, c("x", "y")] <- NA
      y <- pool$e1 + pool$e2
      args <- list(geog = kernel("gaussian", theta = 6e3, max = 1.8e4),
                   env = kernel("gaussian", theta = 0.3, max = 0.9),
                   stat = c("weighted_mean", "sum_weights"), method = "fft",
                   coord_type = "projected")

      untiled <- do.call(analog_impact, c(list(x = focal, pool = pool, y = y), args))
      tiled <- do.call(tiled_analog_search,
                       c(list(x = focal, pool = pool, y = y, n_tiles = 4,
                              fun = analog_impact, progress = FALSE), args))

      expect_identical(tiled$index, untiled$index)
      u <- untiled$weighted_mean
      t <- tiled$weighted_mean
      expect_identical(is.finite(u), is.finite(t))
      ok <- is.finite(u)
      expect_lt(stats::median(abs(t[ok] - u[ok])) / stats::sd(u[ok]), 0.01)
})


# Focals without analogs ------------------------------------------------------

# A patchy pool: geographic holes wider than the kernels' `max`, and no warm
# climates, with focals spread over and beyond it, so many focals have no pool
# point within the limits and many more sit near the edge of the pool's reach.
fft_patchy_data <- function(seed = 7, n = 4000, n_focal = 3000) {
      set.seed(seed)
      pool <- cbind(x = runif(n, 0, 300), y = runif(n, 0, 300),
                    e1 = rnorm(n), e2 = rnorm(n))
      hole <- (pool[, "x"] - 150)^2 + (pool[, "y"] - 150)^2 < 100^2 |
            (pool[, "x"] > 240 & pool[, "y"] < 90) |
            pool[, "e1"] > 0.5
      pool <- pool[!hole, ]
      focal <- cbind(x = runif(n_focal, -30, 330), y = runif(n_focal, -30, 330),
                     e1 = rnorm(n_focal, sd = 1.3), e2 = rnorm(n_focal))
      y <- cbind(a = pool[, "e1"] + rnorm(nrow(pool)),
                 b = rbinom(nrow(pool), 1, 0.3))
      list(pool = pool, focal = focal, y = y,
           covs = cbind(c1 = pool[, "e2"] + rnorm(nrow(pool))))
}

# NA where exact is NA, exact zeros where exact is zero
expect_same_support <- function(exact, fft) {
      for (col in setdiff(names(exact), c("index", "x", "y"))) {
            expect_identical(is.na(fft[[col]]), is.na(exact[[col]]), label = col)
      }
      for (col in intersect(c("count", "sum_weights", "ess"), names(exact))) {
            z <- which(exact[[col]] == 0)
            expect_true(all(fft[[col]][z] == 0), label = col)
      }
}

test_that("fft returns NA and zeros where no pool point is within the limits", {
      d <- fft_patchy_data()
      r <- run_both(d$focal, d$pool, y = d$y, coord_type = "projected",
                    geog = kernel("gaussian", theta = 30, max = 90),
                    env = kernel("gaussian", theta = 0.5, max = 1.5),
                    stat = c("weighted_mean", "mean", "count", "sum_weights", "ess"),
                    bin = 0.25)
      expect_gt(sum(is.na(r$exact$weighted_mean_a)), 100)   # the case is exercised
      expect_same_support(r$exact, r$fft)
      expect_gt(metadata(r$fft)$fft_exact_rows, 0)
})

test_that("fft support matches exact with uniform kernels, box limits, and an annulus", {
      d <- fft_patchy_data()
      r <- run_both(d$focal, d$pool, y = d$y, coord_type = "projected",
                    geog = kernel(max = 90, min = 25),
                    env = kernel("inverse", theta = 0.5, max = c(1, 2)),
                    stat = c("weighted_mean", "count"), bin = 0.25)
      expect_gt(sum(is.na(r$exact$weighted_mean_a)), 100)
      expect_same_support(r$exact, r$fft)
})

test_that("fft support is per response when responses are missing", {
      d <- fft_patchy_data()
      # Response b observed only in the west: focals in the east have analogs
      # for a but none for b
      d$y[d$pool[, "x"] > 180, "b"] <- NA
      r <- run_both(d$focal, d$pool, y = d$y, coord_type = "projected",
                    geog = kernel("gaussian", theta = 30, max = 90),
                    env = kernel("gaussian", theta = 0.5, max = 1.5),
                    covariates = d$covs,
                    stat = c("weighted_mean", "sum_weights", "regression"), bin = 0.25)
      only_a <- !is.na(r$exact$weighted_mean_a) & is.na(r$exact$weighted_mean_b)
      expect_gt(sum(only_a), 100)
      expect_same_support(r$exact, r$fft)
})

test_that("fft support matches exact under leave-one-out", {
      # A sparse pool, so that many points have no other point within the
      # limits once their own is excluded
      set.seed(8)
      pool <- cbind(x = runif(150, 0, 300), y = runif(150, 0, 300), e1 = rnorm(150))
      args <- list(pool, pool, y = pool[, "e1"] + rnorm(150), coord_type = "projected",
                   geog = kernel("gaussian", theta = 7, max = 20),
                   env = kernel("gaussian", theta = 0.5, max = 1.5),
                   stat = c("weighted_mean", "count"), exclude_self = TRUE)
      exact <- do.call(analog_search, args)
      fft <- do.call(analog_search, c(args, method = "fft"))
      expect_gt(sum(is.na(exact$weighted_mean)), 10)
      expect_same_support(exact, fft)
})
