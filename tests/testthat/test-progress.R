# Progress reporting through progressr.

# Evaluate `expr` under a capturing progressr handler; return the result and
# the number of update steps the handler received.
with_captured_progress <- function(expr) {
      steps <- 0L
      h <- progressr::make_progression_handler(
            "capture",
            reporter = list(update = function(config, state, progression, ...) {
                  steps <<- steps + 1L
            }),
            enable = TRUE, interval = 0, times = Inf
      )
      res <- progressr::with_progress(expr, handlers = h, enable = TRUE)
      list(result = res, steps = steps)
}

# Small chunks, so test-sized inputs exercise the chunked exact path
local_small_chunks <- function(min_chunk = 100L, max_chunks = 7L, env = parent.frame()) {
      withr::local_options(analogs.progress_min_chunk = min_chunk,
                           analogs.progress_max_chunks = max_chunks,
                           .local_envir = env)
}

progress_fixture <- function(n = 1000) {
      set.seed(42)
      pool <- data.frame(x = runif(n, 0, 100), y = runif(n, 0, 100),
                         e1 = rnorm(n), e2 = rnorm(n))
      list(pool = pool, focal = pool[sample(n, 650), ], y = rnorm(n),
           args = list(env = kernel("gaussian", theta = 0.5, max = 1.5),
                       geog = kernel("gaussian", theta = 10, max = 30),
                       coord_type = "projected",
                       geog_res_adj = 1, env_res_adj = 1))
}

test_that("exact aggregation reports one step per chunk with identical results", {
      local_small_chunks()
      d <- progress_fixture()
      run <- function(progress) do.call(analog_impact, c(list(
            x = d$focal, pool = d$pool, y = d$y,
            stat = c("weighted_mean", "sum_weights", "count"),
            progress = progress), d$args))
      ref <- run(FALSE)
      got <- with_captured_progress(run(TRUE))
      expect_identical(got$result, ref)
      expect_equal(got$steps, 7L)          # 650 rows: min(7, ceiling(650 / 100))
})

test_that("chunk count respects the minimum chunk size", {
      local_small_chunks(min_chunk = 200L, max_chunks = 100L)
      d <- progress_fixture()
      got <- with_captured_progress(do.call(analog_impact, c(list(
            x = d$focal, pool = d$pool, y = d$y, progress = TRUE), d$args)))
      expect_equal(got$steps, 4L)          # ceiling(650 / 200)
})

test_that("inputs below the minimum chunk size report a single step", {
      d <- progress_fixture()             # default min chunk (2000) > 650 rows
      ref <- do.call(analog_impact, c(list(x = d$focal, pool = d$pool, y = d$y), d$args))
      got <- with_captured_progress(do.call(analog_impact, c(list(
            x = d$focal, pool = d$pool, y = d$y, progress = TRUE), d$args)))
      expect_identical(got$result, ref)
      expect_equal(got$steps, 1L)
})

test_that("chunked pair mode matches the unchunked result", {
      local_small_chunks()
      d <- progress_fixture()
      run <- function(progress) analog_search(
            x = d$focal, pool = d$pool, select = "knn_env", k = 3, stat = "none",
            env = kernel(max = 1.5), geog = kernel(max = 30),
            coord_type = "projected", geog_res_adj = 1, env_res_adj = 1,
            progress = progress)
      ref <- run(FALSE)
      got <- with_captured_progress(run(TRUE))
      expect_identical(got$result, ref)
      expect_equal(got$steps, 7L)
})

test_that("chunked exclude_self matches the unchunked result", {
      local_small_chunks()
      d <- progress_fixture()
      run <- function(progress) do.call(analog_impact, c(list(
            x = d$pool, pool = d$pool, y = d$y, exclude_self = TRUE,
            stat = c("weighted_mean", "count"), progress = progress), d$args))
      ref <- run(FALSE)
      got <- with_captured_progress(run(TRUE))
      expect_identical(got$result, ref)
      expect_equal(got$steps, 7L)
})

test_that("chunked x_cov queries match the unchunked result", {
      local_small_chunks()
      d <- progress_fixture()
      nf <- nrow(d$focal)   # variances, then covariance
      x_cov <- cbind(runif(nf, 0.01, 0.2), runif(nf, 0.01, 0.2), 0.001)
      run <- function(progress) do.call(analog_impact, c(list(
            x = d$focal, pool = d$pool, y = d$y, x_cov = x_cov,
            progress = progress), d$args))
      expect_identical(with_captured_progress(run(TRUE))$result, run(FALSE))
})

test_that("progress = FALSE signals nothing, even inside a handler context", {
      local_small_chunks()
      d <- progress_fixture()
      got <- with_captured_progress(do.call(analog_impact, c(list(
            x = d$focal, pool = d$pool, y = d$y), d$args)))
      expect_equal(got$steps, 0L)
})

test_that("no output without a handler, unless progress is requested", {
      local_small_chunks()
      d <- progress_fixture()
      run <- function(progress) do.call(analog_impact, c(list(
            x = d$focal, pool = d$pool, y = d$y, progress = progress), d$args))
      expect_silent(run(FALSE))

      # progress = TRUE with no handler listening draws a default bar (here
      # forced on; progressr only renders in interactive sessions by default)
      withr::local_options(progressr.enable = TRUE)
      out <- capture.output(res <- run(TRUE), type = "message")
      expect_true(any(grepl("100%", out)))
})

test_that("FFT method reports a single step", {
      d <- progress_fixture()
      args <- d$args
      args$geog_res_adj <- args$env_res_adj <- NULL
      run <- function(progress) do.call(analog_impact, c(list(
            x = d$focal, pool = d$pool, y = d$y, method = "fft",
            stat = c("weighted_mean", "sum_weights"), progress = progress), args))
      got <- with_captured_progress(run(TRUE))
      expect_equal(got$result, run(FALSE))
      expect_equal(got$steps, 1L)
})

test_that("tiled_analog_search reports one step per tile, including empty ones", {
      d <- progress_fixture()
      args <- d$args
      # Focals confined to the left half leave some tiles without focals
      focal <- d$focal[d$focal$x < 50, ]
      # (halo-size warning is expected at this toy scale)
      run <- function(progress) suppressWarnings(do.call(tiled_analog_search, c(list(
            x = focal, pool = d$pool, y = d$y, n_tiles = 4, fun = analog_impact,
            progress = progress), args)))
      # n_tiles is approximate; expect one step per tile actually laid out
      grid <- optimize_tile_grid(4, diff(range(focal$x)) / diff(range(focal$y)))
      got <- with_captured_progress(run(TRUE))
      expect_identical(got$result, run(FALSE))
      expect_equal(got$steps, grid$n_x * grid$n_y)
      expect_equal(with_captured_progress(run(FALSE))$steps, 0L)
})

test_that("analog_cv reports progress per fold (kfold) and per chunk (loo)", {
      local_small_chunks()
      d <- progress_fixture()
      run <- function(progress, ...) do.call(analog_cv, c(list(
            fun = analog_impact, pool = d$pool, y = d$y, progress = progress, ...),
            d$args))

      set.seed(1)
      fold_id <- sample(rep(1:5, length.out = nrow(d$pool)))
      ref <- run(FALSE, cv_method = "kfold", fold_id = fold_id)
      got <- with_captured_progress(run(TRUE, cv_method = "kfold", fold_id = fold_id))
      expect_identical(got$result, ref)
      expect_equal(got$steps, 5L)

      ref <- run(FALSE, cv_method = "loo")
      got <- with_captured_progress(run(TRUE, cv_method = "loo"))
      expect_identical(got$result, ref)
      expect_equal(got$steps, 7L)          # 1000 rows: min(7, ceiling(1000 / 100))

      expect_error(run(NA, cv_method = "loo"), "progress")
})
