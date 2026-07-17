# Regression tests for NA-stripping alignment in analog_cv().
#
# When `pool` is a SpatRaster (or matrix) containing NA cells, .format_data()
# drops those rows and analog_cv() must subset the pool-keyed inputs (y,
# covariates, weight, cell_area_weight) to the same stripped ordering. Prior to
# the fix, y/covariates were validated and passed at the ORIGINAL (unstripped)
# length, causing a spurious "must have exactly N rows/cells" error whenever the
# raster had any NA cells -- even though the identical y worked in a direct
# analog_search()/analog_regression() call.

make_pool_raster_with_na <- function(seed = 1, nc = 12, nr = 12, n_na = 7) {
      skip_if_not_installed("terra")
      set.seed(seed)
      n <- nc * nr
      cwd <- rnorm(n)
      aet <- rnorm(n)
      # Punch NA holes into the env layers (same cells across layers, as in a
      # real masked raster).
      na_idx <- sample(n, n_na)
      cwd[na_idx] <- NA
      aet[na_idx] <- NA
      r <- terra::rast(nrows = nr, ncols = nc,
                       xmin = 0, xmax = nc, ymin = 0, ymax = nr, nlyrs = 2)
      terra::values(r) <- cbind(cwd, aet)
      names(r) <- c("CWD", "AET")
      list(rast = r, na_idx = na_idx, n = n)
}


test_that("analog_cv accepts a raster y against an NA-containing raster pool", {
      skip_if_not_installed("terra")
      f <- make_pool_raster_with_na()
      hist <- f$rast

      # y as a single-layer SpatRaster keyed to the FULL raster (with NAs),
      # exactly as in the vignette (eco_var <- hist$CWD).
      eco_var <- hist$CWD

      # This used to error: "`y` must have exactly <n_complete> rows/cells".
      expect_no_error(
            cv <- analog_cv(
                  fun       = analog_impact,
                  pool      = hist,
                  y         = eco_var,
                  env       = kernel(max = 1, weight = "gaussian", theta = 0.5),
                  geog      = kernel(max = 5, weight = "gaussian", theta = 2),
                  se        = "ess",
                  cv_method = "loo"
            )
      )
})


test_that("CV result aligns with an equivalent direct analog_search call", {
      skip_if_not_installed("terra")
      # Larger, denser raster so a self-excluded weighted mean is finite on
      # plenty of cells in BOTH runs (toy 12x12 grids with tight kernels can
      # leave most focals analog-less, giving no overlap to compare).
      f <- make_pool_raster_with_na(seed = 7, nc = 20, nr = 20, n_na = 15)
      hist <- f$rast
      eco_var <- hist$CWD

      # Geographic neighborhood covering the whole extent, distance-weighted.
      # No env threshold, so every focal has a large analog pool and the
      # self-excluded weighted mean is finite on every non-NA cell -- giving
      # full overlap for the CV-vs-direct comparison.
      k_env  <- NULL
      k_geog <- kernel(max = 5000, weight = "gaussian", theta = 8)

      cv <- analog_cv(
            fun       = analog_search,
            pool      = hist,
            y         = eco_var,
            stat      = "weighted_mean",
            env       = k_env,
            geog      = k_geog,
            cv_method = "loo"
      )

      # Reference: the same self-excluded weighted-mean search, run directly.
      # analog_search() does its own joint NA stripping of pool + y, so passing
      # the full raster here is the correct apples-to-apples comparison.
      ref <- analog_search(
            x            = hist,
            pool         = hist,
            y            = eco_var,
            stat         = "weighted_mean",
            env          = k_env,
            geog         = k_geog,
            exclude_self = TRUE
      )

      cv_wm  <- if (inherits(cv, "SpatRaster"))
            terra::values(cv[["weighted_mean"]])[, 1] else cv$weighted_mean
      ref_wm <- if (inherits(ref, "SpatRaster"))
            terra::values(ref[["weighted_mean"]])[, 1] else ref$weighted_mean

      expect_false(is.null(cv_wm))
      expect_false(is.null(ref_wm))
      expect_equal(length(cv_wm), length(ref_wm))
      ok <- is.finite(cv_wm) & is.finite(ref_wm)
      expect_true(sum(ok) > 0)
      expect_equal(cv_wm[ok], ref_wm[ok], tolerance = 1e-8)
})


test_that("analog_cv handles a full-length vector y with NA raster pool", {
      skip_if_not_installed("terra")
      f <- make_pool_raster_with_na()
      hist <- f$rast

      # y supplied as a plain length-ncell vector (with NAs at the masked
      # cells), which must also be validated at original length and stripped.
      y_vec <- as.vector(terra::values(hist$CWD))
      expect_length(y_vec, f$n)

      expect_no_error(
            cv <- analog_cv(
                  fun       = analog_impact,
                  pool      = hist,
                  y         = y_vec,
                  env       = kernel(max = 1, weight = "gaussian", theta = 0.5),
                  geog      = kernel(max = 5),
                  cv_method = "loo"
            )
      )
})


test_that("covariates raster is aligned under regression CV with NA pool", {
      skip_if_not_installed("terra")
      f <- make_pool_raster_with_na()
      hist <- f$rast
      eco_var <- hist$CWD

      covs <- data.frame(
            aet  = as.vector(terra::values(hist$AET)),
            aet2 = as.vector(terra::values(hist$AET))^2
      )

      expect_no_error(
            cv <- analog_cv(
                  fun        = analog_regression,
                  pool       = hist,
                  y          = eco_var,
                  covariates = covs,
                  select     = "all",
                  geog       = kernel(max = 5, weight = "inverse", theta = 2),
                  lambda     = 0,
                  cv_method  = "loo"
            )
      )
})


test_that("wrong-length y still errors against original pool size", {
      skip_if_not_installed("terra")
      f <- make_pool_raster_with_na()
      hist <- f$rast

      # A y whose length matches neither the original nor stripped pool size
      # must still be rejected, reported against the original size.
      bad_y <- rnorm(f$n - 3)
      expect_error(
            analog_cv(
                  fun       = analog_impact,
                  pool      = hist,
                  y         = bad_y,
                  env       = kernel(max = 1),
                  geog      = kernel(max = 5),
                  cv_method = "loo"
            ),
            "exactly"
      )
})
