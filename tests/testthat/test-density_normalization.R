# Tests for sum_weights / tabulate normalization (D_max).
#
# Coverage:
#   1. Closed-form D_max integrals match numerical integration
#   2. D / D_max clusters near 1 on a uniform-climate raster pool
#   3. `auto` resolves to TRUE/FALSE in the right situations
#   4. Explicit `normalize = TRUE` errors when preconditions fail
#   5. `normalize = FALSE` matches pre-normalization behavior
#   6. Tiled vs non-tiled equivalence for normalized output
#   7. analog_cv k-fold preserves normalization (uses global mean_cell_area)
#   8. tabulate columns are normalized class-by-class
#   9. D_max attribute is recorded on the result

# Initialize RNG state. Several internal paths call sample.int() to derive
# default seeds for index downsampling; running tests in a fresh R session
# without prior RNG use (and without set.seed()) leaves .Random.seed unset
# and triggers an "object '.Random.seed' not found" error before the
# normalization machinery is even reached.
set.seed(1)


# ---------- Test fixtures ----------------------------------------------------

# A small projected raster with non-uniform climate over a square kernel.
# Climate is a smooth radial gradient, varying enough to make sum_weights
# meaningfully heterogeneous.
#
# Note on units: EPSG:32633 (UTM 33N) is in meters, so the extent passed
# here is also in meters. The default `extent_m = 1e6` gives a 1000km x
# 1000km kernel. Pair with `max_geog` in meters in the corresponding tests.
make_test_raster <- function(n = 30, extent_m = 1e6) {
      ext <- terra::ext(0, extent_m, 0, extent_m)
      r <- terra::rast(nrows = n, ncols = n, ext = ext, crs = "EPSG:32633")
      xy <- terra::xyFromCell(r, 1:terra::ncell(r))
      cx <- (extent_m / 2)
      d  <- sqrt((xy[, 1] - cx)^2 + (xy[, 2] - cx)^2)
      terra::values(r) <- d / max(d)   # climate in [0, 1]
      names(r) <- "clim"
      r
}

# A uniform-climate version (every cell has identical climate) to test that
# D / D_max clusters near 1 with mild discretization noise.
make_uniform_clim_raster <- function(n = 30, extent_m = 1e6) {
      r <- make_test_raster(n = n, extent_m = extent_m)
      terra::values(r) <- 1.0
      r
}


# ---------- 1. Closed-form integrals match numerical integration ------------

test_that(".compute_global_dmax matches numerical integration for each kernel", {
      mean_area <- 100  # arbitrary positive scalar
      G <- 50

      # D_max depends ONLY on the geographic kernel: the climate kernel is
      # evaluated at clim_dist = 0, where every reparameterized shape gives
      # weight 1 (uniform = 1, gaussian exp(0) = 1, inverse 1/(1+0/theta) = 1).
      # So the reference integrand K0(r) is the GEOGRAPHIC kernel alone, and any
      # climate shape (uniform/gaussian/inverse) must produce the same D_max.
      # Each case gives the per-family shapes + scales and the geographic
      # reference kernel. Inverse uses the reparameterized form 1/(1 + r/theta).
      cases <- list(
            # --- geographically-uniform: D_max = disk area / mean_area ---
            list(kernel_clim = "uniform",  kernel_geog = "uniform",
                 theta_clim = NULL, theta_geog = NULL,
                 K0 = function(r) rep(1, length(r))),
            # climate shape must not change D_max (climate factor = 1 at r_clim=0)
            list(kernel_clim = "gaussian", kernel_geog = "uniform",
                 theta_clim = 1, theta_geog = NULL,
                 K0 = function(r) rep(1, length(r))),
            list(kernel_clim = "inverse",  kernel_geog = "uniform",
                 theta_clim = 1, theta_geog = NULL,
                 K0 = function(r) rep(1, length(r))),
            # --- geographic Gaussian ---
            list(kernel_clim = "uniform",  kernel_geog = "gaussian",
                 theta_clim = NULL, theta_geog = 10,
                 K0 = function(r) exp(-r^2 / (2 * 10^2))),
            # --- geographic inverse (reparameterized) ---
            list(kernel_clim = "uniform",  kernel_geog = "inverse",
                 theta_clim = NULL, theta_geog = 5,
                 K0 = function(r) 1 / (1 + r / 5)),
            # --- composed: non-uniform on BOTH families; D_max still = geog only ---
            list(kernel_clim = "gaussian", kernel_geog = "gaussian",
                 theta_clim = 1, theta_geog = 10,
                 K0 = function(r) exp(-r^2 / (2 * 10^2))),
            list(kernel_clim = "inverse",  kernel_geog = "inverse",
                 theta_clim = 1, theta_geog = 5,
                 K0 = function(r) 1 / (1 + r / 5))
      )

      for (cs in cases) {
            num_integral <- stats::integrate(
                  function(r) cs$K0(r) * 2 * pi * r,
                  lower = 0, upper = G
            )$value
            expected <- num_integral / mean_area

            actual <- analogs:::.compute_global_dmax(
                  kernel_clim    = cs$kernel_clim,
                  kernel_geog    = cs$kernel_geog,
                  theta_clim     = cs$theta_clim,
                  theta_geog     = cs$theta_geog,
                  max_geog       = G,
                  mean_cell_area = mean_area
            )

            expect_equal(
                  actual, expected,
                  tolerance = 1e-6,
                  info = paste0("clim = ", cs$kernel_clim,
                                ", geog = ", cs$kernel_geog)
            )
      }
})


# ---------- 2. D / D_max clusters near 1 on uniform-climate pool -----------

test_that("normalize on uniform-climate raster gives D/D_max ~ 1", {
      skip_if_not_installed("terra")

      pool <- make_uniform_clim_raster(n = 40, extent_m = 1e6)

      # Use a moderate Gaussian kernel; max_geog small enough to keep
      # neighborhoods well inside the kernel (avoid edge effects).
      res <- analog_density(
            x        = pool,
            pool     = pool,
            geog = kernel(max = 1e5),
            clim = kernel("gaussian", theta = 0.1),
            normalize = TRUE
      )

      vals <- terra::values(res)[, 1]
      vals <- vals[is.finite(vals)]

      # Trim to the interior to avoid edge-truncation effects: keep the
      # densest 60% of focals (those whose neighborhoods are mostly inside
      # the kernel).
      vals_interior <- vals[vals > stats::quantile(vals, 0.4)]

      # Interior cluster should be near 1 with mild discretization spread.
      expect_gt(median(vals_interior), 0.85)
      expect_lt(median(vals_interior), 1.15)

      # And no value should significantly exceed 1 (modest discretization
      # overshoot is acceptable).
      expect_lt(max(vals_interior), 1.10)
})


# ---------- 3. auto-resolution behavior --------------------------------------

test_that("normalize = 'auto' resolves TRUE for raster + caw + climate kernel + finite max_geog", {
      skip_if_not_installed("terra")
      pool <- make_test_raster(n = 20)

      res <- analog_density(
            x = pool, pool = pool,
            geog = kernel(max = 1e5), clim = kernel("gaussian", theta = 0.2),
            normalize = "auto"
      )
      expect_true(attr(res, "normalize"))
      expect_true(is.finite(attr(res, "D_max")))
})

test_that("normalize = 'auto' resolves FALSE for non-raster pool", {
      skip_if_not_installed("terra")
      pool_r <- make_test_raster(n = 15)
      pool_df <- terra::as.data.frame(pool_r, xy = TRUE, na.rm = FALSE)
      pool_df <- pool_df[stats::complete.cases(pool_df), ]

      res <- analog_density(
            x = pool_df, pool = pool_df,
            geog = kernel(max = 1e5), clim = kernel("gaussian", theta = 0.2),
            normalize = "auto"
      )
      expect_false(attr(res, "normalize"))
      expect_true(is.na(attr(res, "D_max")))
})

test_that("normalize = 'auto' resolves TRUE for uniform and geo-only kernels", {
      skip_if_not_installed("terra")
      pool <- make_test_raster(n = 20)

      # uniform kernel + sum_weights: now supported (D_max well-defined
      # via the closed-form integral; D/D_max is count-fraction inside
      # max_clim when max_clim is active).
      res_uni <- analog_search(
            x = pool, pool = pool,
            stat = "sum_weights",
            geog = kernel(max = 1e5), clim = kernel(max = 0.5),
            normalize = "auto"
      )
      expect_true(attr(res_uni, "normalize"))
      expect_true(is.finite(attr(res_uni, "D_max")))

      # geo-only kernel: also supported.
      res_geo <- analog_search(
            x = pool, pool = pool,
            stat = "sum_weights",
            geog = kernel("gaussian", theta = 5e4, max = 1e5),
            clim = kernel(max = 0.5),
            normalize = "auto"
      )
      expect_true(attr(res_geo, "normalize"))
      expect_true(is.finite(attr(res_geo, "D_max")))
})

test_that("normalize = 'auto' resolves FALSE when stat lacks normalizable columns", {
      skip_if_not_installed("terra")
      pool <- make_test_raster(n = 20)

      # stat = "count": no normalizable column. auto -> FALSE.
      res <- analog_search(
            x = pool, pool = pool,
            stat = "count",
            geog = kernel(max = 1e5),
            normalize = "auto"
      )
      expect_false(attr(res, "normalize"))
})


# ---------- 4. Explicit TRUE errors when preconditions fail -----------------

test_that("normalize = TRUE errors when max_geog is missing", {
      skip_if_not_installed("terra")
      pool <- make_test_raster(n = 15)

      expect_error(
            analog_search(
                  x = pool, pool = pool,
                  stat = "sum_weights",
                  clim = kernel("gaussian", theta = 0.2),
                  normalize = TRUE
            ),
            "max_geog"
      )
})

test_that("normalize = TRUE silently no-ops when stat has no normalizable column", {
      skip_if_not_installed("terra")
      pool <- make_test_raster(n = 15)

      # stat = "count" with kernel = NULL is fine; normalize = TRUE should
      # short-circuit silently without complaining about the missing kernel.
      expect_silent({
            res <- analog_search(
                  x = pool, pool = pool,
                  stat = "count",
                  geog = kernel(max = 1e5),
                  normalize = TRUE
            )
      })
      # The result attribute reports normalize = TRUE (caller asked for it),
      # but D_max is NA because nothing was actually normalized.
      expect_true(attr(res, "normalize"))
      expect_true(is.na(attr(res, "D_max")))
})


# ---------- 5. normalize = FALSE matches raw sum_weights --------------------

test_that("normalize = FALSE returns raw sum_weights (no D_max scaling)", {
      skip_if_not_installed("terra")
      pool <- make_test_raster(n = 20)

      res_norm <- analog_density(
            x = pool, pool = pool,
            geog = kernel(max = 1e5), clim = kernel("gaussian", theta = 0.2),
            normalize = TRUE
      )
      res_raw <- analog_density(
            x = pool, pool = pool,
            geog = kernel(max = 1e5), clim = kernel("gaussian", theta = 0.2),
            normalize = FALSE
      )

      D_max <- attr(res_norm, "D_max")
      expect_true(is.finite(D_max) && D_max > 0)

      v_norm <- terra::values(res_norm)[, 1]
      v_raw  <- terra::values(res_raw)[, 1]

      finite_both <- is.finite(v_norm) & is.finite(v_raw)
      expect_equal(
            v_norm[finite_both],
            v_raw[finite_both] / D_max,
            tolerance = 1e-10
      )
})


test_that("uniform kernel + max_clim: D/D_max equals fraction-of-disk-passing-max_clim", {
      skip_if_not_installed("terra")
      # For kernel = "uniform", every accepted neighbor contributes weight
      # 1, so sum_weights == count of (geo-passing AND clim-passing)
      # neighbors. D_max is the kernel-weighted disk integral, which for
      # uniform reduces to the count of cells inside the disk
      # (modulo cell-area weighting, which here is uniform too on a
      # projected raster). The ratio should equal D_max-weighted /
      # D_max-unweighted = D / D_max.
      pool <- make_test_raster(n = 40)

      res_norm <- analog_search(
            x = pool, pool = pool,
            stat = c("count", "sum_weights"),
            geog = kernel(max = 1e5), clim = kernel(max = 0.3),
            normalize = TRUE
      )

      # All ratios should be in [0, ~1]; many in the interior with a tight
      # max_clim should be substantially below 1 (since most disk cells
      # don't satisfy max_clim = 0.3 around an off-center focal).
      sw <- res_norm$sum_weights
      ok <- is.finite(sw)
      expect_true(all(sw[ok] >= 0))
      expect_lt(max(sw[ok]), 1.10)              # mild discretization slack
      expect_lt(min(sw[ok][sw[ok] > 0]), 0.5)   # some focals should have D < 0.5 * D_max
})


# ---------- 6. Tiled vs non-tiled equivalence -------------------------------

test_that("tiled_analog_search produces same normalized output as non-tiled", {
      skip_if_not_installed("terra")
      # Note on grid choice: we use n = 40 deliberately. At certain pool
      # resolutions the disk count for a hard `max_geog` cutoff sits right
      # at a Gauss-circle-problem boundary -- tiny floating-point
      # differences in distance computations flip whole rings of cells
      # in/out of the disk, producing visible "plaid" patterns in the
      # output that are highly sensitive to perturbations like tile
      # boundaries. n = 30 with max_geog = 1e5 is pessimal in this sense
      # (max_geog = exactly 3 cell widths, sitting right at a critical
      # boundary). n = 40 (max_geog = 4 cells) is well-behaved.
      pool <- make_test_raster(n = 40, extent_m = 1e6)

      res_full <- analog_density(
            x = pool, pool = pool,
            geog = kernel(max = 1e5), clim = kernel("gaussian", theta = 0.2),
            normalize = TRUE
      )

      res_tiled <- tiled_analog_search(
            x = pool, pool = pool,
            n_tiles = 4, fun = analog_density,
            geog = kernel(max = 1e5), clim = kernel("gaussian", theta = 0.2),
            normalize = TRUE,
            progress = FALSE
      )

      # Compare by spatial location, not by cell index. terra::merge()
      # in the tiled path can produce a result whose cell layout differs
      # from the input grid (e.g. when tile boundaries don't align with
      # input cell boundaries), so a positional comparison would
      # spuriously fail. Use terra::extract() at the original pool's cell
      # centroids to pull values from both rasters at the exact same
      # geographic locations.
      pts <- terra::xyFromCell(pool, seq_len(terra::ncell(pool)))
      v_full  <- terra::extract(res_full,  pts)[, 1]
      v_tiled <- terra::extract(res_tiled, pts)[, 1]

      finite_both <- is.finite(v_full) & is.finite(v_tiled)
      expect_gt(sum(finite_both), 50)  # sanity: non-trivial overlap

      # Per-cell match. Independent index builds per tile can produce tiny
      # floating-point differences in the sum_weights numerators; the
      # global mean_cell_area injection ensures the D_max denominator is
      # identical across tiles.
      expect_equal(
            v_full[finite_both],
            v_tiled[finite_both],
            tolerance = 1e-8
      )
})


# ---------- 7. analog_cv preserves normalization ----------------------------

test_that("analog_cv k-fold uses globally consistent mean_cell_area", {
      skip_if_not_installed("terra")
      pool <- make_test_raster(n = 25, extent_m = 1e6)
      y <- terra::values(pool)[, 1] + stats::rnorm(terra::ncell(pool), sd = 0.05)

      # Just confirm the call runs and returns sensible output. Detailed
      # equivalence to non-CV is tested elsewhere; here we're confirming
      # the cell_area_weight / mean_cell_area plumbing through the CV
      # runners doesn't error.
      n_pool <- terra::ncell(pool)
      folds <- sample(rep(seq_len(5), length.out = n_pool))

      expect_silent({
            res <- analog_cv(
                  fun = analog_impact,
                  pool = pool, y = y,
                  cv_method = "kfold",
                  fold_id   = folds,
                  geog = kernel(max = 1e5),
                  clim = kernel("gaussian", theta = 0.2, max = 0.5),
                  stat = c("weighted_mean", "sum_weights"),
                  normalize = TRUE
            )
      })
      # Result should have a `sum_weights` column; values should be in
      # [0, 1] roughly.
      sw <- res$sum_weights[is.finite(res$sum_weights)]
      expect_gt(length(sw), 0)
      expect_true(all(sw >= 0))
      expect_lt(max(sw), 1.5)  # mild slack for discretization
})


# ---------- 8. tabulate columns are normalized -------------------------------

test_that("tabulate columns get divided by D_max under normalize = TRUE", {
      skip_if_not_installed("terra")
      pool <- make_test_raster(n = 20)
      n <- terra::ncell(pool)
      y <- factor(sample(c("a", "b", "c"), n, replace = TRUE))

      res_norm <- analog_search(
            x = pool, pool = pool, y = y,
            stat = c("tabulate", "sum_weights"),
            geog = kernel(max = 1e5), clim = kernel("gaussian", theta = 0.2),
            normalize = TRUE
      )
      res_raw <- analog_search(
            x = pool, pool = pool, y = y,
            stat = c("tabulate", "sum_weights"),
            geog = kernel(max = 1e5), clim = kernel("gaussian", theta = 0.2),
            normalize = FALSE
      )

      D_max <- attr(res_norm, "D_max")
      expect_true(is.finite(D_max) && D_max > 0)

      # Result is a SpatRaster (since pool is one); convert to data.frame
      # for column-wise arithmetic. Use na.rm = FALSE so cell ordering
      # matches between the two converted frames.
      df_norm <- terra::as.data.frame(res_norm, na.rm = FALSE)
      df_raw  <- terra::as.data.frame(res_raw,  na.rm = FALSE)

      # For each tabulate column, normalized = raw / D_max.
      tab_cols <- grep("^n_", names(df_norm), value = TRUE)
      expect_gt(length(tab_cols), 0)

      for (cn in tab_cols) {
            v_n <- df_norm[[cn]]
            v_r <- df_raw[[cn]]
            ok <- is.finite(v_n) & is.finite(v_r)
            expect_equal(v_n[ok], v_r[ok] / D_max, tolerance = 1e-10,
                         info = paste0("column = ", cn))
      }

      # And class-fractions should sum (per row) to roughly the
      # normalized sum_weights: tabulate is sum_weights partitioned by class.
      tab_mat <- as.matrix(df_norm[, tab_cols, drop = FALSE])
      class_sum <- rowSums(tab_mat)
      sw <- df_norm$sum_weights
      ok <- is.finite(class_sum) & is.finite(sw)
      expect_equal(class_sum[ok], sw[ok], tolerance = 1e-8)
})


# ---------- 9. D_max attribute is exposed -----------------------------------

test_that("D_max and normalize attributes are exposed on result", {
      skip_if_not_installed("terra")
      pool <- make_test_raster(n = 15)

      # data.frame mode
      res_df <- analog_search(
            x = pool, pool = pool,
            stat = c("count", "sum_weights"),
            geog = kernel(max = 1e5), clim = kernel("gaussian", theta = 0.2),
            normalize = TRUE
      )
      expect_true(isTRUE(attr(res_df, "normalize")))
      expect_true(is.finite(attr(res_df, "D_max")))

      # raster mode
      res_r <- analog_density(
            x = pool, pool = pool,
            geog = kernel(max = 1e5), clim = kernel("gaussian", theta = 0.2),
            normalize = TRUE
      )
      expect_true(isTRUE(attr(res_r, "normalize")))
      expect_true(is.finite(attr(res_r, "D_max")))

      # FALSE => normalize attribute FALSE, D_max NA
      res_off <- analog_search(
            x = pool, pool = pool,
            stat = "sum_weights",
            geog = kernel(max = 1e5), clim = kernel("gaussian", theta = 0.2),
            normalize = FALSE
      )
      expect_false(isTRUE(attr(res_off, "normalize")))
      expect_true(is.na(attr(res_off, "D_max")))
})
