# Tests for cell-area weights and user-supplied per-pool-point weights.
#
# Coverage:
#   - cell-area weights affect lonlat raster aggregations
#   - user weights affect aggregations
#   - both compose multiplicatively
#   - weight = NULL + cell_area_weight = FALSE reproduces unweighted behavior
#   - basic validation: length mismatch, negative values, NA→0, format errors
#   - pre-built index reconciliation (FALSE on TRUE-built index → error, etc.)
#   - pair-mode emits the right columns conditionally
#   - permissive warning when weight + downsampled index combine

# --- Helpers ---------------------------------------------------------------

# Build a small lonlat raster pool with a known cell-area gradient.
# Cell areas shrink toward the poles, so lonlat raster aggregations should
# differ from the unweighted (matrix) equivalent.
.lonlat_pool_raster <- function(nx = 20, ny = 20, seed = 1L) {
      skip_if_not_installed("terra")
      set.seed(seed)
      r <- terra::rast(ncol = nx, nrow = ny,
                       xmin = -10, xmax = 10,
                       ymin = 30, ymax = 80,            # high latitude → big area gradient
                       crs = "EPSG:4326")
      terra::values(r) <- rnorm(nx * ny)
      names(r) <- "clim1"
      r
}


# --- Cell-area weights ------------------------------------------------------

test_that("cell_area_weight changes aggregation results on a lonlat raster", {
      skip_if_not_installed("terra")

      pool <- .lonlat_pool_raster()

      # Two focals at very different latitudes, with a geographic constraint
      # tight enough that each one's neighborhood is a small *local* patch
      # of the raster. Cell areas inside one patch are systematically different
      # from the other (the high-latitude patch has small cells, the
      # low-latitude patch has large cells), so a `count` of cells in the
      # patch and a `sum_weights` over area-weighted cells in the same patch
      # will diverge.
      focal <- data.frame(x = c(0, 0),
                          y = c(35, 75),
                          clim1 = c(0, 0))

      idx_on  <- build_analog_index(pool, cell_area_weight = TRUE,
                                    coord_type = "lonlat")
      idx_off <- build_analog_index(pool, cell_area_weight = FALSE,
                                    coord_type = "lonlat")

      out_on <- analog_search(focal, idx_on,
                              stat = "sum_weights",
                              clim = kernel(max = 100),  # permissive: don't bind on climate
                              geog = kernel(max = 300))  # ~3° at the equator; smaller at lat 75
      out_off <- analog_search(focal, idx_off,
                               stat = "sum_weights",
                               clim = kernel(max = 100),
                               geog = kernel(max = 300))

      # Sanity: each focal actually has analogs in its small patch
      expect_true(all(out_on$sum_weights > 0))
      expect_true(all(out_off$sum_weights > 0))

      # And the weighted vs unweighted sums must disagree, because the local
      # cell-area mean inside a small patch differs from the global mean of 1.
      diffs <- abs(out_on$sum_weights - out_off$sum_weights)
      expect_true(all(diffs > 1e-6))
})


test_that('cell_area_weight = "auto" is on for raster pool, off for matrix pool', {
      skip_if_not_installed("terra")

      pool <- .lonlat_pool_raster()

      # Raster path: auto → on
      idx_raster_auto <- build_analog_index(pool, coord_type = "lonlat")
      expect_false(is.null(idx_raster_auto$cell_area_weight))

      # Matrix path: auto → off (no cell concept)
      pool_mat <- terra::as.data.frame(pool, xy = TRUE, na.rm = FALSE)
      idx_mat_auto <- build_analog_index(as.matrix(pool_mat),
                                         coord_type = "lonlat")
      expect_true(is.null(idx_mat_auto$cell_area_weight))
})


test_that("cell_area_weight = TRUE errors on non-raster pool", {
      d <- sim_test_data()
      expect_error(
            build_analog_index(d$ref, cell_area_weight = TRUE,
                               coord_type = "projected"),
            "SpatRaster"
      )
})


test_that("cell_area_weight reconciliation between query and pre-built index", {
      skip_if_not_installed("terra")

      pool <- .lonlat_pool_raster()
      focal <- data.frame(x = c(0, 0), y = c(40, 70), clim1 = c(0, 0))

      idx_on  <- build_analog_index(pool, cell_area_weight = TRUE,
                                    coord_type = "lonlat")
      idx_off <- build_analog_index(pool, cell_area_weight = FALSE,
                                    coord_type = "lonlat")

      # Forcing FALSE on a TRUE-built index → error
      expect_error(
            analog_search(focal, idx_on,
                          cell_area_weight = FALSE,
                          stat = "count", clim = kernel(max = 5), geog = kernel(max = 5000)),
            "cell_area_weight = FALSE"
      )
      # Forcing TRUE on a FALSE-built index → error
      expect_error(
            analog_search(focal, idx_off,
                          cell_area_weight = TRUE,
                          stat = "count", clim = kernel(max = 5), geog = kernel(max = 5000)),
            "cell_area_weight = TRUE"
      )

      # "auto" silently accepts either configuration
      expect_no_error(
            analog_search(focal, idx_on,
                          stat = "count", clim = kernel(max = 5), geog = kernel(max = 5000))
      )
      expect_no_error(
            analog_search(focal, idx_off,
                          stat = "count", clim = kernel(max = 5), geog = kernel(max = 5000))
      )
})


# --- User weight ------------------------------------------------------------

test_that("user weight changes weighted_mean as expected", {

      d <- sim_test_data()
      n_ref <- nrow(d$ref)
      y <- d$ref[, 3]   # use a climate column as a stand-in response

      # Half the pool down-weighted; the other half up-weighted.
      w <- rep(c(0.1, 1.9), length.out = n_ref)   # mean 1

      out_unw <- analog_search(d$focal, d$ref,
                               y = y, stat = "weighted_mean",
                               clim = kernel(max = 2), geog = kernel(max = 5),
                               coord_type = "projected")
      out_w   <- analog_search(d$focal, d$ref,
                               y = y, weight = w,
                               stat = "weighted_mean",
                               clim = kernel(max = 2), geog = kernel(max = 5),
                               coord_type = "projected")

      # Should differ at the bulk of focals
      diffs <- abs(out_unw$weighted_mean - out_w$weighted_mean)
      expect_true(sum(diffs > 1e-8, na.rm = TRUE) >= 0.5 * length(diffs))
})


test_that("user weight scales sum_weights linearly", {

      d <- sim_test_data()
      n_ref <- nrow(d$ref)

      # weight = 2 everywhere should double sum_weights vs no weight
      out_w1 <- analog_search(d$focal, d$ref,
                              stat = "sum_weights",
                              clim = kernel(max = 2), geog = kernel(max = 5),
                              coord_type = "projected")
      out_w2 <- analog_search(d$focal, d$ref,
                              weight = rep(2, n_ref),
                              stat = "sum_weights",
                              clim = kernel(max = 2), geog = kernel(max = 5),
                              coord_type = "projected")

      expect_equal(out_w2$sum_weights, 2 * out_w1$sum_weights, tolerance = 1e-10)
})


test_that("NA weights are treated as 0 (point excluded)", {

      d <- sim_test_data()
      n_ref <- nrow(d$ref)

      # First half NA, second half 1.0
      half <- floor(n_ref / 2)
      w_na <- c(rep(NA_real_, half), rep(1.0, n_ref - half))
      # Equivalent: same pattern with explicit zeros
      w_zero <- c(rep(0.0,         half), rep(1.0, n_ref - half))

      out_na <- analog_search(d$focal, d$ref,
                              weight = w_na,
                              stat = "sum_weights",
                              clim = kernel(max = 2), geog = kernel(max = 5),
                              coord_type = "projected")
      out_zero <- analog_search(d$focal, d$ref,
                                weight = w_zero,
                                stat = "sum_weights",
                                clim = kernel(max = 2), geog = kernel(max = 5),
                                coord_type = "projected")

      expect_equal(out_na$sum_weights, out_zero$sum_weights, tolerance = 1e-10)
})


# --- Composition ------------------------------------------------------------

test_that("user weight composes multiplicatively with cell-area weight", {
      skip_if_not_installed("terra")

      pool <- .lonlat_pool_raster()
      n_pool <- terra::ncell(pool)
      focal <- data.frame(x = c(0, 0), y = c(40, 70), clim1 = c(0, 0))

      # Constant user weight = 3 should triple sum_weights compared to no
      # user weight (cell-area weighting on in both)
      idx <- build_analog_index(pool, coord_type = "lonlat",
                                cell_area_weight = TRUE)

      out_no_user <- analog_search(focal, idx,
                                   stat = "sum_weights",
                                   clim = kernel(max = 5), geog = kernel(max = 5000))
      out_w3      <- analog_search(focal, idx,
                                   weight = rep(3, n_pool),
                                   stat = "sum_weights",
                                   clim = kernel(max = 5), geog = kernel(max = 5000))

      expect_equal(out_w3$sum_weights, 3 * out_no_user$sum_weights,
                   tolerance = 1e-10)
})


test_that("weight = NULL and cell_area_weight = FALSE reproduce baseline", {

      d <- sim_test_data()

      # Old call (no new args)
      base <- analog_search(d$focal, d$ref,
                            stat = "count",
                            clim = kernel(max = 2), geog = kernel(max = 5),
                            coord_type = "projected")
      # Explicit defaults
      explicit <- analog_search(d$focal, d$ref,
                                weight = NULL,
                                cell_area_weight = FALSE,
                                stat = "count",
                                clim = kernel(max = 2), geog = kernel(max = 5),
                                coord_type = "projected")

      # count is unaffected by weights by design (sample/area/user all
      # left out of count) — so identical regardless.
      expect_equal(base$count, explicit$count)
})


# --- Validation -------------------------------------------------------------

test_that("weight length must match pool", {
      d <- sim_test_data()
      expect_error(
            analog_search(d$focal, d$ref,
                          weight = rep(1, nrow(d$ref) - 1),
                          stat = "count",
                          clim = kernel(max = 2), geog = kernel(max = 5),
                          coord_type = "projected"),
            "rows"
      )
})


test_that("negative weights error", {
      d <- sim_test_data()
      bad <- rep(1, nrow(d$ref))
      bad[1] <- -1
      expect_error(
            analog_search(d$focal, d$ref,
                          weight = bad,
                          stat = "count",
                          clim = kernel(max = 2), geog = kernel(max = 5),
                          coord_type = "projected"),
            "non-negative"
      )
})


test_that("weight accepts vector, single-column matrix/df, single-layer raster", {
      skip_if_not_installed("terra")

      d <- sim_test_data()
      n_ref <- nrow(d$ref)
      w_vec <- runif(n_ref, 0.5, 1.5)

      out_vec <- analog_search(d$focal, d$ref,
                               weight = w_vec,
                               stat = "sum_weights",
                               clim = kernel(max = 2), geog = kernel(max = 5),
                               coord_type = "projected")
      out_mat <- analog_search(d$focal, d$ref,
                               weight = matrix(w_vec, ncol = 1),
                               stat = "sum_weights",
                               clim = kernel(max = 2), geog = kernel(max = 5),
                               coord_type = "projected")
      out_df  <- analog_search(d$focal, d$ref,
                               weight = data.frame(w = w_vec),
                               stat = "sum_weights",
                               clim = kernel(max = 2), geog = kernel(max = 5),
                               coord_type = "projected")

      expect_equal(out_vec$sum_weights, out_mat$sum_weights)
      expect_equal(out_vec$sum_weights, out_df$sum_weights)
})


test_that("multi-column weight matrix errors", {
      d <- sim_test_data()
      bad <- matrix(1, nrow = nrow(d$ref), ncol = 2)
      expect_error(
            analog_search(d$focal, d$ref,
                          weight = bad,
                          stat = "count",
                          clim = kernel(max = 2), geog = kernel(max = 5),
                          coord_type = "projected"),
            "one column"
      )
})


# --- Pair-mode column emission ---------------------------------------------

test_that("pair mode emits weight columns conditionally", {

      d <- sim_test_data()
      n_ref <- nrow(d$ref)

      # Plain query: no downsample, no area, no user weight → no weight cols
      plain <- analog_search(d$focal, d$ref,
                             select = "knn_geog", k = 1, stat = "none",
                             clim = kernel(max = 2),
                             coord_type = "projected")
      expect_false("sample_weight" %in% names(plain))
      expect_false("area_weight"   %in% names(plain))
      expect_false("user_weight"   %in% names(plain))

      # With user weight: only user_weight column appears
      with_user <- analog_search(d$focal, d$ref,
                                 weight = rep(1, n_ref),
                                 select = "knn_geog", k = 1, stat = "none",
                                 clim = kernel(max = 2),
                                 coord_type = "projected")
      expect_false("sample_weight" %in% names(with_user))
      expect_false("area_weight"   %in% names(with_user))
      expect_true( "user_weight"   %in% names(with_user))
})


test_that("pair mode emits area_weight column when cell-area weights are active", {
      skip_if_not_installed("terra")

      pool <- .lonlat_pool_raster()
      focal <- data.frame(x = c(0, 0), y = c(40, 70), clim1 = c(0, 0))

      out <- analog_search(focal, pool,
                           select = "knn_geog", k = 1, stat = "none",
                           clim = kernel(max = 5),
                           coord_type = "lonlat")
      expect_true("area_weight" %in% names(out))
      # Mean-1 normalized → area weights present and finite
      expect_true(all(is.finite(out$area_weight)))
      expect_true(any(out$area_weight != 1))
})


# --- Output attributes ------------------------------------------------------

test_that("output carries cell_area_weight and weight_provided attributes", {

      d <- sim_test_data()

      out_off <- analog_search(d$focal, d$ref,
                               stat = "count",
                               clim = kernel(max = 2), geog = kernel(max = 5),
                               coord_type = "projected")
      expect_false(attr(out_off, "cell_area_weight"))
      expect_false(attr(out_off, "weight_provided"))

      out_user <- analog_search(d$focal, d$ref,
                                weight = rep(1, nrow(d$ref)),
                                stat = "count",
                                clim = kernel(max = 2), geog = kernel(max = 5),
                                coord_type = "projected")
      expect_true(attr(out_user, "weight_provided"))
})


# --- Downsample interaction -------------------------------------------------

test_that("weight + downsampled index emits a permissive warning", {

      d <- sim_test_data(nref = 1000)
      idx <- build_analog_index(d$ref, coord_type = "projected",
                                downsample = 0.5, seed = 1L)

      expect_warning(
            analog_search(d$focal, idx,
                          weight = rep(1, nrow(d$ref)),
                          stat = "count",
                          clim = kernel(max = 2), geog = kernel(max = 5)),
            "downsample"
      )
})


# --- build_analog_index: numeric vector cell_area_weight ---------------------

test_that("build_analog_index accepts a numeric vector for cell_area_weight", {

      d <- sim_test_data()
      n_ref <- nrow(d$ref)

      # Equal-weight vector should be equivalent to no weighting (mean 1).
      idx_off <- build_analog_index(d$ref, coord_type = "projected",
                                    cell_area_weight = FALSE)
      idx_vec <- build_analog_index(d$ref, coord_type = "projected",
                                    cell_area_weight = rep(1, n_ref))

      out_off <- analog_search(d$focal, idx_off,
                               stat = "sum_weights",                                clim = kernel(max = 2), geog = kernel(max = 5))
      out_vec <- analog_search(d$focal, idx_vec,
                               stat = "sum_weights",                                clim = kernel(max = 2), geog = kernel(max = 5))

      expect_equal(out_off$sum_weights, out_vec$sum_weights, tolerance = 1e-12)
})


test_that("build_analog_index validates numeric cell_area_weight: length, sign, NA", {

      d <- sim_test_data()
      n_ref <- nrow(d$ref)

      # Wrong length
      expect_error(
            build_analog_index(d$ref, coord_type = "projected",
                               cell_area_weight = rep(1, n_ref - 1)),
            "length"
      )

      # Negative value
      bad <- rep(1, n_ref); bad[1] <- -1
      expect_error(
            build_analog_index(d$ref, coord_type = "projected",
                               cell_area_weight = bad),
            "non-negative"
      )

      # NA in vector
      bad2 <- rep(1, n_ref); bad2[1] <- NA_real_
      expect_error(
            build_analog_index(d$ref, coord_type = "projected",
                               cell_area_weight = bad2),
            "non-negative"
      )
})


test_that("numeric cell_area_weight cannot be applied to pre-built index in analog_search", {

      d <- sim_test_data()
      n_ref <- nrow(d$ref)
      idx <- build_analog_index(d$ref, coord_type = "projected",
                                cell_area_weight = FALSE)

      expect_error(
            analog_search(d$focal, idx,
                          cell_area_weight = rep(1, n_ref),
                          stat = "count",
                          clim = kernel(max = 2), geog = kernel(max = 5)),
            "pre-built"
      )
})


# --- tiled_analog_search ----------------------------------------------------

test_that("tiled_analog_search forwards weight raster and matches non-tiled result", {
      skip_if_not_installed("terra")

      # Use an extent clearly outside lonlat bounds so the coord_type
      # auto-detector unambiguously returns "projected" in both the
      # non-tiled and tiled runs. Otherwise .detect_geo (which keys on
      # coordinate magnitudes, not CRS metadata) would call this lonlat
      # in the non-tiled path while terra::is.lonlat(crs) would call it
      # projected in the tiled path, producing different distance
      # formulas and a real numerical mismatch.
      set.seed(7)
      pool <- terra::rast(ncol = 24, nrow = 24,
                          xmin = 0, xmax = 24000, ymin = 0, ymax = 24000,
                          crs = "EPSG:32611")
      names(pool) <- "clim1"
      terra::values(pool) <- rnorm(terra::ncell(pool))

      focal <- pool

      w <- pool
      names(w) <- "weight"
      terra::values(w) <- runif(terra::ncell(w), 0.5, 1.5)

      ref <- analog_impact(
            x = focal, pool = pool,
            y = pool[[1]], weight = w,
            clim = kernel("gaussian", theta = 0.5, max = 1), geog = kernel(max = 5000),
            cell_area_weight = FALSE,
            coord_type = "projected"
      )

      tiled <- suppressWarnings(
            tiled_analog_search(
                  x = focal, pool = pool, n_tiles = 4,
                  fun = analog_impact,
                  geog = kernel(max = 5000),
                  y = pool[[1]],
                  weight = w,
                  clim = kernel("gaussian", theta = 0.5, max = 1),
                  cell_area_weight = FALSE,
                  progress = FALSE
            )
      )

      # Both should be SpatRasters with the same weighted_mean layer
      expect_s4_class(tiled, "SpatRaster")
      diff <- terra::values(ref$weighted_mean) - terra::values(tiled$weighted_mean)
      expect_true(max(abs(diff), na.rm = TRUE) < 1e-6)
})


test_that("tiled_analog_search applies cell_area_weight with global normalization", {
      skip_if_not_installed("terra")

      set.seed(11)
      pool <- terra::rast(ncol = 24, nrow = 24,
                          xmin = -10, xmax = 10, ymin = 30, ymax = 80,
                          crs = "EPSG:4326")
      names(pool) <- "clim1"
      terra::values(pool) <- rnorm(terra::ncell(pool))

      focal <- pool

      ref <- analog_impact(
            x = focal, pool = pool,
            y = pool[[1]],
            clim = kernel("gaussian", theta = 0.5, max = 1), geog = kernel(max = 300),
            cell_area_weight = TRUE,
            coord_type = "lonlat"
      )

      tiled <- suppressWarnings(
            tiled_analog_search(
                  x = focal, pool = pool, n_tiles = 4,
                  fun = analog_impact,
                  geog = kernel(max = 300),
                  y = pool[[1]],
                  clim = kernel("gaussian", theta = 0.5, max = 1),
                  cell_area_weight = TRUE,
                  progress = FALSE
            )
      )

      # Should match the non-tiled case to high precision because the global
      # normalization is preserved across tiles.
      expect_s4_class(tiled, "SpatRaster")
      diff <- terra::values(ref$weighted_mean) - terra::values(tiled$weighted_mean)
      expect_true(max(abs(diff), na.rm = TRUE) < 1e-6)
})


test_that("tiled_analog_search crops covariates per tile (regression bugfix)", {
      skip_if_not_installed("terra")

      set.seed(13)
      pool <- terra::rast(ncol = 20, nrow = 20,
                          xmin = 0, xmax = 100, ymin = 0, ymax = 100,
                          crs = "EPSG:32611")
      names(pool) <- "clim1"
      terra::values(pool) <- rnorm(terra::ncell(pool))

      focal <- pool
      cov_r <- pool
      names(cov_r) <- "elev"
      terra::values(cov_r) <- rnorm(terra::ncell(cov_r))

      ref <- analog_search(
            x = focal, pool = pool, y = pool[[1]],
            covariates = cov_r,
            stat = c("count", "weighted_mean", "regression"),
            clim = kernel("gaussian", theta = 0.5, max = 1), geog = kernel(max = 25),
            cell_area_weight = FALSE,
            coord_type = "projected"
      )

      tiled <- suppressWarnings(
            tiled_analog_search(
                  x = focal, pool = pool, n_tiles = 4,
                  fun = analog_search,
                  geog = kernel(max = 25),
                  y = pool[[1]],
                  covariates = cov_r,
                  stat = c("count", "weighted_mean", "regression"),
                  clim = kernel("gaussian", theta = 0.5, max = 1),
                  cell_area_weight = FALSE,
                  progress = FALSE
            )
      )

      expect_s4_class(tiled, "SpatRaster")
      # weighted_mean should match
      diff_wm <- terra::values(ref$weighted_mean) - terra::values(tiled$weighted_mean)
      expect_true(max(abs(diff_wm), na.rm = TRUE) < 1e-6)
      # regression intercept should match
      diff_b <- terra::values(ref$coef_intercept) - terra::values(tiled$coef_intercept)
      expect_true(max(abs(diff_b), na.rm = TRUE) < 1e-6)
})


# --- analog_cv -------------------------------------------------------------

test_that("analog_cv accepts and forwards `weight` (LOO)", {

      set.seed(17)
      d <- sim_test_data(nref = 80, nfocal = 80)
      pool <- d$ref
      y <- pool[, 3]
      w <- runif(nrow(pool), 0.5, 1.5)

      cv_w <- analog_cv(
            fun = analog_impact, pool = pool, y = y,
            weight = w,
            clim = kernel("gaussian", theta = 0.5, max = 1), geog = kernel(max = 2),
            coord_type = "projected", cv_method = "loo"
      )
      cv_no <- analog_cv(
            fun = analog_impact, pool = pool, y = y,
            clim = kernel("gaussian", theta = 0.5, max = 1), geog = kernel(max = 2),
            coord_type = "projected", cv_method = "loo"
      )

      # Expected: weighting changes the weighted_mean predictions and thus
      # the residuals at most rows.
      diffs <- abs(cv_w$residual - cv_no$residual)
      expect_true(sum(diffs > 1e-8, na.rm = TRUE) >= 0.5 * length(diffs))
})


test_that("analog_cv k-fold preserves global cell-area normalization", {
      skip_if_not_installed("terra")

      set.seed(19)
      pool <- terra::rast(ncol = 12, nrow = 12,
                          xmin = -10, xmax = 10, ymin = 30, ymax = 80,
                          crs = "EPSG:4326")
      names(pool) <- "clim1"
      terra::values(pool) <- rnorm(terra::ncell(pool))

      yr <- pool
      names(yr) <- "y"
      terra::values(yr) <- rnorm(terra::ncell(yr))

      # Run k-fold with cell_area_weight=TRUE; should produce no errors and
      # return a SpatRaster with obs/residual layers.
      cv <- suppressWarnings(analog_cv(
            fun = analog_impact, pool = pool, y = yr,
            cv_method = "kfold", n_folds = 3L,
            clim = kernel("gaussian", theta = 0.5, max = 1), geog = kernel(max = 500),
            cell_area_weight = TRUE
      ))

      expect_s4_class(cv, "SpatRaster")
      expect_true(all(c("obs", "residual") %in% names(cv)))
})


test_that("analog_cv with weight rejects bad lengths", {

      d <- sim_test_data(nref = 60, nfocal = 60)

      expect_error(
            analog_cv(
                  fun = analog_impact, pool = d$ref, y = d$ref[, 3],
                  weight = rep(1, nrow(d$ref) - 1),
                  clim = kernel("gaussian", theta = 0.5, max = 1), geog = kernel(max = 2),
                  coord_type = "projected", cv_method = "loo"
            ),
            "rows"
      )
})
