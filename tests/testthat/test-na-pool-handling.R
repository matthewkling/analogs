# Tests for NA handling in pool and focal data.
#
# Layers covered:
#
# 1. .format_data() strips NA rows from any input (pool OR focal) and
#    attaches a row_map attribute mapping kept-row -> original-row indices.
#
# 2. build_analog_index() captures pool_row_map, exposes n_pool (original)
#    and n_pool_used (post-strip), and aligns cell_area_weight to ref_data
#    after stripping.
#
# 3. Pair-mode queries return analog_index values that index correctly into
#    the *original* pool. Aggregation queries against an NA-containing pool
#    produce results identical to the same query against a manually-stripped
#    pool (proving NA pollution has been eliminated).
#
# 4. Defensive C++ guard: build_analog_index_cpp doesn't crash or return NaN
#    metadata when handed a matrix with NA rows directly.
#
# 5. Empty-pool edge case (all rows NA) doesn't crash and produces all-NA
#    or zero-count results.
#
# 6. Focal NA stripping: query_analog_index strips NA focal rows; results
#    are reconstructed to the original focal length for single-row-per-focal
#    queries (so SpatRaster outputs rasterize correctly), and the `index`
#    column in pairs output is remapped back to original-focal indexing in
#    multi-row-per-focal queries.

# Helpers ------------------------------------------------------------------

# Build a small reproducible pool with specified NA pattern.
sim_pool_with_nas <- function(n = 200, seed = 42, na_frac = 0.3) {
      set.seed(seed)
      m <- matrix(rnorm(n * 4), ncol = 4)
      colnames(m) <- c("x", "y", "env1", "env2")
      n_na <- round(n * na_frac)
      if (n_na > 0L) {
            na_rows <- sample.int(n, n_na)
            # Mix of NA-coord and NA-environmental to exercise both patterns.
            half <- ceiling(n_na / 2)
            m[na_rows[seq_len(half)], "env1"] <- NA_real_
            if (length(na_rows) > half) {
                  m[na_rows[(half + 1):length(na_rows)], "x"] <- NA_real_
            }
      }
      m
}

# Manually strip NA rows the same way the package should.
strip_na_pool <- function(pool) {
      keep <- stats::complete.cases(pool)
      list(
            pool = pool[keep, , drop = FALSE],
            row_map = which(keep)
      )
}


# 1. .format_data behavior --------------------------------------------------

# Note: previous behavior had a `purpose` argument distinguishing focal
# (preserve NAs) from pool (strip NAs). Both sides now strip identically.

test_that(".format_data strips NA rows and sets row_map", {
      pool <- sim_pool_with_nas(n = 100, na_frac = 0.4)
      n_complete <- sum(stats::complete.cases(pool))

      out <- analogs:::.format_data(pool)

      expect_equal(nrow(out), n_complete)
      expect_true(!is.null(attr(out, "row_map")))
      expect_equal(length(attr(out, "row_map")), n_complete)
      expect_true(all(stats::complete.cases(out)))

      # row_map should index back into original input correctly
      row_map <- attr(out, "row_map")
      expect_equal(out[, 3], pool[row_map, "env1"])
      expect_equal(out[, 4], pool[row_map, "env2"])
})

test_that(".format_data with no NAs sets row_map = NULL", {
      pool <- sim_pool_with_nas(n = 100, na_frac = 0)

      out <- analogs:::.format_data(pool)

      expect_equal(nrow(out), 100)
      expect_null(attr(out, "row_map"))
})

test_that(".format_data strips fully-NA input to zero rows", {
      pool <- sim_pool_with_nas(n = 50, na_frac = 1.0)

      out <- analogs:::.format_data(pool)

      expect_equal(nrow(out), 0)
      # row_map is integer(0), NOT NULL, when stripping happened (even if to zero)
      expect_true(!is.null(attr(out, "row_map")))
      expect_equal(length(attr(out, "row_map")), 0L)
})

test_that(".format_data works with SpatRaster input", {
      skip_if_not_installed("terra")

      r <- terra::rast(nrows = 10, ncols = 10, nlyrs = 2,
                       xmin = 0, xmax = 10, ymin = 0, ymax = 10)
      set.seed(123)
      terra::values(r) <- matrix(rnorm(200), ncol = 2)
      vals <- terra::values(r)
      vals[c(5, 17, 33, 50, 71), 1] <- NA
      vals[c(8, 25, 60), 2] <- NA
      terra::values(r) <- vals

      out <- analogs:::.format_data(r)

      n_complete <- sum(stats::complete.cases(terra::values(r)))
      expect_equal(nrow(out), n_complete)

      row_map <- attr(out, "row_map")
      expect_equal(length(row_map), n_complete)
      expect_true(all(stats::complete.cases(out)))
})


# 2. build_analog_index integration -----------------------------------------

test_that("build_analog_index stores pool_row_map and n_pool_used", {
      pool <- sim_pool_with_nas(n = 200, na_frac = 0.3)
      n_complete <- sum(stats::complete.cases(pool))

      idx <- build_analog_index(pool, coord_type = "projected")

      expect_equal(idx$n_pool, 200L)              # original size
      expect_equal(idx$n_pool_used, n_complete)   # post-strip
      expect_equal(nrow(idx$ref_data), n_complete)
      expect_true(!is.null(idx$pool_row_map))
      expect_equal(length(idx$pool_row_map), n_complete)
      expect_true(all(stats::complete.cases(idx$ref_data)))
})

test_that("build_analog_index without NAs has pool_row_map NULL and equal n_pool", {
      pool <- sim_pool_with_nas(n = 100, na_frac = 0)

      idx <- build_analog_index(pool, coord_type = "projected")

      expect_null(idx$pool_row_map)
      expect_equal(idx$n_pool, 100L)
      expect_equal(idx$n_pool_used, 100L)
})

test_that("build_analog_index strips cell_area_weight to align with ref_data", {
      skip_if_not_installed("terra")

      r <- terra::rast(nrows = 10, ncols = 10, nlyrs = 2,
                       xmin = -1, xmax = 1, ymin = -1, ymax = 1,
                       crs = "EPSG:4326")
      set.seed(123)
      vals <- matrix(rnorm(200), ncol = 2)
      vals[c(3, 17, 50, 71), 1] <- NA
      terra::values(r) <- vals

      # Provide a length-100 area weight; expect stripped to length n_pool_used
      area_w <- runif(100, 0.5, 1.5)

      idx <- build_analog_index(r, coord_type = "lonlat",
                                cell_area_weight = area_w)

      expect_equal(length(idx$cell_area_weight), idx$n_pool_used)
      # Stored values are subset of input by row_map (no re-normalization
      # is applied to user-supplied weight vectors)
      expect_equal(idx$cell_area_weight, area_w[idx$pool_row_map])
})

test_that("build_analog_index errors on cell_area_weight length mismatch", {
      pool <- sim_pool_with_nas(n = 100, na_frac = 0.3)

      # Wrong length: not equal to original (100)
      expect_error(
            build_analog_index(pool, coord_type = "projected",
                               cell_area_weight = runif(50)),
            "does not match pool size"
      )
})


# 3. Round-trip: pair-mode analog_index references original pool ------------

test_that("analog_index in pair mode references the ORIGINAL (unstripped) pool", {
      set.seed(11)
      pool <- sim_pool_with_nas(n = 200, na_frac = 0.3)
      focal <- matrix(rnorm(20 * 4), ncol = 4)

      res <- analog_search(
            x = focal,
            pool = pool,
            select = "knn_env",
            geog = kernel(max = 5),
            k = 3,
            coord_type = "projected"
      )

      ai <- res$analog_index
      ai_valid <- ai[!is.na(ai)]

      # Every returned analog_index must reference a valid row of ORIGINAL pool
      expect_true(all(ai_valid >= 1 & ai_valid <= nrow(pool)))
      # And must not point at NA-containing rows
      expect_true(all(stats::complete.cases(pool[ai_valid, , drop = FALSE])))

      # analog_x / analog_y must match pool rows at those indices
      expect_equal(res$analog_x[!is.na(ai)],
                   pool[ai_valid, "x"])
      expect_equal(res$analog_y[!is.na(ai)],
                   pool[ai_valid, "y"])
})

test_that("pair-mode results equivalent between NA-containing and pre-stripped pool", {
      set.seed(22)
      pool_full <- sim_pool_with_nas(n = 300, na_frac = 0.3)
      stripped  <- strip_na_pool(pool_full)
      focal <- matrix(rnorm(15 * 4), ncol = 4)

      res_na <- analog_search(
            x = focal, pool = pool_full,
            select = "knn_env", geog = kernel(max = 10), k = 2,
            coord_type = "projected"
      )
      res_clean <- analog_search(
            x = focal, pool = stripped$pool,
            select = "knn_env", geog = kernel(max = 10), k = 2,
            coord_type = "projected"
      )

      # Same number of pairs returned
      expect_equal(nrow(res_na), nrow(res_clean))

      # Same focal-by-focal environmental/geog distances. (Tie-breaking on the
      # C++ side may pick different rows in the rare case of exact ties,
      # but distances will match either way for valid neighbors.)
      expect_equal(res_na$env_dist, res_clean$env_dist, tolerance = 1e-10)
      expect_equal(res_na$geog_dist, res_clean$geog_dist, tolerance = 1e-10)

      # Critical correctness check: every res_na$analog_index must point
      # back to a valid (non-NA) row of the ORIGINAL pool, and the analog
      # x/y must match the pool row at that index.
      ai_full  <- res_na$analog_index
      ai_valid <- ai_full[!is.na(ai_full)]
      expect_true(all(ai_valid >= 1 & ai_valid <= nrow(pool_full)))
      expect_true(all(stats::complete.cases(pool_full[ai_valid, , drop = FALSE])))
      expect_equal(res_na$analog_x[!is.na(ai_full)], pool_full[ai_valid, "x"])
      expect_equal(res_na$analog_y[!is.na(ai_full)], pool_full[ai_valid, "y"])
})


# 4. Aggregation results immune to NA pool pollution ------------------------

test_that("count / sum_weights / mean_weights NOT polluted by NA pool rows", {
      set.seed(33)
      pool_full <- sim_pool_with_nas(n = 300, na_frac = 0.4)
      stripped  <- strip_na_pool(pool_full)$pool
      focal <- matrix(rnorm(20 * 4), ncol = 4)

      args <- list(
            x = focal,
            select = "all",
            stat = c("count", "sum_weights", "mean_weights"),
            env = kernel("gaussian", theta = 1, max = 1.5),
            geog = kernel(max = 5),
            coord_type = "projected"
      )

      res_na    <- do.call(analog_search, c(args, list(pool = pool_full)))
      res_clean <- do.call(analog_search, c(args, list(pool = stripped)))

      expect_equal(res_na$count,        res_clean$count)
      expect_equal(res_na$sum_weights,  res_clean$sum_weights, tolerance = 1e-10)
      expect_equal(res_na$mean_weights, res_clean$mean_weights, tolerance = 1e-10)
})


# 5. y supplied at original-pool size --------------------------------------

test_that("y supplied at original-pool size is correctly aligned to ref_data", {
      set.seed(44)
      pool_full <- sim_pool_with_nas(n = 200, na_frac = 0.3)
      focal <- matrix(rnorm(10 * 4), ncol = 4)

      # User supplies y at ORIGINAL pool size (length 200)
      y_full <- runif(nrow(pool_full), 0, 100)

      res_na <- analog_search(
            x = focal, pool = pool_full,
            select = "all", stat = "weighted_mean",
            env = kernel("gaussian", theta = 1, max = 1.5),
            geog = kernel(max = 5),
            y = y_full,
            coord_type = "projected"
      )

      # Equivalent: stripped pool with manually stripped y
      stripped <- strip_na_pool(pool_full)
      res_clean <- analog_search(
            x = focal, pool = stripped$pool,
            select = "all", stat = "weighted_mean",
            env = kernel("gaussian", theta = 1, max = 1.5),
            geog = kernel(max = 5),
            y = y_full[stripped$row_map],
            coord_type = "projected"
      )

      expect_equal(res_na$weighted_mean, res_clean$weighted_mean,
                   tolerance = 1e-10)
})


# 6. Defensive C++ guard ----------------------------------------------------

test_that("C++ build does not return NaN ranges when given NA rows directly", {
      set.seed(55)
      m <- sim_pool_with_nas(n = 100, na_frac = 0.3)

      # Direct call to the C++ builder, bypassing .format_data. NaN rows
      # should be skipped in the min/max pass and bin assignment.
      res <- analogs:::build_analog_index_cpp(
            ref_mm = m,
            coord_type = "projected",
            geo_target = 10,
            env_target = 10,
            downsample = 1.0,
            seed = 55L
      )

      expect_true(all(is.finite(res$coord_mins)))
      expect_true(all(is.finite(res$coord_maxs)))
      expect_true(all(is.finite(res$env_mins)))
      expect_true(all(is.finite(res$env_maxs)))

      m_clean <- m[stats::complete.cases(m), , drop = FALSE]
      expect_equal(as.numeric(res$coord_mins), unname(apply(m_clean[, 1:2], 2, min)))
      expect_equal(as.numeric(res$coord_maxs), unname(apply(m_clean[, 1:2], 2, max)))
      expect_equal(as.numeric(res$env_mins),  unname(apply(m_clean[, 3:4], 2, min)))
      expect_equal(as.numeric(res$env_maxs),  unname(apply(m_clean[, 3:4], 2, max)))
})

test_that("C++ build with all-NA matrix returns empty index, no crash", {
      m <- sim_pool_with_nas(n = 50, na_frac = 1.0)

      res <- analogs:::build_analog_index_cpp(
            ref_mm = m,
            coord_type = "projected",
            geo_target = 10,
            env_target = 10,
            downsample = 1.0,
            seed = 1L
      )

      expect_equal(as.integer(res$n_bins_nonempty), 0L)
})


# 7. Empty-pool edge case (focal valid, pool entirely NA) -------------------

test_that("query against fully-NA pool returns NA results without crashing", {
      set.seed(66)
      pool_full <- sim_pool_with_nas(n = 50, na_frac = 1.0)
      focal     <- matrix(rnorm(10 * 4), ncol = 4)

      # Aggregation mode
      expect_silent(
            res_agg <- analog_search(
                  x = focal, pool = pool_full,
                  select = "all", stat = "count",
                  env = kernel(max = 1.5), geog = kernel(max = 5),
                  coord_type = "projected"
            )
      )
      expect_equal(nrow(res_agg), 10)
      # Empty pool means zero matches everywhere (count = 0, not NA, since
      # AggWorker initializes the count accumulator and just never finds
      # candidates).
      expect_true(all(res_agg$count == 0 | is.na(res_agg$count)))

      # Pairs mode
      expect_silent(
            res_pair <- analog_search(
                  x = focal, pool = pool_full,
                  select = "knn_env", k = 1, geog = kernel(max = 5),
                  coord_type = "projected"
            )
      )
      expect_equal(nrow(res_pair), 10)
      expect_true(all(is.na(res_pair$analog_index)))
})

# 8. Focal NA stripping -----------------------------------------------------

# Build a small reproducible focal with specified NA pattern (mirror of
# sim_pool_with_nas for clarity in focal-side tests).
sim_focal_with_nas <- function(n = 50, seed = 99, na_frac = 0.3) {
      sim_pool_with_nas(n = n, seed = seed, na_frac = na_frac)
}

test_that("agg-mode results match between NA-containing and pre-stripped focal", {
      set.seed(101)
      pool <- sim_pool_with_nas(n = 300, na_frac = 0)         # clean pool
      focal_full <- sim_focal_with_nas(n = 50, na_frac = 0.3)
      focal_clean <- focal_full[stats::complete.cases(focal_full), , drop = FALSE]
      kept <- which(stats::complete.cases(focal_full))

      args <- list(
            pool = pool,
            select = "all", stat = "count",
            env = kernel(max = 1.5), geog = kernel(max = 5),
            coord_type = "projected"
      )

      res_na    <- do.call(analog_search, c(args, list(x = focal_full)))
      res_clean <- do.call(analog_search, c(args, list(x = focal_clean)))

      # Output for NA-containing focal must have nrow == n_focal_original (50);
      # output for stripped focal has nrow == n_kept.
      expect_equal(nrow(res_na), nrow(focal_full))
      expect_equal(nrow(res_clean), nrow(focal_clean))

      # At kept positions, results should match
      expect_equal(res_na$count[kept], res_clean$count)
      # At stripped positions, results should be NA
      expect_true(all(is.na(res_na$count[-kept])))
      # x / y at kept positions match focal coords; NA at stripped positions
      expect_equal(res_na$x[kept], focal_full[kept, "x"])
      expect_true(all(is.na(res_na$x[-kept])))
})

test_that("pairs-mode k=1 reconstructs to n_focal_original rows", {
      set.seed(102)
      pool <- sim_pool_with_nas(n = 300, na_frac = 0)
      focal_full <- sim_focal_with_nas(n = 30, na_frac = 0.4)
      kept <- which(stats::complete.cases(focal_full))

      res <- analog_search(
            x = focal_full, pool = pool,
            select = "knn_env", k = 1, geog = kernel(max = 5),
            coord_type = "projected"
      )

      # Output should have nrow == n_focal_original
      expect_equal(nrow(res), nrow(focal_full))
      # Stripped focal positions should have NA pair data
      expect_true(all(is.na(res$analog_index[-kept])))
      expect_true(all(is.na(res$env_dist[-kept])))
      # Kept focal positions should have valid pair data
      expect_true(all(!is.na(res$analog_index[kept])))
})

test_that("pairs-mode k>1 remaps `index` column to original focal indexing", {
      set.seed(103)
      pool <- sim_pool_with_nas(n = 300, na_frac = 0)
      focal_full <- sim_focal_with_nas(n = 30, na_frac = 0.4)
      kept <- which(stats::complete.cases(focal_full))

      res <- analog_search(
            x = focal_full, pool = pool,
            select = "knn_env", k = 3, geog = kernel(max = 5),
            coord_type = "projected"
      )

      # `index` column should reference original focal positions.
      # Every value should be in `kept` (no stripped focal positions
      # appear, since stripped focals don't produce pair rows in k>1
      # mode -- they were dropped at ingestion).
      expect_true(all(res$index %in% kept))
      # And no NA focal indices
      expect_true(all(!is.na(res$index)))
})

test_that("focal SpatRaster with NA cells rasterizes correctly post-strip", {
      skip_if_not_installed("terra")

      # Pool: clean projected matrix
      set.seed(104)
      pool <- sim_pool_with_nas(n = 200, na_frac = 0)

      # Focal: SpatRaster with some NA cells
      r <- terra::rast(nrows = 8, ncols = 8, nlyrs = 2,
                       xmin = -2, xmax = 2, ymin = -2, ymax = 2)
      vals <- matrix(rnorm(128), ncol = 2)
      vals[c(3, 7, 12, 25, 40, 55), 1] <- NA
      vals[c(5, 22), 2] <- NA
      terra::values(r) <- vals

      kept <- which(stats::complete.cases(terra::values(r)))

      res <- analog_search(
            x = r, pool = pool,
            select = "all", stat = "count",
            env = kernel(max = 1.5), geog = kernel(max = 5),
            coord_type = "projected"
      )

      # Should rasterize back to a SpatRaster (single-row-per-focal agg mode)
      expect_true(inherits(res, "SpatRaster"))
      # Count layer values: NA at NA-focal cells, finite elsewhere
      cnt <- terra::values(res[["count"]])[, 1]
      expect_true(all(is.na(cnt[-kept])))
      expect_true(all(!is.na(cnt[kept])))
})

test_that("x_cov is correctly stripped to align with NA-stripped focal", {
      skip_if_not_installed("terra")

      # Build a focal raster with NA cells
      r <- terra::rast(nrows = 6, ncols = 6, nlyrs = 2,
                       xmin = -1, xmax = 1, ymin = -1, ymax = 1)
      set.seed(105)
      vals <- matrix(rnorm(72), ncol = 2)
      vals[c(3, 11, 22), 1] <- NA
      terra::values(r) <- vals
      kept <- which(stats::complete.cases(terra::values(r)))

      # x_cov: 3 components per focal cell (n_env = 2, so 2*(2+1)/2 = 3).
      # Layout (per .reconstruct_cov_matrix): diagonals first, then off-
      # diagonals row-major. For 2 environmental vars: [var(c1), var(c2), cov(c1,c2)].
      n_focal <- terra::ncell(r)
      x_cov <- matrix(c(rep(1, n_focal),    # var(c1)
                        rep(1, n_focal),    # var(c2)
                        rep(0, n_focal)),   # cov(c1,c2)
                      ncol = 3)
      # Set x_cov to NA at the same stripped positions to mirror typical
      # SpatRaster x_cov use (it would have NAs at ocean cells too).
      x_cov[c(3, 11, 22), ] <- NA

      pool <- sim_pool_with_nas(n = 200, na_frac = 0)

      # The point of this test: NA in x_cov at stripped focal positions
      # should NOT cause an error. (Ancillary warnings from elsewhere in the
      # query are unrelated to the NA-in-x_cov scenario and out of scope.)
      expect_no_error(
            res <- analog_search(
                  x = r, pool = pool,
                  select = "knn_env", k = 1, geog = kernel(max = 5),
                  x_cov = x_cov,
                  coord_type = "projected"
            )
      )
      expect_true(inherits(res, "SpatRaster"))
})

test_that("focal NAs and pool NAs both stripped correctly in single query", {
      set.seed(106)
      pool_full <- sim_pool_with_nas(n = 200, na_frac = 0.3)
      focal_full <- sim_focal_with_nas(n = 30, na_frac = 0.3)
      pool_clean <- pool_full[stats::complete.cases(pool_full), , drop = FALSE]
      focal_clean <- focal_full[stats::complete.cases(focal_full), , drop = FALSE]
      kept_focal <- which(stats::complete.cases(focal_full))

      args <- list(
            select = "all", stat = c("count", "sum_weights"),
            env = kernel("gaussian", theta = 1, max = 1.5),
            geog = kernel(max = 5),
            coord_type = "projected"
      )

      res_na    <- do.call(analog_search, c(args, list(x = focal_full, pool = pool_full)))
      res_clean <- do.call(analog_search, c(args, list(x = focal_clean, pool = pool_clean)))

      # Same shape: NA-containing focal output is reconstructed to original size
      expect_equal(nrow(res_na), nrow(focal_full))
      expect_equal(nrow(res_clean), nrow(focal_clean))

      # At kept focal positions, results should match
      expect_equal(res_na$count[kept_focal],       res_clean$count)
      expect_equal(res_na$sum_weights[kept_focal], res_clean$sum_weights,
                   tolerance = 1e-10)
      # At stripped focal positions, results should be NA
      expect_true(all(is.na(res_na$count[-kept_focal])))
})

test_that("fully-NA focal produces all-NA output without crashing", {
      set.seed(107)
      pool <- sim_pool_with_nas(n = 100, na_frac = 0)
      focal_all_na <- sim_focal_with_nas(n = 20, na_frac = 1.0)

      expect_silent(
            res <- analog_search(
                  x = focal_all_na, pool = pool,
                  select = "all", stat = "count",
                  env = kernel(max = 1.5), geog = kernel(max = 5),
                  coord_type = "projected"
            )
      )
      expect_equal(nrow(res), nrow(focal_all_na))
      expect_true(all(is.na(res$count)))
})
