# Tests for NA handling in pool data.
#
# Layers covered:
#
# 1. .format_data() with purpose = "pool" strips NA rows and attaches a
#    row_map; with purpose = "focal" preserves all rows.
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

# Helpers ------------------------------------------------------------------

# Build a small reproducible pool with specified NA pattern.
sim_pool_with_nas <- function(n = 200, seed = 42, na_frac = 0.3) {
      set.seed(seed)
      m <- matrix(rnorm(n * 4), ncol = 4)
      colnames(m) <- c("x", "y", "clim1", "clim2")
      n_na <- round(n * na_frac)
      if (n_na > 0L) {
            na_rows <- sample.int(n, n_na)
            # Mix of NA-coord and NA-climate to exercise both patterns.
            half <- ceiling(n_na / 2)
            m[na_rows[seq_len(half)], "clim1"] <- NA_real_
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

test_that(".format_data(purpose = 'focal') preserves NA rows", {
      pool <- sim_pool_with_nas(n = 100, na_frac = 0.4)

      out <- analogs:::.format_data(pool, purpose = "focal")

      expect_equal(nrow(out), 100)
      expect_null(attr(out, "row_map"))
      expect_equal(sum(stats::complete.cases(out)),
                   sum(stats::complete.cases(pool)))
})

test_that(".format_data(purpose = 'pool') strips NA rows and sets row_map", {
      pool <- sim_pool_with_nas(n = 100, na_frac = 0.4)
      n_complete <- sum(stats::complete.cases(pool))

      out <- analogs:::.format_data(pool, purpose = "pool")

      expect_equal(nrow(out), n_complete)
      expect_true(!is.null(attr(out, "row_map")))
      expect_equal(length(attr(out, "row_map")), n_complete)
      expect_true(all(stats::complete.cases(out)))

      # row_map should index back into original pool correctly
      row_map <- attr(out, "row_map")
      expect_equal(out[, 3], pool[row_map, "clim1"])
      expect_equal(out[, 4], pool[row_map, "clim2"])
})

test_that(".format_data(purpose = 'pool') with no NAs sets row_map = NULL", {
      pool <- sim_pool_with_nas(n = 100, na_frac = 0)

      out <- analogs:::.format_data(pool, purpose = "pool")

      expect_equal(nrow(out), 100)
      expect_null(attr(out, "row_map"))
})

test_that(".format_data(purpose = 'pool') strips fully-NA pool to zero rows", {
      pool <- sim_pool_with_nas(n = 50, na_frac = 1.0)

      out <- analogs:::.format_data(pool, purpose = "pool")

      expect_equal(nrow(out), 0)
      # row_map is integer(0), NOT NULL, when stripping happened (even if to zero)
      expect_true(!is.null(attr(out, "row_map")))
      expect_equal(length(attr(out, "row_map")), 0L)
})

test_that(".format_data 'pool' works with SpatRaster input", {
      skip_if_not_installed("terra")

      r <- terra::rast(nrows = 10, ncols = 10, nlyrs = 2,
                       xmin = 0, xmax = 10, ymin = 0, ymax = 10)
      set.seed(123)
      terra::values(r) <- matrix(rnorm(200), ncol = 2)
      vals <- terra::values(r)
      vals[c(5, 17, 33, 50, 71), 1] <- NA
      vals[c(8, 25, 60), 2] <- NA
      terra::values(r) <- vals

      out <- analogs:::.format_data(r, purpose = "pool")

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

      idx <- build_analog_index(pool, coord_type = "projected", index_res = 8)

      expect_equal(idx$n_pool, 200L)              # original size
      expect_equal(idx$n_pool_used, n_complete)   # post-strip
      expect_equal(nrow(idx$ref_data), n_complete)
      expect_true(!is.null(idx$pool_row_map))
      expect_equal(length(idx$pool_row_map), n_complete)
      expect_true(all(stats::complete.cases(idx$ref_data)))
})

test_that("build_analog_index without NAs has pool_row_map NULL and equal n_pool", {
      pool <- sim_pool_with_nas(n = 100, na_frac = 0)

      idx <- build_analog_index(pool, coord_type = "projected", index_res = 8)

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

      idx <- build_analog_index(r, coord_type = "lonlat", index_res = 6,
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
            build_analog_index(pool, coord_type = "projected", index_res = 8,
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
            select = "knn_clim",
            max_geog = 5,
            k = 3,
            coord_type = "projected",
            index_res = 8
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
            select = "knn_clim", max_geog = 10, k = 2,
            coord_type = "projected", index_res = 8
      )
      res_clean <- analog_search(
            x = focal, pool = stripped$pool,
            select = "knn_clim", max_geog = 10, k = 2,
            coord_type = "projected", index_res = 8
      )

      # Same number of pairs returned
      expect_equal(nrow(res_na), nrow(res_clean))

      # Same focal-by-focal climate/geog distances. (Tie-breaking on the
      # C++ side may pick different rows in the rare case of exact ties,
      # but distances will match either way for valid neighbors.)
      expect_equal(res_na$clim_dist, res_clean$clim_dist, tolerance = 1e-10)
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
            kernel = "gaussian_clim", theta = 1,
            max_clim = 1.5, max_geog = 5,
            coord_type = "projected", index_res = 8
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
            kernel = "gaussian_clim", theta = 1,
            max_clim = 1.5, max_geog = 5,
            y = y_full,
            coord_type = "projected", index_res = 8
      )

      # Equivalent: stripped pool with manually stripped y
      stripped <- strip_na_pool(pool_full)
      res_clean <- analog_search(
            x = focal, pool = stripped$pool,
            select = "all", stat = "weighted_mean",
            kernel = "gaussian_clim", theta = 1,
            max_clim = 1.5, max_geog = 5,
            y = y_full[stripped$row_map],
            coord_type = "projected", index_res = 8
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
            index_res = 8L,
            downsample = 1.0,
            seed = 55L
      )

      expect_true(all(is.finite(res$coord_mins)))
      expect_true(all(is.finite(res$coord_maxs)))
      expect_true(all(is.finite(res$clim_mins)))
      expect_true(all(is.finite(res$clim_maxs)))

      m_clean <- m[stats::complete.cases(m), , drop = FALSE]
      expect_equal(as.numeric(res$coord_mins), unname(apply(m_clean[, 1:2], 2, min)))
      expect_equal(as.numeric(res$coord_maxs), unname(apply(m_clean[, 1:2], 2, max)))
      expect_equal(as.numeric(res$clim_mins),  unname(apply(m_clean[, 3:4], 2, min)))
      expect_equal(as.numeric(res$clim_maxs),  unname(apply(m_clean[, 3:4], 2, max)))
})

test_that("C++ build with all-NA matrix returns empty index, no crash", {
      m <- sim_pool_with_nas(n = 50, na_frac = 1.0)

      res <- analogs:::build_analog_index_cpp(
            ref_mm = m,
            coord_type = "projected",
            index_res = 4L,
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
                  max_clim = 1.5, max_geog = 5,
                  coord_type = "projected", index_res = 4
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
                  select = "knn_clim", k = 1, max_geog = 5,
                  coord_type = "projected", index_res = 4
            )
      )
      expect_equal(nrow(res_pair), 10)
      expect_true(all(is.na(res_pair$analog_index)))
})
