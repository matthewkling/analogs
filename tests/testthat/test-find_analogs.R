test_that("find_analogs uses build+query architecture correctly", {

      d <- sim_test_data()

      # Test with explicit index_res (should build index then query)
      res1 <- find_analogs(
            x = d$focal,
            pool = d$ref,
            mode = "knn_geog",
            max_clim = 1,
            k = 1,
            coord_type = "projected",
            index_res = 10
      )

      expect_s3_class(res1, "data.frame")
      expect_equal(nrow(res1), nrow(d$focal))
      expect_true(all(c("focal_index", "analog_index", "geog_dist") %in% names(res1)))
})


test_that("find_analogs auto-tuning works with new architecture", {

      # Create larger dataset for auto-tuning to trigger
      set.seed(456)
      large_focal <- matrix(rnorm(2500 * 4), ncol = 4)
      ref_data <- matrix(rnorm(1000 * 4), ncol = 4)

      # Should use tune_index_res internally when index_res = "auto"
      res <- find_analogs(
            x = large_focal,
            pool = ref_data,
            mode = "count",
            max_clim = 1,
            coord_type = "projected",
            index_res = "auto"
      )

      expect_s3_class(res, "data.frame")
      expect_equal(nrow(res), nrow(large_focal))
      expect_true("value" %in% names(res))
})


test_that("find_analogs with numeric index_res builds and queries correctly", {

      d <- sim_test_data()

      # Different resolutions should all work
      for (res_val in c(8, 12, 16)) {
            result <- find_analogs(
                  x = d$focal,
                  pool = d$ref,
                  mode = "count",
                  max_clim = 1,
                  coord_type = "projected",
                  index_res = res_val
            )

            expect_s3_class(result, "data.frame")
            expect_equal(nrow(result), nrow(d$focal))
      }
})


test_that("find_analogs raw data path matches index path results", {

      d <- sim_test_data()

      # Build index explicitly
      index <- build_analog_index(d$ref, coord_type = "projected", index_res = 12)

      # Query with index
      res_index <- find_analogs(
            x = d$focal,
            pool = index,
            mode = "knn_geog",
            max_clim = 1,
            k = 1
      )

      # Query with raw data at same resolution
      res_raw <- find_analogs(
            x = d$focal,
            pool = d$ref,
            mode = "knn_geog",
            max_clim = 1,
            k = 1,
            coord_type = "projected",
            index_res = 12
      )

      # Should get similar results (may differ slightly in edge cases)
      expect_equal(nrow(res_index), nrow(res_raw))
      expect_equal(res_index$focal_index, res_raw$focal_index)

      # Analog indices should be highly correlated
      expect_true(cor(res_index$analog_index, res_raw$analog_index) > 0.95)
})


test_that("find_analogs preserves diagnostic attributes", {

      d <- sim_test_data()

      res <- find_analogs(
            x = d$focal,
            pool = d$ref,
            mode = "count",
            max_clim = 1,
            coord_type = "projected",
            index_res = 10
      )

      # Should have diagnostic attributes from C++
      expect_true(!is.null(attr(res, "n_ref")))
      expect_true(!is.null(attr(res, "n_clim")))
      expect_true(!is.null(attr(res, "total_bins")))
      expect_true(!is.null(attr(res, "binning_method")))
})


test_that("find_analogs works for all modes with new architecture", {

      d <- sim_test_data()

      # knn_geog
      v <- find_analogs(d$focal, d$ref, mode = "knn_geog", max_clim = 1, k = 1,
                        coord_type = "projected", index_res = 10)
      expect_equal(nrow(v), nrow(d$focal))

      # knn_clim
      i <- find_analogs(d$focal, d$ref, mode = "knn_clim", max_geog = 2, k = 3,
                        coord_type = "projected", index_res = 10)
      expect_true(nrow(i) <= nrow(d$focal) * 3)

      # count
      c <- find_analogs(d$focal, d$ref, mode = "count", max_clim = 1,
                        coord_type = "projected", index_res = 10)
      expect_equal(nrow(c), nrow(d$focal))

      # sum
      s <- find_analogs(d$focal, d$ref, mode = "sum", max_clim = 1,
                        weight = "uniform", coord_type = "projected", index_res = 10)
      expect_equal(nrow(s), nrow(d$focal))

      # all
      a <- find_analogs(d$focal, d$ref, mode = "all", max_clim = 1,
                        coord_type = "projected", index_res = 10)
      expect_true(nrow(a) >= nrow(d$focal))
})


test_that("find_analogs handles lon/lat with new architecture", {

      d <- sim_test_data(lonlat = TRUE)

      res <- find_analogs(
            x = d$focal,
            pool = d$ref,
            mode = "knn_geog",
            max_clim = 1,
            k = 1,
            coord_type = "lonlat",
            index_res = 10
      )

      expect_equal(nrow(res), nrow(d$focal))
      expect_true(all(is.finite(res$geog_dist)))
})


test_that("find_analogs dispatches correctly on analog_index", {

      d <- sim_test_data()

      # Build index
      index <- build_analog_index(d$ref, coord_type = "projected", index_res = 10)

      # Query with index should work
      res_index <- find_analogs(
            x = d$focal,
            pool = index,
            mode = "knn_geog",
            max_clim = 1,
            k = 1
      )

      # Query with raw data should work
      res_raw <- find_analogs(
            x = d$focal,
            pool = d$ref,
            mode = "knn_geog",
            max_clim = 1,
            k = 1,
            coord_type = "projected",
            index_res = 10
      )

      # Results should be similar (may differ slightly due to floating point)
      expect_equal(nrow(res_index), nrow(res_raw))
      expect_equal(res_index$focal_index, res_raw$focal_index)

      # Index-based should have similar analog indices (within reason)
      # Some tolerance since lattice queries can differ in edge cases
      expect_true(cor(res_index$analog_index, res_raw$analog_index) > 0.9)
})


test_that("find_analogs with index works for different modes", {

      d <- sim_test_data()
      index <- build_analog_index(d$ref, coord_type = "projected", index_res = 8)

      # knn_geog
      v <- find_analogs(d$focal, index, mode = "knn_geog", max_clim = 1, k = 1)
      expect_equal(nrow(v), nrow(d$focal))
      expect_true(all(c("focal_index", "analog_index", "geog_dist") %in% names(v)))

      # knn_clim
      i <- find_analogs(d$focal, index, mode = "knn_clim", max_geog = 2, k = 3)
      expect_true(nrow(i) <= nrow(d$focal) * 3)

      # count
      c <- find_analogs(d$focal, index, mode = "count", max_clim = 1, max_geog = 2)
      expect_equal(nrow(c), nrow(d$focal))
      expect_true("value" %in% names(c))
      expect_true(all(c$value >= 0))

      # sum
      s <- find_analogs(d$focal, index, mode = "sum", max_clim = 1, max_geog = 2,
                        weight = "uniform")
      expect_equal(nrow(s), nrow(d$focal))

      # all
      a <- find_analogs(d$focal, index, mode = "all", max_clim = 1, max_geog = 2)
      expect_true(nrow(a) >= nrow(d$focal))
})


test_that("find_analogs with index validates inputs", {

      d <- sim_test_data()
      index <- build_analog_index(d$ref, coord_type = "projected", index_res = 8)

      # Mismatched climate dimensions should error
      bad_focal <- matrix(rnorm(20 * 5), ncol = 5)  # 5 columns instead of 4
      expect_error(
            find_analogs(bad_focal, index, mode = "count"),
            "has 5 columns but index expects 4"
      )
})


test_that("find_analogs preserves backward compatibility with raw data", {

      d <- sim_test_data()

      # Original call pattern should still work
      res <- find_analogs(
            x = d$focal,
            pool = d$ref,
            mode = "knn_geog",
            max_clim = 1,
            k = 1,
            coord_type = "projected"
      )

      expect_s3_class(res, "data.frame")
      expect_true(nrow(res) > 0)
      expect_true(all(c("focal_index", "analog_index") %in% names(res)))
})


test_that("find_analogs with index works for lon/lat data", {

      d <- sim_test_data(lonlat = TRUE)

      # Build lon/lat index
      index <- build_analog_index(d$ref, coord_type = "lonlat", index_res = 8)
      expect_true(index$use_ecef)

      # Query should work
      res <- find_analogs(
            x = d$focal,
            pool = index,
            mode = "knn_geog",
            max_clim = 1,
            k = 1
      )

      expect_equal(nrow(res), nrow(d$focal))
      expect_true(all(is.finite(res$geog_dist)))
})


test_that("find_analogs index path handles all constraint combinations", {

      d <- sim_test_data()
      index <- build_analog_index(d$ref, coord_type = "projected", index_res = 8)

      # Only climate constraint
      r1 <- find_analogs(d$focal, index, mode = "count", max_clim = 1)
      expect_true(all(r1$value >= 0))

      # Only geographic constraint
      r2 <- find_analogs(d$focal, index, mode = "count", max_geog = 2)
      expect_true(all(r2$value >= 0))

      # Both constraints
      r3 <- find_analogs(d$focal, index, mode = "count", max_clim = 1, max_geog = 2)
      expect_true(all(r3$value >= 0))
      expect_true(all(r3$value <= r1$value))  # Combined should be subset
      expect_true(all(r3$value <= r2$value))

      # No constraints (should return all ref points as analogs)
      r4 <- find_analogs(d$focal, index, mode = "count")
      expect_true(all(r4$value > 0))
})
