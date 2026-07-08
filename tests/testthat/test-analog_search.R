test_that("analog_search uses build+query architecture correctly", {

      d <- sim_test_data()

      # Test with explicit per-family resolution (should build index then query)
      res1 <- analog_search(
            x = d$focal,
            pool = d$ref,
            select = "knn_geog",
            clim = kernel(max = 1),
            k = 1,
            coord_type = "projected",
            geog_res_adj = 2, clim_res_adj = 1.25
      )

      expect_s3_class(res1, "data.frame")
      expect_equal(nrow(res1), nrow(d$focal))
      expect_true(all(c("index", "analog_index", "geog_dist") %in% names(res1)))
})


test_that("analog_search auto-tuning works with new architecture", {

      # Create larger dataset for auto-tuning to trigger
      set.seed(456)
      large_focal <- matrix(rnorm(2500 * 4), ncol = 4)
      ref_data <- matrix(rnorm(1000 * 4), ncol = 4)

      # Should use tune_index_res internally when geog_res_adj = "auto"
      res <- analog_search(
            x = large_focal,
            pool = ref_data,
            stat = "count",
            clim = kernel(max = 1),
            coord_type = "projected",
            geog_res_adj = "auto"
      )

      expect_s3_class(res, "data.frame")
      expect_equal(nrow(res), nrow(large_focal))
      expect_true("count" %in% names(res))
})


test_that("analog_search with numeric geog_res_adj builds and queries correctly", {

      d <- sim_test_data()

      # Different resolutions should all work
      for (res_val in c(.5, 1, 2)) {
            result <- analog_search(
                  x = d$focal,
                  pool = d$ref,
                  stat = "count",
                  clim = kernel(max = 1),
                  coord_type = "projected",
                  geog_res_adj = res_val
            )

            expect_s3_class(result, "data.frame")
            expect_equal(nrow(result), nrow(d$focal))
      }
})


test_that("analog_search raw data path matches index path results", {

      d <- sim_test_data()

      # Build index explicitly
      index <- build_analog_index(d$ref, coord_type = "projected")

      # Query with index
      res_index <- analog_search(
            x = d$focal,
            pool = index,
            select = "knn_geog",
            clim = kernel(max = 1),
            k = 1
      )

      # Query with raw data at same resolution
      res_raw <- analog_search(
            x = d$focal,
            pool = d$ref,
            select = "knn_geog",
            clim = kernel(max = 1),
            k = 1,
            coord_type = "projected"
      )

      # Should get similar results (may differ slightly in edge cases)
      expect_equal(nrow(res_index), nrow(res_raw))
      expect_equal(res_index$index, res_raw$index)

      # Analog indices should be highly correlated
      expect_true(cor(res_index$analog_index, res_raw$analog_index) > 0.95)
})


test_that("analog_search preserves diagnostic attributes", {

      d <- sim_test_data()

      res <- analog_search(
            x = d$focal,
            pool = d$ref,
            stat = "count",
            clim = kernel(max = 1),
            coord_type = "projected"
      )

      # Should have diagnostic attributes from C++
      expect_true(!is.null(attr(res, "n_pool")))
      expect_true(!is.null(attr(res, "n_clim")))
      expect_true(!is.null(attr(res, "total_bins")))
      expect_true(!is.null(attr(res, "binning_method")))
})


test_that("analog_search works for all modes with new architecture", {

      d <- sim_test_data()

      # knn_geog
      v <- analog_search(d$focal, d$ref, select = "knn_geog", clim = kernel(max = 1), k = 1,
                         coord_type = "projected")
      expect_equal(nrow(v), nrow(d$focal))

      # knn_clim
      i <- analog_search(d$focal, d$ref, select = "knn_clim", geog = kernel(max = 2), k = 3,
                         coord_type = "projected")
      expect_true(nrow(i) <= nrow(d$focal) * 3)

      # count
      c <- analog_search(d$focal, d$ref, stat = "count", clim = kernel(max = 1),
                         coord_type = "projected")
      expect_equal(nrow(c), nrow(d$focal))

      # sum
      s <- analog_search(d$focal, d$ref, stat = "sum_weights", clim = kernel(max = 1),
                         coord_type = "projected")
      expect_equal(nrow(s), nrow(d$focal))

      # all
      a <- analog_search(d$focal, d$ref, select = "all", clim = kernel(max = 1),
                         coord_type = "projected")
      expect_true(nrow(a) >= nrow(d$focal))
})


test_that("analog_search handles lon/lat with new architecture", {

      d <- sim_test_data(lonlat = TRUE)

      res <- analog_search(
            x = d$focal,
            pool = d$ref,
            select = "knn_geog",
            clim = kernel(max = 1),
            k = 1,
            coord_type = "lonlat"
      )

      expect_equal(nrow(res), nrow(d$focal))
      expect_true(all(is.finite(res$geog_dist)))
})


test_that("analog_search dispatches correctly on analog_index", {

      d <- sim_test_data()

      # Build index
      index <- build_analog_index(d$ref, coord_type = "projected")

      # Query with index should work
      res_index <- analog_search(
            x = d$focal,
            pool = index,
            select = "knn_geog",
            clim = kernel(max = 1),
            k = 1
      )

      # Query with raw data should work
      res_raw <- analog_search(
            x = d$focal,
            pool = d$ref,
            select = "knn_geog",
            clim = kernel(max = 1),
            k = 1,
            coord_type = "projected"
      )

      # Results should be similar (may differ slightly due to floating point)
      expect_equal(nrow(res_index), nrow(res_raw))
      expect_equal(res_index$index, res_raw$index)

      # Index-based should have similar analog indices (within reason)
      # Some tolerance since lattice queries can differ in edge cases
      expect_true(cor(res_index$analog_index, res_raw$analog_index) > 0.9)
})


test_that("analog_search with index works for different modes", {

      d <- sim_test_data()
      index <- build_analog_index(d$ref, coord_type = "projected")

      # knn_geog
      v <- analog_search(d$focal, index, select = "knn_geog", clim = kernel(max = 1), k = 1)
      expect_equal(nrow(v), nrow(d$focal))
      expect_true(all(c("index", "analog_index", "geog_dist") %in% names(v)))

      # knn_clim
      i <- analog_search(d$focal, index, select = "knn_clim", geog = kernel(max = 2), k = 3)
      expect_true(nrow(i) <= nrow(d$focal) * 3)

      # count
      c <- analog_search(d$focal, index, stat = "count", clim = kernel(max = 1), geog = kernel(max = 2))
      expect_equal(nrow(c), nrow(d$focal))
      expect_true("count" %in% names(c))
      expect_true(all(c$value >= 0))

      # sum
      s <- analog_search(d$focal, index, stat = "sum_weights", clim = kernel(max = 1), geog = kernel(max = 2))
      expect_equal(nrow(s), nrow(d$focal))

      # all
      a <- analog_search(d$focal, index, select = "all", clim = kernel(max = 1), geog = kernel(max = 2))
      expect_true(nrow(a) >= nrow(d$focal))
})


test_that("analog_search with index validates inputs", {

      d <- sim_test_data()
      index <- build_analog_index(d$ref, coord_type = "projected")

      # Mismatched climate dimensions should error
      bad_focal <- matrix(rnorm(20 * 5), ncol = 5)  # 5 columns instead of 4
      expect_error(
            analog_search(bad_focal, index, stat = "count"),
            "has 5 columns but index expects 4"
      )
})


test_that("analog_search with index works for lon/lat data", {

      d <- sim_test_data(lonlat = TRUE)

      # Build lon/lat index
      index <- build_analog_index(d$ref, coord_type = "lonlat")
      expect_true(index$use_ecef)

      # Query should work
      res <- analog_search(
            x = d$focal,
            pool = index,
            select = "knn_geog",
            clim = kernel(max = 1),
            k = 1
      )

      expect_equal(nrow(res), nrow(d$focal))
      expect_true(all(is.finite(res$geog_dist)))
})


test_that("analog_search index path handles all constraint combinations", {

      d <- sim_test_data()
      index <- build_analog_index(d$ref, coord_type = "projected")

      # Only climate constraint
      r1 <- analog_search(d$focal, index, stat = "count", clim = kernel(max = 1))
      expect_true(all(r1$value >= 0))

      # Only geographic constraint
      r2 <- analog_search(d$focal, index, stat = "count", geog = kernel(max = 2))
      expect_true(all(r2$value >= 0))

      # Both constraints
      r3 <- analog_search(d$focal, index, stat = "count", clim = kernel(max = 1), geog = kernel(max = 2))
      expect_true(all(r3$value >= 0))
      expect_true(all(r3$value <= r1$value))  # Combined should be subset
      expect_true(all(r3$value <= r2$value))

      # No constraints (should return all ref points as analogs)
      r4 <- analog_search(d$focal, index, stat = "count")
      expect_true(all(r4$value > 0))
})


test_that("analog_search runs error-free with all valid select-stat-kernel combinations", {

      d <- sim_test_data()

      # Per-family kernel specs to sweep. Each entry is a list(clim=, geog=)
      # giving the kernel objects passed to analog_search. Covers uniform,
      # single-family gaussian/inverse on each side, and composed kernels
      # (a non-uniform shape on BOTH families) which the product model now
      # supports. All carry clim = kernel(max = 1) / geog = kernel(max = 2.)
      kernel_specs <- list(
            list(clim = kernel(max = 1),                     geog = kernel(max = 2)),                     # uniform
            list(clim = kernel("gaussian", theta = 1, max = 1), geog = kernel(max = 2)),                  # gaussian_clim
            list(clim = kernel(max = 1),                     geog = kernel("gaussian", theta = 1, max = 2)), # gaussian_geog
            list(clim = kernel("gaussian", theta = 1, max = 1), geog = kernel("gaussian", theta = 1, max = 2)), # gaussian both (was joint)
            list(clim = kernel("inverse", theta = 1, max = 1),  geog = kernel(max = 2)),                  # inverse_clim
            list(clim = kernel(max = 1),                     geog = kernel("inverse", theta = 1, max = 2)),  # inverse_geog
            list(clim = kernel("inverse", theta = 1, max = 1),  geog = kernel("gaussian", theta = 1, max = 2)), # composed
            list(clim = kernel("gaussian", theta = 1, max = 1), geog = kernel("inverse", theta = 1, max = 2))   # composed
      )

      for(s in c("all", "knn_clim", "knn_geog")){
            for(a in c("none", "count", "sum_weights", "mean_weights")){
                  for(ks in kernel_specs){

                        # Non-weighted stats take no kernel: use plain max kernels.
                        if(a %in% c("none", "count")){
                              clim_dom <- kernel(max = 1)
                              geog_dom <- kernel(max = 2)
                        } else {
                              clim_dom <- ks$clim
                              geog_dom <- ks$geog
                        }

                        expect_no_error(
                              analog_search(d$focal, d$ref,
                                            select = s, stat = a,
                                            clim = clim_dom, geog = geog_dom)
                        )
                  }
            }
      }

})


test_that("analog_search handles multiple stats correctly", {

      d <- sim_test_data(lonlat = TRUE)

      stats <- c("count", "sum_weights", "mean_weights")

      results <- analog_search(
            x = d$focal,
            pool = d$ref,
            select = "all",
            stat = stats,
            clim = kernel("gaussian", theta = 0.2, max = 0.5),
            geog = kernel(max = 100)
      )

      expect_equal(tail(colnames(results), 3), stats)
})

test_that("analog_search handles user-supplied values correctly", {

      d <- sim_test_data(lonlat = TRUE)

      y <- matrix(runif(nrow(d$ref)*2), ncol = 2)
      vars <- c("var1", "var2")
      colnames(y) <- vars

      stats <- c("count", "sum", "weighted_mean")

      results <- analog_search(
            x = d$focal,
            pool = d$ref,
            select = "all",
            stat = stats,
            y = y,
            clim = kernel("gaussian", theta = 0.2, max = NULL),
            geog = kernel(max = NULL)
      )

      # result should contain the correct columns
      expect_true(all(c("count", "sum_var1", "sum_var2", "weighted_mean_var1", "weighted_mean_var2")
                      %in% names(results)))

      # with both filters NULL, all focals should match global for sums
      expect_true(all(results$sum_var1 == sum(y[,"var1"])))

})


test_that("analog_search handles raster output correctly", {

      skip_if_not_installed("terra")

      requireNamespace("terra", quietly = TRUE)

      focal <- terra::rast(matrix(runif(100), 10))
      ref <- terra::rast(matrix(runif(100), 10))

      # aggregation mode
      result <- analog_search(
            x = focal,
            pool = ref,
            select = "all",
            stat = c("sum_weights", "count"),
            geog = kernel("inverse"),
            clim = kernel(max = .1)
      )


      expect_true(inherits(result, "SpatRaster"))
      expect_true(all(c("sum_weights", "count") %in% names(result)))

      # pairs mode: knn with k = 1
      result <- analog_search(
            x = focal,
            pool = ref,
            select = "knn_geog",
            k = 1,
            stat = "none",
            clim = kernel(max = .1)
      )
      expect_true(inherits(result, "SpatRaster"))
      expect_true(all(c("clim_dist", "geog_dist", "analog_index", "analog_x", "analog_y") %in%
                            names(result)))

      # pairs mode: knn with k > 1 (should NOT return raster)
      result <- analog_search(
            x = focal,
            pool = ref,
            select = "knn_geog",
            k = 3,
            stat = "none",
            clim = kernel(max = .1)
      )
      expect_true(inherits(result, "data.frame"))
      expect_true(all(c("clim_dist", "geog_dist", "analog_index", "analog_x", "analog_y") %in%
                            names(result)))

})



test_that("chunking/progress bars work correctly", {

      d <- sim_test_data(nref = 100, nfocal = 2000)

      # test both agg and pair modes
      for(stat in c("count", "none")){

            # should print a progress bar
            out <- capture.output(result <- analog_search(
                  x = d$focal,
                  pool = d$ref,
                  select = "knn_geog", k = 5,
                  stat = "count",
                  clim = kernel(max = .1),
                  progress = TRUE
            ))
            expect_true(grepl("100%", out))

            # should not affect results
            result2 <- analog_search(
                  x = d$focal,
                  pool = d$ref,
                  select = "knn_geog", k = 5,
                  stat = "count",
                  clim = kernel(max = .1),
                  progress = FALSE
            )
            expect_equal(result, result2)
      }

})
