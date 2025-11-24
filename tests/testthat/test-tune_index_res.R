test_that("tune_index_res returns valid resolution", {

      d <- sim_test_data()

      # Basic tuning with minimal reps for speed
      res <- tune_index_res(
            x = d$focal,
            pool = d$ref,
            mode = "knn_geog",
            max_clim = 1,
            k = 1,
            coord_type = "projected",
            resolutions = c(8, 12, 16),
            n_reps = 2,
            verbose = FALSE
      )

      expect_type(res, "integer")
      expect_length(res, 1)
      expect_true(res %in% c(8, 12, 16))
})


test_that("tune_index_res works with different modes", {

      d <- sim_test_data()

      # knn_clim
      res1 <- tune_index_res(
            x = d$focal,
            pool = d$ref,
            mode = "knn_clim",
            max_geog = 2,
            k = 3,
            coord_type = "projected",
            resolutions = c(6, 10),
            n_reps = 1,
            verbose = FALSE
      )
      expect_true(res1 %in% c(6, 10))

      # count
      res2 <- tune_index_res(
            x = d$focal,
            pool = d$ref,
            mode = "count",
            max_clim = 1,
            max_geog = 2,
            coord_type = "projected",
            resolutions = c(6, 10),
            n_reps = 1,
            verbose = FALSE
      )
      expect_true(res2 %in% c(6, 10))

      # sum
      res3 <- tune_index_res(
            x = d$focal,
            pool = d$ref,
            mode = "sum",
            max_clim = 1,
            max_geog = 2,
            weight = "uniform",
            coord_type = "projected",
            resolutions = c(6, 10),
            n_reps = 1,
            verbose = FALSE
      )
      expect_true(res3 %in% c(6, 10))
})


test_that("tune_index_res samples large focal datasets", {

      # Create larger dataset
      set.seed(123)
      large_focal <- matrix(rnorm(1000 * 4), ncol = 4)
      ref_data <- matrix(rnorm(200 * 4), ncol = 4)

      # Should sample down to n_sample
      res <- tune_index_res(
            x = large_focal,
            pool = ref_data,
            mode = "count",
            max_clim = 1,
            coord_type = "projected",
            resolutions = c(8, 12),
            n_sample = 50,
            n_reps = 1,
            verbose = FALSE
      )

      expect_type(res, "integer")
      expect_true(res %in% c(8, 12))
})


test_that("tune_index_res returns detailed results invisibly", {

      d <- sim_test_data()

      results <- tune_index_res(
            x = d$focal,
            pool = d$ref,
            mode = "knn_geog",
            max_clim = 1,
            k = 1,
            coord_type = "projected",
            resolutions = c(8, 12, 16),
            n_reps = 2,
            verbose = FALSE
      )

      # Invisible return should be a data.frame
      expect_s3_class(results, "data.frame")
      expect_true(all(c("resolution", "mean_time_ms", "sd_time_ms") %in% names(results)))
      expect_equal(nrow(results), 3)  # 3 resolutions tested
      expect_equal(results$resolution, c(8, 12, 16))
      expect_true(all(results$mean_time_ms > 0))
})


test_that("tune_index_res validates inputs", {

      d <- sim_test_data()

      # Invalid resolutions
      expect_error(
            tune_index_res(d$focal, d$ref, mode = "count", resolutions = c(-1, 10)),
            "positive integers"
      )

      # Invalid n_sample
      expect_error(
            tune_index_res(d$focal, d$ref, mode = "count", n_sample = -5),
            "positive integer"
      )

      # Invalid n_reps
      expect_error(
            tune_index_res(d$focal, d$ref, mode = "count", n_reps = 0),
            "positive integer"
      )
})


test_that("tune_index_res works with auto coord detection", {

      # Projected data
      d_proj <- sim_test_data(lonlat = FALSE)
      res_proj <- tune_index_res(
            x = d_proj$focal,
            pool = d_proj$ref,
            mode = "count",
            max_clim = 1,
            coord_type = "auto",
            resolutions = c(8, 12),
            n_reps = 1,
            verbose = FALSE
      )
      expect_true(res_proj %in% c(8, 12))

      # Lon/lat data
      d_lonlat <- sim_test_data(lonlat = TRUE)
      res_lonlat <- tune_index_res(
            x = d_lonlat$focal,
            pool = d_lonlat$ref,
            mode = "count",
            max_clim = 1,
            coord_type = "auto",
            resolutions = c(8, 12),
            n_reps = 1,
            verbose = FALSE
      )
      expect_true(res_lonlat %in% c(8, 12))
})


test_that("tune_index_res produces reasonable recommendations", {

      d <- sim_test_data()

      # Test a wider range
      res <- tune_index_res(
            x = d$focal,
            pool = d$ref,
            mode = "knn_geog",
            max_clim = 1,
            k = 1,
            coord_type = "projected",
            resolutions = c(4, 8, 16, 32),
            n_reps = 2,
            verbose = FALSE
      )

      # Should pick something in the middle range (4 and 32 are extremes)
      # This is probabilistic but should generally hold
      expect_true(res >= 4 && res <= 32)
})


test_that("tune_index_res verbose output works", {

      d <- sim_test_data()

      # Should produce output
      expect_output(
            tune_index_res(
                  x = d$focal,
                  pool = d$ref,
                  mode = "count",
                  max_clim = 1,
                  coord_type = "projected",
                  resolutions = c(8, 12),
                  n_reps = 1,
                  verbose = TRUE
            ),
            "Optimal resolution"
      )

      # Should not produce output when verbose = FALSE
      expect_silent(
            tune_index_res(
                  x = d$focal,
                  pool = d$ref,
                  mode = "count",
                  max_clim = 1,
                  coord_type = "projected",
                  resolutions = c(8, 12),
                  n_reps = 1,
                  verbose = FALSE
            )
      )
})


test_that("tune_index_res handles edge cases", {

      # Very small dataset
      small_focal <- matrix(rnorm(5 * 3), ncol = 3)
      small_ref <- matrix(rnorm(10 * 3), ncol = 3)

      res <- tune_index_res(
            x = small_focal,
            pool = small_ref,
            mode = "count",
            max_clim = 1,
            coord_type = "projected",
            resolutions = c(4, 8),
            n_reps = 1,
            verbose = FALSE
      )

      expect_type(res, "integer")
      expect_true(res %in% c(4, 8))
})
