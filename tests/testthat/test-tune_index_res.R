test_that("tune_index_res returns valid resolution", {

      d <- sim_test_data()

      # Basic tuning
      res <- tune_index_res(
            x = d$focal,
            pool = d$ref,
            select = "knn_geog",
            max_clim = 1,
            k = 1,
            coord_type = "projected",
            verbose = FALSE
      )

      expect_type(res, "integer")
      expect_length(res, 1)
      expect_true(res > 0)  # Should return a positive integer
})


test_that("tune_index_res works with different modes", {

      d <- sim_test_data()

      # knn_clim
      res1 <- tune_index_res(
            x = d$focal,
            pool = d$ref,
            select = "knn_clim",
            max_geog = 2,
            k = 3,
            coord_type = "projected",
            verbose = FALSE
      )
      expect_type(res1, "integer")
      expect_true(res1 > 0)

      # count
      res2 <- tune_index_res(
            x = d$focal,
            pool = d$ref,
            stat = "count",
            max_clim = 1,
            max_geog = 2,
            coord_type = "projected",
            verbose = FALSE
      )
      expect_type(res2, "integer")
      expect_true(res2 > 0)

      # sum
      res3 <- tune_index_res(
            x = d$focal,
            pool = d$ref,
            stat = "sum_weights",
            max_clim = 1,
            max_geog = 2,
            weight = "uniform",
            coord_type = "projected",
            verbose = FALSE
      )
      expect_type(res3, "integer")
      expect_true(res3 > 0)
})


test_that("tune_index_res uses subsampling for large datasets", {

      # Create larger dataset that triggers tuning
      set.seed(123)
      large_focal <- matrix(rnorm(3000 * 4), ncol = 4)
      ref_data <- matrix(rnorm(500 * 4), ncol = 4)

      # Should perform tuning (n > 2000)
      res <- tune_index_res(
            x = large_focal,
            pool = ref_data,
            stat = "count",
            max_clim = 1,
            coord_type = "projected",
            verbose = FALSE
      )

      expect_type(res, "integer")
      expect_true(res > 0)
})


test_that("tune_index_res returns default for small datasets", {

      # Small dataset (< 2000 rows)
      d <- sim_test_data()

      res <- tune_index_res(
            x = d$focal,
            pool = d$ref,
            select = "knn_geog",
            max_clim = 1,
            k = 1,
            coord_type = "projected",
            default_res = 20L,
            verbose = FALSE
      )

      # Should return default_res since dataset is small
      expect_equal(res, 20L)
})


test_that("tune_index_res respects custom default_res", {

      d <- sim_test_data()

      # Custom default
      res <- tune_index_res(
            x = d$focal,
            pool = d$ref,
            stat = "count",
            max_clim = 1,
            coord_type = "projected",
            default_res = 24L,
            verbose = FALSE
      )

      # Small dataset should return the custom default
      expect_equal(res, 24L)
})


test_that("tune_index_res works with auto coord detection", {

      # Projected data
      d_proj <- sim_test_data(lonlat = FALSE)
      res_proj <- tune_index_res(
            x = d_proj$focal,
            pool = d_proj$ref,
            stat = "count",
            max_clim = 1,
            coord_type = "auto",
            verbose = FALSE
      )
      expect_type(res_proj, "integer")
      expect_true(res_proj > 0)

      # Lon/lat data
      d_lonlat <- sim_test_data(lonlat = TRUE)
      res_lonlat <- tune_index_res(
            x = d_lonlat$focal,
            pool = d_lonlat$ref,
            stat = "count",
            max_clim = 1,
            coord_type = "auto",
            verbose = FALSE
      )
      expect_type(res_lonlat, "integer")
      expect_true(res_lonlat > 0)
})


test_that("tune_index_res adaptive bracketing works", {

      # Create dataset large enough to trigger tuning
      set.seed(456)
      large_focal <- matrix(rnorm(2500 * 4), ncol = 4)
      ref_data <- matrix(rnorm(500 * 4), ncol = 4)

      # Should use adaptive bracketing (3-5 evaluations)
      res <- tune_index_res(
            x = large_focal,
            pool = ref_data,
            select = "knn_geog",
            max_clim = 1,
            k = 1,
            coord_type = "projected",
            default_res = 16L,
            verbose = FALSE
      )

      # Result should be reasonable (between 4 and 32 given default of 16)
      expect_type(res, "integer")
      expect_true(res >= 4 && res <= 32)
})


test_that("tune_index_res verbose output works", {

      # Create dataset large enough to trigger tuning
      set.seed(789)
      large_focal <- matrix(rnorm(2500 * 4), ncol = 4)
      ref_data <- matrix(rnorm(500 * 4), ncol = 4)

      # Should produce output
      expect_message(
            tune_index_res(
                  x = large_focal,
                  pool = ref_data,
                  stat = "count",
                  max_clim = 1,
                  coord_type = "projected",
                  verbose = TRUE
            ),
            "Selected resolution"
      )

      # Should not produce messages when verbose = FALSE
      expect_silent(
            tune_index_res(
                  x = large_focal,
                  pool = ref_data,
                  stat = "count",
                  max_clim = 1,
                  coord_type = "projected",
                  verbose = FALSE
            )
      )
})


test_that("tune_index_res handles edge cases", {

      # Very small dataset - should return default immediately
      small_focal <- matrix(rnorm(5 * 3), ncol = 3)
      small_ref <- matrix(rnorm(10 * 3), ncol = 3)

      res <- tune_index_res(
            x = small_focal,
            pool = small_ref,
            stat = "count",
            max_clim = 1,
            coord_type = "projected",
            default_res = 12L,
            verbose = FALSE
      )

      expect_type(res, "integer")
      expect_equal(res, 12L)  # Should return default for small dataset
})


test_that("tune_index_res works with all weight options", {

      # Create dataset large enough to trigger tuning
      set.seed(999)
      large_focal <- matrix(rnorm(2500 * 4), ncol = 4)
      ref_data <- matrix(rnorm(500 * 4), ncol = 4)

      # Test with inverse_clim weight
      res1 <- tune_index_res(
            x = large_focal,
            pool = ref_data,
            stat = "sum_weights",
            max_clim = 1,
            max_geog = 2,
            weight = "inverse_clim",
            coord_type = "projected",
            verbose = FALSE
      )
      expect_type(res1, "integer")
      expect_true(res1 > 0)

      # Test with inverse_geog weight
      res2 <- tune_index_res(
            x = large_focal,
            pool = ref_data,
            stat = "sum_weights",
            max_clim = 1,
            max_geog = 2,
            weight = "inverse_geog",
            coord_type = "projected",
            verbose = FALSE
      )
      expect_type(res2, "integer")
      expect_true(res2 > 0)
})


test_that("tune_index_res works with lonlat coordinates", {

      # Create lon/lat dataset large enough to trigger tuning
      set.seed(111)
      large_focal <- matrix(nrow = 2500, ncol = 4)
      large_focal[, 1] <- runif(2500, -180, 180)  # lon
      large_focal[, 2] <- runif(2500, -90, 90)    # lat
      large_focal[, 3:4] <- rnorm(2500 * 2)

      ref_data <- matrix(nrow = 500, ncol = 4)
      ref_data[, 1] <- runif(500, -180, 180)
      ref_data[, 2] <- runif(500, -90, 90)
      ref_data[, 3:4] <- rnorm(500 * 2)

      res <- tune_index_res(
            x = large_focal,
            pool = ref_data,
            stat = "count",
            max_clim = 1,
            coord_type = "lonlat",
            verbose = FALSE
      )

      expect_type(res, "integer")
      expect_true(res > 0)
})
