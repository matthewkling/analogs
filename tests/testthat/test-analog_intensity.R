test_that("`analog_intensity` result matches manual calculation", {

      d <- sim_test_data()
      max_clim <- 2
      max_geog <- 2

      # intensity
      nn <- analog_intensity(d$focal, d$ref, max_clim = max_clim, max_geog = max_geog,
                             kernel = "inverse_clim", coord_type = "projected")

      # manual intensity calculation
      dclim <- xdist(d$focal, d$ref, "clim")
      dgeog <- xdist(d$focal, d$ref, "geog")
      dclim[dclim > max_clim] <- Inf
      dclim[dgeog > max_geog] <- Inf
      intens <- as.vector(rowSums(1 / dclim))

      expect_equal(nn$sum_weights, intens)
})

test_that("analog_intensity works with x/pool parameter names", {

      d <- sim_test_data()

      # Test with raw data
      i <- analog_intensity(
            x = d$focal,
            pool = d$ref,
            max_clim = 1,
            max_geog = 2,
            kernel = "inverse_clim",
            coord_type = "projected",
            index_res = 10
      )

      expect_s3_class(i, "data.frame")
      expect_equal(nrow(i), nrow(d$focal))
      expect_true("sum_weights" %in% names(i))
      expect_true(all(is.numeric(i$sum_weights)))
})


test_that("analog_intensity works with analog_index", {

      d <- sim_test_data()

      # Build index
      index <- build_analog_index(d$ref, coord_type = "projected", index_res = 12)

      # Query with index
      i <- analog_intensity(
            x = d$focal,
            pool = index,
            max_clim = 1,
            max_geog = 2,
            kernel = "inverse_clim"
      )

      expect_s3_class(i, "data.frame")
      expect_equal(nrow(i), nrow(d$focal))
      expect_true("sum_weights" %in% names(i))
      expect_true(all(is.numeric(i$sum_weights)))
})

test_that("gaussian_clim kernel works", {

      d <- sim_test_data()

      # Test with gaussian_clim
      result <- analog_intensity(
            d$focal, d$ref,
            max_clim = 2,
            max_geog = 2,
            kernel = "gaussian_clim",
            theta = 0.5,
            coord_type = "projected"
      )

      expect_s3_class(result, "data.frame")
      expect_equal(nrow(result), nrow(d$focal))
      expect_true(all(is.finite(result$sum_weights)))
      expect_true(all(result$sum_weights >= 0))  # Gaussian weights are always positive
})


test_that("gaussian_geog kernel works", {

      d <- sim_test_data()

      # Test with gaussian_geog
      result <- analog_intensity(
            d$focal, d$ref,
            max_clim = 2,
            max_geog = 2,
            kernel = "gaussian_geog",
            theta = 0.5,
            coord_type = "projected"
      )

      expect_s3_class(result, "data.frame")
      expect_equal(nrow(result), nrow(d$focal))
      expect_true(all(is.finite(result$sum_weights)))
      expect_true(all(result$sum_weights >= 0))
})


test_that("gaussian_joint kernel works", {

      d <- sim_test_data()

      # Test with gaussian_joint (requires theta length 2)
      result <- analog_intensity(
            d$focal, d$ref,
            max_clim = 2,
            max_geog = 2,
            kernel = "gaussian_joint",
            theta = c(0.5, 1.0),  # c(theta_clim, theta_geog)
            coord_type = "projected"
      )

      expect_s3_class(result, "data.frame")
      expect_equal(nrow(result), nrow(d$focal))
      expect_true(all(is.finite(result$sum_weights)))
      expect_true(all(result$sum_weights >= 0))
})


test_that("inverse_joint kernel works", {

      d <- sim_test_data()

      # Test with inverse_joint
      result <- analog_intensity(
            d$focal, d$ref,
            max_clim = 2,
            max_geog = 2,
            kernel = "inverse_joint",
            theta = c(1e-6, 1e-3),  # c(eps_clim, eps_geog)
            coord_type = "projected"
      )

      expect_s3_class(result, "data.frame")
      expect_equal(nrow(result), nrow(d$focal))
      expect_true(all(is.finite(result$sum_weights)))
      expect_true(all(result$sum_weights > 0))
})


test_that("gaussian kernels respect theta parameter", {

      d <- sim_test_data()

      # Smaller theta should produce more localized (smaller) kernels
      r1 <- analog_intensity(
            d$focal, d$ref,
            max_clim = 2,
            max_geog = 2,
            kernel = "gaussian_clim",
            theta = 0.2,  # small bandwidth
            coord_type = "projected"
      )

      r2 <- analog_intensity(
            d$focal, d$ref,
            max_clim = 2,
            max_geog = 2,
            kernel = "gaussian_clim",
            theta = 2.0,  # large bandwidth
            coord_type = "projected"
      )

      # Larger bandwidth should generally give larger sums (less decay)
      # Though this isn't guaranteed for every focal point
      mean_ratio <- mean(r2$sum_weights / r1$sum_weights, na.rm = TRUE)
      expect_true(mean_ratio > 1.0)
})


test_that("gaussian_joint requires 2-element theta", {

      d <- sim_test_data()

      # Should error with scalar theta
      expect_error(
            analog_intensity(
                  d$focal, d$ref,
                  max_clim = 2,
                  kernel = "gaussian_joint",
                  theta = 0.5,
                  coord_type = "projected"
            ),
            "length 2"
      )

      # Should error with NULL theta
      expect_error(
            analog_intensity(
                  d$focal, d$ref,
                  max_clim = 2,
                  kernel = "gaussian_joint",
                  theta = NULL,
                  coord_type = "projected"
            ),
            "length 2"
      )
})


test_that("inverse_joint requires 2-element theta", {

      d <- sim_test_data()

      # Should error with scalar theta
      expect_error(
            analog_intensity(
                  d$focal, d$ref,
                  max_clim = 2,
                  kernel = "inverse_joint",
                  theta = 1e-6,
                  coord_type = "projected"
            ),
            "length 2"
      )
})


test_that("new kernels work with analog_index", {

      d <- sim_test_data()

      # Build index
      index <- build_analog_index(d$ref, coord_type = "projected", index_res = 10)

      # Test each new kernel with index
      r1 <- analog_intensity(
            d$focal, index,
            max_clim = 2,
            max_geog = 2,
            kernel = "gaussian_clim",
            theta = 0.5
      )
      expect_true(all(is.finite(r1$sum_weights)))

      r2 <- analog_intensity(
            d$focal, index,
            max_clim = 2,
            max_geog = 2,
            kernel = "gaussian_geog",
            theta = 0.5
      )
      expect_true(all(is.finite(r2$sum_weights)))

      r3 <- analog_intensity(
            d$focal, index,
            max_clim = 2,
            max_geog = 2,
            kernel = "gaussian_joint",
            theta = c(0.5, 1.0)
      )
      expect_true(all(is.finite(r3$sum_weights)))

      r4 <- analog_intensity(
            d$focal, index,
            max_clim = 2,
            max_geog = 2,
            kernel = "inverse_joint",
            theta = c(1e-6, 1e-3)
      )
      expect_true(all(is.finite(r4$sum_weights)))
})


test_that("gaussian and inverse kernels produce different results", {

      d <- sim_test_data()

      # Compare gaussian vs inverse with similar parameterization
      r_gauss <- analog_intensity(
            d$focal, d$ref,
            max_clim = 2,
            max_geog = 2,
            kernel = "gaussian_clim",
            theta = 0.5,
            coord_type = "projected"
      )

      r_inv <- analog_intensity(
            d$focal, d$ref,
            max_clim = 2,
            max_geog = 2,
            kernel = "inverse_clim",
            theta = 1e-6,
            coord_type = "projected"
      )

      # Results should differ (different functional forms)
      expect_false(all(abs(r_gauss$sum_weights - r_inv$sum_weights) < 1e-10))
})


test_that("joint kernels combine both dimensions appropriately", {

      d <- sim_test_data()

      # Joint should differ from single-dimension kernels
      r_joint <- analog_intensity(
            d$focal, d$ref,
            max_clim = 2,
            max_geog = 2,
            kernel = "gaussian_joint",
            theta = c(0.5, 0.5),
            coord_type = "projected"
      )

      r_clim_only <- analog_intensity(
            d$focal, d$ref,
            max_clim = 2,
            max_geog = 2,
            kernel = "gaussian_clim",
            theta = 0.5,
            coord_type = "projected"
      )

      # Results should generally differ since joint considers both dimensions
      expect_false(all(abs(r_joint$sum_weights - r_clim_only$sum_weights) < 1e-10))
})
