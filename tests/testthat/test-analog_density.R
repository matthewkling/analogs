test_that("`analog_density` result matches manual calculation", {

      d <- sim_test_data()
      max_env <- 2
      max_geog <- 2

      # density
      nn <- analog_density(d$focal, d$ref, env = kernel("inverse", max = max_env),
                           geog = kernel(max = max_geog), coord_type = "projected")

      # manual density calculation. The inverse kernel is 1/(1 + d/theta) with
      # theta defaulting to 1, so each in-window analog contributes
      # 1/(1 + env_dist). Out-of-window analogs (env or geog beyond max) get
      # distance Inf, hence weight 1/(1+Inf) = 0.
      denv <- xdist(d$focal, d$ref, "env")
      dgeog <- xdist(d$focal, d$ref, "geog")
      denv[denv > max_env] <- Inf
      denv[dgeog > max_geog] <- Inf
      intens <- as.vector(rowSums(1 / (1 + denv)))

      expect_equal(nn$sum_weights, intens)
})

test_that("analog_density works with x/pool parameter names", {

      d <- sim_test_data()

      # Test with raw data
      i <- analog_density(
            x = d$focal,
            pool = d$ref,
            env = kernel("inverse", max = 1),
            geog = kernel(max = 2),
            coord_type = "projected"
      )

      expect_s3_class(i, "data.frame")
      expect_equal(nrow(i), nrow(d$focal))
      expect_true("sum_weights" %in% names(i))
      expect_true(all(is.numeric(i$sum_weights)))
})


test_that("analog_density works with analog_index", {

      d <- sim_test_data()

      # Build index
      index <- build_analog_index(d$ref, coord_type = "projected")

      # Query with index
      i <- analog_density(
            x = d$focal,
            pool = index,
            env = kernel("inverse", max = 1),
            geog = kernel(max = 2)
      )

      expect_s3_class(i, "data.frame")
      expect_equal(nrow(i), nrow(d$focal))
      expect_true("sum_weights" %in% names(i))
      expect_true(all(is.numeric(i$sum_weights)))
})

test_that("gaussian_env kernel works", {

      d <- sim_test_data()

      # Test with gaussian_env
      result <- analog_density(
            d$focal, d$ref,
            env = kernel("gaussian", theta = 0.5, max = 2),
            geog = kernel(max = 2),
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
      result <- analog_density(
            d$focal, d$ref,
            geog = kernel("gaussian", theta = 0.5, max = 2),
            env = kernel(max = 2),
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
      result <- analog_density(
            d$focal, d$ref,
            env = kernel("gaussian", theta = 0.5, max = 2),
            geog = kernel("gaussian", theta = 1.0, max = 2),  # c(theta_env, theta_geog)
            coord_type = "projected"
      )

      expect_s3_class(result, "data.frame")
      expect_equal(nrow(result), nrow(d$focal))
      expect_true(all(is.finite(result$sum_weights)))
      expect_true(all(result$sum_weights >= 0))
})


test_that("gaussian kernels respect theta parameter", {

      d <- sim_test_data()

      # Smaller theta should produce more localized (smaller) kernels
      r1 <- analog_density(
            d$focal, d$ref,
            env = kernel("gaussian", theta = 0.2, max = 2),
            geog = kernel(max = 2),  # small bandwidth
            coord_type = "projected"
      )

      r2 <- analog_density(
            d$focal, d$ref,
            env = kernel("gaussian", theta = 2.0, max = 2),
            geog = kernel(max = 2),  # large bandwidth
            coord_type = "projected"
      )

      # Larger bandwidth should generally give larger sums (less decay)
      # Though this isn't guaranteed for every focal point
      mean_ratio <- mean(r2$sum_weights / r1$sum_weights, na.rm = TRUE)
      expect_true(mean_ratio > 1.0)
})


test_that("new kernels work with analog_index", {

      d <- sim_test_data()

      # Build index
      index <- build_analog_index(d$ref, coord_type = "projected")

      # Test each new kernel with index
      r1 <- analog_density(
            d$focal, index,
            env = kernel("gaussian", theta = 0.5, max = 2),
            geog = kernel(max = 2)
      )
      expect_true(all(is.finite(r1$sum_weights)))

      r2 <- analog_density(
            d$focal, index,
            geog = kernel("gaussian", theta = 0.5, max = 2),
            env = kernel(max = 2)
      )
      expect_true(all(is.finite(r2$sum_weights)))

      r3 <- analog_density(
            d$focal, index,
            env = kernel("gaussian", theta = 0.5, max = 2),
            geog = kernel("gaussian", theta = 1.0, max = 2)
      )
      expect_true(all(is.finite(r3$sum_weights)))
})


test_that("gaussian and inverse kernels produce different results", {

      d <- sim_test_data()

      # Compare gaussian vs inverse with similar parameterization
      r_gauss <- analog_density(
            d$focal, d$ref,
            env = kernel("gaussian", theta = 0.5, max = 2),
            geog = kernel(max = 2),
            coord_type = "projected"
      )

      r_inv <- analog_density(
            d$focal, d$ref,
            env = kernel("inverse", theta = 0.5, max = 2),
            geog = kernel(max = 2),
            coord_type = "projected"
      )

      # Results should differ (different functional forms)
      expect_false(all(abs(r_gauss$sum_weights - r_inv$sum_weights) < 1e-10))
})


test_that("joint kernels combine both dimensions appropriately", {

      d <- sim_test_data()

      # Joint should differ from single-dimension kernels
      r_joint <- analog_density(
            d$focal, d$ref,
            env = kernel("gaussian", theta = 0.5, max = 2),
            geog = kernel("gaussian", theta = 0.5, max = 2),
            coord_type = "projected"
      )

      r_env_only <- analog_density(
            d$focal, d$ref,
            env = kernel("gaussian", theta = 0.5, max = 2),
            geog = kernel(max = 2),
            coord_type = "projected"
      )

      # Results should generally differ since joint considers both dimensions
      expect_false(all(abs(r_joint$sum_weights - r_env_only$sum_weights) < 1e-10))
})
