test_that("regression stat runs without error and returns correct columns", {

      d <- sim_test_data()
      n_ref <- nrow(d$ref)

      y <- runif(n_ref)
      covariates <- matrix(rnorm(n_ref * 2), ncol = 2,
                           dimnames = list(NULL, c("northness", "eastness")))

      result <- analog_search(
            x = d$focal,
            pool = d$ref,
            select = "all",
            stat = "regression",
            y = y,
            covariates = covariates,
            max_clim = 2,
            max_geog = 2,
            weight = "uniform",
            coord_type = "projected",
            index_res = 10
      )

      expect_s3_class(result, "data.frame")
      expect_equal(nrow(result), nrow(d$focal))
      expect_true(all(c("intercept", "northness", "eastness") %in% names(result)))
})


test_that("regression recovers known coefficients on synthetic data", {

      set.seed(42)

      # Create a pool where y = 2 + 3*x1 - 1*x2 + noise
      n <- 500
      ref <- matrix(0, nrow = n, ncol = 4)
      ref[, 1] <- runif(n, 0, 10)   # x coord
      ref[, 2] <- runif(n, 0, 10)   # y coord
      ref[, 3] <- rnorm(n)          # clim1
      ref[, 4] <- rnorm(n)          # clim2

      x1 <- rnorm(n)
      x2 <- rnorm(n)
      covariates <- cbind(x1 = x1, x2 = x2)

      y <- 2 + 3 * x1 - 1 * x2 + rnorm(n, sd = 0.01)

      # Focal point at the center with matching climate
      focal <- matrix(c(5, 5, 0, 0), nrow = 1)

      result <- analog_search(
            x = focal,
            pool = ref,
            select = "all",
            stat = "regression",
            y = y,
            covariates = covariates,
            max_clim = NULL,  # no climate filter: use all pool points
            max_geog = NULL,  # no geog filter
            weight = "uniform",
            coord_type = "projected",
            index_res = 10
      )

      # With uniform weights and all points used, should closely recover true coefficients
      expect_equal(result$intercept, 2, tolerance = 0.1)
      expect_equal(result$x1, 3, tolerance = 0.1)
      expect_equal(result$x2, -1, tolerance = 0.1)
})


test_that("large lambda makes intercept approach weighted mean", {

      d <- sim_test_data()
      n_ref <- nrow(d$ref)

      y <- rnorm(n_ref)
      covariates <- matrix(rnorm(n_ref * 2), ncol = 2,
                           dimnames = list(NULL, c("c1", "c2")))

      # Run with very large lambda
      result_ridge <- analog_search(
            x = d$focal,
            pool = d$ref,
            select = "all",
            stat = c("weighted_mean", "regression"),
            y = y,
            covariates = covariates,
            max_clim = 2,
            max_geog = 2,
            weight = "gaussian_clim",
            theta = 0.5,
            lambda = 1e8,
            coord_type = "projected",
            index_res = 10
      )

      # Intercept should be very close to weighted_mean
      # (Not exact because covariates aren't centered, but with huge lambda
      # the covariate coefficients are ~0 so intercept ≈ weighted mean)
      has_analogs <- !is.na(result_ridge$weighted_mean) & !is.na(result_ridge$intercept)
      if (any(has_analogs)) {
            expect_equal(
                  result_ridge$intercept[has_analogs],
                  result_ridge$weighted_mean[has_analogs],
                  tolerance = 0.01
            )

            # Covariate coefficients should be near zero
            expect_true(all(abs(result_ridge$c1[has_analogs]) < 0.01))
            expect_true(all(abs(result_ridge$c2[has_analogs]) < 0.01))
      }
})


test_that("zero analogs returns NA for regression coefficients", {

      d <- sim_test_data()
      n_ref <- nrow(d$ref)

      y <- rnorm(n_ref)
      covariates <- matrix(rnorm(n_ref), ncol = 1,
                           dimnames = list(NULL, "cov1"))

      # Use impossibly tight constraints so no analogs are found
      result <- analog_search(
            x = d$focal,
            pool = d$ref,
            select = "all",
            stat = c("count", "regression"),
            y = y,
            covariates = covariates,
            max_clim = 0.0001,
            max_geog = 0.0001,
            weight = "uniform",
            coord_type = "projected",
            index_res = 10
      )

      # Rows with count == 0 should have NA coefficients
      zero_rows <- result$count == 0
      if (any(zero_rows)) {
            expect_true(all(is.na(result$intercept[zero_rows])))
            expect_true(all(is.na(result$cov1[zero_rows])))
      }
})


test_that("regression works with multiple values variables", {

      d <- sim_test_data()
      n_ref <- nrow(d$ref)

      y <- matrix(rnorm(n_ref * 2), ncol = 2,
                       dimnames = list(NULL, c("biomass", "richness")))
      covariates <- matrix(rnorm(n_ref * 2), ncol = 2,
                           dimnames = list(NULL, c("slope", "aspect")))

      result <- analog_search(
            x = d$focal,
            pool = d$ref,
            select = "all",
            stat = "regression",
            y = y,
            covariates = covariates,
            max_clim = 2,
            max_geog = 2,
            weight = "uniform",
            coord_type = "projected",
            index_res = 10
      )

      # Should have 3 columns per y variable: intercept, slope, aspect
      expected_cols <- c("intercept_biomass", "slope_biomass", "aspect_biomass",
                         "intercept_richness", "slope_richness", "aspect_richness")
      expect_true(all(expected_cols %in% names(result)))
})


test_that("regression combines with other stats correctly", {

      d <- sim_test_data()
      n_ref <- nrow(d$ref)

      y <- rnorm(n_ref)
      covariates <- matrix(rnorm(n_ref), ncol = 1,
                           dimnames = list(NULL, "topo"))

      result <- analog_search(
            x = d$focal,
            pool = d$ref,
            select = "all",
            stat = c("count", "ess", "weighted_mean", "regression"),
            y = y,
            covariates = covariates,
            max_clim = 2,
            max_geog = 2,
            weight = "gaussian_clim",
            theta = 0.5,
            coord_type = "projected",
            index_res = 10
      )

      expect_true(all(c("count", "ess", "weighted_mean",
                        "intercept", "topo") %in% names(result)))
      expect_equal(nrow(result), nrow(d$focal))
})


test_that("regression validation: requires y", {

      d <- sim_test_data()
      covariates <- matrix(rnorm(nrow(d$ref)), ncol = 1)

      expect_error(
            analog_search(
                  x = d$focal, pool = d$ref,
                  stat = "regression",
                  covariates = covariates,
                  max_clim = 2,
                  weight = "uniform",
                  coord_type = "projected"
            ),
            "y"
      )
})


test_that("regression validation: requires covariates", {

      d <- sim_test_data()
      y <- rnorm(nrow(d$ref))

      expect_error(
            analog_search(
                  x = d$focal, pool = d$ref,
                  stat = "regression",
                  y = y,
                  max_clim = 2,
                  weight = "uniform",
                  coord_type = "projected"
            ),
            "covariates"
      )
})


test_that("regression validation: requires weight", {

      d <- sim_test_data()
      n_ref <- nrow(d$ref)

      expect_error(
            analog_search(
                  x = d$focal, pool = d$ref,
                  stat = "regression",
                  y = rnorm(n_ref),
                  covariates = matrix(rnorm(n_ref), ncol = 1),
                  max_clim = 2,
                  coord_type = "projected"
            ),
            "weight"
      )
})


test_that("regression validation: lambda must be non-negative", {

      d <- sim_test_data()
      n_ref <- nrow(d$ref)

      expect_error(
            analog_search(
                  x = d$focal, pool = d$ref,
                  stat = "regression",
                  y = rnorm(n_ref),
                  covariates = matrix(rnorm(n_ref), ncol = 1),
                  max_clim = 2,
                  weight = "uniform",
                  lambda = -1,
                  coord_type = "projected"
            ),
            "lambda"
      )
})


test_that("analog_impact passes covariates and lambda through", {

      d <- sim_test_data()
      n_ref <- nrow(d$ref)

      y <- rnorm(n_ref)
      covariates <- matrix(rnorm(n_ref * 2), ncol = 2,
                           dimnames = list(NULL, c("north", "east")))

      result <- analog_impact(
            x = d$focal,
            pool = d$ref,
            y = y,
            covariates = covariates,
            stat = c("count", "weighted_mean", "regression"),
            max_geog = 2,
            max_clim = 2,
            weight = "gaussian_clim",
            theta = 0.5,
            lambda = 0.1,
            coord_type = "projected",
            index_res = 10
      )

      expect_true(all(c("count", "weighted_mean",
                        "intercept", "north", "east") %in% names(result)))
})


test_that("regression with lambda = 0 returns NA for singular systems", {

      set.seed(99)

      # Create scenario where some focals will have very few analogs
      # (fewer than the number of covariates)
      n_ref <- 50
      ref <- matrix(rnorm(n_ref * 4), ncol = 4)
      focal <- matrix(c(100, 100, 0, 0), nrow = 1)  # far from all pool points

      y <- rnorm(n_ref)
      # 3 covariates — need at least 4 analogs for full rank
      covariates <- matrix(rnorm(n_ref * 3), ncol = 3,
                           dimnames = list(NULL, c("a", "b", "c")))

      result <- analog_search(
            x = focal,
            pool = ref,
            select = "all",
            stat = c("count", "regression"),
            y = y,
            covariates = covariates,
            max_clim = 0.01,  # very tight — may get 0 or very few analogs
            max_geog = 0.01,
            weight = "uniform",
            lambda = 0,
            coord_type = "projected",
            index_res = 10
      )

      # If count is 0, all coefficients should be NA
      if (result$count == 0) {
            expect_true(is.na(result$intercept))
            expect_true(is.na(result$a))
            expect_true(is.na(result$b))
            expect_true(is.na(result$c))
      }
})
