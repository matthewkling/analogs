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
            env = kernel(max = 2),
            geog = kernel(max = 2),
            coord_type = "projected"
      )

      expect_s3_class(result, "data.frame")
      expect_equal(nrow(result), nrow(d$focal))
      expect_true(all(c("coef_intercept", "coef_northness", "coef_eastness") %in% names(result)))
})


test_that("regression recovers known coefficients on synthetic data", {

      set.seed(42)

      # Create a pool where y = 2 + 3*x1 - 1*x2 + noise
      n <- 500
      ref <- matrix(0, nrow = n, ncol = 4)
      ref[, 1] <- runif(n, 0, 10)   # x coord
      ref[, 2] <- runif(n, 0, 10)   # y coord
      ref[, 3] <- rnorm(n)          # env1
      ref[, 4] <- rnorm(n)          # env2

      x1 <- rnorm(n)
      x2 <- rnorm(n)
      covariates <- cbind(x1 = x1, x2 = x2)

      y <- 2 + 3 * x1 - 1 * x2 + rnorm(n, sd = 0.01)

      # Focal point at the center with matching environmental
      focal <- matrix(c(5, 5, 0, 0), nrow = 1)

      result <- analog_search(
            x = focal,
            pool = ref,
            select = "all",
            stat = "regression",
            y = y,
            covariates = covariates,
            env = NULL,  # no environmental filter: use all pool points
            geog = NULL,  # no geog filter
            coord_type = "projected"
      )

      # With uniform weights and all points used, should closely recover true coefficients
      expect_equal(result$coef_intercept, 2, tolerance = 0.1)
      expect_equal(result$coef_x1, 3, tolerance = 0.1)
      expect_equal(result$coef_x2, -1, tolerance = 0.1)
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
            env = kernel("gaussian", theta = 0.5, max = 2),
            geog = kernel(max = 2),
            lambda = 1e8,
            coord_type = "projected"
      )

      # Intercept should be very close to weighted_mean
      # (Not exact because covariates aren't centered, but with huge lambda
      # the covariate coefficients are ~0 so intercept ≈ weighted mean)
      has_analogs <- !is.na(result_ridge$weighted_mean) & !is.na(result_ridge$coef_intercept)
      if (any(has_analogs)) {
            expect_equal(
                  result_ridge$coef_intercept[has_analogs],
                  result_ridge$weighted_mean[has_analogs],
                  tolerance = 0.01
            )

            # Covariate coefficients should be near zero
            expect_true(all(abs(result_ridge$coef_c1[has_analogs]) < 0.01))
            expect_true(all(abs(result_ridge$coef_c2[has_analogs]) < 0.01))
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
            env = kernel(max = 0.0001),
            geog = kernel(max = 0.0001),
            coord_type = "projected",
            se = "ess"
      )

      # Rows with count == 0 should have NA coefficients
      zero_rows <- result$count == 0
      if (any(zero_rows)) {
            expect_true(all(is.na(result$coef_intercept[zero_rows])))
            expect_true(all(is.na(result$coef_cov1[zero_rows])))
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
            env = kernel(max = 2),
            geog = kernel(max = 2),
            coord_type = "projected"
      )

      # Should have 3 columns per y variable: intercept, slope, aspect
      expected_cols <- c("coef_intercept_biomass", "coef_slope_biomass", "coef_aspect_biomass",
                         "coef_intercept_richness", "coef_slope_richness", "coef_aspect_richness")
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
            env = kernel("gaussian", theta = 0.5, max = 2),
            geog = kernel(max = 2),
            coord_type = "projected"
      )

      expect_true(all(c("count", "ess", "weighted_mean",
                        "coef_intercept", "coef_topo") %in% names(result)))
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
                  env = kernel(max = 2),
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
                  env = kernel(max = 2),
                  coord_type = "projected"
            ),
            "covariates"
      )
})


test_that("regression runs with uniform weights (ordinary regression)", {

      d <- sim_test_data()
      n_ref <- nrow(d$ref)

      # A regression with no distance weighting (uniform env/geog) is valid:
      # it is ordinary (unweighted) least-squares regression over the selected
      # analogs. A non-uniform kernel would make it kernel-weighted regression.
      expect_no_error(
            analog_search(
                  x = d$focal, pool = d$ref,
                  stat = "regression",
                  y = rnorm(n_ref),
                  covariates = matrix(rnorm(n_ref), ncol = 1),
                  env = kernel(max = 2),
                  coord_type = "projected"
            )
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
                  env = kernel(max = 2),
                  lambda = -1,
                  coord_type = "projected"
            ),
            "lambda"
      )
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
            env = kernel(max = 0.01),  # very tight — may get 0 or very few analogs
            geog = kernel(max = 0.01),
            lambda = 0,
            coord_type = "projected"
      )

      # If count is 0, all coefficients should be NA
      if (result$count == 0) {
            expect_true(is.na(result$coef_intercept))
            expect_true(is.na(result$coef_a))
            expect_true(is.na(result$coef_b))
            expect_true(is.na(result$coef_c))
      }
})



# x_covariates prediction --------------------------------

# New tests for the x_covariates parameter in analog_regression()
# Add these to test-regression.R

test_that("analog_regression without x_covariates returns coefficients only", {
      d <- sim_test_data()
      y <- rnorm(nrow(d$ref))
      covs <- data.frame(a = rnorm(nrow(d$ref)))

      fit <- analog_regression(
            x = d$ref, pool = d$ref, y = y, covariates = covs,
            env = kernel("gaussian", theta = 0.3, max = 1),
            geog = kernel(max = 2),
            coord_type = "projected", select = "all"
      )

      expect_true("coef_intercept" %in% names(fit))
      expect_true("coef_a"         %in% names(fit))
      expect_false("pred"          %in% names(fit))
})


test_that("analog_regression with x_covariates adds pred column (single-y)", {
      d <- sim_test_data()
      y <- rnorm(nrow(d$ref))
      covs <- data.frame(a = rnorm(nrow(d$ref)), b = rnorm(nrow(d$ref)))

      fit <- analog_regression(
            x = d$ref, pool = d$ref, y = y,
            covariates   = covs,
            x_covariates = covs,
            env = kernel("gaussian", theta = 0.3, max = 1),
            geog = kernel(max = 2),
            coord_type = "projected", select = "all"
      )

      expect_true("pred" %in% names(fit))

      # Manual arithmetic check
      ok <- !is.na(fit$coef_intercept)
      manual <- fit$coef_intercept[ok] +
            fit$coef_a[ok] * covs$a[ok] +
            fit$coef_b[ok] * covs$b[ok]
      expect_equal(fit$pred[ok], manual, tolerance = 1e-10)
})


test_that("analog_regression with x_covariates adds pred_{y} columns (multi-y)", {
      d <- sim_test_data()
      y <- cbind(a = rnorm(nrow(d$ref)), b = rnorm(nrow(d$ref)))
      covs <- data.frame(u = rnorm(nrow(d$ref)), v = rnorm(nrow(d$ref)))

      fit <- analog_regression(
            x = d$ref, pool = d$ref, y = y,
            covariates   = covs,
            x_covariates = covs,
            env = kernel("gaussian", theta = 0.3, max = 1),
            geog = kernel(max = 2),
            coord_type = "projected", select = "all"
      )

      expect_true(all(c("pred_a", "pred_b") %in% names(fit)))

      # Manual check for y = "a"
      ok <- !is.na(fit$coef_intercept_a)
      manual_a <- fit$coef_intercept_a[ok] +
            fit$coef_u_a[ok] * covs$u[ok] +
            fit$coef_v_a[ok] * covs$v[ok]
      expect_equal(fit$pred_a[ok], manual_a, tolerance = 1e-10)
})


test_that("analog_regression x_covariates allows distinct focal/pool inputs", {
      # Realistic use case: focal is a different dataset from pool, and the
      # user supplies covariate values at the focal locations separately.
      d1 <- sim_test_data()
      d2 <- sim_test_data()  # distinct draw
      y_pool <- rnorm(nrow(d1$ref))
      cov_pool  <- data.frame(a = rnorm(nrow(d1$ref)))
      cov_focal <- data.frame(a = rnorm(nrow(d2$ref)))

      fit <- analog_regression(
            x = d2$ref, pool = d1$ref, y = y_pool,
            covariates   = cov_pool,
            x_covariates = cov_focal,
            env = kernel("gaussian", theta = 0.3, max = 1),
            geog = kernel(max = 2),
            coord_type = "projected", select = "all"
      )

      expect_true("pred" %in% names(fit))
      # pred uses focal covariates, not pool covariates
      ok <- !is.na(fit$coef_intercept)
      manual <- fit$coef_intercept[ok] + fit$coef_a[ok] * cov_focal$a[ok]
      expect_equal(fit$pred[ok], manual, tolerance = 1e-10)
})


test_that("analog_regression errors when x_covariates has wrong row count", {
      d <- sim_test_data()
      y <- rnorm(nrow(d$ref))
      covs <- data.frame(a = rnorm(nrow(d$ref)))

      expect_error(
            analog_regression(
                  x = d$ref, pool = d$ref, y = y,
                  covariates   = covs,
                  x_covariates = covs[1:3, , drop = FALSE],
                  env = kernel("gaussian", theta = 0.3, max = 1),
                  geog = kernel(max = 2),
                  coord_type = "projected", select = "all"
            ),
            "rows/cells"
      )
})


test_that("analog_regression errors when x_covariates has missing covariate", {
      d <- sim_test_data()
      y <- rnorm(nrow(d$ref))
      covs <- data.frame(a = rnorm(nrow(d$ref)), b = rnorm(nrow(d$ref)))

      # x_covariates is missing 'b'
      expect_error(
            analog_regression(
                  x = d$ref, pool = d$ref, y = y,
                  covariates   = covs,
                  x_covariates = covs[, "a", drop = FALSE],
                  env = kernel("gaussian", theta = 0.3, max = 1),
                  geog = kernel(max = 2),
                  coord_type = "projected", select = "all"
            ),
            "missing columns"
      )
})


test_that("analog_regression x_covariates accepts column-reordered input", {
      # Name-based alignment: column order in x_covariates need not match
      # covariates.
      d <- sim_test_data()
      y <- rnorm(nrow(d$ref))
      covs <- data.frame(a = rnorm(nrow(d$ref)), b = rnorm(nrow(d$ref)))

      fit <- analog_regression(
            x = d$ref, pool = d$ref, y = y,
            covariates   = covs,
            x_covariates = covs[, c("b", "a")],   # reversed
            env = kernel("gaussian", theta = 0.3, max = 1),
            geog = kernel(max = 2),
            coord_type = "projected", select = "all"
      )

      ok <- !is.na(fit$coef_intercept)
      manual <- fit$coef_intercept[ok] +
            fit$coef_a[ok] * covs$a[ok] +
            fit$coef_b[ok] * covs$b[ok]
      expect_equal(fit$pred[ok], manual, tolerance = 1e-10)
})


test_that("analog_regression with x_covariates works on SpatRaster inputs", {
      skip_if_not_installed("terra")

      set.seed(42)
      r <- terra::rast(matrix(runif(100), 10))
      names(r) <- "var1"
      yr <- terra::setValues(r, rnorm(100))
      names(yr) <- "y"
      cov_r <- c(
            setNames(terra::setValues(r, rnorm(100)), "a"),
            setNames(terra::setValues(r, rnorm(100)), "b")
      )

      fit <- analog_regression(
            x = r, pool = r, y = yr,
            covariates   = cov_r,
            x_covariates = cov_r,
            env = kernel("gaussian", theta = 0.3, max = 1),
            geog = kernel(max = 5),
            coord_type = "projected", select = "all"
      )

      expect_s4_class(fit, "SpatRaster")
      expect_true("pred" %in% terra::names(fit))
})


test_that("analog_cv on analog_regression receives pred via x_covariates path", {
      # Residuals should still agree with manual arithmetic -- this verifies
      # the refactored pipeline computes predictions correctly.
      d <- sim_test_data()
      y <- rnorm(nrow(d$ref))
      covs <- data.frame(a = rnorm(nrow(d$ref)), b = rnorm(nrow(d$ref)))

      cv <- analog_cv(
            fun = analog_regression,
            pool = d$ref,
            y = y,
            covariates = covs,
            env = kernel("gaussian", theta = 0.3, max = 1),
            geog = kernel(max = 2),
            coord_type = "projected", select = "all",
            cv_method = "loo"
      )

      expect_true(all(c("pred", "obs", "residual") %in% names(cv)))

      ok <- !is.na(cv$coef_intercept)
      manual_pred <- cv$coef_intercept[ok] +
            cv$coef_a[ok] * covs$a[ok] +
            cv$coef_b[ok] * covs$b[ok]
      expect_equal(cv$pred[ok], manual_pred, tolerance = 1e-10)
      expect_equal(cv$residual[ok], cv$obs[ok] - cv$pred[ok], tolerance = 1e-10)
})
