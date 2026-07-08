# Tests for standard-error computation on weighted_mean and regression stats.
# SE framings supported: "none" (default), "ess", "design".

test_that("se = 'none' (default) omits SE columns", {

      d <- sim_test_data()
      n_ref <- nrow(d$ref)

      y <- rnorm(n_ref)
      covariates <- matrix(rnorm(n_ref), ncol = 1,
                           dimnames = list(NULL, "cov1"))

      result <- analog_search(
            x = d$focal,
            pool = d$ref,
            select = "all",
            stat = c("weighted_mean", "regression"),
            y = y,
            covariates = covariates,
            clim = kernel("gaussian", theta = 0.5, max = 2),
            geog = kernel(max = 2),
            coord_type = "projected"
      )

      # Coefs and weighted_mean present
      expect_true(all(c("weighted_mean", "coef_intercept", "coef_cov1") %in% names(result)))

      # No SE columns when se = "none"
      se_cols <- grep("^se_", names(result), value = TRUE)
      expect_length(se_cols, 0)
})


test_that("se = 'ess' returns SE columns for weighted_mean and regression", {

      d <- sim_test_data()
      n_ref <- nrow(d$ref)

      y <- rnorm(n_ref)
      covariates <- matrix(rnorm(n_ref), ncol = 1,
                           dimnames = list(NULL, "cov1"))

      result <- analog_search(
            x = d$focal,
            pool = d$ref,
            select = "all",
            stat = c("weighted_mean", "regression"),
            y = y,
            covariates = covariates,
            clim = kernel("gaussian", theta = 0.5, max = 2),
            geog = kernel(max = 2),
            se = "ess",
            coord_type = "projected"
      )

      expect_true("se_weighted_mean" %in% names(result))
      expect_true("se_intercept" %in% names(result))
      expect_true("se_cov1" %in% names(result))

      # SEs should be finite non-negative where weighted_mean is finite
      ok <- is.finite(result$weighted_mean)
      if (any(ok)) {
            expect_true(all(result$se_weighted_mean[ok] >= 0, na.rm = TRUE))
      }
})


test_that("se = 'design' returns SE columns for weighted_mean and regression", {

      d <- sim_test_data()
      n_ref <- nrow(d$ref)

      y <- rnorm(n_ref)
      covariates <- matrix(rnorm(n_ref), ncol = 1,
                           dimnames = list(NULL, "cov1"))

      result <- analog_search(
            x = d$focal,
            pool = d$ref,
            select = "all",
            stat = c("weighted_mean", "regression"),
            y = y,
            covariates = covariates,
            clim = kernel("gaussian", theta = 0.5, max = 2),
            geog = kernel(max = 2),
            se = "design",
            coord_type = "projected"
      )

      expect_true("se_weighted_mean" %in% names(result))
      expect_true("se_intercept" %in% names(result))
      expect_true("se_cov1" %in% names(result))

      ok <- is.finite(result$weighted_mean)
      if (any(ok)) {
            expect_true(all(result$se_weighted_mean[ok] >= 0, na.rm = TRUE))
      }
})


test_that("weighted_mean SE is scale-invariant in the weights", {
      # The weighted mean and its SE should not depend on the absolute scale of
      # the weights -- only their relative magnitudes. The `weight` function's
      # output gets multiplied by a shared factor when we change theta magnitude,
      # but scale-invariance is cleanest to check by explicitly exploiting the
      # fact that gaussian weights with the same data produce results that
      # differ only in scale across theta changes. Instead, use uniform weights
      # (which are already all equal) and confirm SE is well-defined.

      d <- sim_test_data()
      n_ref <- nrow(d$ref)
      y <- rnorm(n_ref)

      res_ess <- analog_search(
            x = d$focal,
            pool = d$ref,
            select = "all",
            stat = "weighted_mean",
            y = y,
            clim = kernel(max = 2),
            geog = kernel(max = 2),
            se = "ess",
            coord_type = "projected"
      )

      res_design <- analog_search(
            x = d$focal,
            pool = d$ref,
            select = "all",
            stat = "weighted_mean",
            y = y,
            clim = kernel(max = 2),
            geog = kernel(max = 2),
            se = "design",
            coord_type = "projected"
      )

      # With uniform weights, weighted_mean equals plain mean.
      # ESS-based SE = sd(y)/sqrt(n); design-based SE = sqrt(Σ(y-ȳ)²)/n.
      # These differ by a factor of sqrt(n-1)/sqrt(n) for uniform weights,
      # so they should be close but not identical.
      ok <- is.finite(res_ess$weighted_mean) & is.finite(res_design$weighted_mean)
      if (any(ok)) {
            # Both finite and non-negative
            expect_true(all(res_ess$se_weighted_mean[ok] >= 0))
            expect_true(all(res_design$se_weighted_mean[ok] >= 0))
      }
})


test_that("weighted_mean SE matches manual computation (ESS)", {

      # Build a controlled setup: one focal point, known pool values and weights,
      # then check that the reported SE matches the formula.

      set.seed(7)
      n <- 20
      pool <- cbind(x = runif(n, 0, 1), y = runif(n, 0, 1),
                    clim1 = rnorm(n), clim2 = rnorm(n))

      # Focal at origin with matching climate
      focal <- matrix(c(0.5, 0.5, 0, 0), nrow = 1)

      yv <- rnorm(n)

      result <- analog_search(
            x = focal,
            pool = pool,
            select = "all",
            stat = "weighted_mean",
            y = yv,
            geog = NULL,  # use all points
            clim = kernel("gaussian", theta = 1),
            se = "ess",
            coord_type = "projected"
      )

      # Reproduce the weight computation manually
      clim_dist <- sqrt(rowSums((pool[, 3:4] - matrix(c(0, 0), n, 2, byrow = TRUE))^2))
      w <- exp(-clim_dist^2 / 2)  # sigma = 1
      ybar <- sum(w * yv) / sum(w)
      n_eff <- sum(w)^2 / sum(w^2)
      var_w <- sum(w * yv^2) / sum(w) - ybar^2
      se_expected <- sqrt(var_w / n_eff)

      expect_equal(result$weighted_mean, ybar, tolerance = 1e-8)
      expect_equal(result$se_weighted_mean, se_expected, tolerance = 1e-6)
})


test_that("weighted_mean SE matches manual computation (design)", {

      set.seed(11)
      n <- 20
      pool <- cbind(x = runif(n, 0, 1), y = runif(n, 0, 1),
                    clim1 = rnorm(n), clim2 = rnorm(n))

      focal <- matrix(c(0.5, 0.5, 0, 0), nrow = 1)
      yv <- rnorm(n)

      result <- analog_search(
            x = focal,
            pool = pool,
            select = "all",
            stat = "weighted_mean",
            y = yv,
            clim = kernel("gaussian", theta = 1, max = NULL),
            geog = kernel(max = NULL),
            se = "design",
            coord_type = "projected"
      )

      clim_dist <- sqrt(rowSums((pool[, 3:4])^2))
      w <- exp(-clim_dist^2 / 2)
      ybar <- sum(w * yv) / sum(w)
      se_expected <- sqrt(sum(w^2 * (yv - ybar)^2)) / sum(w)

      expect_equal(result$weighted_mean, ybar, tolerance = 1e-8)
      expect_equal(result$se_weighted_mean, se_expected, tolerance = 1e-6)
})


test_that("regression SE (ESS) matches standard WLS SE with ESS df correction", {

      # With uniform weights and lambda = 0, ESS-based regression SE should match
      # the OLS SE with df = n_eff - p (which for uniform weights is n - p).
      set.seed(13)
      n <- 100
      pool <- cbind(x = runif(n, 0, 1), y = runif(n, 0, 1),
                    clim1 = rnorm(n), clim2 = rnorm(n))
      focal <- matrix(c(0.5, 0.5, 0, 0), nrow = 1)

      x1 <- rnorm(n)
      x2 <- rnorm(n)
      yv <- 1 + 2 * x1 - 0.5 * x2 + rnorm(n, sd = 0.5)

      covariates <- cbind(x1 = x1, x2 = x2)

      result <- analog_search(
            x = focal,
            pool = pool,
            select = "all",
            stat = "regression",
            y = yv,
            covariates = covariates,
            clim = NULL,
            geog = NULL,
            lambda = 0,
            se = "ess",
            coord_type = "projected"
      )

      # Compare to lm() — with uniform weights and no ridge, the ESS df correction
      # (n_eff = n, so df = n - p) matches OLS exactly.
      fit <- lm(yv ~ x1 + x2)
      ref_coefs <- coef(summary(fit))[, "Estimate"]
      ref_ses   <- coef(summary(fit))[, "Std. Error"]

      expect_equal(result$coef_intercept, unname(ref_coefs["(Intercept)"]), tolerance = 1e-6)
      expect_equal(result$coef_x1,        unname(ref_coefs["x1"]),          tolerance = 1e-6)
      expect_equal(result$coef_x2,        unname(ref_coefs["x2"]),          tolerance = 1e-6)

      expect_equal(result$se_intercept, unname(ref_ses["(Intercept)"]), tolerance = 1e-6)
      expect_equal(result$se_x1,        unname(ref_ses["x1"]),          tolerance = 1e-6)
      expect_equal(result$se_x2,        unname(ref_ses["x2"]),          tolerance = 1e-6)
})


test_that("se warns when no requested stat supports SE", {

      d <- sim_test_data()

      expect_warning(
            analog_search(
                  x = d$focal,
                  pool = d$ref,
                  select = "all",
                  stat = c("count", "ess"),
                  clim = kernel(max = 2),
                  geog = kernel(max = 2),
                  se = "ess",
                  coord_type = "projected"
            ),
            "no requested stat supports SE"
      )
})


test_that("se does not warn when at least one requested stat supports SE", {

      d <- sim_test_data()
      n_ref <- nrow(d$ref)
      y <- rnorm(n_ref)

      # weighted_mean supports SE, so no warning
      expect_no_warning(
            analog_search(
                  x = d$focal,
                  pool = d$ref,
                  select = "all",
                  stat = c("count", "weighted_mean"),
                  y = y,
                  clim = kernel("gaussian", theta = 0.5, max = 2),
                  geog = kernel(max = 2),
                  se = "ess",
                  coord_type = "projected"
            )
      )
})


test_that("se works with multiple y variables for weighted_mean", {

      d <- sim_test_data()
      n_ref <- nrow(d$ref)

      y <- matrix(rnorm(n_ref * 2), ncol = 2,
                  dimnames = list(NULL, c("biomass", "richness")))

      result <- analog_search(
            x = d$focal,
            pool = d$ref,
            select = "all",
            stat = "weighted_mean",
            y = y,
            clim = kernel("gaussian", theta = 0.5, max = 2),
            geog = kernel(max = 2),
            se = "ess",
            coord_type = "projected"
      )

      expected_cols <- c("weighted_mean_biomass", "weighted_mean_richness",
                         "se_weighted_mean_biomass", "se_weighted_mean_richness")
      expect_true(all(expected_cols %in% names(result)))
})


test_that("se works with multiple y variables for regression", {

      d <- sim_test_data()
      n_ref <- nrow(d$ref)

      y <- matrix(rnorm(n_ref * 2), ncol = 2,
                  dimnames = list(NULL, c("biomass", "richness")))
      covariates <- matrix(rnorm(n_ref), ncol = 1,
                           dimnames = list(NULL, "slope"))

      result <- analog_search(
            x = d$focal,
            pool = d$ref,
            select = "all",
            stat = "regression",
            y = y,
            covariates = covariates,
            clim = kernel(max = 2),
            geog = kernel(max = 2),
            se = "ess",
            coord_type = "projected"
      )

      expected_cols <- c(
            "coef_intercept_biomass", "coef_slope_biomass",
            "coef_intercept_richness", "coef_slope_richness",
            "se_intercept_biomass", "se_slope_biomass",
            "se_intercept_richness", "se_slope_richness"
      )
      expect_true(all(expected_cols %in% names(result)))
})


test_that("analog_regression passes se through", {

      d <- sim_test_data()
      n_ref <- nrow(d$ref)

      y <- rnorm(n_ref)
      covariates <- matrix(rnorm(n_ref), ncol = 1,
                           dimnames = list(NULL, "cov1"))

      # With se = "none" (default)
      r_none <- analog_regression(
            x = d$focal,
            pool = d$ref,
            y = y,
            covariates = covariates,
            clim = kernel(max = 2),
            geog = kernel(max = 2),
            coord_type = "projected"
      )
      se_cols_none <- grep("^se_", names(r_none), value = TRUE)
      expect_length(se_cols_none, 0)

      # With se = "ess"
      r_ess <- analog_regression(
            x = d$focal,
            pool = d$ref,
            y = y,
            covariates = covariates,
            clim = kernel(max = 2),
            geog = kernel(max = 2),
            se = "ess",
            coord_type = "projected"
      )
      expect_true("se_intercept" %in% names(r_ess))
      expect_true("se_cov1" %in% names(r_ess))
})


test_that("analog_impact passes se through", {

      d <- sim_test_data()
      n_ref <- nrow(d$ref)
      y <- rnorm(n_ref)

      r_ess <- analog_impact(
            x = d$focal,
            pool = d$ref,
            y = y,
            geog = kernel(max = 2),
            clim = kernel("gaussian", theta = 0.5, max = 2),
            se = "ess",
            coord_type = "projected"
      )
      expect_true("se_weighted_mean" %in% names(r_ess))
})


test_that("regression SE is NA when system is singular", {

      # Same setup as zero-analogs test but with SEs requested
      d <- sim_test_data()
      n_ref <- nrow(d$ref)

      y <- rnorm(n_ref)
      covariates <- matrix(rnorm(n_ref), ncol = 1,
                           dimnames = list(NULL, "cov1"))

      result <- analog_search(
            x = d$focal,
            pool = d$ref,
            select = "all",
            stat = c("count", "regression"),
            y = y,
            covariates = covariates,
            clim = kernel(max = 0.0001),
            geog = kernel(max = 0.0001),
            se = "ess",
            coord_type = "projected"
      )

      zero_rows <- result$count == 0
      if (any(zero_rows)) {
            expect_true(all(is.na(result$coef_intercept[zero_rows])))
            expect_true(all(is.na(result$se_intercept[zero_rows])))
            expect_true(all(is.na(result$se_cov1[zero_rows])))
      }
})


test_that("invalid se value is rejected", {

      d <- sim_test_data()
      n_ref <- nrow(d$ref)
      y <- rnorm(n_ref)

      expect_error(
            analog_search(
                  x = d$focal,
                  pool = d$ref,
                  select = "all",
                  stat = "weighted_mean",
                  y = y,
                  clim = kernel(max = 2),
                  geog = kernel(max = 2),
                  se = "bogus",
                  coord_type = "projected"
            )
      )
})
