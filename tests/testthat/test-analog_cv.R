
# analog_search(exclude_self) ----------------------------------

test_that("exclude_self = TRUE validates inputs", {
      d <- sim_test_data()

      # Must be same object
      expect_error(
            analog_search(
                  x = d$focal,
                  pool = d$ref,
                  stat = "count",
                  max_clim = 1,
                  coord_type = "projected",
                  exclude_self = TRUE
            ),
            "identical"
      )

      # Pre-built index disallowed
      idx <- build_analog_index(d$ref, coord_type = "projected", index_res = 8)
      expect_error(
            analog_search(
                  x = idx,
                  pool = idx,
                  stat = "count",
                  max_clim = 1,
                  exclude_self = TRUE
            )
      )

      # Progress disallowed
      expect_error(
            analog_search(
                  x = d$ref,
                  pool = d$ref,
                  stat = "count",
                  max_clim = 1,
                  coord_type = "projected",
                  exclude_self = TRUE,
                  progress = TRUE
            ),
            "progress"
      )

      # Downsample != 1 disallowed
      expect_error(
            analog_search(
                  x = d$ref,
                  pool = d$ref,
                  stat = "count",
                  max_clim = 1,
                  coord_type = "projected",
                  exclude_self = TRUE,
                  downsample = 0.5
            ),
            "downsample"
      )
})


test_that("exclude_self reduces count by exactly 1 within-pool when focal row passes its own filters", {
      d <- sim_test_data()

      cnt_with <- analog_search(
            x = d$ref, pool = d$ref,
            stat = "count", max_clim = Inf,
            coord_type = "projected", index_res = 8,
            exclude_self = FALSE
      )
      cnt_without <- analog_search(
            x = d$ref, pool = d$ref,
            stat = "count", max_clim = Inf,
            coord_type = "projected", index_res = 8,
            exclude_self = TRUE
      )

      expect_equal(cnt_with$count - cnt_without$count, rep(1, nrow(d$ref)))
})



# analog_cv() ----------------------------------

test_that("analog_cv(fun = analog_impact, cv_method = 'loo') runs and returns expected columns", {
      d <- sim_test_data()
      y <- rnorm(nrow(d$ref))

      res <- analog_cv(
            fun = analog_impact,
            pool = d$ref,
            y = y,
            max_clim = 1,
            max_geog = 2,
            kernel = "gaussian_clim",
            theta = 0.3,
            coord_type = "projected",
            cv_method = "loo"
      )

      expect_s3_class(res, "data.frame")
      expect_equal(nrow(res), nrow(d$ref))
      expect_true("weighted_mean" %in% names(res))
      expect_true("obs" %in% names(res))
      expect_true("residual" %in% names(res))
      expect_equal(res$obs, y)

      # Residual = obs - weighted_mean
      expect_equal(res$residual, res$obs - res$weighted_mean)
})


test_that("analog_cv kfold produces correct number of rows and fold column", {
      d <- sim_test_data()
      y <- rnorm(nrow(d$ref))

      res <- analog_cv(
            fun = analog_impact,
            pool = d$ref,
            y = y,
            max_clim = 1,
            max_geog = 2,
            kernel = "gaussian_clim",
            theta = 0.3,
            coord_type = "projected",
            cv_method = "kfold",
            n_folds = 5
      )

      expect_equal(nrow(res), nrow(d$ref))
      expect_true("fold" %in% names(res))
      expect_true(all(res$fold %in% 1:5))
      expect_true(all(table(res$fold) > 0))

      # Output row order matches pool row order
      expect_equal(res$obs, y)
})


test_that("analog_cv supports analog_regression with covariates", {
      d <- sim_test_data()
      y <- rnorm(nrow(d$ref))
      covs <- data.frame(a = rnorm(nrow(d$ref)), b = rnorm(nrow(d$ref)))

      res <- analog_cv(
            fun = analog_regression,
            pool = d$ref,
            y = y,
            covariates = covs,
            max_clim = 1,
            max_geog = 2,
            kernel = "gaussian_clim",
            theta = 0.3,
            coord_type = "projected",
            select = "all",
            cv_method = "loo"
      )

      expect_true(all(c("coef_intercept", "coef_a", "coef_b", "obs", "residual")
                      %in% names(res)))

      # Residual should equal obs - (intercept + a*cov_a + b*cov_b),
      # for rows where the regression succeeded.
      ok <- !is.na(res$coef_intercept)
      pred_manual <- res$coef_intercept[ok] +
            res$coef_a[ok] * covs$a[ok] +
            res$coef_b[ok] * covs$b[ok]
      expect_equal(res$residual[ok], y[ok] - pred_manual, tolerance = 1e-10)
})


test_that("analog_cv rejects unsupported fun and invalid configs", {
      d <- sim_test_data()
      y <- rnorm(nrow(d$ref))

      # Wrong function
      expect_error(
            analog_cv(
                  fun = analog_velocity,
                  pool = d$ref, y = y,
                  max_clim = 1,
                  cv_method = "loo"
            ),
            "analog_impact.*analog_regression.*analog_search"
      )

      # Missing y
      expect_error(
            analog_cv(
                  fun = analog_impact,
                  pool = d$ref,
                  max_clim = 1,
                  cv_method = "loo"
            ),
            "y"
      )

      # kfold without n_folds or fold_id
      expect_error(
            analog_cv(
                  fun = analog_impact,
                  pool = d$ref, y = y,
                  max_clim = 1,
                  cv_method = "kfold"
            ),
            "n_folds|fold_id"
      )

      # loo with n_folds
      expect_error(
            analog_cv(
                  fun = analog_impact,
                  pool = d$ref, y = y,
                  max_clim = 1,
                  cv_method = "loo",
                  n_folds = 5
            ),
            "NULL"
      )

      # analog_regression without covariates
      expect_error(
            analog_cv(
                  fun = analog_regression,
                  pool = d$ref, y = y,
                  max_clim = 1,
                  cv_method = "loo"
            ),
            "covariates"
      )

      # progress = TRUE not allowed
      expect_error(
            analog_cv(
                  fun = analog_impact,
                  pool = d$ref, y = y,
                  max_clim = 1,
                  cv_method = "loo",
                  progress = TRUE
            ),
            "progress"
      )
})


test_that("analog_cv handles multi-y correctly", {
      d <- sim_test_data()
      y <- cbind(a = rnorm(nrow(d$ref)), b = rnorm(nrow(d$ref)))

      res <- analog_cv(
            fun = analog_impact,
            pool = d$ref,
            y = y,
            max_clim = 1,
            max_geog = 2,
            kernel = "gaussian_clim",
            theta = 0.3,
            coord_type = "projected",
            cv_method = "loo"
      )

      expect_true(all(c("weighted_mean_a", "weighted_mean_b",
                        "obs_a", "obs_b",
                        "residual_a", "residual_b") %in% names(res)))
      expect_equal(res$obs_a, y[, "a"])
      expect_equal(res$obs_b, y[, "b"])
})


test_that("analog_cv with fold_id respects user-supplied folds", {
      d <- sim_test_data()
      y <- rnorm(nrow(d$ref))
      fold_id <- rep(1:3, length.out = nrow(d$ref))

      res <- analog_cv(
            fun = analog_impact,
            pool = d$ref,
            y = y,
            max_clim = 1,
            max_geog = 2,
            kernel = "gaussian_clim",
            theta = 0.3,
            coord_type = "projected",
            cv_method = "kfold",
            fold_id = fold_id
      )

      expect_equal(res$fold, fold_id)
})


test_that("analog_cv(fun = analog_search) works with weighted_mean stat", {
      d <- sim_test_data()
      y <- rnorm(nrow(d$ref))

      res <- analog_cv(
            fun        = analog_search,
            pool       = d$ref,
            y          = y,
            select     = "knn_clim",
            k          = 5,
            stat       = c("count", "weighted_mean"),
            kernel     = "gaussian_clim",
            theta      = 0.3,
            max_geog   = 2,
            coord_type = "projected",
            cv_method  = "loo"
      )

      expect_true(all(c("weighted_mean", "obs", "residual") %in% names(res)))
      expect_equal(res$obs, y)
      expect_equal(res$residual, res$obs - res$weighted_mean)
})


test_that("analog_cv(fun = analog_search) errors on ambiguous stat", {
      d <- sim_test_data()
      y <- rnorm(nrow(d$ref))
      covs <- data.frame(a = rnorm(nrow(d$ref)))

      expect_error(
            analog_cv(
                  fun        = analog_search,
                  pool       = d$ref,
                  y          = y,
                  covariates = covs,
                  select     = "all",
                  stat       = c("weighted_mean", "regression"),
                  kernel     = "gaussian_clim",
                  theta      = 0.3,
                  max_clim   = 1,
                  max_geog   = 2,
                  coord_type = "projected",
                  cv_method  = "loo"
            ),
            "ambiguous"
      )
})


test_that("analog_cv(fun = analog_search) skips residuals when stat lacks prediction target", {
      d <- sim_test_data()
      y <- rnorm(nrow(d$ref))

      expect_message(
            res <- analog_cv(
                  fun        = analog_search,
                  pool       = d$ref,
                  y          = y,
                  select     = "all",
                  stat       = "count",
                  max_clim   = 1,
                  max_geog   = 2,
                  coord_type = "projected",
                  cv_method  = "loo"
            ),
            "skipping obs/residual"
      )

      expect_true("count" %in% names(res))
      expect_false("obs" %in% names(res))
      expect_false("residual" %in% names(res))
})

