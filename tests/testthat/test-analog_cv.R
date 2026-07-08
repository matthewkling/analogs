# analog_search(exclude_self) ----------------------------------

test_that("exclude_self = TRUE validates inputs", {
      d <- sim_test_data()

      # Must be same object
      expect_error(
            analog_search(
                  x = d$focal,
                  pool = d$ref,
                  stat = "count",
                  env = kernel(max = 1),
                  coord_type = "projected",
                  exclude_self = TRUE
            ),
            "identical"
      )

      # Pre-built index disallowed
      idx <- build_analog_index(d$ref, coord_type = "projected")
      expect_error(
            analog_search(
                  x = idx,
                  pool = idx,
                  stat = "count",
                  env = kernel(max = 1),
                  exclude_self = TRUE
            )
      )

      # Progress disallowed
      expect_error(
            analog_search(
                  x = d$ref,
                  pool = d$ref,
                  stat = "count",
                  env = kernel(max = 1),
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
                  env = kernel(max = 1),
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
            stat = "count", env = kernel(max = Inf),
            coord_type = "projected",
            exclude_self = FALSE
      )
      cnt_without <- analog_search(
            x = d$ref, pool = d$ref,
            stat = "count", env = kernel(max = Inf),
            coord_type = "projected",
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
            env = kernel("gaussian", theta = 0.3, max = 1),
            geog = kernel(max = 2),
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
            env = kernel("gaussian", theta = 0.3, max = 1),
            geog = kernel(max = 2),
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
            env = kernel("gaussian", theta = 0.3, max = 1),
            geog = kernel(max = 2),
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
                  env = kernel(max = 1),
                  cv_method = "loo"
            ),
            "analog_impact.*analog_regression.*analog_search"
      )

      # Missing y
      expect_error(
            analog_cv(
                  fun = analog_impact,
                  pool = d$ref,
                  env = kernel(max = 1),
                  cv_method = "loo"
            ),
            "y"
      )

      # kfold without n_folds or fold_id
      expect_error(
            analog_cv(
                  fun = analog_impact,
                  pool = d$ref, y = y,
                  env = kernel(max = 1),
                  cv_method = "kfold"
            ),
            "n_folds|fold_id"
      )

      # loo with n_folds
      expect_error(
            analog_cv(
                  fun = analog_impact,
                  pool = d$ref, y = y,
                  env = kernel(max = 1),
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
                  env = kernel(max = 1),
                  cv_method = "loo"
            ),
            "covariates"
      )

      # progress = TRUE not allowed
      expect_error(
            analog_cv(
                  fun = analog_impact,
                  pool = d$ref, y = y,
                  env = kernel(max = 1),
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
            env = kernel("gaussian", theta = 0.3, max = 1),
            geog = kernel(max = 2),
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
            env = kernel("gaussian", theta = 0.3, max = 1),
            geog = kernel(max = 2),
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
            select     = "knn_env",
            k          = 5,
            stat       = c("count", "weighted_mean"),
            env       = kernel("gaussian", theta = 0.3),
            geog       = kernel(max = 2),
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
                  env       = kernel("gaussian", theta = 0.3, max = 1),
                  geog       = kernel(max = 2),
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
                  env       = kernel(max = 1),
                  geog       = kernel(max = 2),
                  coord_type = "projected",
                  cv_method  = "loo"
            ),
            "skipping obs/residual"
      )

      expect_true("count" %in% names(res))
      expect_false("obs" %in% names(res))
      expect_false("residual" %in% names(res))
})



# test stat = "tablulate" -----------------------------------------

# Reuse the local synthetic dataset helper from test-tabulate.R
.tabcv_sim <- function(n = 80, n_env = 2, seed = 1L) {
      set.seed(seed)
      x <- runif(n, 0, 10); y <- runif(n, 0, 10)
      env <- matrix(rnorm(n * n_env), ncol = n_env,
                     dimnames = list(NULL, paste0("c", seq_len(n_env))))
      pool <- cbind(x = x, y = y, env)
      list(pool = pool)
}


test_that("LOO tabulate CV adds obs / primary / brier columns", {
      d <- .tabcv_sim()
      veg <- factor(rep(c("forest", "grass", "shrub"),
                        length.out = nrow(d$pool)))

      cv <- analog_cv(
            fun        = analog_impact,
            pool       = d$pool,
            y          = veg,
            cv_method  = "loo",
            stat       = c("count", "sum_weights", "tabulate"),
            env       = kernel("gaussian", theta = 0.5, max = 5),
            geog       = kernel(max = 100),
            coord_type = "projected"
      )

      # All three CV-specific tabulate columns should be present
      expect_true(all(c("obs", "primary", "brier") %in% names(cv)))

      # obs and primary are character/character-like
      expect_type(cv$obs, "character")
      expect_type(cv$primary, "character")

      # Brier in 0..2
      expect_true(all(is.na(cv$brier) | (cv$brier >= 0 & cv$brier <= 2)))

      # The n_<level> columns should still be present
      expect_true(all(c("n_forest", "n_grass", "n_shrub") %in% names(cv)))
})


test_that("LOO tabulate CV: brier and primary are consistent with vote columns", {
      d <- .tabcv_sim()
      veg <- factor(rep(c("a", "b", "c"), length.out = nrow(d$pool)))

      cv <- analog_cv(
            fun = analog_impact,
            pool = d$pool, y = veg,
            stat = c("sum_weights", "tabulate"),
            env = kernel("gaussian", theta = 0.5, max = 5),
            geog = kernel(max = 100), coord_type = "projected"
      )

      # primary should match argmax of n_<level>
      vote_cols <- cv[, c("n_a", "n_b", "n_c")]
      expected_primary <- c("a", "b", "c")[max.col(as.matrix(vote_cols),
                                                   ties.method = "first")]
      # NA where row sums are 0
      row_sums <- rowSums(as.matrix(vote_cols))
      expected_primary[row_sums == 0 | !is.finite(row_sums)] <- NA_character_

      expect_equal(cv$primary, expected_primary)

      # Brier hand-computation on first non-NA row
      finite_idx <- which(!is.na(cv$brier))
      expect_true(length(finite_idx) > 0)
      i <- finite_idx[1]
      probs <- as.numeric(vote_cols[i, ]) / row_sums[i]
      ind <- as.numeric(c("a", "b", "c") == cv$obs[i])
      expected_brier <- sum((probs - ind)^2)
      expect_equal(cv$brier[i], expected_brier, tolerance = 1e-10)
})


test_that("kfold tabulate CV produces consistent columns across folds", {
      d <- .tabcv_sim(n = 100)
      veg <- factor(rep(c("a", "b"), length.out = nrow(d$pool)))

      cv <- analog_cv(
            fun = analog_impact,
            pool = d$pool, y = veg,
            cv_method = "kfold", n_folds = 5,
            stat = c("sum_weights", "tabulate"),
            env = kernel("gaussian", theta = 0.5, max = 5),
            geog = kernel(max = 100), coord_type = "projected"
      )

      expect_true(all(c("obs", "primary", "brier", "n_a", "n_b", "fold")
                      %in% names(cv)))
      expect_equal(nrow(cv), nrow(d$pool))
      expect_equal(length(unique(cv$fold)), 5L)
})


test_that("tabulate CV handles NA in y by setting CV cols to NA", {
      d <- .tabcv_sim()
      veg <- factor(rep(c("a", "b", "c"), length.out = nrow(d$pool)))
      veg[c(3, 11, 27)] <- NA

      cv <- analog_cv(
            fun = analog_impact,
            pool = d$pool, y = veg,
            stat = c("sum_weights", "tabulate"),
            env = kernel("gaussian", theta = 0.5, max = 5),
            geog = kernel(max = 100), coord_type = "projected"
      )

      # For NA-observed focals, obs and brier should be NA
      expect_true(all(is.na(cv$obs[c(3, 11, 27)])))
      expect_true(all(is.na(cv$brier[c(3, 11, 27)])))
      # primary may still be defined (analogs may have non-NA y), but obs
      # being NA means we can't score correctness for those focals
})


test_that("tabulate CV with multi-y produces var-suffixed CV columns", {
      d <- .tabcv_sim()
      yA <- factor(rep(c("x", "y"), length.out = nrow(d$pool)))
      yB <- factor(rep(c("p", "q", "r"), length.out = nrow(d$pool)))
      Y <- data.frame(habitat = yA, soil = yB)

      cv <- analog_cv(
            fun = analog_impact,
            pool = d$pool, y = Y,
            stat = "tabulate",
            env = kernel("gaussian", theta = 0.5, max = 5),
            geog = kernel(max = 100), coord_type = "projected"
      )

      # Per-y CV columns
      expect_true(all(c("obs_habitat", "primary_habitat", "brier_habitat")
                      %in% names(cv)))
      expect_true(all(c("obs_soil", "primary_soil", "brier_soil")
                      %in% names(cv)))
      # And the underlying vote columns
      expect_true(all(c("habitat_n_x", "habitat_n_y",
                        "soil_n_p", "soil_n_q", "soil_n_r")
                      %in% names(cv)))
})


test_that("incompatible tabulate + weighted_mean errors in CV", {
      d <- .tabcv_sim()
      veg <- factor(rep(c("a", "b"), length.out = nrow(d$pool)))

      # The .validate_query_params block in utils.R catches this before
      # analog_cv's own check, but either way an error is expected.
      expect_error(
            analog_cv(
                  fun = analog_impact,
                  pool = d$pool, y = veg,
                  stat = c("weighted_mean", "tabulate"),
                  env = kernel("gaussian", theta = 0.5, max = 5),
                  geog = kernel(max = 100), coord_type = "projected"
            )
      )
})


test_that("all-NA tabulate y errors clearly under CV", {
      d <- .tabcv_sim()
      veg <- factor(rep(NA_character_, nrow(d$pool)))

      expect_error(
            analog_cv(
                  fun = analog_impact,
                  pool = d$pool, y = veg,
                  stat = "tabulate",
                  env = kernel("gaussian", theta = 0.5, max = 5),
                  geog = kernel(max = 100), coord_type = "projected"
            ),
            "no non-NA"
      )
})


test_that("character / integer y is accepted in tabulate CV", {
      d <- .tabcv_sim()
      veg_chr <- rep(c("a", "b", "c"), length.out = nrow(d$pool))
      veg_int <- rep(1:3, length.out = nrow(d$pool))

      cv_chr <- analog_cv(
            fun = analog_impact,
            pool = d$pool, y = veg_chr,
            stat = "tabulate",
            env = kernel("gaussian", theta = 0.5, max = 5),
            geog = kernel(max = 100), coord_type = "projected"
      )
      expect_type(cv_chr$obs, "character")
      expect_true(all(cv_chr$obs %in% c("a", "b", "c", NA_character_)))

      cv_int <- analog_cv(
            fun = analog_impact,
            pool = d$pool, y = veg_int,
            stat = "tabulate",
            env = kernel("gaussian", theta = 0.5, max = 5),
            geog = kernel(max = 100), coord_type = "projected"
      )
      # Levels are "1", "2", "3" (factor coercion of integer)
      expect_type(cv_int$obs, "character")
      expect_true(all(cv_int$obs %in% c("1", "2", "3", NA_character_)))
})


test_that("CV metadata attributes set correctly for tabulate", {
      d <- .tabcv_sim()
      veg <- factor(rep(c("a", "b"), length.out = nrow(d$pool)))

      cv <- analog_cv(
            fun = analog_impact,
            pool = d$pool, y = veg,
            cv_method = "kfold", n_folds = 4,
            stat = "tabulate",
            env = kernel("gaussian", theta = 0.5, max = 5),
            geog = kernel(max = 100), coord_type = "projected"
      )

      m <- metadata(cv)
      expect_identical(m$cv_method, "kfold")
      expect_identical(m$cv_pred_target, "tabulate")
      expect_identical(m$cv_n_folds, 4L)
      expect_identical(m$cv_fun, "analog_impact")
})
