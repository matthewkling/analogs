test_that("cv_performance returns long-format data.frame for continuous single-y", {
      d <- sim_test_data()
      y <- rnorm(nrow(d$ref))

      cv <- analog_cv(
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

      perf <- cv_performance(cv)
      expect_s3_class(perf, "data.frame")
      expect_equal(names(perf), c("variable", "type", "metric", "value"))
      expect_true(all(perf$variable == "y"))
      expect_true(all(perf$type == "continuous"))
      expect_equal(perf$metric, c("n", "rmse", "mae", "bias", "r2"))

      # Manually verify against residuals
      ok <- is.finite(cv$obs) & is.finite(cv$residual)
      r <- cv$residual[ok]; o <- cv$obs[ok]

      expect_equal(perf$value[perf$metric == "n"],    length(o))
      expect_equal(perf$value[perf$metric == "rmse"], sqrt(mean(r^2)), tolerance = 1e-12)
      expect_equal(perf$value[perf$metric == "mae"],  mean(abs(r)),     tolerance = 1e-12)
      expect_equal(perf$value[perf$metric == "bias"], mean(r),          tolerance = 1e-12)
      expect_equal(perf$value[perf$metric == "r2"],
                   1 - sum(r^2) / sum((o - mean(o))^2), tolerance = 1e-12)
})


test_that("cv_performance handles multi-y layout", {
      d <- sim_test_data()
      y <- cbind(a = rnorm(nrow(d$ref)), b = rnorm(nrow(d$ref)))

      cv <- analog_cv(
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

      perf <- cv_performance(cv)
      expect_equal(sort(unique(perf$variable)), c("a", "b"))
      # 5 continuous metrics per variable
      expect_equal(nrow(perf), 10L)
      expect_true(all(perf$type == "continuous"))
})


test_that("cv_performance auto-detects binary outcomes and returns binary metrics", {
      d <- sim_test_data()
      set.seed(1)
      # Construct binary y that's more predictable than noise so metrics
      # are stable and interpretable.
      y <- as.integer(d$ref[, 3] > median(d$ref[, 3]))

      cv <- analog_cv(
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

      perf <- cv_performance(cv)
      expect_true(all(perf$type == "binary"))
      expect_equal(perf$metric, c("n", "auc", "tss", "tss_threshold", "brier"))

      auc   <- perf$value[perf$metric == "auc"]
      tss   <- perf$value[perf$metric == "tss"]
      brier <- perf$value[perf$metric == "brier"]

      expect_true(auc > 0 && auc <= 1)
      # With a predictive signal, AUC should beat chance
      expect_true(auc > 0.5)
      expect_true(tss >= 0 && tss <= 1)
      expect_true(brier >= 0 && brier <= 1)
})


test_that("cv_performance AUC matches manual Mann-Whitney U on a small example", {
      # Handcraft a small case where AUC has a known value.
      # obs: 1,1,0,0,0
      # pred: 0.9, 0.8, 0.7, 0.3, 0.1
      # Pairs (pos, neg): (0.9,0.7),(0.9,0.3),(0.9,0.1), (0.8,0.7),(0.8,0.3),(0.8,0.1)
      # All 6 pairs are positive > negative, so AUC = 6/6 = 1.
      obs  <- c(1, 1, 0, 0, 0)
      pred <- c(0.9, 0.8, 0.7, 0.3, 0.1)
      w    <- rep(1, length(obs))
      auc1 <- analogs:::.auc_weighted(obs, pred, w)
      expect_equal(auc1, 1.0)

      # Shift one pair to create a tie
      pred2 <- c(0.9, 0.7, 0.7, 0.3, 0.1)  # pos=0.7 ties with one neg=0.7
      # Pairs: (0.9,0.7): win, (0.9,0.3): win, (0.9,0.1): win,
      #        (0.7,0.7): tie=0.5, (0.7,0.3): win, (0.7,0.1): win
      # Sum = 5.5 / 6 = 0.916...
      auc2 <- analogs:::.auc_weighted(obs, pred2, w)
      expect_equal(auc2, 5.5 / 6, tolerance = 1e-12)

      # All-same class → NA
      auc_na <- analogs:::.auc_weighted(c(1, 1, 1), c(0.1, 0.2, 0.3), rep(1, 3))
      expect_true(is.na(auc_na))
})


test_that("cv_performance TSS optimization finds a sane threshold", {
      # Perfectly separable case: positives at pred > 0.5, negatives at pred < 0.5
      obs  <- c(0, 0, 0, 1, 1, 1)
      pred <- c(0.1, 0.2, 0.3, 0.7, 0.8, 0.9)
      w    <- rep(1, length(obs))

      res <- analogs:::.tss_optim_weighted(obs, pred, w)
      expect_equal(res$tss, 1.0)
      # Threshold should be between the highest negative (0.3) and
      # the lowest positive (0.7), i.e. somewhere in (0.3, 0.7].
      expect_true(res$threshold > 0.3 && res$threshold <= 0.7)

      # All-same class → NA
      res_na <- analogs:::.tss_optim_weighted(c(1, 1, 1), c(0.1, 0.2, 0.3), rep(1, 3))
      expect_true(is.na(res_na$tss))
      expect_true(is.na(res_na$threshold))
})


test_that("cv_performance outcome_type override works", {
      d <- sim_test_data()
      y <- as.integer(d$ref[, 3] > median(d$ref[, 3]))

      cv <- analog_cv(
            fun = analog_impact, pool = d$ref, y = y,
            max_clim = 1, max_geog = 2,
            kernel = "gaussian_clim", theta = 0.3,
            coord_type = "projected", cv_method = "loo"
      )

      # Force continuous treatment even though obs is 0/1
      perf <- cv_performance(cv, outcome_type = "continuous")
      expect_true(all(perf$type == "continuous"))
      expect_equal(perf$metric, c("n", "rmse", "mae", "bias", "r2"))

      # Force binary
      perf2 <- cv_performance(cv, outcome_type = "binary")
      expect_true(all(perf2$type == "binary"))
      expect_equal(perf2$metric, c("n", "auc", "tss", "tss_threshold", "brier"))

      # Per-variable named vector (single y case -- named by "y")
      perf3 <- cv_performance(cv, outcome_type = c(y = "binary"))
      expect_true(all(perf3$type == "binary"))
})


test_that("cv_performance handles mixed continuous and binary variables", {
      d <- sim_test_data()
      y <- cbind(
            biomass  = rnorm(nrow(d$ref)),
            presence = as.integer(d$ref[, 3] > median(d$ref[, 3]))
      )

      cv <- analog_cv(
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

      perf <- cv_performance(cv)
      expect_equal(perf$type[perf$variable == "biomass"][1], "continuous")
      expect_equal(perf$type[perf$variable == "presence"][1], "binary")

      cont_metrics <- perf$metric[perf$variable == "biomass"]
      bin_metrics  <- perf$metric[perf$variable == "presence"]
      expect_equal(cont_metrics, c("n", "rmse", "mae", "bias", "r2"))
      expect_equal(bin_metrics,  c("n", "auc", "tss", "tss_threshold", "brier"))
})


test_that("cv_performance rejects invalid outcome_type args", {
      d <- sim_test_data()
      y <- rnorm(nrow(d$ref))
      cv <- analog_cv(
            fun = analog_impact, pool = d$ref, y = y,
            max_clim = 1, max_geog = 2,
            kernel = "gaussian_clim", theta = 0.3,
            coord_type = "projected", cv_method = "loo"
      )

      expect_error(cv_performance(cv, outcome_type = "nonsense"), "outcome_type")
      expect_error(cv_performance(cv, outcome_type = c("binary", "binary")),
                   "named character vector")
      expect_error(cv_performance(cv, outcome_type = c(wrong_name = "binary")),
                   "missing entries")
      expect_error(cv_performance(cv, outcome_type = c(y = "bogus")),
                   "must be one of")
})


test_that("cv_performance accepts weights and computes weighted versions", {
      d <- sim_test_data()
      y <- rnorm(nrow(d$ref))

      cv <- analog_cv(
            fun = analog_impact, pool = d$ref, y = y,
            max_clim = 1, max_geog = 2,
            kernel = "gaussian_clim", theta = 0.3,
            coord_type = "projected", cv_method = "loo"
      )

      # Uniform weights match unweighted
      w <- rep(1, nrow(cv))
      p_unw <- cv_performance(cv)
      p_wgt <- cv_performance(cv, weights = w)
      expect_equal(p_unw$value, p_wgt$value, tolerance = 1e-12)

      # Non-uniform weights: manual check on rmse
      set.seed(1)
      w2 <- runif(nrow(cv))
      p_w2 <- cv_performance(cv, weights = w2)

      ok <- is.finite(cv$obs) & is.finite(cv$residual) & is.finite(w2) & w2 > 0
      o <- cv$obs[ok]; r <- cv$residual[ok]; ww <- w2[ok]
      sum_w <- sum(ww)
      expected_rmse <- sqrt(sum(ww * r^2) / sum_w)
      expect_equal(p_w2$value[p_w2$metric == "rmse"], expected_rmse, tolerance = 1e-12)
})


test_that("cv_performance rejects invalid weights", {
      d <- sim_test_data()
      y <- rnorm(nrow(d$ref))

      cv <- analog_cv(
            fun = analog_impact, pool = d$ref, y = y,
            max_clim = 1, max_geog = 2,
            kernel = "gaussian_clim", theta = 0.3,
            coord_type = "projected", cv_method = "loo"
      )

      expect_error(cv_performance(cv, weights = rep(1, 3)), "length")
      expect_error(cv_performance(cv, weights = rep(-1, nrow(cv))), "non-negative")
})


test_that("cv_performance errors when obs/residual columns are missing", {
      df <- data.frame(x = 1:3, y = 1:3, count = c(5, 6, 7))
      expect_error(cv_performance(df), "obs/residual")
})


test_that("cv_performance works on SpatRaster CV output", {
      skip_if_not_installed("terra")

      set.seed(42)
      r <- terra::rast(matrix(runif(100), 10))
      names(r) <- "var1"
      yr <- terra::setValues(r, rnorm(100))
      names(yr) <- "y"

      cv_rast <- analog_cv(
            fun = analog_impact,
            pool = r,
            y = yr,
            max_clim = 1,
            max_geog = 5,           # permissive to ensure non-NA predictions
            kernel = "gaussian_clim",
            theta = 0.3,
            coord_type = "projected",
            cv_method = "loo"
      )

      expect_s4_class(cv_rast, "SpatRaster")
      expect_true(all(c("obs", "residual") %in% terra::names(cv_rast)))

      perf <- cv_performance(cv_rast)
      expect_s3_class(perf, "data.frame")
      expect_equal(perf$variable[1], "y")
      expect_equal(perf$type[1], "continuous")
      expect_true(perf$value[perf$metric == "n"] > 0)
})
