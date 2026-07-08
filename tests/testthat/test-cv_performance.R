test_that("cv_performance returns long-format data.frame for continuous single-y", {
      d <- sim_test_data()
      y <- rnorm(nrow(d$ref))

      cv <- analog_cv(
            fun = analog_impact,
            pool = d$ref,
            y = y,
            env = kernel("gaussian", theta = 0.3, max = 1),
            geog = kernel(max = 2),
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
            env = kernel("gaussian", theta = 0.3, max = 1),
            geog = kernel(max = 2),
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
            env = kernel("gaussian", theta = 0.3, max = 1),
            geog = kernel(max = 2),
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
            env = kernel("gaussian", theta = 0.3, max = 1),
            geog = kernel(max = 2),
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
            env = kernel("gaussian", theta = 0.3, max = 1),
            geog = kernel(max = 2),
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
            env = kernel("gaussian", theta = 0.3, max = 1),
            geog = kernel(max = 2),
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
            env = kernel("gaussian", theta = 0.3, max = 1),
            geog = kernel(max = 2),
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
            env = kernel("gaussian", theta = 0.3, max = 1),
            geog = kernel(max = 2),
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
            env = kernel("gaussian", theta = 0.3, max = 1),
            geog = kernel(max = 5),           # permissive to ensure non-NA predictions
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


# categorical ----------------------

# Tests for cv_performance() with categorical CV results.
#
# These tests exercise the categorical path end-to-end (analog_cv with
# stat = "tabulate" -> cv_performance) plus hand-computable unit tests on
# the metric helpers.

# Helper: build a small synthetic categorical CV setup and run analog_cv.
# Returns the CV result data.frame with obs/primary/brier columns and the
# n_<level> vote columns.
.cat_cv_fixture <- function(n = 60, seed = 1L) {
      set.seed(seed)
      x <- runif(n, 0, 10); y <- runif(n, 0, 10)
      env <- matrix(rnorm(n * 2), ncol = 2,
                     dimnames = list(NULL, c("c1", "c2")))
      pool <- cbind(x = x, y = y, env)
      veg <- factor(rep(c("forest", "grass", "shrub"), length.out = n))

      analog_cv(
            fun = analog_impact,
            pool = pool,
            y = veg,
            stat = c("count", "sum_weights", "tabulate"),
            env = kernel("gaussian", theta = 0.5, max = 5),
            geog = kernel(max = 100),
            coord_type = "projected",
            cv_method = "loo"
      )
}


test_that("categorical CV detection triggers categorical metrics", {
      cv <- .cat_cv_fixture()
      perf <- cv_performance(cv)

      # Should be categorical
      expect_true(all(perf$type == "categorical"))
      expect_true("accuracy" %in% perf$metric)
      expect_true("brier" %in% perf$metric)
      expect_true("n_classes" %in% perf$metric)

      # Confusion rows present, K^2 of them
      conf_rows <- perf[startsWith(perf$metric, "confusion["), ]
      n_classes <- perf$value[perf$metric == "n_classes"]
      expect_equal(nrow(conf_rows), n_classes^2)
})


test_that("accuracy and confusion sum to n", {
      cv <- .cat_cv_fixture()
      perf <- cv_performance(cv)

      # Sum of confusion cells equals n (count of valid focals)
      conf_rows <- perf[startsWith(perf$metric, "confusion["), ]
      total_conf <- sum(conf_rows$value)
      n_val <- perf$value[perf$metric == "n"]
      expect_equal(total_conf, n_val)

      # Accuracy = sum of diagonal / total
      # Parse confusion rows to extract obs and pred labels
      m <- regmatches(conf_rows$metric,
                      regexec("^confusion\\[([^|]+)\\|([^]]+)\\]$",
                              conf_rows$metric))
      obs_lbl <- vapply(m, `[`, character(1), 2L)
      pred_lbl <- vapply(m, `[`, character(1), 3L)
      diag_sum <- sum(conf_rows$value[obs_lbl == pred_lbl])
      expected_accuracy <- if (total_conf > 0) diag_sum / total_conf else NA_real_
      observed_accuracy <- perf$value[perf$metric == "accuracy"]
      expect_equal(observed_accuracy, expected_accuracy)
})


test_that("brier is in [0, 2] and is finite when there are valid focals", {
      cv <- .cat_cv_fixture()
      perf <- cv_performance(cv)
      brier_val <- perf$value[perf$metric == "brier"]
      expect_true(is.finite(brier_val))
      expect_gte(brier_val, 0)
      expect_lte(brier_val, 2)
})


test_that("outcome_type errors when overridden on categorical CV", {
      cv <- .cat_cv_fixture()
      expect_error(cv_performance(cv, outcome_type = "binary"),
                   "does not apply to categorical")
      expect_error(cv_performance(cv, outcome_type = "continuous"),
                   "does not apply to categorical")
})


test_that(".cv_metrics_categorical: hand-computable confusion and accuracy", {
      # Construct labels directly with known counts
      obs <-     c("a", "a", "a", "b", "b", "c", "c", "c")
      primary <- c("a", "a", "b", "b", "c", "c", "a", "c")
      brier <-   rep(0.5, 8)  # arbitrary; not testing brier here
      levels <- c("a", "b", "c")

      m <- .cv_metrics_categorical(obs, primary, brier, levels, w = NULL)

      # n is 8, all valid
      expect_equal(m$value[m$metric == "n"], 8)

      # By hand: confusion matrix with rows = obs, cols = pred
      #          a b c
      #     a   2 1 0
      #     b   0 1 1
      #     c   1 0 2
      # accuracy = (2 + 1 + 2) / 8 = 0.625
      expect_equal(m$value[m$metric == "accuracy"], 0.625)

      # Diagonal cells
      expect_equal(m$value[m$metric == "confusion[a|a]"], 2)
      expect_equal(m$value[m$metric == "confusion[b|b]"], 1)
      expect_equal(m$value[m$metric == "confusion[c|c]"], 2)

      # Off-diagonal
      expect_equal(m$value[m$metric == "confusion[a|b]"], 1)
      expect_equal(m$value[m$metric == "confusion[a|c]"], 0)
      expect_equal(m$value[m$metric == "confusion[b|a]"], 0)
      expect_equal(m$value[m$metric == "confusion[b|c]"], 1)
      expect_equal(m$value[m$metric == "confusion[c|a]"], 1)
      expect_equal(m$value[m$metric == "confusion[c|b]"], 0)

      # n_classes
      expect_equal(m$value[m$metric == "n_classes"], 3)
})


test_that(".cv_metrics_categorical: hand-computable brier", {
      # Two focals, both observed class "a", probs given by row-normalized
      # votes. brier_vec is computed upstream so we pass it in directly.
      obs <- c("a", "b")
      primary <- c("a", "b")
      # For obs == "a" with p = (0.7, 0.2, 0.1), brier = (0.7-1)^2 + 0.2^2 + 0.1^2
      #                                                  = 0.09 + 0.04 + 0.01 = 0.14
      # For obs == "b" with p = (0.1, 0.8, 0.1), brier = 0.1^2 + (0.8-1)^2 + 0.1^2
      #                                                  = 0.01 + 0.04 + 0.01 = 0.06
      brier <- c(0.14, 0.06)
      levels <- c("a", "b", "c")

      m <- .cv_metrics_categorical(obs, primary, brier, levels, w = NULL)

      # mean brier
      expect_equal(m$value[m$metric == "brier"], mean(brier))
})


test_that(".cv_metrics_categorical: NA in obs and primary are handled", {
      obs <-     c("a", NA, "b", "c", NA)
      primary <- c("a", "b", NA, "c", NA)
      brier <-   c(0.1, 0.2, 0.3, 0.4, NA)
      levels <- c("a", "b", "c")

      m <- .cv_metrics_categorical(obs, primary, brier, levels, w = NULL)

      # Only rows 1 and 4 have both obs and primary non-NA
      expect_equal(m$value[m$metric == "n"], 2)
      expect_equal(m$value[m$metric == "accuracy"], 1)  # both correct

      # Brier mean uses only rows where obs, primary, and brier are all
      # non-NA (rows 1 and 4): mean(0.1, 0.4) = 0.25
      expect_equal(m$value[m$metric == "brier"], 0.25)
})


test_that(".cv_metrics_categorical: weights affect accuracy and brier", {
      obs     <- c("a", "a", "b", "b")
      primary <- c("a", "b", "a", "b")
      brier   <- c(0.1, 0.2, 0.3, 0.4)
      levels  <- c("a", "b")
      # weights upweight the wrong predictions
      w <- c(1, 4, 4, 1)

      m <- .cv_metrics_categorical(obs, primary, brier, levels, w = w)

      # weighted accuracy: correct mass / total mass
      # row 1: obs=a, pred=a, correct, w=1
      # row 2: obs=a, pred=b, wrong, w=4
      # row 3: obs=b, pred=a, wrong, w=4
      # row 4: obs=b, pred=b, correct, w=1
      # accuracy = (1 + 1) / 10 = 0.2
      expect_equal(m$value[m$metric == "accuracy"], 0.2)

      # weighted brier mean = (1*0.1 + 4*0.2 + 4*0.3 + 1*0.4) / 10 = 2.5/10 = 0.25
      expect_equal(m$value[m$metric == "brier"], 0.25)
})


test_that("multi-y categorical CV produces per-variable metrics", {
      n <- 60
      set.seed(7)
      x <- runif(n, 0, 10); y <- runif(n, 0, 10)
      env <- matrix(rnorm(n * 2), ncol = 2,
                     dimnames = list(NULL, c("c1", "c2")))
      pool <- cbind(x = x, y = y, env)

      yA <- factor(rep(c("p", "q"), length.out = n))
      yB <- factor(rep(c("x", "y", "z"), length.out = n))
      Y <- data.frame(habitat = yA, soil = yB)

      cv <- analog_cv(
            fun = analog_impact,
            pool = pool,
            y = Y,
            stat = "tabulate",
            env = kernel("gaussian", theta = 0.5, max = 5),
            geog = kernel(max = 100),
            coord_type = "projected",
            cv_method = "loo"
      )

      perf <- cv_performance(cv)

      # Both variables present
      expect_setequal(unique(perf$variable), c("habitat", "soil"))

      # habitat has K=2 -> 4 confusion rows; soil has K=3 -> 9 confusion rows
      h_conf <- perf[perf$variable == "habitat" &
                           startsWith(perf$metric, "confusion["), ]
      s_conf <- perf[perf$variable == "soil" &
                           startsWith(perf$metric, "confusion["), ]
      expect_equal(nrow(h_conf), 4)
      expect_equal(nrow(s_conf), 9)

      # n_classes correct
      expect_equal(perf$value[perf$variable == "habitat" &
                                    perf$metric == "n_classes"], 2)
      expect_equal(perf$value[perf$variable == "soil" &
                                    perf$metric == "n_classes"], 3)
})


test_that("categorical detection: empty CV result errors clearly", {
      # data.frame missing all categorical columns
      df <- data.frame(index = 1:5, x = 1:5, y = 1:5,
                       count = c(2, 3, 1, 4, 0))
      # This should fall through to the continuous path's error since
      # there are no obs_/residual_ columns either
      expect_error(cv_performance(df), "No obs/residual columns")
})


test_that("categorical CV from raster input errors with a hint", {
      skip_if_not_installed("terra")
      # Build a SpatRaster whose layer names match the single-y categorical
      # pattern (obs, primary, brier). This triggers categorical detection,
      # at which point cv_performance should error with a hint about
      # rasters dropping the character class labels.
      r <- terra::rast(nrows = 4, ncols = 4,
                       xmin = 0, xmax = 4, ymin = 0, ymax = 4,
                       nlyrs = 3)
      terra::values(r) <- matrix(rnorm(48), ncol = 3)
      names(r) <- c("obs", "primary", "brier")

      expect_error(cv_performance(r),
                   "Categorical CV metrics require observed and predicted")
})
