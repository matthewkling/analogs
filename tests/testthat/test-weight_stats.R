# Tests for the analog-support stats: sum_weights, ess, effective_weight, and
# weighted_sd, against hand-computed kernel weights.

support_test_data <- function(seed = 1, n = 400) {
      set.seed(seed)
      list(pool = data.frame(x = runif(n, 0, 100), y = runif(n, 0, 100), t = rnorm(n)),
           y = data.frame(a = rnorm(n, 5, 2), b = rbinom(n, 1, 0.3)),
           # the third focal is out of geographic range of every pool point
           focal = data.frame(x = c(50, 10, 500), y = c(50, 90, 500), t = c(0, 1, 0)),
           geog = kernel("gaussian", theta = 10, max = 30),
           env = kernel("gaussian", theta = 1, max = 3))
}

# Kernel weights of every pool point for focal i, computed directly
hand_weights <- function(d, i) {
      dg <- sqrt((d$pool$x - d$focal$x[i])^2 + (d$pool$y - d$focal$y[i])^2)
      de <- abs(d$pool$t - d$focal$t[i])
      w <- exp(-dg^2 / (2 * 10^2)) * exp(-de^2 / 2)
      w[dg > 30 | de > 3] <- 0
      w
}

hand_sd <- function(w, v) {
      keep <- !is.na(v)
      w <- w[keep]
      v <- v[keep]
      m <- sum(w * v) / sum(w)
      sqrt(sum(w * (v - m)^2) / sum(w))
}


test_that("effective_weight and weighted_sd match manual computation", {
      d <- support_test_data()
      r <- analog_impact(d$focal, d$pool, y = d$y, env = d$env, geog = d$geog,
                         stat = c("sum_weights", "ess", "effective_weight",
                                  "weighted_mean", "weighted_sd"),
                         coord_type = "projected")

      expect_identical(names(r), c("index", "x", "y", "sum_weights", "ess",
                                   "effective_weight", "weighted_mean_a", "weighted_sd_a",
                                   "weighted_mean_b", "weighted_sd_b"))
      for (i in 1:2) {
            w <- hand_weights(d, i)
            expect_equal(r$sum_weights[i], sum(w))
            expect_equal(r$effective_weight[i], sum(w^2) / sum(w))
            expect_equal(r$weighted_sd_a[i], hand_sd(w, d$y$a))
            expect_equal(r$weighted_sd_b[i], hand_sd(w, d$y$b))
      }
      # sum_weights = ess x effective_weight
      expect_equal(r$sum_weights[1:2], r$ess[1:2] * r$effective_weight[1:2])
      # For a 0/1 response the weighted SD is determined by the weighted mean
      p <- r$weighted_mean_b[1:2]
      expect_equal(r$weighted_sd_b[1:2], sqrt(p * (1 - p)))
})


test_that("a focal with no analogs gets NA closeness and spread", {
      d <- support_test_data()
      r <- analog_impact(d$focal, d$pool, y = d$y$a, env = d$env, geog = d$geog,
                         stat = c("sum_weights", "effective_weight", "weighted_sd"),
                         coord_type = "projected")
      expect_identical(r$sum_weights[3], 0)
      expect_true(is.na(r$effective_weight[3]))
      expect_true(is.na(r$weighted_sd[3]))
})


test_that("weighted_sd needs neither weighted_mean nor se, and works beside them", {
      d <- support_test_data()
      args <- list(d$focal, d$pool, y = d$y$a, env = d$env, geog = d$geog,
                   coord_type = "projected")
      alone <- do.call(analog_impact, c(args, list(stat = "weighted_sd")))
      both  <- do.call(analog_impact, c(args, list(stat = c("weighted_mean", "weighted_sd"),
                                                   se = "ess")))
      # single y: bare names, SE directly after weighted_mean
      expect_identical(names(both), c("index", "x", "y", "weighted_mean",
                                      "se_weighted_mean", "weighted_sd"))
      expect_equal(alone$weighted_sd, both$weighted_sd)
      # se = "ess" is the weighted SD over the square root of the ESS
      ess <- do.call(analog_impact, c(args, list(stat = "ess")))$ess
      expect_equal(both$se_weighted_mean[1:2], both$weighted_sd[1:2] / sqrt(ess[1:2]))
})


test_that("weighted_sd skips missing responses", {
      d <- support_test_data()
      y <- d$y
      y$a[seq(1, nrow(y), by = 7)] <- NA
      r <- analog_impact(d$focal, d$pool, y = y, env = d$env, geog = d$geog,
                         stat = c("weighted_mean", "weighted_sd"),
                         coord_type = "projected")
      for (i in 1:2) {
            w <- hand_weights(d, i)
            expect_equal(r$weighted_sd_a[i], hand_sd(w, y$a))
            expect_equal(r$weighted_sd_b[i], hand_sd(w, y$b))
      }
})


test_that("with uniform weights the stats reduce to their unweighted forms", {
      d <- support_test_data()
      r <- analog_search(d$focal, d$pool, y = d$y$a, select = "all",
                         geog = kernel(max = 30), env = kernel(max = 3),
                         stat = c("count", "ess", "effective_weight", "weighted_sd"),
                         coord_type = "projected")
      expect_equal(r$ess[1:2], r$count[1:2])
      expect_equal(r$effective_weight[1:2], c(1, 1))
      for (i in 1:2) {
            w <- as.numeric(hand_weights(d, i) > 0)
            expect_equal(r$weighted_sd[i], hand_sd(w, d$y$a))   # population SD
      }
})


test_that("the new stats are validated like their neighbors", {
      d <- support_test_data()
      expect_error(analog_search(d$focal, d$pool, select = "all", env = d$env, geog = d$geog,
                                 stat = "weighted_sd", coord_type = "projected"),
                   "require 'y'")
      expect_error(analog_impact(d$focal, d$pool, y = sample(c("p", "q"), 400, TRUE),
                                 env = d$env, geog = d$geog,
                                 stat = c("tabulate", "weighted_sd"),
                                 coord_type = "projected"),
                   "cannot be combined")
      r <- analog_density(d$focal, d$pool, env = d$env, geog = d$geog,
                          stat = c("sum_weights", "effective_weight"),
                          coord_type = "projected")
      expect_true("effective_weight" %in% names(r))
})


test_that("effective_weight is the mean kernel weight, in [0, 1] whatever the other weights", {
      d <- support_test_data()
      u <- runif(nrow(d$pool), 0.2, 3)
      ew <- function(...) analog_search(d$focal, d$pool, select = "all", env = d$env, geog = d$geog,
                                        stat = "effective_weight", coord_type = "projected", ...)$effective_weight
      plain  <- ew()
      scaled <- ew(weight = rep(10, nrow(d$pool)))
      uneven <- ew(weight = u)
      expect_equal(scaled, plain)
      for (i in 1:2) {
            k <- hand_weights(d, i)
            expect_equal(uneven[i], sum(u * k * k) / sum(u * k))
      }
      down <- ew(downsample = 0.5, seed = 1, geog_res_adj = 1, env_res_adj = 1)
      expect_true(all(c(plain, uneven, down)[is.finite(c(plain, uneven, down))] <= 1))
})
