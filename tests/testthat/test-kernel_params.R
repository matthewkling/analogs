test_that("Gaussian theta achieves target expected weight at random MVN cell", {
      for (d in c(1, 2, 4, 6)) {
            for (fraction in c(0.01, 0.05, 0.1, 0.3)) {
                  theta <- kernel_params(fraction = fraction, d = d,
                                         kernel = "gaussian")$theta
                  chi_density <- function(r) {
                        (r^(d - 1) * exp(-r^2 / 2)) /
                              (2^(d / 2 - 1) * gamma(d / 2))
                  }
                  expected_weight <- integrate(
                        function(r) chi_density(r) * exp(-r^2 / (2 * theta^2)),
                        lower = 0, upper = Inf
                  )$value
                  expect_equal(expected_weight, fraction, tolerance = 1e-4,
                               info = sprintf("d=%d, fraction=%g", d, fraction))
            }
      }
})

test_that("Uniform theta achieves target expected weight at random MVN cell", {
      for (d in c(1, 2, 4, 6)) {
            for (fraction in c(0.01, 0.05, 0.1, 0.3)) {
                  theta <- kernel_params(fraction = fraction, d = d,
                                         kernel = "uniform")$max
                  # Expected weight under uniform kernel = P(chi(d) < theta)
                  expected_weight <- pchisq(theta^2, df = d)
                  expect_equal(expected_weight, fraction, tolerance = 1e-6,
                               info = sprintf("d=%d, fraction=%g", d, fraction))
            }
      }
})

test_that("Inverse-distance theta achieves target expected weight at random MVN cell", {
      for (d in c(2, 4, 6)) {
            for (fraction in c(0.01, 0.05, 0.1, 0.3)) {
                  theta <- kernel_params(fraction = fraction, d = d,
                                         kernel = "inverse_distance")$theta
                  chi_density <- function(r) {
                        (r^(d - 1) * exp(-r^2 / 2)) /
                              (2^(d / 2 - 1) * gamma(d / 2))
                  }
                  expected_weight <- integrate(
                        function(r) chi_density(r) / (r + theta),
                        lower = 0, upper = Inf
                  )$value
                  expect_equal(expected_weight, fraction, tolerance = 1e-4,
                               info = sprintf("d=%d, fraction=%g", d, fraction))
            }
      }
})

test_that("Expected weight is equal across kernel shapes at matched fraction", {
      # The key PDF-equivalence property: switching kernels at fixed
      # fraction preserves total expected weight delivered to MVN cells.
      for (d in c(2, 4)) {
            for (fraction in c(0.05, 0.2)) {
                  th_gauss <- kernel_params(fraction = fraction, d = d,
                                            kernel = "gaussian")$theta
                  th_unif  <- kernel_params(fraction = fraction, d = d,
                                            kernel = "uniform")$max
                  th_inv   <- kernel_params(fraction = fraction, d = d,
                                            kernel = "inverse_distance")$theta

                  chi_density <- function(r) {
                        (r^(d - 1) * exp(-r^2 / 2)) /
                              (2^(d / 2 - 1) * gamma(d / 2))
                  }

                  ew_gauss <- integrate(
                        function(r) chi_density(r) *
                              exp(-r^2 / (2 * th_gauss^2)),
                        0, Inf)$value
                  ew_unif <- pchisq(th_unif^2, df = d)
                  ew_inv <- integrate(
                        function(r) chi_density(r) / (r + th_inv),
                        0, Inf)$value

                  expect_equal(ew_gauss, fraction, tolerance = 1e-4)
                  expect_equal(ew_unif, fraction, tolerance = 1e-6)
                  expect_equal(ew_inv, fraction, tolerance = 1e-4)
                  expect_equal(ew_gauss, ew_unif, tolerance = 1e-4)
                  expect_equal(ew_unif, ew_inv, tolerance = 1e-4)
            }
      }
})

test_that("Gaussian max truncates at target weight loss under MVN data", {
      for (d in c(1, 2, 4)) {
            for (theta in c(0.3, 1, 3)) {
                  for (loss in c(0.01, 0.001)) {
                        max_val <- kernel_params(theta = theta, d = d,
                                                 loss = loss,
                                                 kernel = "gaussian")$max
                        chi_density <- function(r) {
                              (r^(d - 1) * exp(-r^2 / 2)) /
                                    (2^(d / 2 - 1) * gamma(d / 2))
                        }
                        weight_beyond <- integrate(
                              function(r) chi_density(r) *
                                    exp(-r^2 / (2 * theta^2)),
                              lower = max_val, upper = Inf
                        )$value
                        weight_total <- integrate(
                              function(r) chi_density(r) *
                                    exp(-r^2 / (2 * theta^2)),
                              lower = 0, upper = Inf
                        )$value
                        actual_loss <- weight_beyond / weight_total
                        expect_equal(actual_loss, loss, tolerance = 1e-4,
                                     info = sprintf("d=%d, theta=%g, loss=%g",
                                                    d, theta, loss))
                  }
            }
      }
})

test_that("Gaussian max truncates at target weight loss under uniform data", {
      # Under uniform data, cell density at radius r is proportional to r^(d-1)
      # without the MVN exp(-r^2/2) factor.
      for (d in c(1, 2, 4)) {
            for (theta in c(0.5, 2, 10)) {
                  for (loss in c(0.01, 0.001)) {
                        max_val <- kernel_params(theta = theta, d = d,
                                                 loss = loss,
                                                 kernel = "gaussian",
                                                 data_dist = "uniform")$max
                        # For uniform data: aggregate weight at r ∝ r^(d-1) * kernel(r)
                        weight_beyond <- integrate(
                              function(r) r^(d - 1) *
                                    exp(-r^2 / (2 * theta^2)),
                              lower = max_val, upper = Inf
                        )$value
                        weight_total <- integrate(
                              function(r) r^(d - 1) *
                                    exp(-r^2 / (2 * theta^2)),
                              lower = 0, upper = Inf
                        )$value
                        actual_loss <- weight_beyond / weight_total
                        expect_equal(actual_loss, loss, tolerance = 1e-4,
                                     info = sprintf("d=%d, theta=%g, loss=%g",
                                                    d, theta, loss))
                  }
            }
      }
})

test_that("Inverse-distance max truncates at target weight loss under MVN data", {
      for (d in c(2, 4)) {
            for (theta in c(0.1, 1, 5)) {
                  for (loss in c(0.01, 0.001)) {
                        max_val <- kernel_params(theta = theta, d = d,
                                                 loss = loss,
                                                 kernel = "inverse_distance")$max
                        chi_density <- function(r) {
                              (r^(d - 1) * exp(-r^2 / 2)) /
                                    (2^(d / 2 - 1) * gamma(d / 2))
                        }
                        weight_beyond <- integrate(
                              function(r) chi_density(r) / (r + theta),
                              lower = max_val, upper = Inf
                        )$value
                        weight_total <- integrate(
                              function(r) chi_density(r) / (r + theta),
                              lower = 0, upper = Inf
                        )$value
                        actual_loss <- weight_beyond / weight_total
                        expect_equal(actual_loss, loss, tolerance = 1e-4,
                                     info = sprintf("d=%d, theta=%g, loss=%g",
                                                    d, theta, loss))
                  }
            }
      }
})

test_that("Uniform kernel returns only `max` (matching package interface)", {
      for (d in c(1, 2, 4)) {
            for (fraction in c(0.05, 0.2)) {
                  result <- kernel_params(fraction = fraction, d = d,
                                          kernel = "uniform")
                  expect_named(result, "max")
                  expect_equal(result$max, sqrt(qchisq(fraction, df = d)))
            }

            # `loss` argument is ignored for uniform kernels (single parameter)
            result_with_loss <- kernel_params(fraction = 0.05, d = d,
                                              loss = 0.01, kernel = "uniform")
            expect_named(result_with_loss, "max")
      }
})
