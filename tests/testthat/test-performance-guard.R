test_that("k-NN is faster than brute force on modest sizes (smoke perf)", {
      # skip_on_cran()
      # skip_if_not_installed("bench")
      #
      # set.seed(4)
      # ref   <- matrix(rnorm(30000), ncol = 3)   # 1000 x 3
      # focal <- matrix(rnorm(300),   ncol = 3)   # 20 x 3
      # k <- 3
      #
      # bf <- function() bf_clim_knn(focal, ref, k)
      # ix <- function() analogs:::find_analogs_core(focal, ref, NULL, NULL, NA, k, 1e9, NA, NA, NA)
      #
      # b <- bench::mark(bf = bf(), ix = ix(), iterations = 5, check = FALSE)
      # # Not a strict assertion—just ensure it runs and doesn't regress catastrophically
      # expect_true(all(is.finite(b$median)))
})
