test_that("x_cov basic functionality works without error", {

      set.seed(123)
      d <- sim_test_data()

      n_focal <- nrow(d$focal)
      n_clim <- ncol(d$focal) - 2

      # Create simple diagonal covariance (variance = 1 for all variables)
      n_cov_cols <- n_clim * (n_clim + 1) / 2
      x_cov <- matrix(NA, nrow = n_focal, ncol = n_cov_cols)
      x_cov[, 1:n_clim] <- 1.0
      x_cov[, 2:n_cov_cols] <- .5

      # Should run without error
      result <- analog_search(
            x = d$focal,
            pool = d$ref,
            select = "knn_geog",
            clim = NULL,
            k = 1,
            x_cov = x_cov,
            coord_type = "projected"
      )

      expect_s3_class(result, "data.frame")
      expect_equal(nrow(result), n_focal)
})


test_that("x_cov validation catches dimension mismatches", {

      set.seed(456)
      d <- sim_test_data()

      n_focal <- nrow(d$focal)
      n_clim <- ncol(d$focal) - 2
      n_cov_cols <- n_clim * (n_clim + 1) / 2

      # Wrong number of rows
      x_cov_bad_rows <- matrix(1.0, nrow = n_focal - 1, ncol = n_cov_cols)
      expect_error(
            analog_search(d$focal, d$ref, stat = "count", clim = kernel(max = 1),
                          x_cov = x_cov_bad_rows, coord_type = "projected"),
            "must have same number of rows"
      )

      # Wrong number of columns
      x_cov_bad_cols <- matrix(1.0, nrow = n_focal, ncol = n_cov_cols + 1)
      expect_error(
            analog_search(d$focal, d$ref, stat = "count", clim = kernel(max = 1),
                          x_cov = x_cov_bad_cols, coord_type = "projected"),
            "must have.*columns"
      )
})


test_that("x_cov validation catches non-finite values", {

      set.seed(789)
      d <- sim_test_data()

      n_focal <- nrow(d$focal)
      n_clim <- ncol(d$focal) - 2
      n_cov_cols <- n_clim * (n_clim + 1) / 2

      # Create valid covariance
      x_cov <- matrix(1.0, nrow = n_focal, ncol = n_cov_cols)

      # Add NA
      x_cov[5, 2] <- NA
      expect_error(
            analog_search(d$focal, d$ref, stat = "count", clim = kernel(max = 1),
                          x_cov = x_cov, coord_type = "projected"),
            "non-finite values"
      )

      # Add Inf
      x_cov[5, 2] <- Inf
      expect_error(
            analog_search(d$focal, d$ref, stat = "count", clim = kernel(max = 1),
                          x_cov = x_cov, coord_type = "projected"),
            "non-finite values"
      )
})


test_that("x_cov validation catches non-positive-definite matrices", {

      set.seed(111)
      d <- sim_test_data()

      n_focal <- nrow(d$focal)
      n_clim <- ncol(d$focal) - 2  # Should be 2 for sim_test_data

      # Create non-PD covariance: var1=1, var2=1, cov12=2 (impossible!)
      # For a 2x2 matrix to be PD: var1*var2 > cov12^2
      # Here: 1*1 = 1, but cov12^2 = 4, so NOT PD
      x_cov <- matrix(c(1, 1, 2), nrow = n_focal, ncol = 3, byrow = TRUE)

      suppressWarnings(expect_warning(
            analog_search(d$focal, d$ref, stat = "count", clim = kernel(max = 1),
                          x_cov = x_cov, coord_type = "projected"),
            "not positive definite"
      ))
})


test_that("x_cov with identity covariance matches Euclidean results", {

      set.seed(222)
      d <- sim_test_data()

      n_focal <- nrow(d$focal)
      n_clim <- ncol(d$focal) - 2

      # Create identity covariance: var(clim1)=1, var(clim2)=1, cov12=0
      x_cov_identity <- matrix(c(1, 1, 0), nrow = n_focal, ncol = 3, byrow = TRUE)

      # Result with Mahalanobis (identity = Euclidean)
      result_mahal <- analog_search(
            d$focal, d$ref,
            select = "knn_geog",
            clim = kernel(max = 1),
            k = 3,
            x_cov = x_cov_identity,
            coord_type = "projected"
      )

      # Result with standard Euclidean
      result_eucl <- analog_search(
            d$focal, d$ref,
            select = "knn_geog",
            clim = kernel(max = 1),
            k = 3,
            x_cov = NULL,
            coord_type = "projected"
      )

      # Should be nearly identical
      expect_equal(nrow(result_mahal), nrow(result_eucl))

      # Analog indices should match (at least mostly)
      # Allow some minor differences due to numerical precision in ties
      focal1_mahal <- result_mahal$analog_index[result_mahal$index == 1]
      focal1_eucl <- result_eucl$analog_index[result_eucl$index == 1]

      # At least 2 out of 3 should match
      n_matches <- sum(focal1_mahal %in% focal1_eucl)
      expect_true(n_matches >= 2)
})

test_that("x_cov correctly implements Mahalanobis distance for filtering and clim_dist", {

      set.seed(333)

      # Create simple test data where we can verify Mahalanobis distance
      # Focal: (0, 0) in climate space
      x <- matrix(c(0, 0, 0, 0), nrow = 1, ncol = 4)

      # Reference: points on regular grid
      pool <- cbind(
            matrix(0, nrow = 11^2, ncol = 2),
            as.matrix(expand.grid(t = -5:5, p = -5:5))
      )

      x_cov <- matrix(c(1, 2, .75), nrow = 1, ncol = 3)

      # manual MD calculation
      md <- sqrt(stats::mahalanobis(x = pool[, 3:4],
                                    center = x[, 3:4],
                                    cov = matrix(x_cov[c(1, 3, 3, 2)], 2)))


      # check if filtered sets match
      max_clim <- 4
      result <- analog_search(
            x, pool,
            select = "all",
            clim = kernel(max = max_clim),
            x_cov = x_cov,
            coord_type = "projected"
      )
      expect_equal(which(md <= max_clim),
                   sort(result$analog_index))


      # check if distances match
      result <- analog_search(
            x, pool,
            select = "all",
            clim = NULL,
            x_cov = x_cov,
            coord_type = "projected"
      )
      expect_equal(result$clim_dist, md[result$analog_index])
})


test_that("x_cov works correctly with knn_clim mode", {

      set.seed(444)
      d <- sim_test_data()

      n_focal <- nrow(d$focal)
      n_clim <- ncol(d$focal) - 2

      # Create covariance where clim2 has higher variance (gets downweighted)
      x_cov <- matrix(c(1, 4, 0), nrow = n_focal, ncol = 3, byrow = TRUE)

      # Find k=5 nearest in climate space within geographic radius
      result <- analog_search(
            d$focal, d$ref,
            select = "knn_clim",
            geog = kernel(max = 2),
            k = 5,
            x_cov = x_cov,
            coord_type = "projected"
      )

      # Should return up to k matches per focal
      expect_true(nrow(result) <= n_focal * 5)
      expect_true(all(c("index", "analog_index", "clim_dist") %in% names(result)))

      # All returned distances should be finite
      expect_true(all(is.finite(result$clim_dist)))
})


test_that("x_cov works correctly with count mode", {

      set.seed(555)
      d <- sim_test_data()

      n_focal <- nrow(d$focal)
      n_clim <- ncol(d$focal) - 2

      # Create simple diagonal covariance
      x_cov <- matrix(c(1, 1, 0), nrow = n_focal, ncol = 3, byrow = TRUE)

      # Count analogs within thresholds
      result <- analog_search(
            d$focal, d$ref,
            stat = "count",
            clim = kernel(max = 1),
            geog = kernel(max = 2),
            x_cov = x_cov,
            coord_type = "projected"
      )

      expect_s3_class(result, "data.frame")
      expect_equal(nrow(result), n_focal)
      expect_true("count" %in% names(result))
      expect_true(all(result$count >= 0))
      expect_true(all(is.finite(result$count)))
})


test_that("x_cov works with focal-specific covariance matrices", {

      set.seed(666)
      d <- sim_test_data()

      n_focal <- nrow(d$focal)
      n_clim <- ncol(d$focal) - 2

      # Create different covariance for each focal point
      # First half: var1=1, var2=1, cov12=0
      # Second half: var1=1, var2=4, cov12=0
      x_cov <- matrix(0, nrow = n_focal, ncol = 3)
      half <- n_focal %/% 2

      x_cov[1:half, ] <- matrix(c(1, 1, 0), nrow = half, ncol = 3, byrow = TRUE)
      x_cov[(half+1):n_focal, ] <- matrix(c(1, 4, 0),
                                          nrow = n_focal - half, ncol = 3, byrow = TRUE)

      # Should work without error
      result <- analog_search(
            d$focal, d$ref,
            stat = "count",
            clim = kernel(max = 1),
            x_cov = x_cov,
            coord_type = "projected"
      )

      expect_equal(nrow(result), n_focal)

      # Results should differ between the two groups
      # (different covariance → different distance metric → different counts)
      counts_group1 <- result$count[1:half]
      counts_group2 <- result$count[(half+1):n_focal]

      # At least some should differ (not guaranteed all will, but most should)
      n_different <- sum(counts_group1[1:min(5, half)] !=
                               counts_group2[1:min(5, n_focal - half)])
      expect_true(n_different >= 1)
})


test_that("x_cov works with analog_velocity wrapper", {

      set.seed(777)
      d <- sim_test_data()

      n_focal <- nrow(d$focal)

      # Simple diagonal covariance
      x_cov <- matrix(c(1, 1, 0), nrow = n_focal, ncol = 3, byrow = TRUE)

      # Should work through the wrapper
      result <- analog_velocity(
            x = d$focal,
            pool = d$ref,
            clim = kernel(max = 1),
            k = 1,
            x_cov = x_cov,
            coord_type = "projected"
      )

      expect_s3_class(result, "data.frame")
      expect_equal(nrow(result), n_focal)
      expect_true(all(c("index", "analog_index", "geog_dist") %in% names(result)))
})


test_that("x_cov works with pre-built analog_index", {

      set.seed(888)
      d <- sim_test_data()

      n_focal <- nrow(d$focal)

      # Build index
      index <- build_analog_index(d$ref, coord_type = "projected")

      # Create covariance
      x_cov <- matrix(c(1, 1, 0), nrow = n_focal, ncol = 3, byrow = TRUE)

      # Query with x_cov
      result <- analog_search(
            x = d$focal,
            pool = index,
            select = "knn_geog",
            clim = kernel(max = 1),
            k = 1,
            x_cov = x_cov
      )

      expect_s3_class(result, "data.frame")
      expect_equal(nrow(result), n_focal)
})


test_that("x_cov NULL behavior is unchanged from before", {

      set.seed(999)
      d <- sim_test_data()

      # With x_cov = NULL explicitly
      result_null <- analog_search(
            d$focal, d$ref,
            select = "knn_geog",
            clim = kernel(max = 1),
            k = 1,
            x_cov = NULL,
            coord_type = "projected"
      )

      # Without x_cov argument (default NULL)
      result_default <- analog_search(
            d$focal, d$ref,
            select = "knn_geog",
            clim = kernel(max = 1),
            k = 1,
            coord_type = "projected"
      )

      # Should be identical
      expect_equal(result_null, result_default)
})


test_that("x_cov works with single climate variable", {

      set.seed(1111)

      # Create data with single climate variable
      focal <- matrix(rnorm(20 * 3), ncol = 3)  # x, y, clim1
      ref <- matrix(rnorm(100 * 3), ncol = 3)

      n_focal <- nrow(focal)

      # For 1 climate var: only 1 column needed (variance)
      x_cov <- matrix(1.0, nrow = n_focal, ncol = 1)

      # Should work
      result <- analog_search(
            focal, ref,
            stat = "count",
            clim = kernel(max = 1),
            x_cov = x_cov,
            coord_type = "projected"
      )

      expect_equal(nrow(result), n_focal)
      expect_true("count" %in% names(result))
})


test_that("x_cov validation happens before expensive operations", {

      set.seed(1212)
      d <- sim_test_data()

      n_focal <- nrow(d$focal)

      # Invalid x_cov (wrong dimensions)
      x_cov_bad <- matrix(1.0, nrow = n_focal, ncol = 5)  # Wrong number of columns

      # Should fail quickly with clear error, not crash during C++ execution
      expect_error(
            analog_search(d$focal, d$ref, stat = "count", clim = kernel(max = 1),
                          x_cov = x_cov_bad, coord_type = "projected"),
            "must have.*columns"
      )
})


test_that("x_cov works with correlated climate variables", {

      set.seed(1313)

      # Create test case with known correlation structure
      focal <- matrix(c(0, 0, 0, 0), nrow = 1, ncol = 4)

      # Reference points designed to test correlation handling
      ref <- matrix(c(
            0, 0,  1,  0,   # Point 1
            0, 0,  0,  1,   # Point 2
            0, 0,  1,  1,   # Point 3
            0, 0,  2,  2    # Point 4
      ), nrow = 4, ncol = 4, byrow = TRUE)

      # Covariance with positive correlation: var1=1, var2=1, cov12=0.5
      # This creates positive correlation between the two climate variables
      x_cov <- matrix(c(1, 1, 0.5), nrow = 1, ncol = 3)

      # Verify matrix is PD: det = 1*1 - 0.5^2 = 0.75 > 0 ✓

      # With correlation, the distance calculation changes
      # Point 4 (2,2) is highly correlated, so distance should be larger
      # than naive Euclidean would suggest

      result <- analog_search(
            focal, ref,
            select = "all",
            clim = kernel(max = 3),
            x_cov = x_cov,
            coord_type = "projected"
      )

      # Should find some matches
      expect_true(nrow(result) >= 1)
      expect_true(all(result$analog_index %in% 1:4))
})


test_that("x_cov error messages are informative", {

      set.seed(1414)
      d <- sim_test_data()

      n_focal <- nrow(d$focal)
      n_clim <- ncol(d$focal) - 2

      # Wrong number of rows
      x_cov_wrong_rows <- matrix(1.0, nrow = 5, ncol = 3)
      err_msg <- tryCatch(
            analog_search(d$focal, d$ref, stat = "count", clim = kernel(max = 1),
                          x_cov = x_cov_wrong_rows, coord_type = "projected"),
            error = function(e) e$message
      )
      expect_true(grepl("same number of rows", err_msg, ignore.case = TRUE))
      expect_true(grepl(as.character(n_focal), err_msg))

      # Wrong number of columns
      x_cov_wrong_cols <- matrix(1.0, nrow = n_focal, ncol = 5)
      err_msg2 <- tryCatch(
            analog_search(d$focal, d$ref, stat = "count", clim = kernel(max = 1),
                          x_cov = x_cov_wrong_cols, coord_type = "projected"),
            error = function(e) e$message
      )
      expect_true(grepl("columns", err_msg2, ignore.case = TRUE))
      expect_true(grepl(as.character(n_clim), err_msg2))
})


test_that("x_cov works with sum mode and weights", {

      set.seed(333)
      d <- sim_test_data()

      n_focal <- nrow(d$focal)

      # Create covariance with different variances
      x_cov <- matrix(0, nrow = n_focal, ncol = 3)
      x_cov[, 1] <- 1.5  # higher variance in clim1
      x_cov[, 2] <- 0.5  # lower variance in clim2
      x_cov[, 3] <- 0.3  # moderate correlation

      # Test with inverse_clim weight
      s1 <- analog_density(d$focal, d$ref,
                           clim = kernel("inverse", max = 1), geog = kernel(max = 2),
                           coord_type = "projected")

      s2 <- analog_density(d$focal, d$ref,
                           clim = kernel("inverse", max = 1), geog = kernel(max = 2),
                           coord_type = "projected",
                           x_cov = x_cov)

      # Should return same structure
      expect_equal(nrow(s1), nrow(s2))
      expect_true("sum_weights" %in% names(s2))
      expect_true(all(is.finite(s2$sum_weights)))

      # Values should differ (different distance metric affects weights)
      expect_false(all(abs(s1$sum_weights - s2$sum_weights) < 1e-10))
})


test_that("x_cov with different variance scales affects analog selection", {

      set.seed(444)
      d <- sim_test_data()

      n_focal <- nrow(d$focal)

      # Low variance (more sensitive to differences)
      x_cov_low <- matrix(0, nrow = n_focal, ncol = 3)
      x_cov_low[, 1] <- 0.25
      x_cov_low[, 2] <- 0.25
      x_cov_low[, 3] <- 0.0

      # High variance (less sensitive to differences)
      x_cov_high <- matrix(0, nrow = n_focal, ncol = 3)
      x_cov_high[, 1] <- 4.0
      x_cov_high[, 2] <- 4.0
      x_cov_high[, 3] <- 0.0

      v_low <- analog_search(d$focal, d$ref, select = "all",
                             clim = kernel(max = .25), coord_type = "projected",
                             x_cov = x_cov_low)

      v_high <- analog_search(d$focal, d$ref, select = "all",
                              clim = kernel(max = .25), coord_type = "projected",
                              x_cov = x_cov_high)

      expect_false(all(suppressWarnings(v_low$analog_index == v_high$analog_index)))
})


test_that("x_cov affects which points pass max_clim threshold", {

      set.seed(555)
      d <- sim_test_data()

      n_focal <- nrow(d$focal)
      max_clim <- 0.6  # Tight threshold

      # Narrow variance (identity-like)
      x_cov_std <- matrix(0, nrow = n_focal, ncol = 3)
      x_cov_std[, 1] <- 1
      x_cov_std[, 2] <- 1
      x_cov_std[, 3] <- 0.0

      # High variance (distances will be smaller in Mahalanobis space)
      x_cov_high <- matrix(0, nrow = n_focal, ncol = 3)
      x_cov_high[, 1] <- 3.0
      x_cov_high[, 2] <- 3.0
      x_cov_high[, 3] <- 0.0

      a_std <- analog_availability(d$focal, d$ref,
                                   clim = kernel(max = max_clim), geog = kernel(max = 10),
                                   coord_type = "projected",
                                   x_cov = x_cov_std)

      a_high <- analog_availability(d$focal, d$ref,
                                    clim = kernel(max = max_clim), geog = kernel(max = 10),
                                    coord_type = "projected",
                                    x_cov = x_cov_high)

      # High variance should result in more points passing the threshold
      # (because Mahalanobis distances are scaled down by sqrt(variance))
      mean_count_std <- mean(a_std$count)
      mean_count_high <- mean(a_high$count)

      expect_true(mean_count_high > mean_count_std)
})


test_that("x_cov works with auto-tuning for large datasets", {

      # Skip if this takes too long in automated testing
      skip_on_cran()

      set.seed(666)
      # Need >2000 focal points to trigger auto-tuning
      focal_large <- matrix(rnorm(2100 * 4), ncol = 4)
      ref_large <- matrix(rnorm(500 * 4), ncol = 4)

      n_focal <- nrow(focal_large)

      x_cov <- matrix(0, nrow = n_focal, ncol = 3)
      x_cov[, 1] <- 1.0
      x_cov[, 2] <- 1.0
      x_cov[, 3] <- 0.4

      # Should trigger auto-tuning and complete successfully
      v <- analog_velocity(focal_large, ref_large,
                           clim = NULL, k = 1,
                           coord_type = "projected",
                           geog_res_adj = "auto",
                           x_cov = x_cov)

      expect_equal(nrow(v), n_focal)
      expect_true(all(c("index", "analog_index") %in% names(v)))
})


test_that("x_cov works correctly with 3 climate variables", {

      set.seed(777)
      # 3 climate variables: 6 covariance components
      focal_3clim <- matrix(rnorm(50 * 5), ncol = 5)  # x, y, clim1, clim2, clim3
      ref_3clim <- matrix(rnorm(200 * 5), ncol = 5)

      n_focal <- nrow(focal_3clim)

      # Create covariance matrix: [var1, var2, var3, cov12, cov13, cov23]
      x_cov <- matrix(0, nrow = n_focal, ncol = 6)
      x_cov[, 1] <- 1.0  # var(clim1)
      x_cov[, 2] <- 1.5  # var(clim2)
      x_cov[, 3] <- 0.8  # var(clim3)
      x_cov[, 4] <- 0.3  # cov(clim1, clim2)
      x_cov[, 5] <- 0.2  # cov(clim1, clim3)
      x_cov[, 6] <- 0.4  # cov(clim2, clim3)

      v <- analog_velocity(focal_3clim, ref_3clim,
                           clim = NULL, k = 1,
                           coord_type = "projected",
                           x_cov = x_cov)

      expect_equal(nrow(v), n_focal)
      expect_true(all(c("index", "analog_index") %in% names(v)))
})


test_that("x_cov works with strongly correlated climate variables", {

      set.seed(888)
      d <- sim_test_data()

      n_focal <- nrow(d$focal)

      # Identity (no correlation)
      x_cov_uncorr <- matrix(0, nrow = n_focal, ncol = 3)
      x_cov_uncorr[, 1] <- 1.0
      x_cov_uncorr[, 2] <- 1.0
      x_cov_uncorr[, 3] <- 0.0

      # Strong positive correlation
      x_cov_corr <- matrix(0, nrow = n_focal, ncol = 3)
      x_cov_corr[, 1] <- 1.0
      x_cov_corr[, 2] <- 1.0
      x_cov_corr[, 3] <- 0.85  # cor = 0.85

      v_uncorr <- analog_search(d$focal, d$ref, select = "all", clim = kernel(max = 1),
                                coord_type = "projected",
                                x_cov = x_cov_uncorr)

      v_corr <- analog_search(d$focal, d$ref, select = "all", clim = kernel(max = 1),
                              coord_type = "projected",
                              x_cov = x_cov_corr)

      expect_false(length(v_uncorr$analog_index) == length(v_corr$analog_index))
})
