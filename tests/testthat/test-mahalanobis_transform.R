test_that("mahalanobis_transform works with single matrix input", {

      set.seed(123)
      # Create correlated environmental data
      n <- 100
      x <- matrix(rnorm(n * 2), ncol = 2)
      # Add correlation
      x[, 2] <- x[, 2] + 0.5 * x[, 1]

      # Add coordinates
      coords <- matrix(rnorm(n * 2), ncol = 2)
      data <- cbind(coords, x)

      # Transform
      result <- mahalanobis_transform(data)

      expect_true(is.matrix(result))
      expect_equal(nrow(result), n)
      expect_equal(ncol(result), 4)  # 2 coords + 2 environmental

      # Check that transformed data has near-zero correlation
      transformed_env <- result[, 3:4]
      cor_mat <- cor(transformed_env)
      expect_true(abs(cor_mat[1, 2]) < 0.1)  # Should be uncorrelated

      # Check that variances are near 1
      vars <- apply(transformed_env, 2, var)
      expect_true(all(abs(vars - 1) < 0.1))
})


test_that("mahalanobis_transform works with paired inputs", {

      set.seed(456)
      n_x <- 50
      n_pool <- 100

      # Create correlated data
      x_env <- matrix(rnorm(n_x * 2), ncol = 2)
      x_env[, 2] <- x_env[, 2] + 0.7 * x_env[, 1]
      x_data <- cbind(matrix(rnorm(n_x * 2), ncol = 2), x_env)

      pool_env <- matrix(rnorm(n_pool * 2), ncol = 2)
      pool_env[, 2] <- pool_env[, 2] + 0.7 * pool_env[, 1]
      pool_data <- cbind(matrix(rnorm(n_pool * 2), ncol = 2), pool_env)

      # Transform jointly
      result <- mahalanobis_transform(x_data, pool_data)

      expect_type(result, "list")
      expect_true("x" %in% names(result))
      expect_true("pool" %in% names(result))

      expect_equal(nrow(result$x), n_x)
      expect_equal(nrow(result$pool), n_pool)

      # Combined transformed data should be decorrelated
      combined_transformed <- rbind(result$x[, 3:4], result$pool[, 3:4])
      cor_mat <- cor(combined_transformed)
      expect_true(abs(cor_mat[1, 2]) < 0.1)
})


test_that("mahalanobis_transform respects center and scale parameters", {

      set.seed(789)
      n <- 100

      # Data with different scales
      x_env <- matrix(c(rnorm(n, mean = 10, sd = 2),
                         rnorm(n, mean = 1000, sd = 500)), ncol = 2)
      x_data <- cbind(matrix(rnorm(n * 2), ncol = 2), x_env)

      # With scaling (default)
      result_scaled <- mahalanobis_transform(x_data, center = TRUE, scale = TRUE)

      # Without scaling
      result_unscaled <- mahalanobis_transform(x_data, center = TRUE, scale = FALSE)

      # Results should differ
      expect_false(all(abs(result_scaled[, 3:4] - result_unscaled[, 3:4]) < 1e-10))

      # Both should be centered
      means_scaled <- colMeans(result_scaled[, 3:4])
      means_unscaled <- colMeans(result_unscaled[, 3:4])
      expect_true(all(abs(means_scaled) < 0.1))
      expect_true(all(abs(means_unscaled) < 0.1))
})


test_that("mahalanobis_transform validates input dimensions", {

      set.seed(111)
      x_data <- matrix(rnorm(50 * 4), ncol = 4)
      pool_data <- matrix(rnorm(100 * 5), ncol = 5)  # Different number of columns

      expect_error(
            mahalanobis_transform(x_data, pool_data),
            "same number of environmental variables"
      )
})


test_that("mahalanobis_transform handles data.frame input", {

      set.seed(222)
      n <- 50
      x_env <- matrix(rnorm(n * 2), ncol = 2)
      x_df <- data.frame(
            x = rnorm(n),
            y = rnorm(n),
            temp = x_env[, 1],
            precip = x_env[, 2]
      )

      result <- mahalanobis_transform(x_df)

      expect_s3_class(result, "data.frame")
      expect_equal(nrow(result), n)
      expect_equal(ncol(result), 4)
})


test_that("mahalanobis_transform preserves coordinate columns", {

      set.seed(333)
      n <- 50
      coords <- matrix(c(seq(1, 50), seq(101, 150)), ncol = 2)
      x_env <- matrix(rnorm(n * 2), ncol = 2)
      x_data <- cbind(coords, x_env)

      result <- mahalanobis_transform(x_data)

      # Coordinates should be unchanged (compare values, not attributes)
      expect_equal(as.numeric(result[, 1:2]), as.numeric(coords))
})


test_that("transformation actually enables Mahalanobis distance", {

      set.seed(444)
      n <- 100

      # Create data with known covariance structure
      # Var(X1) = 4, Var(X2) = 1, Cov(X1,X2) = 0
      x1 <- rnorm(n, sd = 2)
      x2 <- rnorm(n, sd = 1)
      x_env <- cbind(x1, x2)
      x_data <- cbind(matrix(rnorm(n * 2), ncol = 2), x_env)

      # Transform
      result <- mahalanobis_transform(x_data, scale = FALSE)  # Covariance-based
      transformed_env <- result[, 3:4]

      # In transformed space, covariance should be identity
      cov_transformed <- cov(transformed_env)
      expect_true(all(abs(cov_transformed - diag(2)) < 0.2))

      # Distance in transformed space = Mahalanobis in original space
      # Pick two points
      p1_orig <- x_env[1, ]
      p2_orig <- x_env[2, ]
      p1_transformed <- transformed_env[1, ]
      p2_transformed <- transformed_env[2, ]

      # Euclidean distance in transformed space
      eucl_transformed <- sqrt(sum((p1_transformed - p2_transformed)^2))

      # Mahalanobis distance in original space
      cov_orig <- cov(x_env)
      diff_orig <- p1_orig - p2_orig
      mahal_orig <- sqrt(t(diff_orig) %*% solve(cov_orig) %*% diff_orig)

      expect_equal(eucl_transformed, as.numeric(mahal_orig), tolerance = 0.1)
})


test_that("mahalanobis_transform detects zero variance", {

      set.seed(555)
      n <- 50

      # Second environmental variable has zero variance
      x_env <- cbind(rnorm(n), rep(5, n))
      x_data <- cbind(matrix(rnorm(n * 2), ncol = 2), x_env)

      expect_error(
            mahalanobis_transform(x_data),
            "zero or near-zero variance"
      )
})
