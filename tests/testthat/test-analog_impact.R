test_that("`analog_impact` runs without error", {

      d <- sim_test_data(lonlat = TRUE)

      expect_no_error(
            impact <- analog_impact(
                  x = d$focal,
                  pool = d$ref,
                  y = matrix(rnorm(2*nrow(d$ref)), ncol = 2),
                  max_geog = 1000,
                  max_clim = 1,
                  theta = .25
            )
      )

      expect_true(all(c("count", "sum_weights",
                        "weighted_mean_y1", "weighted_mean_y2") %in%
                            names(impact)))

})


test_that("analog_impact passes covariates and lambda through, in regression mode", {

      d <- sim_test_data()
      n_ref <- nrow(d$ref)

      y <- rnorm(n_ref)
      covariates <- matrix(rnorm(n_ref * 2), ncol = 2,
                           dimnames = list(NULL, c("north", "east")))

      result <- analog_impact(
            x = d$focal,
            pool = d$ref,
            y = y,
            covariates = covariates,
            stat = c("count", "weighted_mean", "regression"),
            max_geog = 2,
            max_clim = 2,
            kernel = "gaussian_clim",
            theta = 0.5,
            lambda = 0.1,
            coord_type = "projected",
            index_res = 10
      )

      expect_true(all(c("count", "weighted_mean",
                        "coef_intercept", "coef_north", "coef_east") %in% names(result)))
})
