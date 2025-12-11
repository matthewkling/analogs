test_that("`analog_impact` runs without error", {

      d <- sim_test_data(lonlat = TRUE)

      expect_no_error(
            impact <- analog_impact(
                  x = d$focal,
                  pool = d$ref,
                  values = matrix(rnorm(2*nrow(d$ref)), ncol = 2),
                  max_geog = 1000,
                  max_clim = 1,
                  theta = .25
            )
      )

      expect_true(all(c("count", "sum_weights",
                        "weighted_mean_value_1", "weighted_mean_value_2") %in%
                            names(impact)))

})


