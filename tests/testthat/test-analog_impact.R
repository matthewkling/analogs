test_that("`analog_impact` runs without error", {

      d <- sim_test_data(lonlat = TRUE)

      expect_no_error(
            impact <- analog_impact(
                  x = d$focal,
                  pool = d$ref,
                  values = rnorm(nrow(d$ref)),
                  max_geog = 1000,   # 1000 km dispersal range
                  max_clim = 1,     # Climate analog threshold
                  theta = .25
            )
      )

      expect_true(all(c("count", "sum_weights", "weighted_mean_value_1") %in%
                            names(impact)))

})


