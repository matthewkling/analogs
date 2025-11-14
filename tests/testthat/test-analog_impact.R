test_that("`analog_impact` result matches manual calculation", {

      d <- sim_test_data()
      max_geog <- 1

      # impact
      nn <- analog_impact(d$focal, d$ref, max_geog = max_geog, k = 1, coord_type = "projected")

      # manual impact calculation
      dclim <- xdist(d$focal, d$ref, "clim")
      dgeog <- xdist(d$focal, d$ref, "geog")
      dclim[dgeog > max_geog] <- Inf
      nn_idx <- as.vector(apply(dclim, 1, which.min))
      nn_dst <- as.vector(apply(dclim, 1, min))

      expect_equal(nn$analog_index, nn_idx)
      expect_equal(nn$clim_dist, nn_dst)
})

