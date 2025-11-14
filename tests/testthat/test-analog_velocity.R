test_that("`analog_velocity` result matches manual calculation", {

      d <- sim_test_data()
      max_clim <- 1

      # velocity
      nn <- analog_velocity(d$focal, d$ref, max_clim = max_clim, k = 1, coord_type = "projected")

      # manual velocity calculation
      dclim <- xdist(d$focal, d$ref, "clim")
      dgeog <- xdist(d$focal, d$ref, "geog")
      dgeog[dclim > max_clim] <- Inf
      nn_idx <- as.vector(apply(dgeog, 1, which.min))
      nn_dst <- as.vector(apply(dgeog, 1, min))

      expect_equal(nn$analog_index, nn_idx)
      expect_equal(nn$geog_dist, nn_dst)
})
