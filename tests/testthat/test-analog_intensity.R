
test_that("`analog_intensity` result matches manual calculation", {

      d <- sim_test_data()
      max_clim <- 2
      max_geog <- 2

      # intensity
      nn <- analog_intensity(d$focal, d$ref, max_clim = max_clim, max_geog = max_geog,
                             weight = "inverse_clim", coord_type = "projected")

      # manual intensity calculation
      dclim <- xdist(d$focal, d$ref, "clim")
      dgeog <- xdist(d$focal, d$ref, "geog")
      dclim[dclim > max_clim] <- Inf
      dclim[dgeog > max_geog] <- Inf
      intens <- as.vector(rowSums(1 / dclim))

      expect_equal(nn$value, intens)
})
