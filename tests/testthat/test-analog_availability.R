
test_that("`analog_availability` result matches manual calculation", {

      d <- sim_test_data()
      max_clim <- 2
      max_geog <- 2

      # availability
      nn <- analog_availability(d$focal, d$ref, max_clim = max_clim, max_geog = max_geog, coord_type = "projected")

      # manual availability calculation
      dclim <- xdist(d$focal, d$ref, "clim")
      dgeog <- xdist(d$focal, d$ref, "geog")
      avail <- as.vector(rowSums(dclim < max_clim & dgeog < max_geog))

      expect_equal(nn$value, avail)
})
