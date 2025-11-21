
test_that("`analog_availability` result matches manual calculation for planar coords", {

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


test_that("`analog_availability` result matches manual calculation for LON/LAT coords", {

      d <- sim_test_data(lonlat = TRUE)

      max_clim <- 2
      max_geog <- 5000  # km, large enough that some analogs will typically pass

      # availability using lon/lat mode (AggWorker path)
      nn <- analog_availability(
            d$focal, d$ref,
            max_clim = max_clim,
            max_geog = max_geog,
            coord_type = "lonlat"
      )

      # manual availability calculation (clim + Haversine geog)
      dclim <- xdist(d$focal, d$ref, "clim")
      dgeog <- xdist(d$focal, d$ref, "lonlat")
      avail <- as.vector(rowSums(dclim < max_clim & dgeog < max_geog))

      expect_equal(nn$value, avail)
})
