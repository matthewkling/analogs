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


test_that("analog_availability works with x/pool parameter names", {

      d <- sim_test_data()

      # Test with raw data
      a <- analog_availability(
            x = d$focal,
            pool = d$ref,
            max_clim = 1,
            max_geog = 2,
            coord_type = "projected",
            index_res = 10
      )

      expect_s3_class(a, "data.frame")
      expect_equal(nrow(a), nrow(d$focal))
      expect_true("value" %in% names(a))
      expect_true(all(is.numeric(a$value)))
})


test_that("analog_availability works with analog_index", {

      d <- sim_test_data()

      # Build index
      index <- build_analog_index(d$ref, coord_type = "projected", index_res = 12)

      # Query with index
      a <- analog_availability(
            x = d$focal,
            pool = index,
            max_clim = 1,
            max_geog = 2
      )

      expect_s3_class(a, "data.frame")
      expect_equal(nrow(a), nrow(d$focal))
      expect_true("value" %in% names(a))
      expect_true(all(is.numeric(a$value)))
})
