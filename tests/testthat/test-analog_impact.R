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

test_that("`analog_impact` result matches manual calculation for LON/LAT coords", {

      d <- sim_test_data(lonlat = TRUE)
      max_geog <- 5000  # km, large enough that some analogs will typically pass

      # impact (lon/lat)
      nn <- analog_impact(
            d$focal, d$ref,
            max_geog = max_geog,
            k = 1,
            coord_type = "lonlat"
      )

      # --- manual lon/lat impact calculation ---

      # (1) climate distances (Euclidean in climate space)
      dclim <- xdist(d$focal, d$ref, "clim")

      # (2) great-circle distances (km)
      dgeog <- xdist(d$focal, d$ref, "lonlat")

      # (3) apply geographic filter (those beyond max_geog cannot be analogs)
      dclim[dgeog > max_geog] <- Inf

      # (4) pick minimum climate distance among those passing the geo filter
      nn_idx <- as.vector(apply(dclim, 1, which.min))
      nn_dst <- as.vector(apply(dclim, 1, min))

      expect_equal(nn$analog_index, nn_idx)
      expect_equal(nn$clim_dist, nn_dst, tolerance = 1e-8)
})

test_that("analog_impact works with x/pool parameter names", {

      d <- sim_test_data()

      # Test with raw data
      i <- analog_impact(
            x = d$focal,
            pool = d$ref,
            max_geog = 2,
            k = 3,
            coord_type = "projected",
            index_res = 10
      )

      expect_s3_class(i, "data.frame")
      expect_true(nrow(i) <= nrow(d$focal) * 3)
      expect_true(all(c("focal_index", "analog_index", "clim_dist") %in% names(i)))
})


test_that("analog_impact works with analog_index", {

      d <- sim_test_data()

      # Build index
      index <- build_analog_index(d$ref, coord_type = "projected", index_res = 12)

      # Query with index
      i <- analog_impact(
            x = d$focal,
            pool = index,
            max_geog = 2,
            k = 3
      )

      expect_s3_class(i, "data.frame")
      expect_true(nrow(i) <= nrow(d$focal) * 3)
      expect_true(all(c("focal_index", "analog_index", "clim_dist") %in% names(i)))
})

