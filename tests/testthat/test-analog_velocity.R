test_that("`analog_velocity` result matches manual calculation for planar coords", {

      d <- sim_test_data(lonlat = TRUE)
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


test_that("`analog_velocity` result matches manual calculation for LON/LAT coords", {

      d <- sim_test_data(lonlat = TRUE)
      max_clim <- 1

      # velocity (ECEF chord-space internal; output arc-length km)
      nn <- analog_velocity(d$focal, d$ref,
                            max_clim = max_clim,
                            k = 1,
                            coord_type = "lonlat")

      # --- manual lon/lat nearest-analog calculation ---

      # (1) compute climatic distance matrix (Euclidean)
      dclim <- xdist(d$focal, d$ref, "clim")

      # (2) compute true great-circle geog distances (km)
      dgeog <- xdist(d$focal, d$ref, "lonlat")

      # (3) apply climatic filter
      dgeog[dclim > max_clim] <- Inf

      # (4) find nearest analog index + distance
      nn_idx <- apply(dgeog, 1, which.min)
      nn_dst <- apply(dgeog, 1, min)

      # ---- expectations ----
      expect_equal(nn$analog_index, nn_idx)
      expect_equal(nn$geog_dist, nn_dst, tolerance = 1e-8)
})


test_that("analog_velocity works with analog_index", {

      d <- sim_test_data()

      # Build index
      index <- build_analog_index(d$ref, coord_type = "projected", index_res = 12)

      # Query with index
      v <- analog_velocity(
            x = d$focal,
            pool = index,
            max_clim = 1,
            k = 1
      )

      expect_s3_class(v, "data.frame")
      expect_equal(nrow(v), nrow(d$focal))
      expect_true(all(c("index", "analog_index", "geog_dist") %in% names(v)))
})
