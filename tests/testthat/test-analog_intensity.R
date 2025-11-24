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

test_that("analog_intensity works with x/pool parameter names", {

      d <- sim_test_data()

      # Test with raw data
      i <- analog_intensity(
            x = d$focal,
            pool = d$ref,
            max_clim = 1,
            max_geog = 2,
            weight = "inverse_clim",
            coord_type = "projected",
            index_res = 10
      )

      expect_s3_class(i, "data.frame")
      expect_equal(nrow(i), nrow(d$focal))
      expect_true("value" %in% names(i))
      expect_true(all(is.numeric(i$value)))
})


test_that("analog_intensity works with analog_index", {

      d <- sim_test_data()

      # Build index
      index <- build_analog_index(d$ref, coord_type = "projected", index_res = 12)

      # Query with index
      i <- analog_intensity(
            x = d$focal,
            pool = index,
            max_clim = 1,
            max_geog = 2,
            weight = "inverse_clim"
      )

      expect_s3_class(i, "data.frame")
      expect_equal(nrow(i), nrow(d$focal))
      expect_true("value" %in% names(i))
      expect_true(all(is.numeric(i$value)))
})

