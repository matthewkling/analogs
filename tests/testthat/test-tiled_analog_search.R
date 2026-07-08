test_that("tiled_analog_search requires SpatRaster inputs", {
      # Create non-raster data
      focal_df <- data.frame(x = 1:10, y = 1:10, var1 = rnorm(10))
      ref_df <- data.frame(x = 1:100, y = 1:100, var1 = rnorm(100))

      expect_error(
            tiled_analog_search(focal_df, ref_df, n_tiles = 4, fun = analog_velocity, geog = kernel(max = 100)),
            "x must be a SpatRaster"
      )
})

test_that("tiled_analog_search requires matching CRS", {
      # Create rasters with different CRS
      focal <- terra::rast(ncol = 10, nrow = 10, xmin = 0, xmax = 10, ymin = 0, ymax = 10)
      terra::crs(focal) <- "EPSG:4326"
      terra::values(focal) <- rnorm(100)

      ref <- terra::rast(ncol = 20, nrow = 20, xmin = 0, xmax = 20, ymin = 0, ymax = 20)
      terra::crs(ref) <- "EPSG:3857"
      terra::values(ref) <- rnorm(400)

      expect_error(
            tiled_analog_search(focal, ref, n_tiles = 4, fun = analog_velocity, geog = kernel(max = 5)),
            "x and pool must have the same CRS"
      )
})

test_that("tiled_analog_search produces identical results to non-tiled (projected)", {
      # Create projected test data
      focal <- terra::rast(ncol = 20, nrow = 20, xmin = 0, xmax = 100, ymin = 0, ymax = 100)
      terra::crs(focal) <- "EPSG:32611"  # UTM zone 11N
      names(focal) <- "temp"
      terra::values(focal) <- rnorm(400, mean = 10, sd = 2)

      ref <- terra::rast(ncol = 40, nrow = 40, xmin = -50, xmax = 150, ymin = -50, ymax = 150)
      terra::crs(ref) <- "EPSG:32611"
      names(ref) <- "temp"
      terra::values(ref) <- rnorm(1600, mean = 10, sd = 2)

      # Run both versions
      result_full <- analog_velocity(focal, ref, env = kernel(max = 2), geog = kernel(max = 50),
                                     k = 1, progress = FALSE)
      result_tiled <- suppressWarnings(
            tiled_analog_search(focal, ref, n_tiles = 4, fun = analog_velocity,
                                env = kernel(max = 2), geog = kernel(max = 50), k = 1, progress = FALSE))
      result_full <- terra::subset(result_full, setdiff(names(result_full), "analog_index"))

      # Compare results
      expect_s4_class(result_tiled, "SpatRaster")
      expect_equal(terra::nlyr(result_tiled), terra::nlyr(result_full))
      expect_equal(names(result_tiled), names(result_full))

      # Values should be very close (allowing for numerical precision)
      diff <- terra::values(result_full$geog_dist) - terra::values(result_tiled$geog_dist)
      expect_true(max(abs(diff), na.rm = TRUE) < 1e-6)
})


test_that("tiled_analog_search produces identical results to non-tiled (lonlat)", {
      # Create lonlat test data
      focal <- terra::rast(ncol = 20, nrow = 20, xmin = -120, xmax = -100, ymin = 30, ymax = 50)
      terra::crs(focal) <- "EPSG:4326"
      names(focal) <- "temp"
      terra::values(focal) <- rnorm(400, mean = 10, sd = 2)

      ref <- terra::rast(ncol = 40, nrow = 40, xmin = -125, xmax = -95, ymin = 25, ymax = 55)
      terra::crs(ref) <- "EPSG:4326"
      names(ref) <- "temp"
      terra::values(ref) <- rnorm(1600, mean = 10, sd = 2)

      # Run both versions
      result_full <- analog_velocity(focal, ref, env = kernel(max = 2), geog = kernel(max = 200),
                                     k = 1, progress = FALSE)
      result_tiled <- tiled_analog_search(focal, ref, n_tiles = 4, fun = analog_velocity,
                                          env = kernel(max = 2), geog = kernel(max = 200), k = 1, progress = FALSE)
      result_full <- terra::subset(result_full, setdiff(names(result_full), "analog_index"))

      # Compare results
      diff <- terra::values(result_full$geog_dist) - terra::values(result_tiled$geog_dist)
      expect_true(max(abs(diff), na.rm = TRUE) < 1e-6)
})

test_that("tiled_analog_search works with different tile counts", {
      focal <- terra::rast(ncol = 20, nrow = 20, xmin = 0, xmax = 100, ymin = 0, ymax = 100)
      terra::crs(focal) <- "EPSG:32611"
      names(focal) <- "temp"
      terra::values(focal) <- rnorm(400, mean = 10, sd = 2)

      ref <- terra::rast(ncol = 40, nrow = 40, xmin = -50, xmax = 150, ymin = -50, ymax = 150)
      terra::crs(ref) <- "EPSG:32611"
      names(ref) <- "temp"
      terra::values(ref) <- rnorm(1600, mean = 10, sd = 2)

      result_full <- analog_velocity(focal, ref, env = kernel(max = 2), geog = kernel(max = 50),
                                     k = 1, progress = FALSE)
      result_full <- terra::subset(result_full, setdiff(names(result_full), "analog_index"))

      # Test different tile counts
      for (n_tiles in c(4, 9, 16, 25)) {
            result_tiled <- suppressWarnings(
                  tiled_analog_search(focal, ref, n_tiles = n_tiles, fun = analog_velocity,
                                      env = kernel(max = 2), geog = kernel(max = 50), k = 1, progress = FALSE))
            diff <- terra::values(result_full$geog_dist) - terra::values(result_tiled$geog_dist)
            expect_true(max(abs(diff), na.rm = TRUE) < 1e-6,
                        info = paste("Failed for n_tiles =", n_tiles))
      }
})


test_that("tiled_analog_search works with disk-based output", {
      focal <- terra::rast(ncol = 20, nrow = 20, xmin = 0, xmax = 100, ymin = 0, ymax = 100)
      terra::crs(focal) <- "EPSG:32611"
      names(focal) <- "temp"
      terra::values(focal) <- rnorm(400, mean = 10, sd = 2)

      ref <- terra::rast(ncol = 40, nrow = 40, xmin = -50, xmax = 150, ymin = -50, ymax = 150)
      terra::crs(ref) <- "EPSG:32611"
      names(ref) <- "temp"
      terra::values(ref) <- rnorm(1600, mean = 10, sd = 2)

      output_file <- tempfile(fileext = ".tif")

      result_disk <- suppressWarnings(
            tiled_analog_search(focal, ref, n_tiles = 4, fun = analog_velocity,
                                env = kernel(max = 2), geog = kernel(max = 50), k = 1,
                                output_file = output_file, progress = FALSE))
      result_mem <- suppressWarnings(
            tiled_analog_search(focal, ref, n_tiles = 4, fun = analog_velocity,
                                env = kernel(max = 2), geog = kernel(max = 50), k = 1, progress = FALSE))

      # Check file exists
      expect_true(file.exists(output_file))

      # Check results match
      diff <- terra::values(result_disk$geog_dist) - terra::values(result_mem$geog_dist)
      expect_true(max(abs(diff), na.rm = TRUE) < 1e-3)

      # Clean up
      unlink(output_file)
})

test_that("tiled_analog_search works with non-square kernels", {
      # Create a rectangular kernel (2:1 aspect ratio)
      focal <- terra::rast(ncol = 40, nrow = 20, xmin = 0, xmax = 200, ymin = 0, ymax = 100)
      terra::crs(focal) <- "EPSG:32611"
      names(focal) <- "temp"
      terra::values(focal) <- rnorm(800, mean = 10, sd = 2)

      ref <- terra::rast(ncol = 60, nrow = 40, xmin = -50, xmax = 250, ymin = -50, ymax = 150)
      terra::crs(ref) <- "EPSG:32611"
      names(ref) <- "temp"
      terra::values(ref) <- rnorm(2400, mean = 10, sd = 2)

      result_full <- analog_velocity(focal, ref, env = kernel(max = 2), geog = kernel(max = 50),
                                     k = 1, progress = FALSE)
      result_tiled <- tiled_analog_search(focal, ref, n_tiles = 16, fun = analog_velocity,
                                          env = kernel(max = 2), geog = kernel(max = 50), k = 1, progress = FALSE)

      # Compare results
      diff <- terra::values(result_full$geog_dist) - terra::values(result_tiled$geog_dist)
      expect_true(max(abs(diff), na.rm = TRUE) < 1e-6)
})

test_that("tiled_analog_search works with different analog_* functions", {
      focal <- terra::rast(ncol = 20, nrow = 20, xmin = 0, xmax = 100, ymin = 0, ymax = 100)
      terra::crs(focal) <- "EPSG:32611"
      names(focal) <- "temp"
      terra::values(focal) <- rnorm(400, mean = 10, sd = 2)

      ref <- terra::rast(ncol = 40, nrow = 40, xmin = -50, xmax = 150, ymin = -50, ymax = 150)
      terra::crs(ref) <- "EPSG:32611"
      names(ref) <- "temp"
      terra::values(ref) <- rnorm(1600, mean = 10, sd = 2)

      # Test analog_availability
      result_avail <- suppressWarnings(
            tiled_analog_search(focal, ref, n_tiles = 4, fun = analog_availability,
                                env = kernel(max = 2), geog = kernel(max = 50), progress = FALSE))
      expect_s4_class(result_avail, "SpatRaster")
      expect_true("count" %in% names(result_avail))
})

test_that("tiled_analog_search warns for ineffective tiling", {
      focal <- terra::rast(ncol = 10, nrow = 10, xmin = 0, xmax = 10, ymin = 0, ymax = 10)
      terra::crs(focal) <- "EPSG:32611"
      names(focal) <- "temp"
      terra::values(focal) <- rnorm(100, mean = 10, sd = 2)

      # Small reference kernel with large max_geog
      ref <- terra::rast(ncol = 15, nrow = 15, xmin = -5, xmax = 15, ymin = -5, ymax = 15)
      terra::crs(ref) <- "EPSG:32611"
      names(ref) <- "temp"
      terra::values(ref) <- rnorm(225, mean = 10, sd = 2)

      expect_warning(
            tiled_analog_search(focal, ref, n_tiles = 4, fun = analog_velocity,
                                env = kernel(max = 2), geog = kernel(max = 10), k = 1, progress = FALSE),
            "max_geog is large relative to reference kernel"
      )
})

test_that("tiled_analog_search works with x_cov and y", {

      focal <- terra::rast(ncol = 20, nrow = 20, xmin = 0, xmax = 100, ymin = 0, ymax = 100)
      terra::crs(focal) <- "EPSG:32611"  # UTM zone 11N
      x_cov <- c(focal, focal, focal)
      focal <- c(focal, focal)
      names(focal) <- c("temp", "prec")
      terra::values(focal) <- rnorm(800, mean = 10, sd = 2)
      terra::values(x_cov[[1:2]]) <- runif(800, 1, 2)
      terra::values(x_cov[[3]]) <- runif(400, 0, 1)

      ref <- terra::rast(ncol = 40, nrow = 40, xmin = -50, xmax = 150, ymin = -50, ymax = 150)
      terra::crs(ref) <- "EPSG:32611"
      vals <- c(ref, ref)
      ref <- c(ref, ref)
      names(ref) <- c("temp", "prec")
      names(vals) <- c("var1", "var2")
      terra::values(ref) <- rnorm(3200, mean = 10, sd = 2)
      terra::values(vals) <- rnorm(3200, mean = 10, sd = 2)

      expect_no_error(result <- suppressWarnings(
            tiled_analog_search(focal, ref, n_tiles = 4, fun = analog_impact,
                                y = vals, x_cov = x_cov,
                                env = kernel(max = 1), geog = kernel(max = 50),
                                progress = FALSE)))
      expect_equal(sort(names(result)),
                   sort(c("count", "sum_weights", "weighted_mean_var1", "weighted_mean_var2")))
})

test_that("optimize_tile_grid creates square-ish tiles", {
      # Test on square kernel
      grid <- optimize_tile_grid(n_tiles = 16, kernel_aspect = 1)
      expect_equal(grid$n_x, 4)
      expect_equal(grid$n_y, 4)

      # Test on 2:1 rectangular kernel
      grid <- optimize_tile_grid(n_tiles = 16, kernel_aspect = 2)
      expect_true(grid$n_x > grid$n_y)  # More tiles in x direction

      # Test that n_tiles close to target
      expect_true(grid$n_x * grid$n_y >= 16 * 0.8)
      expect_true(grid$n_x * grid$n_y <= 16 * 1.2)
})
