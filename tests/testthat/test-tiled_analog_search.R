test_that("tiled_analog_search rejects inputs that are neither rasters nor point tables", {
      ref_df <- data.frame(x = 1:100, y = 1:100, var1 = rnorm(100))
      expect_error(
            tiled_analog_search(list(1), ref_df, n_tiles = 4, fun = analog_velocity, geog = kernel(max = 100)),
            "x must be a SpatRaster, matrix, or data.frame"
      )
      expect_error(
            tiled_analog_search(ref_df, build_analog_index(ref_df, coord_type = "projected"),
                                n_tiles = 4, fun = analog_velocity, geog = kernel(max = 100)),
            "pool must be a SpatRaster, matrix, or data.frame"
      )
})


# --- point-table inputs -------------------------------------------------------

# Drop per-call index diagnostics (n_pool, bins_per_axis, ...), which describe
# one tile's index rather than the result
strip_meta <- function(d) {
      attributes(d) <- attributes(d)[c("names", "row.names", "class")]
      d
}

point_fixture <- function(n_pool = 3000, n_focal = 800) {
      set.seed(11)
      pool <- data.frame(x = runif(n_pool, 0, 1000), y = runif(n_pool, 0, 600))
      pool$e1 <- pool$x / 300 + rnorm(n_pool, sd = 0.3)
      pool$e2 <- pool$y / 200 + rnorm(n_pool, sd = 0.3)
      focal <- data.frame(x = runif(n_focal, 0, 1000), y = runif(n_focal, 0, 600))
      focal$e1 <- focal$x / 300 + 0.2
      focal$e2 <- focal$y / 200 + 0.2
      # Missing coordinates and missing environment, both of which must come
      # back as NA rows in their original positions
      focal[5, c("x", "y")] <- NA
      focal[9, "e1"] <- NA
      y <- data.frame(a = pool$e1 + rnorm(n_pool), b = rbinom(n_pool, 1, 0.3))
      list(pool = pool, focal = focal, y = y,
           args = list(env = kernel("gaussian", theta = 0.3, max = 0.9),
                       geog = kernel("gaussian", theta = 40, max = 120),
                       stat = c("weighted_mean", "sum_weights")))
}

test_that("point x and point pool match the non-tiled result exactly", {
      d <- point_fixture()
      untiled <- do.call(analog_impact, c(list(x = d$focal, pool = d$pool, y = d$y,
                                               coord_type = "projected"), d$args))
      tiled <- do.call(tiled_analog_search,
                       c(list(x = d$focal, pool = d$pool, y = d$y, n_tiles = 6,
                              fun = analog_impact, coord_type = "projected",
                              progress = FALSE), d$args))
      expect_equal(strip_meta(tiled), strip_meta(untiled))
      expect_true(all(is.na(tiled[c(5, 9), "weighted_mean_a"])))
})

test_that("point x with point pool forwards and subsets weight", {
      d <- point_fixture()
      w <- runif(nrow(d$pool))
      untiled <- do.call(analog_impact, c(list(x = d$focal, pool = d$pool, y = d$y, weight = w,
                                               coord_type = "projected"), d$args))
      tiled <- do.call(tiled_analog_search,
                       c(list(x = d$focal, pool = d$pool, y = d$y, weight = w, n_tiles = 6,
                              fun = analog_impact, coord_type = "projected",
                              progress = FALSE), d$args))
      expect_equal(strip_meta(tiled), strip_meta(untiled))
})

test_that("point x returns pair results in x's row order with x's indices", {
      d <- point_fixture(n_focal = 60)
      args <- list(x = d$focal, pool = d$pool, coord_type = "projected",
                   env = kernel(max = 0.5), geog = kernel(max = 60), select = "all", stat = NULL)
      untiled <- do.call(analog_search, args)
      tiled <- do.call(tiled_analog_search,
                       c(args, list(n_tiles = 6, fun = analog_search, progress = FALSE)))
      key <- function(r) r[order(r$index, r$analog_index), c("index", "analog_index", "geog_dist")]
      expect_false(is.unsorted(tiled$index, na.rm = TRUE))
      expect_equal(key(tiled), key(untiled), ignore_attr = TRUE)
})

test_that("raster x with a lonlat point pool matches the non-tiled result", {
      set.seed(12)
      x <- terra::rast(nrows = 30, ncols = 40, xmin = -110, xmax = -106, ymin = 38, ymax = 41,
                       crs = "EPSG:4326")
      xy <- terra::xyFromCell(x, seq_len(terra::ncell(x)))
      terra::values(x) <- xy[, 1] / 2 + xy[, 2] / 3
      names(x) <- "e"
      pool <- data.frame(x = runif(4000, -111, -105), y = runif(4000, 37, 42))
      pool$e <- pool$x / 2 + pool$y / 3 + rnorm(4000, sd = 0.2)
      y <- pool$e + rnorm(4000)
      args <- list(env = kernel(max = 0.5), geog = kernel(max = 60), stat = "weighted_mean")
      untiled <- do.call(analog_impact, c(list(x = x, pool = pool, y = y), args))
      tiled <- do.call(tiled_analog_search,
                       c(list(x = x, pool = pool, y = y, n_tiles = 4, fun = analog_impact,
                              coord_type = "lonlat", progress = FALSE), args))
      expect_equal(terra::values(tiled[["weighted_mean"]]), terra::values(untiled[["weighted_mean"]]))
})

test_that("point pool companions must be row-aligned, and cell_area_weight needs a raster", {
      d <- point_fixture()
      base <- c(list(x = d$focal, pool = d$pool, n_tiles = 4, fun = analog_impact,
                     coord_type = "projected", progress = FALSE), d$args)
      expect_error(do.call(tiled_analog_search, c(base, list(y = d$y[-1, ]))),
                   "y must be a vector, matrix, or data.frame with one row per row of pool")
      expect_error(do.call(tiled_analog_search, c(base, list(y = d$y, cell_area_weight = TRUE))),
                   "requires a raster pool")
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

test_that("data.frame y works with an empty pool, as in a tile with no pool in its halo", {
      d <- point_fixture()
      far <- d$focal[1:3, ]
      far$x <- far$x + 5000                       # beyond every pool point's reach
      res <- do.call(analog_impact, c(list(x = far, pool = d$pool[0, ], y = d$y[0, ],
                                           coord_type = "projected"), d$args))
      expect_equal(nrow(res), 3)
      expect_true(all(is.na(res$weighted_mean_a)))

      # Through tiling: focals spanning a gap wider than the halo, so some
      # tiles have focals but no pool points
      focal <- rbind(d$focal, far)
      untiled <- do.call(analog_impact, c(list(x = focal, pool = d$pool, y = d$y,
                                               coord_type = "projected"), d$args))
      tiled <- do.call(tiled_analog_search,
                       c(list(x = focal, pool = d$pool, y = d$y, n_tiles = 30,
                              fun = analog_impact, coord_type = "projected",
                              progress = FALSE), d$args))
      expect_equal(strip_meta(tiled), strip_meta(untiled))
})
