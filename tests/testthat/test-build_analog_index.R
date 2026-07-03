test_that("build_analog_index creates valid index object", {

      set.seed(123)
      ref_data <- matrix(rnorm(200 * 4), ncol = 4)
      colnames(ref_data) <- c("x", "y", "clim1", "clim2")

      # Test building index with projected coords
      index <- build_analog_index(ref_data, coord_type = "projected")

      # Check class
      expect_s3_class(index, "analog_index")
      expect_true(is_analog_index(index))

      # Check components exist
      expect_true(!is.null(index$lattice_xptr))
      expect_true(!is.null(index$ref_data))
      expect_equal(index$coord_type, "projected")
      expect_equal(index$n_pool, 200)
      expect_equal(index$n_clim, 2)

      # Check ranges are sensible
      expect_length(index$coord_mins, 2)
      expect_length(index$coord_maxs, 2)
      expect_length(index$clim_mins, 2)
      expect_length(index$clim_maxs, 2)
      expect_true(all(index$coord_maxs > index$coord_mins))
      expect_true(all(index$clim_maxs > index$clim_mins))

      # Check diagnostics
      expect_true(index$total_bins > 0)
      expect_true(index$n_bins_nonempty > 0)
      expect_true(index$n_bins_nonempty <= index$total_bins)
})


test_that("build_analog_index works with auto coord detection", {

      # Projected data
      set.seed(456)
      ref_proj <- matrix(rnorm(100 * 3), ncol = 3)
      ref_proj[, 1] <- ref_proj[, 1] * 1000  # large x values
      ref_proj[, 2] <- ref_proj[, 2] * 1000  # large y values

      index_proj <- build_analog_index(ref_proj, coord_type = "auto")
      expect_equal(index_proj$coord_type, "projected")

      # Lon/lat data
      ref_lonlat <- matrix(rnorm(100 * 3), ncol = 3)
      ref_lonlat[, 1] <- runif(100, -180, 180)  # lon
      ref_lonlat[, 2] <- runif(100, -90, 90)    # lat

      index_lonlat <- build_analog_index(ref_lonlat, coord_type = "auto")
      expect_equal(index_lonlat$coord_type, "lonlat")
      expect_true(index_lonlat$use_ecef)
})


test_that("build_analog_index handles lonlat coords correctly", {

      set.seed(789)
      ref_lonlat <- matrix(rnorm(150 * 3), ncol = 3)
      ref_lonlat[, 1] <- runif(150, -180, 180)  # lon
      ref_lonlat[, 2] <- runif(150, -90, 90)    # lat

      index <- build_analog_index(ref_lonlat, coord_type = "lonlat")

      expect_equal(index$coord_type, "lonlat")
      expect_true(index$use_ecef)
      expect_true(!is.null(index$ecef_xptr))  # Should have ECEF data

      # Coordinate ranges should be in lon/lat space (not ECEF)
      expect_true(index$coord_mins[1] >= -180 && index$coord_maxs[1] <= 180)
      expect_true(index$coord_mins[2] >= -90 && index$coord_maxs[2] <= 90)
})


test_that("build_analog_index validates inputs", {

      ref_data <- matrix(rnorm(100 * 3), ncol = 3)

      # Invalid geog_res_adj
      expect_error(
            build_analog_index(ref_data, geog_res_adj = -5),
            "non-negative"
      )
      expect_error(
            build_analog_index(ref_data, geog_res_adj = c(1, 2)),
            "non-negative"
      )

      # Invalid coord_type
      expect_error(
            build_analog_index(ref_data, coord_type = "invalid"),
            "should be one of"
      )
})


test_that("is_analog_index correctly identifies objects", {

      ref_data <- matrix(rnorm(100 * 3), ncol = 3)
      index <- build_analog_index(ref_data)

      expect_true(is_analog_index(index))
      expect_false(is_analog_index(ref_data))
      expect_false(is_analog_index(list(a = 1, b = 2)))
      expect_false(is_analog_index(NULL))
})


test_that(".validate_analog_index catches invalid indices", {

      ref_data <- matrix(rnorm(100 * 3), ncol = 3)
      index <- build_analog_index(ref_data)

      # Valid index should pass
      expect_invisible(.validate_analog_index(index))

      # Non-index should fail
      expect_error(
            .validate_analog_index(list(a = 1)),
            "not an analog_index"
      )

      # Index with missing components should fail
      bad_index <- index
      bad_index$lattice_xptr <- NULL
      expect_error(
            .validate_analog_index(bad_index),
            "missing components"
      )
})


test_that(".validate_analog_index validates query data", {

      ref_data <- matrix(rnorm(100 * 3), ncol = 3)
      index <- build_analog_index(ref_data)

      # Matching query data should pass
      query_data <- matrix(rnorm(20 * 3), ncol = 3)
      expect_invisible(.validate_analog_index(index, query_data))

      # Mismatched dimensions should fail
      bad_query <- matrix(rnorm(20 * 4), ncol = 4)  # Too many columns
      expect_error(
            .validate_analog_index(index, bad_query),
            "has 4 columns but index expects 3"
      )
})


test_that(".validate_analog_index warns about out-of-bounds queries", {

      set.seed(999)
      ref_data <- matrix(rnorm(100 * 3), ncol = 3)
      index <- build_analog_index(ref_data)

      # Query data far outside bounds
      query_out <- matrix(0, nrow = 10, ncol = 3)
      query_out[, 1] <- 1000  # Way outside x range
      query_out[, 2] <- 1000  # Way outside y range

      expect_warning(
            .validate_analog_index(index, query_out, validate_ranges = TRUE),
            "fall outside index.*bounds"
      )
})


test_that("print.analog_index produces output", {

      ref_data <- matrix(rnorm(50 * 3), ncol = 3)
      index <- build_analog_index(ref_data)

      # Should not error
      expect_output(print(index), "Analog Index")
      expect_output(print(index), "Reference data")
      expect_output(print(index), "50 locations")
      expect_output(print(index), "1 climate variable")
      expect_output(print(index), "Index structure")
})


test_that("build_analog_index works with different resolutions", {

      ref_data <- matrix(rnorm(100 * 3), ncol = 3)

      # Low resolution
      index_low <- build_analog_index(ref_data, geog_res_adj = 1/2)

      # High resolution
      index_high <- build_analog_index(ref_data, geog_res_adj = 2)

      # Higher resolution setting should be recorded as higher
      expect_true(index_high$geog_res_adj > index_low$geog_res_adj)
})

test_that("build_analog_index works with downsampling", {

      d <- sim_test_data(nref = 1000000)

      # Build full index
      index_full <- build_analog_index(
            d$ref,
            coord_type = "projected",
            downsample = 1.0
      )

      # Build downsampled index
      index_down <- build_analog_index(
            d$ref,
            coord_type = "projected",
            downsample = 0.1,
            seed = 456
      )

      expect_equal(index_down$downsample_target, 0.1)
      expect_equal(index_down$downsample_actual, index_down$downsample_target, tolerance = .02)


      # Test reproducibility
      index_down2 <- build_analog_index(
            d$ref,
            coord_type = "projected",
            downsample = 0.1,
            seed = 456  # Same seed
      )
      expect_equal(index_down$downsample_actual, index_down2$downsample_actual)
})


# test coord_type = "auto" detection ----------------------------

# The detector should prefer CRS metadata (when available via SpatRaster /
# SpatVector) over coordinate-magnitude heuristics. The magnitude heuristic
# is the fallback for matrix / data.frame inputs that have no CRS, and for
# rasters whose CRS is unknown.

test_that(".detect_geo prefers SpatRaster CRS over coordinate magnitudes", {
      skip_if_not_installed("terra")

      # Pool with projected CRS but coordinate values that fall inside lonlat
      # bounds (xmin/xmax/ymin/ymax all in [-180, 180] / [-90, 90]). Without
      # the CRS-aware detector, magnitude alone would say lonlat. With it,
      # the EPSG:32611 (UTM) metadata wins and we get projected.
      r <- terra::rast(ncol = 5, nrow = 5,
                       xmin = 0, xmax = 24, ymin = 0, ymax = 24,
                       crs = "EPSG:32611")
      terra::values(r) <- rnorm(25)

      xy <- terra::crds(r)
      expect_identical(.detect_geo(xy, r), "projected")

      # Sanity: same xy without the SpatRaster falls back to magnitude
      # heuristic, which (incorrectly but unavoidably) returns lonlat.
      expect_identical(.detect_geo(xy), "lonlat")
})


test_that(".detect_geo recognizes a lonlat SpatRaster", {
      skip_if_not_installed("terra")

      r <- terra::rast(ncol = 5, nrow = 5,
                       xmin = -10, xmax = 10, ymin = 30, ymax = 50,
                       crs = "EPSG:4326")
      terra::values(r) <- rnorm(25)

      expect_identical(.detect_geo(terra::crds(r), r), "lonlat")
})


test_that(".detect_geo falls back to magnitude when SpatRaster CRS is unknown", {
      skip_if_not_installed("terra")

      r <- terra::rast(ncol = 5, nrow = 5,
                       xmin = 0, xmax = 24, ymin = 0, ymax = 24,
                       crs = "")                # empty CRS
      terra::values(r) <- rnorm(25)

      # is.lonlat() returns NA when CRS is unknown (terra also emits a
      # warning, which our detector suppresses — assert no warning leaks).
      # Magnitude fallback then says lonlat since values fit in lonlat bounds.
      expect_no_warning(
            res <- .detect_geo(terra::crds(r), r)
      )
      expect_identical(res, "lonlat")
})


test_that(".detect_geo magnitude check works for matrix/data.frame inputs", {
      # Lon/lat-shaped values
      m_ll <- cbind(runif(10, -180, 180), runif(10, -90, 90))
      expect_identical(.detect_geo(m_ll), "lonlat")

      # Out-of-lonlat values (e.g. UTM)
      m_proj <- cbind(runif(10, 4e5, 5e5), runif(10, 4e6, 5e6))
      expect_identical(.detect_geo(m_proj), "projected")
})


test_that("build_analog_index auto-detect uses CRS when pool is a SpatRaster", {
      skip_if_not_installed("terra")

      # Same setup as the .detect_geo test — coords inside lonlat bounds, but
      # CRS is projected. Build the index and confirm the resolved coord_type.
      r <- terra::rast(ncol = 5, nrow = 5,
                       xmin = 0, xmax = 24, ymin = 0, ymax = 24,
                       crs = "EPSG:32611")
      terra::values(r) <- rnorm(25)

      idx <- build_analog_index(r)         # coord_type defaults to "auto"
      expect_identical(idx$coord_type, "projected")
})
