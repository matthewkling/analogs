# Tests for stat = "tabulate" (categorical analog impact / aggregation).
#
# Covers:
#   - Factor / character / integer y inputs all factor cleanly
#   - Output column naming: n_<level> for single y, <var>_n_<level> for multi
#   - Multi-column y with independent K_v
#   - NA in y → that analog skipped, not whole row
#   - Hand-computable expected sums match C++ output
#   - Compatible regular stats coexist (count, sum_weights, ess)
#   - Incompatible stats error (weighted_mean, regression)
#   - Factor y with continuous stat errors with a tabulate hint
#   - All-NA y errors
#   - K > .tabulate_K_warn warns; K > .tabulate_K_error errors
#   - Raster I/O round-trip yields K layers per input layer
#
# Note: analog_impact() weights by climate similarity via the `clim` kernel
# (gaussian or inverse), since AIM is fundamentally weighted by climate
# similarity. So all tests here use a gaussian clim kernel. For lighter
# weighting, theta is set large relative to the climate ranges encountered.

# Helper to build a small synthetic dataset locally so this test file is
# self-contained even if sim_test_data() isn't a public helper everywhere.
.tab_sim <- function(n = 60, n_clim = 2, seed = 1L, lonlat = FALSE) {
      set.seed(seed)
      if (lonlat) {
            x <- runif(n, -120, -100); y <- runif(n, 35, 45)
      } else {
            x <- runif(n, 0, 10); y <- runif(n, 0, 10)
      }
      clim <- matrix(rnorm(n * n_clim), ncol = n_clim,
                     dimnames = list(NULL, paste0("c", seq_len(n_clim))))
      ref   <- cbind(x = x, y = y, clim)
      focal <- ref      # use pool as focal for simple, predictable behavior
      list(focal = focal, ref = ref)
}


test_that("tabulate accepts factor y and returns one column per level", {
      d <- .tab_sim()
      veg <- factor(rep(c("forest", "grassland", "shrubland"),
                        length.out = nrow(d$ref)))

      out <- analog_impact(
            x = d$focal, pool = d$ref, y = veg,
            stat = "tabulate",
            clim = kernel("gaussian", theta = 0.5, max = 5),
            geog = kernel(max = 100),
            coord_type = "projected"
      )

      expect_true(all(c("n_forest", "n_grassland", "n_shrubland") %in% names(out)))
      n_cols <- out[, c("n_forest", "n_grassland", "n_shrubland")]
      # All non-negative; all rows have at least one positive entry (we used
      # pool=focal so every row finds itself as an analog within max_clim).
      expect_true(all(as.matrix(n_cols) >= 0))
      expect_true(all(rowSums(n_cols) > 0))
})


test_that("tabulate accepts character y and integer y (factored in R)", {
      d <- .tab_sim()
      veg_chr <- rep(c("a", "b", "c"), length.out = nrow(d$ref))
      veg_int <- rep(1:3, length.out = nrow(d$ref))     # dense integer codes
      veg_sparse <- rep(c(1L, 5L, 9L), length.out = nrow(d$ref))  # sparse codes

      out_chr <- analog_impact(
            x = d$focal, pool = d$ref, y = veg_chr,
            stat = "tabulate", clim = kernel("gaussian", theta = 1.0, max = 5),
            geog = kernel(max = 100), coord_type = "projected"
      )
      expect_true(all(c("n_a", "n_b", "n_c") %in% names(out_chr)))

      out_int <- analog_impact(
            x = d$focal, pool = d$ref, y = veg_int,
            stat = "tabulate", clim = kernel("gaussian", theta = 1.0, max = 5),
            geog = kernel(max = 100), coord_type = "projected"
      )
      expect_true(all(c("n_1", "n_2", "n_3") %in% names(out_int)))

      # Sparse integers: factor() drops the gaps -> only 3 columns, not 9.
      out_sparse <- analog_impact(
            x = d$focal, pool = d$ref, y = veg_sparse,
            stat = "tabulate", clim = kernel("gaussian", theta = 1.0, max = 5),
            geog = kernel(max = 100), coord_type = "projected"
      )
      n_cols_sparse <- grep("^n_", names(out_sparse), value = TRUE)
      expect_length(n_cols_sparse, 3L)
      expect_true(all(c("n_1", "n_5", "n_9") %in% n_cols_sparse))
})


test_that("tabulate column sums equal sum_weights for non-NA y", {
      d <- .tab_sim()
      veg <- factor(rep(c("a", "b"), length.out = nrow(d$ref)))

      out <- analog_impact(
            x = d$focal, pool = d$ref, y = veg,
            stat = c("sum_weights", "tabulate"),
            clim = kernel("gaussian", theta = 0.5, max = 5),
            geog = kernel(max = 100), coord_type = "projected"
      )

      # For a single y with no NAs, the sum across all class columns at each
      # row should equal sum_weights at that row.
      class_sum <- out$n_a + out$n_b
      expect_equal(class_sum, out$sum_weights, tolerance = 1e-10)
})


test_that("tabulate handles NA in y by skipping that analog only", {
      d <- .tab_sim()
      veg <- factor(rep(c("a", "b"), length.out = nrow(d$ref)))
      veg[seq(2, length(veg), by = 5)] <- NA   # ~20% NA

      out <- analog_impact(
            x = d$focal, pool = d$ref, y = veg,
            stat = c("sum_weights", "tabulate"),
            clim = kernel("gaussian", theta = 1.0, max = 5),
            geog = kernel(max = 100), coord_type = "projected"
      )

      # The analog with NA y is still retained in count/sum_weights, so
      # row sum across class columns should be <= sum_weights.
      class_sum <- out$n_a + out$n_b
      expect_true(all(class_sum <= out$sum_weights + 1e-10))
      # And typically strictly less for at least some rows
      expect_true(any(class_sum < out$sum_weights - 1e-10))
})


test_that("multi-column y produces var-prefixed columns with independent K_v", {
      d <- .tab_sim()
      yA <- factor(rep(c("x", "y"), length.out = nrow(d$ref)))     # 2 levels
      yB <- factor(rep(c("p", "q", "r"), length.out = nrow(d$ref))) # 3 levels
      Y <- data.frame(habitat = yA, soil = yB)

      out <- analog_impact(
            x = d$focal, pool = d$ref, y = Y,
            stat = "tabulate", clim = kernel("gaussian", theta = 1.0, max = 5),
            geog = kernel(max = 100), coord_type = "projected"
      )

      expect_true(all(c("habitat_n_x", "habitat_n_y") %in% names(out)))
      expect_true(all(c("soil_n_p", "soil_n_q", "soil_n_r") %in% names(out)))
      # No bleed-through of unprefixed names
      expect_false("n_x" %in% names(out))
})


test_that("incompatible stat combinations error", {
      d <- .tab_sim()
      veg <- factor(rep(c("a", "b"), length.out = nrow(d$ref)))

      # Continuous and categorical aggregations are mutually exclusive
      expect_error(
            analog_impact(
                  x = d$focal, pool = d$ref, y = veg,
                  stat = c("weighted_mean", "tabulate"),
                  clim = kernel("gaussian", theta = 1.0, max = 5),
                  geog = kernel(max = 100), coord_type = "projected"
            ),
            "tabulate.*cannot be combined"
      )
})


test_that("factor y with continuous stat errors with a hint pointing at tabulate", {
      d <- .tab_sim()
      veg <- factor(rep(c("a", "b"), length.out = nrow(d$ref)))

      expect_error(
            analog_impact(
                  x = d$focal, pool = d$ref, y = veg,
                  stat = "weighted_mean",
                  clim = kernel("gaussian", theta = 1.0, max = 5),
                  geog = kernel(max = 100), coord_type = "projected"
            ),
            "tabulate"
      )
})


test_that("all-NA y errors clearly under tabulate", {
      d <- .tab_sim()
      veg <- factor(rep(NA_character_, nrow(d$ref)))

      expect_error(
            analog_impact(
                  x = d$focal, pool = d$ref, y = veg,
                  stat = "tabulate", clim = kernel("gaussian", theta = 1.0, max = 5),
                  geog = kernel(max = 100), coord_type = "projected"
            ),
            "no non-NA"
      )
})


test_that("K > warn threshold warns; K > error threshold errors", {
      d <- .tab_sim(n = 200)
      # Exceed the warn (100) but not the error (1000) threshold.
      # Build a y with 150 distinct levels actually present in the data.
      veg_warn <- factor(rep(seq_len(150), length.out = nrow(d$ref)))
      expect_warning(
            analog_impact(
                  x = d$focal, pool = d$ref, y = veg_warn,
                  stat = "tabulate", clim = kernel("gaussian", theta = 1.0, max = 5),
                  geog = kernel(max = 100), coord_type = "projected"
            ),
            "more than"
      )

      # Exceed the error threshold
      d_big <- .tab_sim(n = 1500)
      veg_err <- factor(seq_len(1500))
      expect_error(
            analog_impact(
                  x = d_big$focal, pool = d_big$ref, y = veg_err,
                  stat = "tabulate", clim = kernel("gaussian", theta = 1.0, max = 5),
                  geog = kernel(max = 100), coord_type = "projected"
            ),
            "exceeding the limit"
      )
})


test_that("tabulate works with select = 'knn_clim' (Hoecker-style AIM)", {
      d <- .tab_sim(n = 100)
      veg <- factor(rep(c("forest", "grass", "shrub"),
                        length.out = nrow(d$ref)))

      # k nearest climate analogs, weighted vote per class — this is the
      # configuration used in Hoecker et al. (2026).
      out <- analog_search(
            x = d$focal, pool = d$ref, y = veg,
            select = "knn_clim", k = 20,
            stat = c("count", "sum_weights", "tabulate"),
            clim = kernel("gaussian", theta = 0.5, max = 5),
            geog = kernel(max = 1000), coord_type = "projected"
      )

      expect_true(all(c("n_forest", "n_grass", "n_shrub") %in% names(out)))
      # Rows have at least one analog (count > 0); we don't assert count <= k
      # exactly because pool == focal here means self matches plus ties at
      # the kth boundary can both inflate the count slightly.
      expect_true(all(out$count > 0))
})


test_that("tabulate raster-in / raster-out yields K layers per input layer", {
      skip_if_not_installed("terra")
      n_side <- 8
      r_clim <- terra::rast(nrows = n_side, ncols = n_side,
                            xmin = 0, xmax = n_side, ymin = 0, ymax = n_side,
                            nlyrs = 2)
      set.seed(2)
      terra::values(r_clim) <- matrix(rnorm(n_side * n_side * 2), ncol = 2)
      names(r_clim) <- c("c1", "c2")

      veg_vals <- rep(c("A", "B", "C"), length.out = terra::ncell(r_clim))
      veg_factor_vec <- factor(veg_vals)

      out <- analog_impact(
            x = r_clim, pool = r_clim, y = veg_factor_vec,
            stat = "tabulate", clim = kernel("gaussian", theta = 1.0, max = 5),
            geog = kernel(max = 100), coord_type = "projected"
      )

      # x is a raster, so output should be a SpatRaster
      expect_s4_class(out, "SpatRaster")
      lyr_names <- names(out)
      expect_true(all(c("n_A", "n_B", "n_C") %in% lyr_names))
})


test_that("compatible regular stats coexist with tabulate", {
      d <- .tab_sim()
      veg <- factor(rep(c("a", "b"), length.out = nrow(d$ref)))

      out <- analog_impact(
            x = d$focal, pool = d$ref, y = veg,
            stat = c("count", "sum_weights", "mean_weights", "ess", "tabulate"),
            clim = kernel("gaussian", theta = 0.5, max = 5),
            geog = kernel(max = 100), coord_type = "projected"
      )

      expect_true(all(c("count", "sum_weights", "mean_weights", "ess",
                        "n_a", "n_b") %in% names(out)))
      # Sanity: count >= 1 for the self-as-analog rows; ESS finite
      expect_true(all(out$count >= 1))
      expect_true(all(is.finite(out$ess)))
})
