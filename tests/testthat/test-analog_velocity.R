test_that("`analog_velocity` result matches manual calculation for planar coords", {

      d <- sim_test_data(lonlat = TRUE)
      max_env <- 1

      # velocity
      nn <- analog_velocity(d$focal, d$ref, env = kernel(max = max_env), k = 1, coord_type = "projected")

      # manual velocity calculation
      denv <- xdist(d$focal, d$ref, "env")
      dgeog <- xdist(d$focal, d$ref, "geog")
      dgeog[denv > max_env] <- Inf
      nn_idx <- as.vector(apply(dgeog, 1, which.min))
      nn_dst <- as.vector(apply(dgeog, 1, min))

      expect_equal(nn$analog_index, nn_idx)
      expect_equal(nn$geog_dist, nn_dst)
})


test_that("`analog_velocity` result matches manual calculation for LON/LAT coords", {

      d <- sim_test_data(lonlat = TRUE)
      max_env <- 1

      # velocity (ECEF chord-space internal; output arc-length km)
      nn <- analog_velocity(d$focal, d$ref,
                            env = kernel(max = max_env),
                            k = 1,
                            coord_type = "lonlat")

      # --- manual lon/lat nearest-analog calculation ---

      # (1) compute environmental distance matrix (Euclidean)
      denv <- xdist(d$focal, d$ref, "env")

      # (2) compute true great-circle geog distances (km)
      dgeog <- xdist(d$focal, d$ref, "lonlat")

      # (3) apply environmental filter
      dgeog[denv > max_env] <- Inf

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
      index <- build_analog_index(d$ref, coord_type = "projected")

      # Query with index
      v <- analog_velocity(
            x = d$focal,
            pool = index,
            env = kernel(max = 1),
            k = 1
      )

      expect_s3_class(v, "data.frame")
      expect_equal(nrow(v), nrow(d$focal))
      expect_true(all(c("index", "analog_index", "geog_dist") %in% names(v)))
})


# Lon/lat kNN queries use the lattice's expanding search in ECEF chord space.
# These tests compare against brute-force great-circle distances on a regional
# pool large enough that the lattice has several bins per axis. Data are
# continuous, so there are no exact ties at the climate threshold.
lonlat_regional_data <- function(seed = 7, nref = 3000, nfocal = 150) {
      set.seed(seed)
      mk <- function(n) {
            m <- cbind(runif(n, -115, -95), runif(n, 30, 50), rnorm(n), rnorm(n))
            m[, 3] <- m[, 3] + (m[, 2] - 40) / 5   # climate gradient with latitude
            m
      }
      list(ref = mk(nref), focal = mk(nfocal))
}

test_that("lon/lat knn_geog (k > 1, geographic annulus) matches brute force", {
      d <- lonlat_regional_data()
      max_env <- 0.6
      k <- 3
      nn <- analog_search(d$focal, d$ref, select = "knn_geog", k = k,
                          env = kernel(max = max_env),
                          geog = kernel(max = 800, min = 50),
                          coord_type = "lonlat")

      denv <- xdist(d$focal, d$ref, "env")
      dgeog <- xdist(d$focal, d$ref, "lonlat")
      dgeog[denv > max_env | dgeog > 800 | dgeog < 50] <- Inf
      for (i in seq_len(nrow(d$focal))) {
            want <- unname(sort(dgeog[i, ])[1:k])
            want <- want[is.finite(want)]
            got <- sort(nn$geog_dist[nn$index == i])
            expect_equal(got, want, tolerance = 1e-8)
      }
})

test_that("lon/lat knn_env within a geographic radius matches brute force", {
      d <- lonlat_regional_data()
      k <- 4
      nn <- analog_search(d$focal, d$ref, select = "knn_env", k = k,
                          geog = kernel(max = 300), coord_type = "lonlat")

      denv <- xdist(d$focal, d$ref, "env")
      dgeog <- xdist(d$focal, d$ref, "lonlat")
      denv[dgeog > 300] <- Inf
      for (i in seq_len(nrow(d$focal))) {
            want <- unname(sort(denv[i, ])[1:k])
            want <- want[is.finite(want)]
            got <- sort(nn$env_dist[nn$index == i])
            expect_equal(got, want, tolerance = 1e-8)
      }
})
