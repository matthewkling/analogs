# Tests for the geographic `min` kernel parameter (annulus inner radius).
#
# `min` on a geographic kernel excludes candidates strictly closer than `min`,
# so retained candidates satisfy min <= geog_dist <= max. It is geographic-only
# (env min is an error) and its main use is a spatial buffer for cross-
# validation. These tests validate the annulus geometry against brute force,
# the interaction with max / exclude_self / per-variable env max, the lonlat
# path, error handling, and metadata propagation.

# Brute-force geographic annulus membership (projected / Euclidean geo dist).
# Returns, per focal, the set of 1-based ref indices with min <= d <= max.
bf_geog_annulus <- function(focal, ref, min = 0, max = Inf) {
      gd <- xdist(focal, ref, type = "geog")  # nfocal x nref Euclidean geo dist
      out <- vector("list", nrow(focal))
      for (i in seq_len(nrow(focal))) {
            d <- gd[i, ]
            out[[i]] <- unname(which(d >= min & d <= max))
      }
      out
}


test_that("kernel() stores and validates min", {
      k <- kernel(max = 100, min = 5)
      expect_equal(k$min, 5)
      expect_equal(k$max, 100)

      # min defaults to NULL
      expect_null(kernel(max = 100)$min)

      # invalid mins
      expect_error(kernel(min = -1), "positive finite")
      expect_error(kernel(min = 0), "positive finite")
      expect_error(kernel(min = Inf), "positive finite")
      expect_error(kernel(min = c(1, 2)), "single positive finite")
      expect_error(kernel(min = NA_real_), "positive finite")
})


test_that("geographic annulus matches brute force (select = all, pairs)", {
      d <- sim_test_data(seed = 11, nref = 300, nfocal = 25)

      min_r <- 0.4
      max_r <- 1.2

      res <- analog_search(
            x = d$focal,
            pool = d$ref,
            select = "all",
            geog = kernel(max = max_r, min = min_r),
            coord_type = "projected"
      )

      # Expected annulus membership per focal
      exp_sets <- bf_geog_annulus(d$focal, d$ref, min = min_r, max = max_r)

      # Compare per-focal analog index sets
      for (i in seq_len(nrow(d$focal))) {
            got <- res$analog_index[res$index == i]
            got <- unname(got[!is.na(got)])
            expect_setequal(got, exp_sets[[i]])
      }

      # Every reported geog_dist must lie within the annulus
      gd <- res$geog_dist[!is.na(res$geog_dist)]
      expect_true(all(gd >= min_r - 1e-9))
      expect_true(all(gd <= max_r + 1e-9))
})


test_that("min with no max gives an outer-unbounded annulus", {
      d <- sim_test_data(seed = 12, nref = 250, nfocal = 20)

      min_r <- 0.6
      res <- analog_search(
            x = d$focal,
            pool = d$ref,
            select = "all",
            geog = kernel(min = min_r),   # no max
            coord_type = "projected"
      )

      exp_sets <- bf_geog_annulus(d$focal, d$ref, min = min_r, max = Inf)
      for (i in seq_len(nrow(d$focal))) {
            got <- res$analog_index[res$index == i]
            got <- unname(got[!is.na(got)])
            expect_setequal(got, exp_sets[[i]])
      }
      expect_true(all(res$geog_dist[!is.na(res$geog_dist)] >= min_r - 1e-9))
})


test_that("min excludes near duplicates under knn_geog", {
      d <- sim_test_data(seed = 13, nref = 400, nfocal = 15)

      min_r <- 0.5
      max_r <- 3
      k <- 3
      res <- analog_search(
            x = d$focal,
            pool = d$ref,
            select = "knn_geog",
            geog = kernel(max = max_r, min = min_r),
            k = k,
            coord_type = "projected"
      )

      # Every returned neighbor must lie in the annulus [min_r, max_r]. This
      # specifically guards the lattice knn fast path, which ranks by geographic
      # distance and must honor the inner radius.
      gd_reported <- res$geog_dist[!is.na(res$analog_index)]
      expect_true(all(gd_reported >= min_r - 1e-9 & gd_reported <= max_r + 1e-9))

      # And the neighbors must be exactly the k closest annulus-eligible
      # candidates. Compare sorted DISTANCE vectors (robust to ties and to any
      # pool index remapping).
      gd_all <- unname(xdist(d$focal, d$ref, type = "geog"))
      for (i in seq_len(nrow(d$focal))) {
            got_d <- sort(res$geog_dist[res$index == i & !is.na(res$analog_index)])
            elig  <- gd_all[i, ][gd_all[i, ] >= min_r & gd_all[i, ] <= max_r]
            exp_d <- sort(elig)[seq_len(min(k, length(elig)))]
            expect_equal(got_d, unname(exp_d), tolerance = 1e-6,
                         info = paste("focal", i))
      }
})


test_that("min composes with exclude_self (buffered LOO)", {
      # x == pool; each focal's own row is dropped AND everything within the
      # buffer is dropped. The union should equal: annulus membership excluding
      # the self row (which is at distance 0, already outside the annulus for
      # any positive min, so exclude_self is largely redundant here -- but the
      # combination must not error and must remain within the annulus).
      d <- sim_test_data(seed = 14, nref = 200, nfocal = 200)
      pool <- d$ref  # use ref as both focal and pool

      min_r <- 0.7
      res <- analog_search(
            x = pool,
            pool = pool,
            select = "all",
            geog = kernel(max = 2, min = min_r),
            exclude_self = TRUE,
            coord_type = "projected"
      )

      gd_all <- xdist(pool, pool, type = "geog")
      for (i in seq_len(nrow(pool))) {
            got <- res$analog_index[res$index == i]
            got <- got[!is.na(got)]
            # self excluded
            expect_false(i %in% got, info = paste("focal", i, "self present"))
            # all within annulus
            if (length(got)) {
                  expect_true(all(gd_all[i, got] >= min_r - 1e-9 &
                                        gd_all[i, got] <= 2 + 1e-9),
                              info = paste("focal", i))
            }
      }
})


test_that("scalar geog min composes with per-variable env max", {
      # env uses a per-variable box (max vector); geog uses a scalar annulus.
      # Membership is the intersection of the env box and the geo annulus.
      d <- sim_test_data(seed = 15, nref = 300, nfocal = 20)

      env_band <- c(0.8, 1.1)   # length n_env (= 2)
      min_r <- 0.4
      max_r <- 1.5

      res <- analog_search(
            x = d$focal,
            pool = d$ref,
            select = "all",
            env = kernel(max = env_band),
            geog = kernel(max = max_r, min = min_r),
            coord_type = "projected"
      )

      env_focal <- d$focal[, 3:4, drop = FALSE]
      env_ref   <- d$ref[,  3:4, drop = FALSE]
      env_sets  <- bf_env_band(env_focal, env_ref, env_band)
      geo_sets  <- bf_geog_annulus(d$focal, d$ref, min = min_r, max = max_r)

      for (i in seq_len(nrow(d$focal))) {
            got <- res$analog_index[res$index == i]
            got <- unname(got[!is.na(got)])
            expect_setequal(got, unname(intersect(env_sets[[i]], geo_sets[[i]])))
      }
})


test_that("geographic annulus works for lonlat coordinates", {
      d <- sim_test_data(seed = 16, nref = 300, nfocal = 20, lonlat = TRUE)

      min_km <- 500
      max_km <- 3000
      res <- analog_search(
            x = d$focal,
            pool = d$ref,
            select = "all",
            geog = kernel(max = max_km, min = min_km),
            coord_type = "lonlat"
      )

      # NOTE: the engine filters on ECEF chord distance (converted from km),
      # while xdist() here uses haversine (great-circle) km. The two metrics
      # diverge by a few percent at continental distances, so candidates within
      # a small band around either boundary may legitimately differ. We
      # therefore check (a) every returned neighbor is within the annulus up to
      # a chord/arc tolerance, and (b) points comfortably inside the annulus are
      # all returned. Points near the boundaries are exempt from (b).
      gd <- xdist(d$focal, d$ref, type = "lonlat")  # haversine km
      tol <- 0.03 * max_km  # ~3% covers chord-vs-arc divergence at these ranges
      for (i in seq_len(nrow(d$focal))) {
            got <- res$analog_index[res$index == i]
            got <- unname(got[!is.na(got)])

            # (a) returned neighbors respect the annulus within tolerance
            if (length(got)) {
                  expect_true(all(gd[i, got] >= min_km - tol &
                                        gd[i, got] <= max_km + tol),
                              info = paste("focal", i, "annulus"))
            }

            # (b) points well inside the annulus must be present
            core <- unname(which(gd[i, ] >= min_km + tol & gd[i, ] <= max_km - tol))
            expect_true(all(core %in% got), info = paste("focal", i, "core"))
      }
})


test_that("min on the environmental kernel is rejected", {
      d <- sim_test_data(seed = 17, nref = 100, nfocal = 10)
      expect_error(
            analog_search(
                  x = d$focal,
                  pool = d$ref,
                  select = "all",
                  env = kernel(max = 1, min = 0.2),
                  coord_type = "projected"
            ),
            "geographic family"
      )
})


test_that("empty annulus (min >= max) is rejected", {
      d <- sim_test_data(seed = 18, nref = 100, nfocal = 10)
      expect_error(
            analog_search(
                  x = d$focal,
                  pool = d$ref,
                  select = "all",
                  geog = kernel(max = 1, min = 1),
                  coord_type = "projected"
            ),
            "annulus is empty"
      )
      expect_error(
            analog_search(
                  x = d$focal,
                  pool = d$ref,
                  select = "all",
                  geog = kernel(max = 1, min = 2),
                  coord_type = "projected"
            ),
            "annulus is empty"
      )
})


test_that("min_geog is recorded in metadata", {
      d <- sim_test_data(seed = 19, nref = 150, nfocal = 12)

      res <- analog_search(
            x = d$focal,
            pool = d$ref,
            select = "all",
            geog = kernel(max = 2, min = 0.5),
            coord_type = "projected"
      )
      expect_equal(attr(res, "min_geog"), 0.5)
      expect_equal(metadata(res)$min_geog, 0.5)

      # Absent when not set
      res2 <- analog_search(
            x = d$focal,
            pool = d$ref,
            select = "all",
            geog = kernel(max = 2),
            coord_type = "projected"
      )
      expect_null(attr(res2, "min_geog"))
})
