test_that("euclid distance from C++ works", {
      expect_equal(analogs:::analogs_euclid_cpp(c(0,0), c(3,4)), 5)
})

test_that("haversine distance from C++ is sensible", {
      # Using mean Earth radius R = 6371.0088 km in the C++ implementation
      # => ~111.195 km per degree of central angle.
      d1 <- analogs:::analogs_haversine_cpp(c(0,0), c(0,1))  # 1° latitude
      expect_equal(d1, 111.195, tolerance = 1e-3)

      d2 <- analogs:::analogs_haversine_cpp(c(0,0), c(1,0))  # 1° longitude at equator
      expect_equal(d2, 111.195, tolerance = 1e-3)
})

test_that("simple scalar radius works (1D example)", {
      focal <- matrix(c(0, 1), ncol = 1)
      ref   <- matrix(0:10, ncol = 1)

      got <- analogs:::find_analogs_core(focal, ref,
                                         NULL, NULL,
                                         NA, 0,          # k=0 -> radius mode
                                         1,              # max_clim radius
                                         NA, NA, NA)

      expect_equal(got[[1]], 1:2)
      expect_equal(got[[2]], 1:3)
      expect_identical(attr(got, "n_ref"), nrow(ref))
      expect_identical(attr(got, "n_vars"), ncol(ref))
})

test_that("radius vs brute-force agrees (random small, multi-dim)", {
      set.seed(1)
      ref   <- matrix(rnorm(200), ncol = 4)  # 50 x 4
      focal <- matrix(rnorm(40),  ncol = 4)  # 10 x 4
      R <- 1.25

      bf <- bf_clim_radius(focal, ref, R)
      got <- analogs:::find_analogs_core(focal, ref,
                                         NULL, NULL, NA,
                                         0,          # k=0 radius
                                         R, NA, NA, NA)

      # Compare as sets for each focal (order not guaranteed in radius mode)
      for (i in seq_along(bf)) {
            expect_setequal(got[[i]], bf[[i]])
      }
})

test_that("k-NN matches brute-force (indices only, no ties)", {
      set.seed(2)
      ref   <- matrix(rnorm(300), ncol = 3)   # 100 x 3
      focal <- matrix(rnorm(30),  ncol = 3)   # 10 x 3
      k <- 5

      bf <- bf_clim_knn(focal, ref, k)
      got <- analogs:::find_analogs_core(focal, ref,
                                         NULL, NULL, NA,
                                         k,      # k>0 -> best-first top-k
                                         1e9,    # large radius to avoid pruning; ranking drives result
                                         NA, NA, NA)

      # Our core returns ascending distance; bf uses stable order(d, idx)
      # To reduce tie risk, we used continuous rnorm; still compare exact equality
      for (i in seq_along(bf)) {
            expect_equal(got[[i]], bf[[i]])
      }
})

test_that("per-var band filter works like brute-force band", {
      set.seed(3)
      ref   <- matrix(rnorm(240), ncol = 4)   # 60 x 4
      focal <- matrix(rnorm(20),  ncol = 4)   # 5 x 4
      band  <- c(0.5, 0.75, 1.0, 0.25)

      bf <- bf_clim_band(focal, ref, band)
      got <- analogs:::find_analogs_core(focal, ref,
                                         NULL, NULL, NA,
                                         0,        # radius mode, but with band vector
                                         band,     # pass band vector into climate_band_
                                         NA, NA, NA)

      for (i in seq_along(bf)) {
            expect_setequal(got[[i]], bf[[i]])
      }
})

test_that("returns 1-based indices and handles empty results", {
      ref   <- matrix(c(10, 20, 30), ncol = 1)
      focal <- matrix(0, ncol = 1)
      R <- 1e-6

      got <- analogs:::find_analogs_core(focal, ref, NULL, NULL, NA, 0, R, NA, NA, NA)
      expect_length(got[[1]], 0)
      # If we artificially set R large, it should return 1..n_ref
      got2 <- analogs:::find_analogs_core(focal, ref, NULL, NULL, NA, 0, 1e9, NA, NA, NA)
      expect_equal(got2[[1]], seq_len(nrow(ref)))
})
