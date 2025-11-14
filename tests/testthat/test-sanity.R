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
