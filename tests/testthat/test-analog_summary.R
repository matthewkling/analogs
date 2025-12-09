test_that("analog_summary doesn't error", {
      d <- sim_test_data()
      result <- analog_velocity(d$focal, d$ref, k = 5, max_clim = .5)
      expect_no_error(capture.output(analog_summary(result)))
})
