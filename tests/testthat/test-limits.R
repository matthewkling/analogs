# Distance limits are inclusive, and every engine path must apply them the
# same way. Points lying exactly on a limit are common on regular grids and
# with rounded data; floating-point rounding must not include them in one path
# (e.g. select = "all") and exclude them in another (the kNN fast path, which
# compares squared distances).

test_that("climate limit ties are treated the same by all-mode and kNN queries", {
      set.seed(3)
      n <- 400
      pool <- cbind(x = runif(n, 0, 100), y = runif(n, 0, 100),
                    e1 = round(rnorm(n), 2), e2 = round(rnorm(n), 2))
      focal <- pool[1:40, ]
      # one analog per focal at climate distance exactly 0.5 (a 0.3-0.4-0.5 triangle)
      ties <- cbind(x = focal[, 1] + 1, y = focal[, 2],
                    e1 = focal[, 3] - 0.4, e2 = focal[, 4] + 0.3)
      pool <- rbind(pool, ties)

      cnt <- analog_search(focal, pool, stat = "count", env = kernel(max = 0.5),
                           coord_type = "projected")$count
      nn <- analog_search(focal, pool, select = "knn_geog", k = nrow(pool),
                          env = kernel(max = 0.5), coord_type = "projected")
      expect_identical(as.integer(cnt), tabulate(nn$index, nbins = nrow(focal)))

      # and the tied analogs are included (limits are inclusive)
      tie_rows <- n + seq_len(nrow(focal))
      expect_true(all(tie_rows %in% nn$analog_index))
})

test_that("geographic limit ties are treated the same by all-mode and kNN queries", {
      # 0.1-spaced grid: cells at offsets (0.3, 0.4) lie exactly 0.5 away, and
      # at 0.2 (inner radius) along the axes
      g <- as.matrix(expand.grid(x = seq(0, 4, by = 0.1), y = seq(0, 4, by = 0.1)))
      pool <- cbind(g, e1 = 0, e2 = 0)
      set.seed(5)
      focal <- pool[sample(nrow(pool), 60), ]

      for (gk in list(kernel(max = 0.5), kernel(max = 0.5, min = 0.2))) {
            cnt <- analog_search(focal, pool, stat = "count", geog = gk,
                                 coord_type = "projected")$count
            nn <- analog_search(focal, pool, select = "knn_geog", k = nrow(pool),
                                geog = gk, coord_type = "projected")
            expect_identical(as.integer(cnt), tabulate(nn$index, nbins = nrow(focal)))
      }

      # Interior focals see every grid point in the closed disk of radius 5
      # cells (81 points); the inner radius of 2 cells removes the 9 points
      # strictly closer than 2 cells and keeps those exactly 2 cells away
      inner <- focal[, 1] >= 0.5 & focal[, 1] <= 3.5 & focal[, 2] >= 0.5 & focal[, 2] <= 3.5
      cnt_disk <- analog_search(focal, pool, stat = "count", geog = kernel(max = 0.5),
                                coord_type = "projected")$count
      cnt_ann <- analog_search(focal, pool, stat = "count", geog = kernel(max = 0.5, min = 0.2),
                               coord_type = "projected")$count
      expect_true(all(cnt_disk[inner] == 81))
      expect_true(all(cnt_ann[inner] == 81 - 9))
})
