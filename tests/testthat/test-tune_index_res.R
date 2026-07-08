# Tests for tune_index_res(), the per-family lattice resolution tuner.
#
# Contract (post per-family rework):
#   - Returns a list with numeric elements `geog_res_adj` and `env_res_adj`.
#   - Only tunes when there are > 2000 focal points AND the family is active
#     (adjustment > 0). Otherwise inputs are returned unchanged.
#   - A deactivated family (adjustment 0) is passed through untouched.
#   - Tuned adjustments are bounded to [1/32, 32].

ADJ_MIN <- 1 / 32
ADJ_MAX <- 32

# Helper: assert the return has the expected list shape.
expect_adj_list <- function(res) {
      expect_type(res, "list")
      expect_named(res, c("geog_res_adj", "env_res_adj"), ignore.order = TRUE)
      expect_true(is.numeric(res$geog_res_adj) && length(res$geog_res_adj) == 1L)
      expect_true(is.numeric(res$env_res_adj) && length(res$env_res_adj) == 1L)
}

# Large dataset (> 2000 focals) that actually triggers tuning.
make_large <- function(seed = 123, nfocal = 2500, nref = 500, lonlat = FALSE) {
      set.seed(seed)
      focal <- matrix(rnorm(nfocal * 4), ncol = 4)
      ref   <- matrix(rnorm(nref * 4), ncol = 4)
      if (lonlat) {
            focal[, 1] <- runif(nfocal, -180, 180); focal[, 2] <- runif(nfocal, -90, 90)
            ref[, 1]   <- runif(nref, -180, 180);   ref[, 2]   <- runif(nref, -90, 90)
      }
      list(focal = focal, ref = ref)
}


# ---- Pass-through path (small datasets: n_focal <= 2000) -------------------

test_that("tune_index_res returns a valid per-family list", {
      d <- sim_test_data()   # nfocal = 20 -> pass-through

      res <- tune_index_res(
            x = d$focal, pool = d$ref,
            select = "knn_geog", env = kernel(max = 1), k = 1,
            coord_type = "projected", verbose = FALSE
      )

      expect_adj_list(res)
      expect_true(res$geog_res_adj > 0)
      expect_true(res$env_res_adj > 0)
})


test_that("tune_index_res passes inputs through unchanged on small data", {
      d <- sim_test_data()   # nfocal = 20 -> pass-through

      res <- tune_index_res(
            x = d$focal, pool = d$ref,
            select = "all", stat = "count", env = kernel(max = 1), geog = kernel(max = 2),
            coord_type = "projected",
            geog_res_adj = 1.5, env_res_adj = 0.5,
            verbose = FALSE
      )

      # Below the tuning threshold, inputs are returned unchanged.
      expect_equal(res$geog_res_adj, 1.5)
      expect_equal(res$env_res_adj, 0.5)
})


test_that("tune_index_res works with different modes (pass-through)", {
      d <- sim_test_data()

      # knn_env
      res1 <- tune_index_res(
            x = d$focal, pool = d$ref,
            select = "knn_env", geog = kernel(max = 2), k = 3,
            coord_type = "projected", verbose = FALSE
      )
      expect_adj_list(res1)

      # count aggregation
      res2 <- tune_index_res(
            x = d$focal, pool = d$ref,
            stat = "count", env = kernel(max = 1), geog = kernel(max = 2),
            coord_type = "projected", verbose = FALSE
      )
      expect_adj_list(res2)

      # sum_weights aggregation
      res3 <- tune_index_res(
            x = d$focal, pool = d$ref,
            stat = "sum_weights", env = kernel(max = 1), geog = kernel(max = 2),
            coord_type = "projected", verbose = FALSE
      )
      expect_adj_list(res3)
})


test_that("tune_index_res works with auto coord detection (pass-through)", {
      d_proj <- sim_test_data(lonlat = FALSE)
      res_proj <- tune_index_res(
            x = d_proj$focal, pool = d_proj$ref,
            stat = "count", env = kernel(max = 1),
            coord_type = "auto", verbose = FALSE
      )
      expect_adj_list(res_proj)

      d_lonlat <- sim_test_data(lonlat = TRUE)
      res_lonlat <- tune_index_res(
            x = d_lonlat$focal, pool = d_lonlat$ref,
            stat = "count", env = kernel(max = 1),
            coord_type = "auto", verbose = FALSE
      )
      expect_adj_list(res_lonlat)
})


test_that("tune_index_res handles tiny datasets (pass-through)", {
      small_focal <- matrix(rnorm(5 * 3), ncol = 3)
      small_ref   <- matrix(rnorm(10 * 3), ncol = 3)

      res <- tune_index_res(
            x = small_focal, pool = small_ref,
            stat = "count", env = kernel(max = 1),
            coord_type = "projected",
            geog_res_adj = 2, env_res_adj = 2,
            verbose = FALSE
      )

      expect_adj_list(res)
      expect_equal(res$geog_res_adj, 2)
      expect_equal(res$env_res_adj, 2)
})


# ---- Tuning path (large datasets: n_focal > 2000) -------------------------

test_that("tune_index_res tunes large datasets within bounds", {
      d <- make_large(seed = 123)

      res <- tune_index_res(
            x = d$focal, pool = d$ref,
            stat = "count", env = kernel(max = 1), geog = kernel(max = 2),
            coord_type = "projected", verbose = FALSE
      )

      expect_adj_list(res)
      expect_true(res$geog_res_adj >= ADJ_MIN && res$geog_res_adj <= ADJ_MAX)
      expect_true(res$env_res_adj >= ADJ_MIN && res$env_res_adj <= ADJ_MAX)
})


test_that("tune_index_res keeps tuned adjustments in [1/32, 32]", {
      d <- make_large(seed = 456)

      res <- tune_index_res(
            x = d$focal, pool = d$ref,
            select = "knn_geog", env = kernel(max = 1), k = 1,
            coord_type = "projected", verbose = FALSE
      )

      expect_adj_list(res)
      # knn_geog constrains environmental (filter) and sorts geographically, so both
      # families are active and tuned; both must stay within bounds.
      expect_true(res$geog_res_adj >= ADJ_MIN && res$geog_res_adj <= ADJ_MAX)
      expect_true(res$env_res_adj >= ADJ_MIN && res$env_res_adj <= ADJ_MAX)
})


test_that("tune_index_res passes a deactivated family through untouched", {
      d <- make_large(seed = 321)

      # environmental deactivated (adj 0): should come back exactly 0, geo tuned.
      res <- tune_index_res(
            x = d$focal, pool = d$ref,
            stat = "count", geog = kernel(max = 2),
            coord_type = "projected",
            geog_res_adj = 1, env_res_adj = 0,
            verbose = FALSE
      )
      expect_adj_list(res)
      expect_equal(res$env_res_adj, 0)                 # untouched
      expect_true(res$geog_res_adj >= ADJ_MIN &&
                        res$geog_res_adj <= ADJ_MAX)          # tuned, in bounds

      # Geo deactivated (adj 0): should come back exactly 0, environmental tuned.
      res2 <- tune_index_res(
            x = d$focal, pool = d$ref,
            stat = "count", env = kernel(max = 1),
            coord_type = "projected",
            geog_res_adj = 0, env_res_adj = 1,
            verbose = FALSE
      )
      expect_adj_list(res2)
      expect_equal(res2$geog_res_adj, 0)                # untouched
      expect_true(res2$env_res_adj >= ADJ_MIN &&
                        res2$env_res_adj <= ADJ_MAX)         # tuned, in bounds
})


test_that("tune_index_res returns inputs when both families deactivated", {
      d <- make_large(seed = 654)

      res <- tune_index_res(
            x = d$focal, pool = d$ref,
            stat = "count",
            coord_type = "projected",
            geog_res_adj = 0, env_res_adj = 0,
            verbose = FALSE
      )
      expect_adj_list(res)
      expect_equal(res$geog_res_adj, 0)
      expect_equal(res$env_res_adj, 0)
})


test_that("tune_index_res verbose output works", {
      d <- make_large(seed = 789)

      # Verbose tuning should emit sweep / convergence progress.
      expect_message(
            tune_index_res(
                  x = d$focal, pool = d$ref,
                  stat = "count", env = kernel(max = 1), geog = kernel(max = 2),
                  coord_type = "projected", verbose = TRUE
            ),
            "Sweep|Converged"
      )

      # Quiet by default.
      expect_silent(
            tune_index_res(
                  x = d$focal, pool = d$ref,
                  stat = "count", env = kernel(max = 1), geog = kernel(max = 2),
                  coord_type = "projected", verbose = FALSE
            )
      )
})


test_that("tune_index_res works with kernel weight options", {
      d <- make_large(seed = 999)

      res1 <- tune_index_res(
            x = d$focal, pool = d$ref,
            stat = "sum_weights", env = kernel("inverse", max = 1), geog = kernel(max = 2),
            coord_type = "projected", verbose = FALSE
      )
      expect_adj_list(res1)

      res2 <- tune_index_res(
            x = d$focal, pool = d$ref,
            stat = "sum_weights", env = kernel(max = 1), geog = kernel("inverse", max = 2),
            coord_type = "projected", verbose = FALSE
      )
      expect_adj_list(res2)
})


test_that("tune_index_res works with lonlat coordinates", {
      d <- make_large(seed = 111, lonlat = TRUE)

      res <- tune_index_res(
            x = d$focal, pool = d$ref,
            stat = "count", env = kernel(max = 1),
            coord_type = "lonlat", verbose = FALSE
      )

      expect_adj_list(res)
      expect_true(res$geog_res_adj >= ADJ_MIN && res$geog_res_adj <= ADJ_MAX)
      expect_true(res$env_res_adj >= ADJ_MIN && res$env_res_adj <= ADJ_MAX)
})
