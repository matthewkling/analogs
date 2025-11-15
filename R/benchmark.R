## internal dev functions for benchmarking -- will not be included in final production package ##

# simulate reasonably realistic data with spatiotemporal structure -------------

simulate_climate_rasters <- function(d = 500) {
      v1 <- matrix(rep(1:d, each = d), d) %>% "/"(40) %>% sin()
      v1 <- matrix(rep(1:d, d), d) %>% "/"(20) %>% cos() %>% "+"(v1) %>% rast()
      v2 <- matrix(rep(1:d, d), d) %>% rast()
      clim1 <- c(v1, v2) %>% setNames(c("t", "p"))
      clim1[] <- scale(clim1[])

      clim2 <- clim1
      clim2[[1]] <- clim2[[1]] + .5
      clim2[[2]] <- clim2[[2]] * .5

      list(clim1 = clim1, clim2 = clim2)

      # focal <- as.data.frame(clim1, xy = TRUE) %>% sample_n(100)
      # colnames(focal) <- c("x", "y", "t", "p")
}


# load example climate data ---------------------------------

load_climate_rasters <- function() {
      require(terra)
      require(tidyverse)

      scale_rast <- function(x) {
            v <- values(x, mat = FALSE) # Returns a vector, not matrix
            v <- v[is.finite(v)] # Remove NAs before calculating stats
            (x - mean(v)) / sd(v)
      }

      # prep climate data
      bb <- ext(-114, -108, 43, 49)

      t1 <- rast(
            "/Volumes/T7/CHELSA/v2/raw/CHELSA_bio1_1981-2010_V.2.1.tif"
      ) %>%
            crop(bb)
      t2 <- list.files("/Volumes/T7/CHELSA/v2/cmip/", full.names = T)
      t2 <- t2[grepl("2071", t2) & grepl("126_tas_", t2)] %>%
            rast() %>%
            crop(bb) %>%
            terra::mean()
      t <- c(t1, t2) %>% scale_rast()

      p1 <- rast(
            "/Volumes/T7/CHELSA/v2/raw/CHELSA_bio12_1981-2010_V.2.1.tif"
      ) %>%
            crop(bb)
      p2 <- list.files("/Volumes/T7/CHELSA/v2/cmip/", full.names = T)
      p2 <- p2[grepl("2071", p2) & grepl("126_pr_", p2)] %>%
            rast() %>%
            crop(bb) %>%
            terra::mean()
      p <- c(p1, p2) %>% log() %>% scale_rast()

      clim1 <- c(t[[1]], p[[1]])
      clim2 <- c(t[[2]], p[[2]])

      return(list(clim1 = clim1, clim2 = clim2))
}


# benchmarks ----------------------------------------

benchmark_velocity <- function(clim, n_focal = 100) {
      # climate velocity: few focals, many refs
      focal <- as.data.frame(clim$clim1, xy = TRUE) %>% sample_n(n_focal)
      ref <- clim$clim2
      st <- system.time({
            a <- analog_velocity(focal, ref, max_clim = .5, k = 1)
      })
      return(st[["elapsed"]])
      # cat("velocity: ", st[["elapsed"]], "\n")
}

benchmark_impact <- function(clim, n_focal = 100) {
      # analog impact: dist and geo constraints; few focals, many refs
      focal <- as.data.frame(clim$clim1, xy = TRUE) %>% sample_n(n_focal)
      ref <- clim$clim2
      st <- system.time({
            a <- analog_impact(focal, ref, max_clim = NULL, max_geog = 3, k = 20)
      })
      return(st[["elapsed"]])
      # cat("impact: ", st[["elapsed"]], "\n")
}

benchmark_availability <- function(clim, n_focal = 100) {
      # single-era analog availability: few focals, many refs
      focal <- as.data.frame(clim$clim1, xy = TRUE) %>% sample_n(n_focal)
      ref <- clim$clim1
      st <- system.time({
            a <- analog_availability(focal, ref, max_geog = 3, max_clim = .5)
      })
      return(st[["elapsed"]])
      # cat("availability: ", st[["elapsed"]], "\n")
}

benchmark_intensity <- function(clim, n_focal = 100) {
      # analog intensity: few focals, many refs
      focal <- as.data.frame(clim$clim1, xy = TRUE) %>% sample_n(n_focal)
      ref <- clim$clim2
      st <- system.time({
            a <- analog_intensity(focal, ref, max_geog = 3, max_clim = .5, weight = "inverse_clim")
      })
      return(st[["elapsed"]])
      # cat("intensity: ", st[["elapsed"]], "\n")
}

# benchmark_dissimlarity <- function(clim) {
#       # wall-to-wall climate dissimilarity for a single focal site
#       focal <- as.data.frame(clim$clim1, xy = TRUE) %>% sample_n(1)
#       ref <- clim$clim2
#       st <- system.time({
#             a <- find_analogs(
#                   focal = focal,
#                   ref = ref,
#                   max_dist = NULL,
#                   max_clim = NULL,
#                   weight = "inverse_clim",
#                   fun = "all"
#             )
#       })
#       cat("dissimilarity: ", st[["elapsed"]], "\n")
# }

run_benchmarks <- function(
            clim = simulate_climate_rasters(1000)
) {

      n_focal <- c(100, 300, 1000, 3000)

      bm <- function(n){
            data.frame(
                  n = n,
                  velocity = benchmark_velocity(clim, n),
                  impact = benchmark_impact(clim, n),
                  availability = benchmark_availability(clim, n),
                  intensity = benchmark_intensity(clim, n)
            )
      }

      do.call(rbind, lapply(n_focal, bm))
}
