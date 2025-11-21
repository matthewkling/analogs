## internal dev functions for benchmarking -- will not be included in final production package ##

# simulate reasonably realistic data with spatiotemporal structure -------------

simulate_climate_rasters <- function(d = 500) {
      set.seed(123)
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

      p1 <- rast("/Volumes/T7/CHELSA/v2/raw/CHELSA_bio12_1981-2010_V.2.1.tif") %>%
            crop(bb)
      p2 <- list.files("/Volumes/T7/CHELSA/v2/cmip/", full.names = T)
      p2 <- p2[grepl("2071", p2) & grepl("126_pr_", p2)] %>%
            rast() %>%
            crop(bb) %>%
            terra::mean() %>% # average across months and GCMs
            "*"(12) # convert monthly to annual to match baseline
      p <- c(p1, p2) %>% log() %>% scale_rast()

      clim1 <- c(t[[1]], p[[1]])
      clim2 <- c(t[[2]], p[[2]])

      names(clim1) <- names(clim2) <- c("t", "p")

      return(list(clim1 = clim1, clim2 = clim2))
}


# benchmarks ----------------------------------------

time <- function(expr) system.time({ expr })[["elapsed"]]

prof <- function(expr){
      Rprof("out", line.profiling = TRUE)
      expr
      Rprof(NULL)
      summaryRprof("out")$by.total
}

benchmark <- function(clim, fun = time, n_focal = 1e5, ...){

      set.seed(123)
      focal <- as.data.frame(clim$clim1, xy = TRUE) %>% sample_n(n_focal)
      ref <- clim$clim2

      max_geog <- 3 # degrees
      dots <- list(...)
      if("coord_type" %in% names(dots) && dots$coord_type == "lonlat") max_geog <- max_geog * 100 # km

      list(
            velocity = fun(analog_velocity(
                  focal, ref,
                  max_clim = .5, k = 1,
                  ...)),
            impact = fun(analog_impact(
                  focal, ref,
                  max_clim = NULL, max_geog = 3, k = 20,
                  ...)),
            availability = fun(analog_availability(
                  focal, ref,
                  max_geog = 3, max_clim = .5,
                  ...)),
            intensity = fun(analog_intensity(
                  focal, ref,
                  max_geog = 3, max_clim = .5, weight = "inverse_clim",
                  ...))
      )
}

