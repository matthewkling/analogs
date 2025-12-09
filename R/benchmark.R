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
      # data constructed in .sandbox/clim_prep.R
      require(terra)
      clim <- rast(".sandbox/clim.tif")
      clim <- list(clim1 = clim[[1:2]], clim2 = clim[[3:4]])
      clim
}


# benchmarks ----------------------------------------

time <- function(expr) system.time({ expr })[["elapsed"]]

prof <- function(expr){
      Rprof("out", line.profiling = TRUE)
      expr
      Rprof(NULL)
      summaryRprof("out")$by.total
}

benchmark <- function(clim, fun = time, n_focal = 1e5,
                      max_geog = 3, # degrees
                      ...){

      set.seed(123)
      focal <- as.data.frame(clim$clim1, xy = TRUE) %>% sample_n(n_focal)
      ref <- clim$clim2

      # approx deg -> km conversion if applicable (for comparable projected vs lonlat benchmarks)
      dots <- list(...)
      if("coord_type" %in% names(dots) && dots$coord_type == "lonlat") max_geog <- max_geog * 100

      list(
            velocity = fun(analog_velocity(
                  focal, ref,
                  max_clim = .5, k = 1,
                  ...)),

            impact = fun(analog_similarity(
                  focal, ref,
                  max_clim = NULL, max_geog = max_geog, k = 20,
                  ...)),

            availability = fun(analog_availability(
                  focal, ref,
                  max_geog = max_geog, max_clim = .5,
                  ...)),

            intensity = fun(analog_intensity(
                  focal, ref,
                  max_geog = max_geog, max_clim = .5, weight = "inverse_clim",
                  ...))
      )
}

