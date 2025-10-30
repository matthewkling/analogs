## internal dev functions for benchmarking -- will not be included in final production package ##


# simulate reasonably realistic data with spatiotemporal structure -------------

simulate_climate_rasters <- function(d = 500){
      v1 <- matrix(rep(1:d, each = d), d) %>% "/"(40) %>% sin()
      v1 <- matrix(rep(1:d, d), d) %>% "/"(20) %>% cos() %>% "+"(v1) %>% rast()
      v2 <- matrix(rep(1:d, d), d) %>% rast()
      clim1 <- c(v1, v2) %>% setNames(c("t", "p"))
      clim1[] <- scale(clim1[])

      clim2 <- clim1
      clim2[[1]] <- clim2[[1]] + .5
      clim2[[2]] <- clim2[[2]] * .5

      focal <- as.data.frame(clim1, xy = TRUE) %>% sample_n(100)
      colnames(focal) <- c("x", "y", "t", "p")
}



# load example climate data ---------------------------------

load_climate_rasters <- function(){

      require(terra)
      require(tidyverse)

      scale_rast <- function(x) (x - mean(values(x), na.rm = T)) / sd(values(x), na.rm = T)

      # prep climate data

      bb <- ext(-114, -108, 43, 49)

      t1 <- rast("/Volumes/T7/CHELSA/v2/raw/CHELSA_bio1_1981-2010_V.2.1.tif") %>% crop(bb)
      t2 <- list.files("/Volumes/T7/CHELSA/v2/cmip/", full.names = T)
      t2 <- t2[grepl("2071", t2) & grepl("126_tas_", t2)] %>% rast() %>% crop(bb) %>% mean()
      t <- c(t1, t2) %>% scale_rast()

      p1 <- rast("/Volumes/T7/CHELSA/v2/raw/CHELSA_bio12_1981-2010_V.2.1.tif") %>% crop(bb)
      p2 <- list.files("/Volumes/T7/CHELSA/v2/cmip/", full.names = T)
      p2 <- p2[grepl("2071", p2) & grepl("126_pr_", p2)] %>% rast() %>% crop(bb) %>% mean()
      p <- c(p1, p2) %>% log() %>% scale_rast()

      clim1 <- c(t[[1]], p[[1]])
      clim2 <- c(t[[2]], p[[2]])

      return(list(clim1 = clim1, clim2 = clim2))
}


# benchmarks ----------------------------------------

benchmark_velocity <- function(clim){

      # climate velocity: few focals, many refs
      focal <- as.data.frame(clim$clim1, xy = TRUE) %>% sample_n(100)
      ref <- clim$clim2
      st <- system.time({
            a <- find_analogs(focal = focal, ref = ref,
                              max_dist = NULL, max_clim = .5, n = 1,
                              weight = "inverse_dist", fun = "nmax")
      })
      cat("velocity: ", st[["elapsed"]], "\n")
}

benchmark_impact <- function(clim){

      # analog impact: dist and geo constraints; few focals, many refs
      focal <- as.data.frame(clim$clim1, xy = TRUE) %>% sample_n(100)
      ref <- clim$clim2
      st <- system.time({
            a <- find_analogs(focal = focal, ref = ref,
                              max_dist = 100, max_clim = .5, n = 20,
                              weight = "inverse_clim", fun = "nmax")
      })
      cat("impact: ", st[["elapsed"]], "\n")
}

benchmark_availability <- function(clim){

      # single-era analog availability: few focals, many refs
      focal <- as.data.frame(clim$clim1, xy = TRUE) %>% sample_n(100)
      ref <- clim$clim1
      st <- system.time({
            a <- find_analogs(focal = focal, ref = ref,
                              # max_dist = 50, max_clim = 1,
                              max_dist = 100, max_clim = .5,
                              weight = "uniform", fun = "count")
      })
      cat("availability: ", st[["elapsed"]], "\n")
}

benchmark_dissimlarity <- function(clim){

      # wall-to-wall climate dissimilarity for a single focal site
      focal <- as.data.frame(clim$clim1, xy = TRUE) %>% sample_n(1)
      ref <- clim$clim2
      st <- system.time({
            a <- find_analogs(focal = focal, ref = ref,
                              max_dist = NULL, max_clim = NULL,
                              weight = "inverse_clim", fun = "all")
      })
      cat("dissimilarity: ", st[["elapsed"]], "\n")
}

run_benchmarks <- function(
            clim # e.g. from load_climate_rasters()
            ){

      benchmark_velocity(clim)
      benchmark_impact(clim)
      benchmark_availability(clim)
      benchmark_dissimlarity(clim)
}
