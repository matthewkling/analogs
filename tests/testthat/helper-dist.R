xdist <- function(focal, ref, type = "clim"){
      i <- switch(type, geog = 1:2, clim = 3:ncol(focal))
      focal <- focal[, i]
      ref <- ref[, i]
      as.matrix(dist(rbind(focal, ref)))[1:nrow(focal), nrow(focal) + (1:nrow(ref))]
}

sim_test_data <- function(seed = 123){
      set.seed(seed)
      list(ref = matrix(rnorm(200), ncol = 4),   # 50 x 4 (2 geog, 2 clim)
           focal = matrix(rnorm(20),  ncol = 4))   # 5 x 4
}
