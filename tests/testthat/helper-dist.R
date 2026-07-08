xdist <- function(focal, ref, type = "env") {

      # variable selection
      i <- switch(type,
                  geog = 1:2,
                  env = 3:ncol(focal),
                  lonlat = 1:2)

      f <- focal[, i, drop=FALSE]
      r <- ref[, i, drop=FALSE]

      if (type == "lonlat") {
            # compute haversine distance matrix in km
            R <- 6371.0088
            # lon/lat in degrees
            to_rad <- pi / 180
            f_lon <- f[,1] * to_rad
            f_lat <- f[,2] * to_rad
            r_lon <- r[,1] * to_rad
            r_lat <- r[,2] * to_rad

            out <- matrix(NA_real_, nrow(f), nrow(r))
            for (i in seq_len(nrow(f))) {
                  dlon <- r_lon - f_lon[i]
                  dlat <- r_lat - f_lat[i]
                  a <- sin(dlat/2)^2 +
                        cos(f_lat[i]) * cos(r_lat) * sin(dlon/2)^2
                  out[i,] <- 2 * R * asin(pmin(1, sqrt(a)))
            }
            return(out)
      }

      # DEFAULT: Euclidean for projected or environmental dims
      return(as.matrix(dist(rbind(f, r)))[1:nrow(f), nrow(f) + (1:nrow(r))])
}


sim_test_data <- function(seed = 123, nref = 200, nfocal = 20, lonlat = FALSE){
      set.seed(seed)

      d <- list(ref = matrix(rnorm(nref*4), ncol = 4),
                focal = matrix(rnorm(nfocal*4),  ncol = 4))

      if(lonlat){
            d$ref[, 1] <- runif(nref, -180, 180)
            d$focal[, 1] <- runif(nfocal, -180, 180)
            d$ref[, 2] <- runif(nref, -90, 90)
            d$focal[, 2] <- runif(nfocal, -90, 90)
      }

      d
}
