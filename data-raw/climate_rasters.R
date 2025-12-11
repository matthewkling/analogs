library(raster)
library(hydro)
library(tidyverse)
library(terra)

bb <- extent(-115, -95, 30, 50)


# derive variables -------------------

# historic
f <- list.files("/Volumes/T7/CHELSA/v2/raw/", full.names = T)
f <- f[grepl("tas|pr", f)]
hst <- f %>%
      stack() %>%
      crop(bb) %>%
      hydro(already_latlong = TRUE, ncores = 8)
layers <- tolower(names(hst))

# future
f <- list.files("/Volumes/T7/CHELSA/v2/cmip/", full.names = T)
f <- f[grepl("tas|pr", f) &
             grepl("2071", f) &
             grepl("ssp585", f) &
             grepl("gfdl-esm4", f)]
fut <- f %>%
      stack() %>%
      crop(bb) %>%
      hydro(already_latlong = TRUE, ncores = 8)




# munge -----------------

# select focal variables
v <- c("CWD", "AET")
clim1 <- rast(hst[[v]])
clim2 <- rast(fut[[v]])

# scale data
scales <- list()
for(var in names(clim1)){
      x <- values(clim1[[var]])

      m <- mean(x, na.rm = T)
      s <- sd(x, na.rm = T)

      clim1[[var]] <- (clim1[[var]] - m) / s
      clim2[[var]] <- (clim2[[var]] - m) / s

      scales[[var]]$mean <- m
      scales[[var]]$sd <- s
}

# export
writeRaster(
      clim1,
      "inst/extdata/historic_climate.tif",
      datatype = "INT2S",
      scale = 0.01,
      gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2"),
      overwrite = TRUE
)
writeRaster(
      clim2,
      "inst/extdata/future_climate.tif",
      datatype = "INT2S",
      scale = 0.01,
      gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2"),
      overwrite = TRUE
)



