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

# mask water
land <- rast("/Volumes/T7/CHELSA/v2/raw//CHELSA_ai_1981-2010_V.2.1.tif") %>%
      crop(clim1)
clim1 <- mask(clim1, land)
clim2 <- mask(clim2, land)

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

clim1 <- aggregate(clim1, 5)
clim2 <- aggregate(clim2, 5)

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



