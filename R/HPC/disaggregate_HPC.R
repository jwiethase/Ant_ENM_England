rm(list = ls())
library(tidyverse)
library(tidyterra)
library(terra)
setwd('/users/jhw538/scratch/ant_modelling')
source('source/misc_functions.R')

# Load data ----------------------------------------------------
climate_stack_1km <- rast("data/climate_stack_1km.tif")
topo_stack <- rast("data/topo_stack_30m.tif")

## Climate+Topo PCA ----------------------------------------------------
climate_stack_1km$median_total_rain_coldest_log <- log(climate_stack_1km$median_total_rain_coldest+1)
climate_stack_1km$median_total_rain_hottest_log <- log(climate_stack_1km$median_total_rain_hottest+1)
topo_stack$slope_log <- log(topo_stack$slope+1)

clim_finer <- climate_stack_1km %>%
      terra::disagg(fact = c(res(climate_stack_1km)/res(topo_stack))[1]) 

topo_resampled <- topo_stack %>% 
      tidyterra::select(-slope) %>% 
      terra::resample(clim_finer)

clim_topo <- c(clim_finer %>% 
                     tidyterra::select(-lat_raster, -median_total_rain_coldest, -median_total_rain_hottest), 
               topo_resampled)
print("climtopo stack done")
terra::writeRaster(clim_finer, file = "data/clim_finer_30m.tif", steps = 20, overwrite=TRUE)
terra::writeRaster(clim_topo, file = "data/clim_topo_30m.tif", steps = 20, overwrite=TRUE)
