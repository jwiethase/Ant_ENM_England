rm(list = ls())
library(tidyverse)
library(tidyterra)
library(terra)
source('source/misc_functions.R')

raster_for_projection <- rast('/Users/joriswiethase/Downloads/CEDA_climate_2010-2013_1km/rainfall_hadukgrid_uk_1km_day_20100101-20100131.nc') 
ROI <- vect('spatial_other/ROI_buffered_27700.shp')

list <- list.files("/Users/joriswiethase/Downloads/CEDA_climate_2010-2013_1km", full.names = T)

crop_and_export <- function(filepath){
      img <- rast(filepath) %>% 
            crop(ROI) %>% 
            mask(ROI) %>% 
            terra::project(crs(km_proj))

      writeRaster(img,
                  paste0("/Users/joriswiethase/Library/CloudStorage/GoogleDrive-joris.wiethase@york.ac.uk/My Drive/Work/Ant modelling/Ant_ENM_England/covariates/raw/clim_daily/",
                         basename(filepath),
                         ".tif"),
                  overwrite = T)
      
}

lapply(list, crop_and_export)
