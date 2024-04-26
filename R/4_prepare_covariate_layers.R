# HEADER --------------------------------------------
#
# Author: Joris Wiethase
# Email: j.wiethase@gmail.com
# 
# Script Description:  
# Prepare the covariate layers needed for the suitability models

rm(list = ls())
library(tidyverse)
library(tidyterra)
library(terra)
source('source/misc_functions.R')

ROI_27700 <- vect('spatial_other/ROI_outline_27700.shp') 
ROI <- ROI_27700 %>% 
      terra::project(crs(km_proj))

# Climate ----------------------------------------------------
## Soil climate layers from Lembrechts et al. 2022  ----------------------------------------------------
soil_temp <- rast('covariates/raw/soil_temp_full.tif') %>% 
      terra::project(crs(km_proj)) %>% 
      crop(ROI) %>% 
      mask(ROI) %>% 
      tidyterra::select(-SBIO8_Mean_Temperature_of_Wettest_Quarter, -SBIO9_Mean_Temperature_of_Driest_Quarter)

## CEDA UK climate  ----------------------------------------------------
## Climate downloaded from: https://data.ceda.ac.uk/badc/ukmo-hadobs/data/insitu/MOHC/HadOBS/HadUK-Grid/v1.2.0.ceda/1km
## Total rain of hottest and coldest quarter, maximum dryspell of hottest quarter of each of 4 years, 
## aggregated using 90th percentile
## Median temperature of hottest and coldest quarter, maximum dryspell of hottest quarter of each of 4 years, 
## aggregated using 90th percentile
CEDA_clim_stack <- rast("covariates/raw/CEDA_clim_stack.tif") %>% 
      terra::project(crs(km_proj)) %>% 
      crop(ROI) %>% 
      mask(ROI)

## Latitude raster ----------------------------------------------------
## Simple raster of latitude, helps with model fitting
lat_raster <- rast(ext = ext(ROI), crs = crs(km_proj), res = 1)
lat_raster[] <- rep(terra::yFromRow(lat_raster, 1:nrow(lat_raster)), each = ncol(lat_raster))

lat_raster <- lat_raster %>% 
      mask(ROI) %>% 
      terra::resample(CEDA_clim_stack)
names(lat_raster) <- 'lat_raster'

## Final climate stack ----------------------------------------------------
climate_stack_1km <- c(soil_temp, CEDA_clim_stack, lat_raster)

# Forest ----------------------------------------------------
## Forest characteristics ----------------------------------------------------
## Aggregated from 1m VOM Lidar data layer, Environment Agency
cover_VOM_30m <- rast('covariates/raw/cover_VOM_30m_27700.tif') %>%    
      terra::project(crs(km_proj)) %>%
      mask(ROI)

perc09_height_VOM_30m <- rast('covariates/raw/perc09_VOM_30m_27700.tif') %>% 
      terra::project(crs(km_proj)) %>%
      mask(ROI)

sd_height_VOM_30m <- rast('covariates/raw/sd_VOM_30m_27700.tif') %>%    
      terra::project(crs(km_proj)) %>%
      mask(ROI)

mean_height_VOM_30m <- rast('covariates/raw/mean_VOM_30m_27700.tif') %>%    
      terra::project(crs(km_proj)) %>%
      mask(ROI)

## Distance to ancient woodland  ----------------------------------------------------
## Raster created from Natural England data layer
ancient_30m <- rast('covariates/raw/ancient_raster_30m.tif')
distance_ancient_30m <- ancient_30m %>%
      mask(ROI) %>%  
      terra::distance(target = 1, unit = "m") %>% 
      +1 %>%
      log() 

## Forest dummy layer ----------------------------------------------------
## Create a mask layer of forest patches. Buffer to include ants at the edges
forest_mask <- cover_VOM_30m
values(forest_mask) <- ifelse(values(forest_mask) < 0.3, NA, 1)

## Create a distance to forest patch layer
distance_forest <- forest_mask %>% 
      terra::distance(target = 1, unit = "m") 

## How far away were most ant records found?
sporadic <- vect('species_data/processed_shp/sporadic_combined.shp')
vals <- terra::extract(distance_forest, sporadic, ID = F)
distance_upper <- quantile(vals$cover_VOM, probs = 0.95)

## 95% of nests are within 67.1 meters of the forest edge. Round up and use as buffer
forest_mask_buff <- buffer(forest_mask, width = 0.07, background = 0)
values(forest_mask_buff) <- ifelse(values(forest_mask_buff) == TRUE, 1, 0)
forest_mask_buff <- mask(forest_mask_buff, ROI)

# How many nests are outside this buffer?
nests_captured <- terra::extract(forest_mask_buff, sporadic, ID = F)
table(nests_captured) # 23 fall outside the buffered zone, 525 are inside

## Final forest stack ----------------------------------------------------
forest_stack_30m <- c(distance_ancient_30m, 
                      cover_VOM_30m, perc09_height_VOM_30m, sd_height_VOM_30m, mean_height_VOM_30m)
names(forest_stack_30m) <- c('distance_ancient', 
                           'cover_VOM', 'perc09_height_VOM', 'sd_height_VOM', 'mean_height_VOM')

# Make coarser stacks for local testing
forest_stack_300m <- forest_stack_30m %>% terra::aggregate(fact = 10, fun = 'median', threads = T)

forest_mask_300m <- ifel(forest_stack_300m$cover_VOM > 0.3, 1, NA)
forest_mask_buff_300m <- buffer(forest_mask_300m, width = 0.4, background = 0)
values(forest_mask_buff_300m) <- ifelse(values(forest_mask_buff_300m) == TRUE, 1, 0)
names(forest_mask_buff_300m) <- "forest_mask_buff"
nests_captured_300m <- terra::extract(forest_mask_buff_300m, sporadic, ID = F)
table(nests_captured_300m) # 22 fall outside buffer

# Topography ----------------------------------------------------
# Derived from NASA DEM layer
NASA_dem_30m <- rast('covariates/raw/nasa_dem_30m.tif') %>% 
      crop(ROI_27700) %>% # Use original projection, as elevation units have to be same as map units
      mask(ROI_27700) 
      
aspect_30m <- terrain(NASA_dem_30m, v="aspect")
northness_30m <- cos(aspect_30m * pi / 180)
eastness_30m <- sin(aspect_30m * pi / 180)

slope_30m <- terrain(NASA_dem_30m, v="slope")
hillshade_30m <- shade(slope_30m, aspect_30m, angle=45)

topo_stack_30m <- c(northness_30m, eastness_30m, hillshade_30m, slope_30m) %>%    
      terra::project(crs(km_proj)) 
      
names(topo_stack_30m) <- c('northness', 'eastness', 'hillshade', 'slope')

topo_stack_300m <- topo_stack_30m %>% terra::aggregate(fact = 10, fun = 'median', threads = T)

# Effort raster ----------------------------------------------------
effort_lgcp_vector <- vect("spatial_other/effort_lgcp_10km.shp") 
effort_integrated_vector <- vect("spatial_other/effort_integrated_10km.shp") 

rast_OS_grid <- rast(ext = ext(effort_lgcp_vector), crs = crs(effort_lgcp_vector), res = 10, vals = 1)
    
effort_rast_lgcp_10km <- rasterize(effort_lgcp_vector, rast_OS_grid, field = "days_sampl") %>% 
      crop(ROI) %>% 
      mask(ROI) %>% 
      +1 %>% 
      log() 

effort_rast_integrated_10km <- rasterize(effort_integrated_vector, rast_OS_grid, field = "days_sampl") %>% 
      crop(ROI) %>% 
      mask(ROI) %>% 
      +1 %>% 
      log()

# Export ----------------------------------------------------
writeVector(ROI, "spatial_other/ROI_kmproj.shp", overwrite=TRUE)

writeRaster(climate_stack_1km, "covariates/processed/climate_stack_1km.tif", overwrite=TRUE)

writeRaster(topo_stack_30m, "covariates/processed/topo_stack_30m.tif", overwrite=TRUE)
writeRaster(topo_stack_300m, "covariates/processed/topo_stack_300m.tif", overwrite=TRUE)

writeRaster(forest_stack_30m, "covariates/processed/forest_stack_30m.tif", overwrite=TRUE)
writeRaster(forest_stack_300m, "covariates/processed/forest_stack_300m.tif", overwrite=TRUE)

writeRaster(forest_mask_buff, "covariates/processed/forest_mask_buff_30m.tif", overwrite=TRUE)
writeRaster(forest_mask_buff_300m, "covariates/processed/forest_mask_buff_300m.tif", overwrite=TRUE)

writeRaster(effort_rast_lgcp_10km, "covariates/processed/effort_rast_lgcp_10km.tif", overwrite=TRUE)
writeRaster(effort_rast_integrated_10km, "covariates/processed/effort_rast_integrated_10km.tif", overwrite=TRUE)

writeRaster(distance_forest, "covariates/processed/distance_forest_30m.tif", overwrite=TRUE)

