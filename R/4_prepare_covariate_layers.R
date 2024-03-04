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
writeRaster(lat_raster, 'test2.tif')
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
forest_mask_buff <- buffer(forest_mask, width = 0.15, background = 0)
values(forest_mask_buff) <- ifelse(values(forest_mask_buff) == TRUE, 1, 0)

# Check if this captures most ant presence records
sporadic <- read.csv('species_data/processed_csv/sporadic_combined.csv') %>% 
      dplyr::select(x, y, species)
exhaustive <- read.csv('species_data/processed_csv/exhaustive_combined.csv') %>% 
      dplyr::select(x, y, species)
all_spatial <- vect(rbind(exhaustive, sporadic), crs = crs(km_proj), geom=c("x", "y"))

all_spatial$val <- terra::extract(forest_mask_buff, all_spatial)[, 2]
table(all_spatial$val) # Most ants are within the forest mask

sporadic$val <- terra::extract(forest_mask_buff, sporadic %>% 
                                     vect(crs = crs(km_proj), geom=c("x", "y")))[, 2]
table(sporadic$val) 

## Final forest stack ----------------------------------------------------
forest_stack_30m <- c(distance_ancient_30m, 
                      cover_VOM_30m, perc09_height_VOM_30m, sd_height_VOM_30m, mean_height_VOM_30m,
                      forest_mask_buff)
names(forest_stack_30m) <- c('distance_ancient', 
                           'cover_VOM', 'perc09_height_VOM', 'sd_height_VOM', 'mean_height_VOM',
                           'forest_mask_buff')

# Make coarser stack for local testing
forest_stack_300m <- forest_stack_30m %>% terra::aggregate(fact = 10, fun = 'median', threads = T)
forest_stack_300m$forest_mask_buff <- subst(forest_stack_300m$forest_mask_buff, 0.5, 1) # Aggregation introduced 0.5, convert

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
effort_vector <- vect("spatial_other/effort_10km.shp") 
rast_OS_grid <- rast(ext = ext(effort_vector), crs = crs(effort_vector), res = 10, vals = 1)
    
effort_rast_10km <- rasterize(effort_vector, rast_OS_grid, field = "days_sampl") %>% 
      crop(ROI) %>% 
      mask(ROI) %>% 
      +1 %>% 
      log() %>% 
      normalise_raster() 

# Export ----------------------------------------------------
writeVector(ROI, "spatial_other/ROI_kmproj.shp", overwrite=TRUE)

writeRaster(climate_stack_1km, "covariates/processed/climate_stack_1km.tif", overwrite=TRUE)

writeRaster(topo_stack_30m, "covariates/processed/topo_stack_30m.tif", overwrite=TRUE)
writeRaster(topo_stack_300m, "covariates/processed/topo_stack_300m.tif", overwrite=TRUE)

writeRaster(forest_stack_30m, "covariates/processed/forest_stack_30m.tif", overwrite=TRUE)
writeRaster(forest_stack_300m, "covariates/processed/forest_stack_300m.tif", overwrite=TRUE)

writeRaster(effort_rast_10km, "covariates/processed/effort_rast_10km.tif", overwrite=TRUE)


