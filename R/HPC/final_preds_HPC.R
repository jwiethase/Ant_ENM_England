library(inlabru)
library(terra)
library(tidyterra)
library(tidyverse)
library(sf)
library(sp)
library(data.table)
library(landscapemetrics)
setwd('/users/jhw538/scratch/ant_modelling')
source('source/misc_functions.R')

args <- commandArgs(trailingOnly=TRUE)
job <- as.integer(args[1])

# SET PARAMETERS ---------------------------------------
predictions_resolution <- 0.03
species_choices = c('Formica rufa', "Formica lugubris")
slice_ids <- 1:20

mult_combs <- crossing(species_choices, slice_ids)

comb_values <- mult_combs[job, ]

print(comb_values)

species_choice = comb_values$species_choices
slice_id <- comb_values$slice_ids

print(paste('Starting', species_choice))

# DATA FILES ------------------------------------------
if(species_choice == 'Formica lugubris'){
      model_name = '13564.52_all30m_3k_tp_Formica_lugubris_E11_thin0_r250_0.5_s1_0.01'
      load(paste0('model_out/Formica_lugubris/lgcp/', model_name, '.RData'))
}

if(species_choice == 'Formica rufa'){
      model_name = '8478334.22_all30m_3k_tp_Formica_rufa_E5_thin100_r200_0.1_s0.1_0.01'
      load(paste0('model_out/Formica_rufa/lgcp/', model_name, '.RData'))
}

ROI <- vect('data/ROI_kmproj.shp')
forest_stack <- rast("data/forest_stack_30m.tif")
distance_forest <- rast("data/distance_forest_30m.tif")
names(distance_forest) <- "distance_forest"
FE_managed <- vect('data/Forestry_England_managed_forest.shp') %>% 
      terra::project(crs(km_proj))

clim_topo_covariates <- rast(paste0("data/6clim_topo_tp_3k.tif"))
clim_topo_covariates$lat_raster <- terra::scale(clim_topo_covariates$lat_raster)

forest_covariates <- rast(paste0("data/2forest_30m_tp_3k.tif")) 

effort_rast_10km <- rast('data/effort_rast_lgcp_10km.tif')  %>% 
      terra::scale()

forestdummy <- rast("data/forest_mask_buff_30m.tif") %>%
      mask(ROI)
names(forestdummy) <- "forest_mask_buff"

print('Covariate prep done')

# PREDICT -----------------------------------------
forest_buffer_mask <- ifel(forestdummy > 0, 1, NA)

grid_points <- rast(res = predictions_resolution, ext = ext(forestdummy), crs = crs(forestdummy), vals = 1) %>% 
      mask(resample(forest_buffer_mask, .)) %>% 
      as.data.table(xy = T) %>% 
      SpatialPoints(proj4string = CRS(km_proj))

lat_min <- ext(forestdummy)$ymin
lat_max <- ext(forestdummy)$ymax
slice_height <- (lat_max - lat_min) / length(slice_ids)

# Generate slices based on latitude and select the one for this job
slices <- list()
for (i in slice_ids) {
      lat_lower <- lat_min + (i - 1) * slice_height
      lat_upper <- lat_min + i * slice_height
      current_slice <- grid_points[grid_points@coords[, 2] >= lat_lower & grid_points@coords[, 2] < lat_upper, ]
      slices[[i]] <- current_slice
}

current_slice_points <- slices[[slice_id]]

print('Grid points done')
gc()
fixed_only <- predict(object = model, 
                      newdata = current_slice_points,
                      formula = as.formula(paste0('~ ', fixed_effects_effort)),
                      n.samples = 250,
                      num.threads = 64,
                      seed = 42)

fixed_only_raster <- as.data.frame(fixed_only) %>% 
      dplyr::select(x, y, q0.025, q0.5, q0.975, sd) %>% 
      rast(type = 'xyz', crs = crs(km_proj))

writeRaster(fixed_only_raster, 
            filename = paste0('model_out/', gsub(' ', '_', species_choice),
                              '/lgcp/', gsub(' ', '_', species_choice), '_final_LGCP_fixedOnly_all', 
                              predictions_resolution*1000, 'm_', 
                              "slice_", slice_id, '.tif'), overwrite=TRUE)

print('fixed_only_raster done')