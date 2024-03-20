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
rsq <- function(x, y) summary(lm(y~x))$r.squared

args <- commandArgs(trailingOnly=TRUE)
job <- as.integer(args[1])

# SET PARAMETERS ---------------------------------------
predictions_resolution_list <- c(0.1)
species_choices = c('Formica rufa', "Formica lugubris")
model_list <- c("allFixed", "fixedRandom")

mult_combs <- crossing(predictions_resolution_list, species_choices, model_list) # 8
comb_values <- mult_combs[job, ]

print(comb_values)

species_choice = comb_values$species_choices
predictions_resolution = comb_values$predictions_resolution_list
model_choice = comb_values$model_list

print(paste('Starting', species_choice))

# DATA FILES ------------------------------------------
if(species_choice == 'Formica lugubris'){
      model_name = '29411.81_all30m_3k_tp_Formica_lugubris_E9_thin100_r50_0.1_s1_0.01'
      load(paste0('model_out/Formica_lugubris/lgcp/', model_name, '.RData'))
}

if(species_choice == 'Formica rufa'){
      model_name = '3867505.72_all30m_3k_tp_Formica_rufa_E7_thin0_r250_0.1_s0.1_0.01'
      load(paste0('model_out/Formica_rufa/lgcp/', model_name, '.RData'))
}

ROI <- vect('data/ROI_kmproj.shp')
forest_stack <- rast('data/forest_stack_30m.tif')

clim_topo_covariates <- rast(paste0('data/6clim_topo_tp_3k.tif'))
clim_topo_covariates$lat_raster <- terra::scale(clim_topo_covariates$lat_raster)

forest_covariates <- rast(paste0('data/2forest_tp_3k.tif')) 

distance_ancient <- forest_stack %>% 
      tidyterra::select(distance_ancient) %>% 
      terra::scale()

forestdummy <- forest_stack %>% 
      tidyterra::select(forest_mask_buff) %>% 
      terra::scale()

effort_rast_10km <- rast('data/effort_rast_lgcp_10km.tif')  %>% 
      terra::scale()

print('Covariate prep done')

# PREDICT -----------------------------------------
forest_buffer_mask <- ifel(forestdummy > 0, 1, NA)

grid_points <- rast(res = predictions_resolution, ext = ext(forestdummy), crs = crs(forestdummy), vals = 1) %>% 
      mask(resample(forest_buffer_mask, .)) %>% 
      as.data.table(xy = T) %>% 
      SpatialPoints(proj4string = CRS(km_proj))

print('Grid points done')

if(model_choice == "fixedRandom"){
      all_preds <- predict(object = model,
                           newdata = grid_points,
                           formula = as.formula(paste0('~ mySPDE + ', fixed_effects_effort)),
                           n.samples = 500,
                           num.threads = 64,
                           seed = 42)
      
      all_preds_raster <- as.data.frame(all_preds) %>%
            dplyr::select(x, y, q0.025, q0.5, q0.975, sd) %>%
            rast(type = 'xyz', crs = crs(km_proj))
      
      writeRaster(all_preds_raster,
                  filename = paste0('model_out/', gsub(' ', '_', species_choice),
                                    '/lgcp/final_allPreds_all', predictions_resolution*1000, 'm_',
                                    model_name, '.tif'), overwrite=TRUE)
      
      print('all_preds_raster done')
}

if(model_choice == "allFixed"){
      fixed_only <- predict(object = model, 
                            newdata = grid_points,
                            formula = as.formula(paste0('~ ', fixed_effects_effort)),
                            n.samples = 500,
                            num.threads = 64,
                            seed = 42)
      
      fixed_only_raster <- as.data.frame(fixed_only) %>% 
            dplyr::select(x, y, q0.025, q0.5, q0.975, sd) %>% 
            rast(type = 'xyz', crs = crs(km_proj))
      
      writeRaster(fixed_only_raster, 
                  filename = paste0('model_out/', gsub(' ', '_', species_choice),
                                    '/lgcp/final_fixedOnly_all', predictions_resolution*1000, 'm_', 
                                    model_name, '.tif'), overwrite=TRUE)
      
      print('fixed_only_raster done')
}

