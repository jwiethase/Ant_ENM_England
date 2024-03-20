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
species_choices = c('Formica rufa', "Formica lugubris")

species_choice = species_choices[job]
predictions_resolution = 0.5

print(paste("Starting", species_choice))

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
forest_stack <- rast("data/forest_stack_30m.tif")

clim_topo_covariates <- rast(paste0("data/6clim_topo_tp_3k.tif"))
clim_topo_covariates$lat_raster <- terra::scale(clim_topo_covariates$lat_raster)

forest_covariates <- rast(paste0("data/2forest_tp_3k.tif")) 

distance_ancient <- forest_stack %>% 
      tidyterra::select(distance_ancient) %>% 
      terra::scale()

forestdummy <- forest_stack %>% 
      tidyterra::select(forest_mask_buff) %>% 
      terra::scale()

effort_rast_10km <- rast('data/effort_rast_lgcp_10km.tif')  %>% 
      terra::scale()

print("Covariate prep done")

# VARIABLE IMPORTANCE -----------------------------------------
fixed_effects <- paste(model[["names.fixed"]], collapse = " + ")
fixed_effects_effort <- gsub("effort", "quantile(effort, probs = 0.95, na.rm=T)", fixed_effects)

forest_buffer_mask <- ifel(forestdummy > 0, 1, NA)

grid_points <- rast(res = predictions_resolution, ext = ext(forestdummy), crs = crs(forestdummy), vals = 1) %>% 
      mask(resample(forest_buffer_mask, .)) %>% 
      as.data.table(xy = T) %>% 
      SpatialPoints(proj4string = CRS(km_proj))

all_preds <- predict(object = model,
                     newdata = grid_points,
                     formula = as.formula(paste0('~ ', fixed_effects_effort)),
                     n.samples = 500,
                     num.threads = 64,
                     seed = 42)

all_preds_df <- as.data.table(all_preds) %>%
      dplyr::select(x, y, q0.025, q0.5, q0.975, sd) 

no_forest <- predict(object = model,
                     newdata = grid_points,
                     formula = ~ lat_raster + forestdummy + distance_ancient + 
                           quantile(effort, probs = 0.95, na.rm=T) +
                           clim_topo_PC1_spline1 + clim_topo_PC1_spline2 +
                           clim_topo_PC2_spline1 + clim_topo_PC2_spline2 +
                           clim_topo_PC3_spline1 + clim_topo_PC3_spline2 +
                           clim_topo_PC4_spline1 + clim_topo_PC4_spline2 +
                           clim_topo_PC5_spline1 + clim_topo_PC5_spline2 +
                           clim_topo_PC6_spline1 + clim_topo_PC6_spline2,
                     n.samples = 500,
                     num.threads = 64,
                     seed = 42)

no_forest_df <- as.data.table(no_forest) %>%
      dplyr::select(x, y, q0.025, q0.5, q0.975, sd) 

no_clim_topo <- predict(object = model,
                        newdata = grid_points,
                        formula = ~ lat_raster + forestdummy + distance_ancient + 
                              quantile(effort, probs = 0.95, na.rm=T) +
                              forest_PC1_spline1 + forest_PC1_spline2 +
                              forest_PC2_spline1 + forest_PC2_spline2,
                        n.samples = 500,
                        num.threads = 64,
                        seed = 42)

no_clim_topo_df <- as.data.table(no_clim_topo) %>%
      dplyr::select(x, y, q0.025, q0.5, q0.975, sd) 

forest_structure_contribution <- 1 - rsq(all_preds_df$q0.5, no_forest_df$q0.5)
clim_topo_contribution <- 1 - rsq(all_preds_df$q0.5, no_clim_topo_df$q0.5)

print(forest_structure_contribution)
print(clim_topo_contribution)

print(paste("Forest variable contribution: ", forest_structure_contribution))
print(paste("Climate and topography variable contribution: ", clim_topo_contribution))


