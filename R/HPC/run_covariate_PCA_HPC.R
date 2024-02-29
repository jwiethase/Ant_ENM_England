rm(list = ls())
library(tidyverse)
library(tidyterra)
library(terra)
library(INLA)
library(mgcv)
library(data.table)

setwd('/users/jhw538/scratch/ant_modelling')
source('source/misc_functions.R')

args <- commandArgs(trailingOnly=TRUE)
job <- as.integer(args[1])

# Load data ----------------------------------------------------
ROI <- vect('data/ROI_kmproj.shp')
clim_finer <- rast("data/clim_finer_30m.tif")
clim_topo <- rast("data/clim_topo_30m.tif")
forest_stack <- rast("data/forest_stack_30m.tif") %>% 
      tidyterra::select(-mean_height_VOM)

eigen_threshold <- 0.85

forest_stack$cover_VOM_log <- log(forest_stack$cover_VOM+1)
forest_stack$perc09_height_VOM_log <- log(forest_stack$perc09_height_VOM+1)
forest_stack$sd_height_VOM_log <- log(forest_stack$sd_height_VOM+1)

# PCA ----------------------------------------------------
## Climate+Topo PCA ----------------------------------------------------
# Run the PCA
rpc_clim_topo <- terra::prcomp(clim_topo, center = T, scale. = T, maxcell = 2e5)
summary(rpc_clim_topo)

nclim_topo_included <- sum(rpc_clim_topo$sdev >= eigen_threshold)

clim_topo_PC_stack <- predict(clim_topo, rpc_clim_topo) %>% 
      subset(1:nclim_topo_included) %>% 
      c(., clim_finer$lat_raster)
names(clim_topo_PC_stack) <- c(paste0("clim_topo_PC", 1:nclim_topo_included), "lat_raster")

print("climtopo PCA done")

## Forest PCA ----------------------------------------------------
forest_stack_sub <- forest_stack %>%
      tidyterra::select(cover_VOM_log, perc09_height_VOM_log, sd_height_VOM_log, distance_ancient) 

forest_stack_PCA <- forest_stack_sub %>% 
      tidyterra::select(cover_VOM_log, perc09_height_VOM_log, sd_height_VOM_log) %>% 
      mask(subst(forest_stack$forest_mask_buff, 0, NA)) %>% 
      c(forest_stack_sub$distance_ancient)

rpc_forest <- terra::prcomp(forest_stack_PCA, center = T, scale. = T, maxcell = 2e5)
summary(rpc_forest)
forest_PC_stack <- predict(forest_stack_sub, rpc_forest) %>% 
      c(., forest_stack$forest_mask_buff)
names(forest_PC_stack) <- gsub("PC", "forest_PC", names(forest_PC_stack))

print("forest PCA done")

# Export ----------------------------------------------------
writeRaster(clim_topo_PC_stack, paste0("data/clim_topo_PCA_30m.tif"), overwrite=TRUE)
writeRaster(forest_PC_stack, paste0("data/forest_PCA_30m.tif"), overwrite=TRUE)


