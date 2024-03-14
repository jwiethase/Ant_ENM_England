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

forest_stack$cover_VOM_sqrt <- sqrt(forest_stack$cover_VOM)
forest_stack$perc09_height_VOM_sqrt <- sqrt(forest_stack$perc09_height_VOM)
forest_stack$sd_height_VOM_sqrt <- sqrt(forest_stack$sd_height_VOM)

# pdf("figures/clim_topo_densities_30m.pdf", width = 12, height = 12)
# par(mfrow = c(5, 4))
# lapply(names(clim_topo), function(x){hist(values(clim_topo[[x]], na.rm = T), main = x)})
# dev.off()
# 
# pdf("figures/forest_densities_30m.pdf", width = 12, height = 7)
# par(mfrow = c(3, 3))
# lapply(names(forest_stack), function(x){hist(values(forest_stack[[x]] , na.rm = T), main = x)})
# dev.off()
# 
# par(mfrow = c(1, 1))

# PCA ----------------------------------------------------
## Climate+Topo PCA ----------------------------------------------------
# Run the PCA
rpc_clim_topo <- terra::prcomp(clim_topo, center = T, scale. = T, maxcell = 10e6)
summary(rpc_clim_topo)

nclim_topo_included <- sum(rpc_clim_topo$sdev >= eigen_threshold)

clim_topo_cumulative <- rpc_clim_topo$sdev^2 / sum(rpc_clim_topo$sdev^2)
clim_topo_cumulative_df <- data.frame(cumsum(clim_topo_cumulative)) %>% 
      mutate(PCA_comp = paste0("PC", rownames(.)),
             prop_var = clim_topo_cumulative) %>% 
      rename(cum_var = cumsum.clim_topo_cumulative.)

clim_topo_scores_df <- data.frame(rpc_clim_topo$x)

clim_topo_loadings_df <- data.frame(rpc_clim_topo$rotation) %>% 
      mutate(PCA_comp = rownames(.)) %>% 
      dplyr::select(PCA_comp, everything())

clim_topo_PC_stack <- predict(clim_topo, rpc_clim_topo) %>% 
      subset(1:nclim_topo_included) %>% 
      c(., clim_finer$lat_raster)
names(clim_topo_PC_stack) <- c(paste0("clim_topo_PC", 1:nclim_topo_included), "lat_raster")

print("climtopo PCA done")

## Forest PCA ----------------------------------------------------
forest_stack_sub <- forest_stack %>%
      tidyterra::select(cover_VOM_sqrt, perc09_height_VOM_sqrt, sd_height_VOM_sqrt, distance_ancient) 

forest_stack_PCA <- forest_stack_sub %>% 
      tidyterra::select(cover_VOM_sqrt, perc09_height_VOM_sqrt, sd_height_VOM_sqrt) %>% 
      mask(subst(forest_stack$forest_mask_buff, 0, NA)) %>% 
      c(forest_stack_sub$distance_ancient)

rpc_forest <- terra::prcomp(forest_stack_PCA, center = T, scale. = T, maxcell = 10e6)
summary(rpc_forest)

forest_cumulative <- rpc_forest$sdev^2 / sum(rpc_forest$sdev^2)
forest_cumulative_df <- data.frame(cumsum(forest_cumulative)) %>% 
      mutate(PCA_comp = paste0("PC", rownames(.)),
             prop_var = forest_cumulative) %>% 
      rename(cum_var = cumsum.forest_cumulative.)

forest_scores_df <- data.frame(rpc_forest$x)

forest_loadings_df <- data.frame(rpc_forest$rotation) %>% 
      mutate(PCA_comp = rownames(.)) %>% 
      dplyr::select(PCA_comp, everything())

forest_PC_stack <- predict(forest_stack_sub, rpc_forest) %>% 
      c(., forest_stack$forest_mask_buff)
names(forest_PC_stack) <- gsub("PC", "forest_PC", names(forest_PC_stack))

print("forest PCA done")

# Export ----------------------------------------------------
write.csv(clim_topo_cumulative_df, "data/HPC_clim_topo_cumulative_df.csv")
fwrite(clim_topo_scores_df, "data/HPC_clim_topo_scores_df.csv")
write.csv(clim_topo_loadings_df, "data/HPC_clim_topo_loadings_df.csv")

write.csv(forest_cumulative_df, "data/HPC_forest_cumulative_df.csv")
fwrite(forest_scores_df, "data/HPC_forest_scores_df.csv")
write.csv(forest_loadings_df, "data/HPC_forest_loadings_df.csv")

writeRaster(clim_topo_PC_stack, paste0("data/clim_topo_PCA_30m.tif"), overwrite=TRUE)
writeRaster(forest_PC_stack, paste0("data/forest_PCA_30m.tif"), overwrite=TRUE)


