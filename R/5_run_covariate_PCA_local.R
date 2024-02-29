rm(list = ls())
library(tidyverse)
library(tidyterra)
library(terra)
library(factoextra)
library(data.table)

source('source/misc_functions.R')

# Load data ----------------------------------------------------
climate_stack_1km <- rast("covariates/processed/climate_stack_1km.tif")
forest_stack <- rast("covariates/processed/forest_stack_300m.tif") %>% 
      tidyterra::select(-mean_height_VOM)
topo_stack <- rast("covariates/processed/topo_stack_300m.tif")

eigen_threshold <- 0.85

# Rainfall makes probably more sense on log scale (e.g. 0-50mm more meaningful than 1500-1550mm)
climate_stack_1km$median_total_rain_coldest_log <- log(climate_stack_1km$median_total_rain_coldest)
climate_stack_1km$median_total_rain_hottest_log <- log(climate_stack_1km$median_total_rain_hottest)

# VOM layers are heavily skewed; run transformation
forest_stack$cover_VOM_log <- log(forest_stack$cover_VOM+1)
forest_stack$perc09_height_VOM_log <- log(forest_stack$perc09_height_VOM+1)
forest_stack$sd_height_VOM_log <- log(forest_stack$sd_height_VOM+1)

# Slope makes more sense on log scale
topo_stack$slope_log <- log(topo_stack$slope+1)

# Check value distribution of covariates
pdf("PCA_output/figures/climate_densities.pdf", width = 12, height = 7)
par(mfrow = c(3, 6))
lapply(names(climate_stack_1km), function(x){hist(values(climate_stack_1km[[x]], na.rm = T), main = x)})
dev.off()

pdf("PCA_output/figures/forest_densities.pdf", width = 12, height = 7)
par(mfrow = c(3, 3))
lapply(names(forest_stack), function(x){hist(values(forest_stack[[x]], na.rm = T), main = x)})
dev.off()

pdf("PCA_output/figures/topo_densities.pdf", width = 12, height = 7)
par(mfrow = c(3, 2))
lapply(names(topo_stack), function(x){hist(values(topo_stack[[x]], na.rm = T), main = x)})
dev.off()

par(mfrow = c(1, 1))

# PCA ----------------------------------------------------
# Run PCA on covariates that might be correlated, to create components that 
# maximize variance explained. This is to reduce multi-colinearity, and 
# reduce the number of fixed effects that go into the model, for better
# stability.
# Retain groupings for interpretability.

## Climate+Topo PCA ----------------------------------------------------
## Disaggregate climate and stack
clim_finer <- climate_stack_1km %>%
      terra::disagg(fact = c(res(climate_stack_1km)/res(topo_stack))[1]) 

topo_resampled <- topo_stack %>% 
      tidyterra::select(-slope) %>% 
      terra::resample(clim_finer)

clim_topo <- c(clim_finer %>% 
                     tidyterra::select(-lat_raster, -median_total_rain_coldest, -median_total_rain_hottest), 
               topo_resampled)

# Run the PCA
rpc_clim_topo <- terra::prcomp(clim_topo, center = T, scale. = T)
summary(rpc_clim_topo)

pdf("PCA_output/figures/climTopo_PCA_scree.pdf", width = 7, height = 7)
fviz_eig(rpc_clim_topo, addlabels = TRUE)
dev.off()

pdf("PCA_output/figures/climTopo_PCA_variables.pdf", width = 7, height = 7)
fviz_pca_var(rpc_clim_topo, col.var = "black")
dev.off()

# Extract results
clim_topo_cumulative <- rpc_clim_topo$sdev^2 / sum(rpc_clim_topo$sdev^2)
clim_topo_cumulative_df <- data.frame(cumsum(clim_topo_cumulative)) %>% 
      mutate(PCA_comp = paste0("PC", rownames(.)),
             prop_var = clim_topo_cumulative) %>% 
      rename(cum_var = cumsum.clim_topo_cumulative.)

clim_topo_scores_df <- data.frame(rpc_clim_topo$x)

clim_topo_loadings_df <- data.frame(rpc_clim_topo$rotation) %>% 
      mutate(PCA_comp = rownames(.)) %>% 
      dplyr::select(PCA_comp, everything())

nclim_topo_included <- sum(rpc_clim_topo$sdev >= eigen_threshold)

clim_topo_PC_stack <- predict(clim_topo, rpc_clim_topo) %>% 
      subset(1:nclim_topo_included) %>% 
      c(., clim_finer$lat_raster)
names(clim_topo_PC_stack) <- c(paste0("clim_topo_PC", 1:nclim_topo_included), "lat_raster")

## Forest PCA ----------------------------------------------------
## Run the PCA only on non-zero values (there are many of those). These will be mostly be ignored
## later anyways due to the forest/non-forest dummy layers
forest_stack_sub <- forest_stack %>%
      tidyterra::select(cover_VOM_log, perc09_height_VOM_log, sd_height_VOM_log, distance_ancient) 

forest_stack_PCA <- forest_stack_sub %>% 
      tidyterra::select(cover_VOM_log, perc09_height_VOM_log, sd_height_VOM_log) %>% 
      mask(subst(forest_stack$forest_mask_buff, 0, NA)) %>% 
      c(forest_stack_sub$distance_ancient)

rpc_forest <- terra::prcomp(forest_stack_PCA, center = T, scale. = T)
summary(rpc_forest)

pdf("PCA_output/figures/forest_PCA_scree.pdf", width = 7, height = 7)
fviz_eig(rpc_forest, addlabels = TRUE)
dev.off()

pdf("PCA_output/figures/forest_PCA_variables.pdf", width = 7, height = 7)
fviz_pca_var(rpc_forest, col.var = "black")
dev.off()

forest_cumulative <- rpc_forest$sdev^2 / sum(rpc_forest$sdev^2)
forest_cumulative_df <- data.frame(cumsum(forest_cumulative)) %>% 
      mutate(PCA_comp = paste0("PC", rownames(.)),
             prop_var = forest_cumulative) %>% 
      rename(cum_var = cumsum.forest_cumulative.)

forest_scores_df <- data.frame(rpc_forest$x)

forest_loadings_df <- data.frame(rpc_forest$rotation) %>% 
      mutate(PCA_comp = rownames(.)) %>% 
      dplyr::select(PCA_comp, everything())

# PC1 is basically just forest/non-forest.
forest_PC_stack <- predict(forest_stack_sub, rpc_forest) %>% 
      c(., forest_stack$forest_mask_buff)
names(forest_PC_stack) <- gsub("PC", "forest_PC", names(forest_PC_stack))

## Export PCA loadings ----------------------------------------------------
write.csv(clim_topo_cumulative_df, "PCA_output/csv/clim_topo_cumulative_df.csv")
fwrite(clim_topo_scores_df, "PCA_output/csv/clim_topo_scores_df.csv")
write.csv(clim_topo_loadings_df, "PCA_output/csv/clim_topo_loadings_df.csv")

write.csv(forest_cumulative_df, "PCA_output/csv/forest_cumulative_df.csv")
fwrite(forest_scores_df, "PCA_output/csv/forest_scores_df.csv")
write.csv(forest_loadings_df, "PCA_output/csv/forest_loadings_df.csv")

writeRaster(clim_topo_PC_stack, "covariates/processed/clim_topo_PC_stack.tif")
writeRaster(forest_PC_stack, "covariates/processed/forest_PC_stack.tif")

