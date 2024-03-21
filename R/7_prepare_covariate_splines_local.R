# HEADER --------------------------------------------
#
# Author: Joris Wiethase
# Email: j.wiethase@gmail.com
# 
# Script Description:  
# Prepares spline base functions for the PCA covariates, to be included in the final SDMs  

rm(list = ls())
library(tidyverse)
library(tidyterra)
library(terra)
library(mgcv)
library(data.table)

source('source/misc_functions.R')

# Load data ----------------------------------------------------
clim_topo_PC_stack <- rast("covariates/processed/clim_topo_PC_stack.tif")
forest_PC_stack <- rast("covariates/processed/forest_PC_stack_30m.tif")
forest_mask_buff <- rast('covariates/processed/forest_mask_buff_30m.tif') %>% 
      subst(0, NA)
forest_PC_stack <- forest_PC_stack %>% 
      mask(forest_mask_buff)

smoother = 'tp'
n_knots = 3

# Make splines ----------------------------------------------------
# ## Climate splines ---------------------------------------------------- 
# clim_topo_df <- as.data.table(clim_topo_PC_stack, xy = TRUE) %>%
#       drop_na()
# 
# for(i in 1:sum(stringr::str_detect(colnames(clim_topo_df), "PC"))){
#       clim_topo_df <- prepareMGCVsplines(clim_topo_df, paste0("clim_topo_PC", i), n_knots = n_knots, smooth = smoother)
# }
# 
# clim_topo_covariates <- rast(clim_topo_df, type = 'xyz', crs = crs(km_proj))

## Forest splines ----------------------------------------------------
forest_df <- as.data.table(forest_PC_stack, xy = TRUE) 

# Create dataset for splines, this is forest only
forest_df_splines <- forest_df %>% 
      dplyr::select(-forest_mask) %>% 
      drop_na()

for(i in 1:sum(stringr::str_detect(colnames(forest_df), "PC"))){
      forest_df_splines <- prepareMGCVsplines(forest_df_splines, paste0("forest_PC", i), n_knots = n_knots, smooth = smoother)
}

# Add splines to full data set. Should be okay to leave the non-forest values as NA
forest_df_splines_sub <- forest_df_splines %>%
      dplyr::select(contains(c("x", "y", "spline")))

forest_df_merged <- merge(forest_df, forest_df_splines_sub, by = c("x", "y"), all.x = T)

forest_df_merged_sub <- forest_df_merged %>%
      dplyr::select(contains(c("x", "y", "spline")))

forest_covariates <- rast(forest_df_merged_sub, type = 'xyz', crs = crs(km_proj))
# Export ----------------------------------------------------
# writeRaster(clim_topo_covariates, paste0("covariates/processed/", sum(stringr::str_detect(names(clim_topo_PC_stack), "PC")), "clim_topo_300m_", smoother, "_", n_knots, "k.tif"), overwrite=TRUE)
writeRaster(forest_covariates, paste0("covariates/processed/", sum(stringr::str_detect(names(forest_PC_stack), "PC")), "forest_30m_", smoother, "_", n_knots, "k.tif"), overwrite=TRUE)

