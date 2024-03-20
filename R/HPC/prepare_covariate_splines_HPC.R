rm(list = ls())
library(tidyverse)
library(tidyterra)
library(terra)
library(INLA)
library(mgcv)
library(data.table)

setwd('/users/jhw538/scratch/ant_modelling')
source('source/misc_functions.R')

# args <- commandArgs(trailingOnly=TRUE)
# job <- as.integer(args[1])

# Load data ----------------------------------------------------
ROI <- vect('data/ROI_kmproj.shp')
forest_mask_buff <- rast('data/forest_mask_buff_30m.tif') %>% 
      subst(0, NA)
forest_PC_stack <- rast("data/forest_PCA_30m.tif") %>% 
      mask(forest_mask_buff)
# clim_topo_PC_stack <- rast("data/clim_topo_PCA_30m.tif")

# smoother_list = c('tp')
n_knots = 3

smoother = 'tp'
print(n_knots)

# Make splines ----------------------------------------------------
# ## Climate splines ---------------------------------------------------- 
# clim_topo_df <- as.data.table(clim_topo_PC_stack, xy = TRUE) %>%
#       drop_na()
# 
# for(i in 1:6){
#       clim_topo_df <- prepareMGCVsplines(clim_topo_df, paste0("clim_topo_PC", i), n_knots = n_knots, smooth = smoother)
# }
# 
# print("climtopo splines done")
# 
# clim_topo_df <- clim_topo_df %>% dplyr::select(contains(c("x", "y", "spline", "lat_raster")))
# clim_topo_covariates <- rast(clim_topo_df, type = 'xyz', crs = crs(km_proj))
# 
# writeRaster(clim_topo_covariates, paste0("data/6clim_topo_", smoother, "_", n_knots, "k.tif"), overwrite=TRUE)
# print("climtopo spline stack done")

## Forest splines ----------------------------------------------------
forest_df <- as.data.table(forest_PC_stack, xy = TRUE)

summary(forest_df)

# Create dataset for splines, this is forest only
forest_df_splines <- forest_df %>% 
      dplyr::select(-forest_mask) %>% 
      drop_na()

summary(forest_df_splines)

for(i in 1:2){
      forest_df_splines <- prepareMGCVsplines(forest_df_splines, paste0("forest_PC", i), n_knots = n_knots, smooth = smoother)
}

print("forest splines done")

# Add splines to full data set. Leave the non-forest values as NA
forest_df_splines_sub <- forest_df_splines %>% dplyr::select(contains(c("x", "y", "spline")))

forest_df_merged <- merge(forest_df,
                          forest_df_splines_sub,
                          by = c("x", "y"), all.x = T)
print("forest merge done")
forest_df_merged_sub <- forest_df_merged %>%
                                dplyr::select(contains(c("x", "y", "spline")))

forest_covariates <- rast(forest_df_merged_sub, type = 'xyz', crs = crs(km_proj))

print("forest spline stack done")

writeRaster(forest_covariates, paste0("data/2forest_", smoother, "_", n_knots, "k.tif"), overwrite=TRUE)

