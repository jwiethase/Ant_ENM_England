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
forest_PC_stack <- rast("data/forest_PCA_30m.tif")
clim_topo_PC_stack <- rast("data/clim_topo_PCA_30m.tif")

smoother_list = c('tp', 'ps', 'cr')
n_knots_list = c(3, 4)
mult_combs <- crossing(smoother_list, n_knots_list) %>% filter(!(smoother_list == 'ps' & n_knots_list == 3))
comb_values <- mult_combs[job, ]
print(comb_values)
n_knots = comb_values$n_knots_list
smoother = comb_values$smoother_list

# Make splines ----------------------------------------------------
## Climate splines ---------------------------------------------------- 
clim_topo_df <- as.data.frame(clim_topo_PC_stack, xy = TRUE) %>%
      drop_na()

for(i in 1:6){
      clim_topo_df <- prepareMGCVsplines(clim_topo_df, paste0("clim_topo_PC", i), n_knots = n_knots, smooth = smoother)
}

print("climtopo splines done")

# Scale covariates not in splines
clim_topo_df$lat_raster <- scale(clim_topo_df$lat_raster)
clim_topo_df <- clim_topo_df %>%
                                dplyr::select(contains(c("x", "y", "spline", "lat_raster")))
clim_topo_covariates <- rast(clim_topo_df, type = 'xyz', crs = crs(km_proj))

writeRaster(clim_topo_covariates, paste0("data/6clim_topo_", smoother, "_", n_knots, "k.tif"), overwrite=TRUE)
print("climtopo spline stack done")

## Forest splines ----------------------------------------------------
forest_df <- as.data.frame(forest_PC_stack, xy = TRUE) %>%
      drop_na()

# Create dataset for splines, this is forest only
forest_df_splines <- forest_df %>% 
      filter(forest_mask_buff == 1)

for(i in 2:3){
      forest_df_splines <- prepareMGCVsplines(forest_df_splines, paste0("forest_PC", i), n_knots = n_knots, smooth = smoother)
}
#write.csv(forest_df_splines, paste0("data/", sum(stringr::str_detect(colnames(forest_df), "PC")), "forestDFsplines_", smoother, "_", n_knots, "k.csv"))
print("forest splines done")
# Add splines to full data set. Should be okay to leave the non-forest values as NA
forest_df_splines_sub <- forest_df_splines %>%
                                dplyr::select(contains(c("x", "y", "spline")))

forest_df_merged <- merge(forest_df,
                          forest_df_splines_sub,
                          by = c("x", "y"), all.x = T)
print("forest merge done")
forest_df_merged_sub <- forest_df_merged %>%
                                dplyr::select(contains(c("x", "y", "PC2_spline", "PC3_spline", "forest_mask_buff")))

forest_covariates <- rast(forest_df_merged_sub, type = 'xyz', crs = crs(km_proj))

print("forest spline stack done")

# Export ----------------------------------------------------
writeRaster(forest_covariates, paste0("data/2forest_", smoother, "_", n_knots, "k.tif"), overwrite=TRUE)

