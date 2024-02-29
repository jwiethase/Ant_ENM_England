rm(list = ls())
library(tidyverse)
library(tidyterra)
library(terra)
library(mgcv)

source('source/misc_functions.R')

# Load data ----------------------------------------------------
clim_topo_PC_stack <- rast("covariates/processed/clim_topo_PC_stack.tif")
forest_PC_stack <- rast("covariates/processed/forest_PC_stack.tif")

smoother = 'tp'
n_knots = 4

# Make splines ----------------------------------------------------
## Climate splines ---------------------------------------------------- 
clim_topo_df <- as.data.table(clim_topo_PC_stack, xy = TRUE) %>%
      drop_na()

for(i in 1:sum(stringr::str_detect(colnames(clim_topo_df), "PC"))){
      clim_topo_df <- prepareMGCVsplines(clim_topo_df, paste0("clim_topo_PC", i), n_knots = n_knots, smooth = smoother)
}

# Scale covariates not in splines
clim_topo_df$lat_raster <- scale(clim_topo_df$lat_raster)
clim_topo_covariates <- rast(clim_topo_df, type = 'xyz', crs = crs(km_proj))

## Forest splines ----------------------------------------------------
forest_df <- as.data.table(forest_PC_stack, xy = TRUE) %>%
      drop_na()

# Create dataset for splines, this is forest only
forest_df_splines <- forest_df %>% 
      filter(forest_mask_buff == 1)

for(i in 1:sum(stringr::str_detect(colnames(forest_df), "PC"))){
      forest_df_splines <- prepareMGCVsplines(forest_df_splines, paste0("forest_PC", i), n_knots = n_knots, smooth = smoother)
}

# Add splines to full data set. Should be okay to leave the non-forest values as NA
forest_df_splines_sub <- forest_df_splines %>%
      dplyr::select(contains(c("x", "y", "spline")))

forest_df_merged <- merge(forest_df, forest_df_splines_sub, by = c("x", "y"), all.x = T)

forest_df_merged_sub <- forest_df_merged %>%
      dplyr::select(contains(c("x", "y", "PC2_spline", "PC3_spline", "PC4_spline", "forest_mask_buff")))

forest_covariates <- rast(forest_df_merged, type = 'xyz', crs = crs(km_proj))

# Export ----------------------------------------------------
writeRaster(clim_topo_covariates, paste0("covariates/processed/", sum(stringr::str_detect(names(clim_topo_PC_stack), "PC")), "clim_topo_300m_", smoother, "_", n_knots, "k.tif"), overwrite=TRUE)
writeRaster(forest_covariates, paste0("covariates/processed/", sum(stringr::str_detect(names(forest_PC_stack), "PC")), "forest_300m_", smoother, "_", n_knots, "k.tif"), overwrite=TRUE)

