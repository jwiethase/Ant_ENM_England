rm(list = ls())
library(terra)
library(tidyverse)
library(multcomp)
source('source/misc_functions.R')

NFI <- vect("covariates/raw/NFI_subset_example.shp")
mean_VOM_25m <- rast("covariates/processed/mean_raster_25m.tif")
std_VOM_25m <- rast("covariates/processed/std_raster_25m.tif")
cover_VOM_25m <- rast("covariates/processed/canopy_cover_raster_25m.tif")
TCD_25m <- rast("covariates/processed/TCD_raster_25m.tif") %>% 
      project(cover_VOM_25m)

# Check if LiDAR maps can distinguish land cover classes
NFI_df <- as.data.frame(NFI)
VOM_mean <- terra::extract(mean_VOM_25m, NFI, fun=median, ID=F)
VOM_std <- terra::extract(std_VOM_25m, NFI, fun=median, ID=F)
VOM_cover <- terra::extract(cover_VOM_25m, NFI, fun=median, ID=F)

NFI_df$VOM_mean <- VOM_mean$mean_raster_25m
NFI_df$VOM_std <- VOM_std$std_raster_25m
NFI_df$VOM_cover <- VOM_cover$canopy_cover_raster_25m

NFI_df <- NFI_df %>% 
      filter(!IFT_IOA %in% c("Other vegetation", "Assumed woodland", "Ground prep", "Low density")) %>% 
      pivot_longer(cols = c("VOM_mean", "VOM_std", "VOM_cover"), names_to = "VOM_var", values_to = "VOM_val")

ggplot(NFI_df, aes(x = VOM_val, y = IFT_IOA)) +
      geom_boxplot() +
      theme_bw() +
      facet_wrap(.~VOM_var, scales = "free") +
      xlab("Mean LiDAR vegetation height")
 
# Check agreement between TCD mapo and LiDAR derived canopy cover
cover_df <- c(cover_VOM_25m, TCD_25m) %>% 
      as.data.frame() %>% 
      drop_na()

plot(cover_df$canopy_cover_raster_25m, cover_df$TCD_raster_25m)
rsq(cover_df$canopy_cover_raster_25m, cover_df$TCD_raster_25m)


