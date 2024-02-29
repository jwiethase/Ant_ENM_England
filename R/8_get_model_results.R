# LOAD PACKAGES -----------------------------------
rm(list = ls())
library(INLA) 
library(inlabru)
library(terra)
library(tidyverse)
library(tidyterra)
library(viridis)
library(patchwork)
library(ggpmisc)
source('source/misc_functions.R')

species_choice <- "Formica rufa"
file_list <- list.files("model_output/lgcp_rufa_V5/", full.names = T); file_list
filechoice <- file_list[1]

predictions_resolution <- 0.25

# LOAD FILES  ---------------------------------------
ROI <- vect('spatial_other/ROI_outline_27700.shp') %>% 
      terra::project(crs(km_proj))

load(filechoice)

obs_spatial <- read.csv('species_data/processed_csv/sporadic_combined.csv') %>% 
      filter(species == species_choice) %>% 
      vect(geom = c("x", "y"), crs = crs(km_proj), keepgeom = TRUE)

# GET MODEL ESTIMATES -----------------------------------------
boundary <- st_as_sf(ROI) %>% as("Spatial")  # sp format
max.edge <- 11
mesh <- inla.mesh.2d(boundary = boundary,
                     loc = unique(as.data.frame(obs_spatial)[, c("x", "y")]),
                     max.edge = c(1 ,5) * max.edge,
                     offset = c(1, 2) * max.edge, 
                     cutoff = max.edge/3,
                     crs = fm_CRS(km_proj),
                     min.angle = 26)

required_nx <- round((max(mesh$loc[,1]) - min(mesh$loc[,1])) / predictions_resolution)
required_ny <- round((max(mesh$loc[,2]) - min(mesh$loc[,2])) / predictions_resolution)

maskborder <- st_as_sf(ROI)
maskborder = st_transform(maskborder, mesh$crs)

proj_grid <- fm_pixels(mesh, 
                       dims = c(required_nx, required_ny), 
                       format = 'sp')
fixed_effects <- paste(model[["names.fixed"]], collapse = " + ")
fixed_effects_effort <- gsub("effort_rast_10km", "max(effort_rast_10km, na.rm=T)", fixed_effects)
all_preds <- predict(object = model, 
                     newdata = proj_grid,
                     mask = maskborder,
                     formula = as.formula(paste0("~ (mySPDE + ", fixed_effects_effort, ")")))

all_preds_raster <- as.data.frame(all_preds) %>% 
      rename(x = coords.x1, y = coords.x2) %>% 
      dplyr::select(x, y, q0.025, q0.5, q0.975, sd) %>% 
      rast(type = "xyz") %>% 
      mask(ROI)

random_field_raster <- as.data.frame(random_field) %>% 
      rename(x = coords.x1, y = coords.x2) %>% 
      dplyr::select(x, y, q0.025, q0.5, q0.975, sd) %>% 
      rast(type = "xyz") %>% 
      mask(ROI)

suitability_raster <- as.data.frame(suitability) %>% 
      rename(x = coords.x1, y = coords.x2) %>% 
      dplyr::select(x, y, q0.025, q0.5, q0.975, sd) %>% 
      rast(type = "xyz") %>% 
      mask(ROI)

median_plot <- ggplot() +
      geom_spatraster(data = exp(all_preds_raster$q0.5)) +
      geom_spatvector(data = ROI, fill = NA, lwd = 1) +
      theme_minimal() +
      ggtitle("Median estimated point density", subtitle = "(Response Scale)") +
      scale_fill_viridis(na.value = "transparent", name = "Median", option = 'H')

sd_plot <- ggplot() +
      geom_spatraster(data = (all_preds_raster$sd)) +
      geom_spatvector(data = ROI, fill = NA, lwd = 1) +
      theme_minimal() +
      ggtitle("Standard deviation", subtitle = "(Linear Scale)") +
      scale_fill_viridis(na.value = "transparent", name = "SD")

field_plot <- ggplot() +
      geom_spatraster(data = (random_field_raster$q0.5)) +
      geom_spatvector(data = ROI, fill = NA, lwd = 1) +
      geom_spatvector(data = obs_spatial, cex = 0.4, col = 'red', alpha = 0.2) +
      theme_minimal() +
      ggtitle("Spatial random field", subtitle = "(Linear Scale)") +
      scale_fill_viridis(na.value = "transparent", name = "Median")

suitability_plot <- ggplot() +
      geom_spatraster(data = exp(suitability_raster$q0.5)) +
      geom_spatvector(data = ROI, fill = NA, lwd = 1) +
      theme_minimal() +
      ggtitle("Suitability", subtitle = "(Response Scale)") +
      scale_fill_viridis(na.value = "transparent", name = "Median")
quantile(values(exp(suitability_raster$q0.5)), na.rm = T, probs = 0.05)
combined_plot <- (median_plot | suitability_plot) / (field_plot | sd_plot)
combined_plot

# EXPORT ESTIMATES -----------------------------------------
comb_raster <- c(all_preds_raster$q0.5, all_preds_raster$sd, 
                 suitability_raster$q0.5, random_field_raster$q0.5,
                 exp(all_preds_raster$q0.5), exp(suitability_raster$q0.5)) 
crs(comb_raster) <- crs(km_proj)
names(comb_raster) <- c("all_linear_median", "all_linear_sd",
                        "suitability_linear_median", "random_field",
                        "all_response_median", "suitability_response_median")

writeRaster(comb_raster, paste0("model_output/rasters/", gsub(".RData", ".tif", basename(filechoice))), overwrite=TRUE)




