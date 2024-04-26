library(inlabru)
library(terra)
library(tidyterra)
library(sf)
library(sp)
source('source/misc_functions.R')

# DATA FILES ------------------------------------------
load('model_out/Formica_lugubris/lgcp/2575498.22_all30m_3k_tp_Formica_lugubris_E7_r150_0.1_s1_0.01.RData')

clim_topo_covariates <- rast("covariates/processed/6clim_topo_tp_3k.tif")
clim_topo_covariates$lat_raster <- scale(clim_topo_covariates$lat_raster)

forest_covariates <- rast("covariates/processed/3forest_tp_3k.tif")

effort_rast_10km <- rast('covariates/processed/effort_rast_lgcp_10km.tif') %>%
      scale()

ROI <- vect('spatial_other/ROI_kmproj.shp')
FE_managed <- vect('spatial_other/Forestry_England_managed_forest.shp') %>% 
      terra::project(crs(km_proj))
cropton <- FE_managed %>% 
      filter(extent == "Cropton") %>% 
      fillHoles

# SET PARAMETERS ---------------------------------------
predictions_resolution <- 2
species_choice = "Formica lugubris"
smoother = 'tp'
n_knots = 3

pred_area <- cropton

# PREDICT -----------------------------------------
# Create grid prediction pixels
x_range <- seq(from = ext(pred_area)$xmin, to = ext(pred_area)$xmax, by = predictions_resolution)
y_range <- seq(from = ext(pred_area)$ymin, to = ext(pred_area)$ymax, by = predictions_resolution)

# Create a grid of points
grid_points <- expand.grid(x = x_range, y = y_range)





coordinates(grid_points) <- ~x+y
grid_pixels <- SpatialPixels(grid_points)

maskborder <- st_as_sf(pred_area) %>% st_transform(., crs = mesh$crs)

fixed_effects <- paste(model[["names.fixed"]], collapse = " + ")
fixed_effects_effort <- gsub("effort", "max(effort)", fixed_effects)

suitability <- predict(object = model, 
                       newdata = grid_pixels,
                       # mask = maskborder,
                       formula = as.formula(paste0("~ Intercept + mySPDE + ", fixed_effects_effort)),
                       n.samples = 100)

suitability_raster <- as.data.frame(suitability) %>% 
      rename(x = coords.x1, y = coords.x2) %>% 
      dplyr::select(x, y, q0.025, q0.5, q0.975, sd) %>% 
      rast(type = "xyz") %>% 
      mask(pred_area) %>% 
      crop(pred_area)






