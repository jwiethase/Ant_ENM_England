rm(list = ls())
library(dismo)
library(terra)
library(tidyterra)
library(raster)
source('source/misc_functions.R')

smoother = 'cr'
n_knots = 3

species_choice <- "Formica rufa"

# DATA FILES ------------------------------------------
ROI <- vect('spatial_other/ROI_outline_27700.shp') %>% 
      terra::project(crs(km_proj))
sporadic_filtered <- read.csv('species_data/processed_csv/sporadic_combined.csv') %>% 
      filter(species == species_choice) %>% 
      dplyr::select(x, y)
exhaustive <- read.csv('species_data/processed_csv/exhaustive_combined.csv') %>% 
      filter(species == species_choice) %>% 
      dplyr::select(x, y)

combined_presences <- sporadic_filtered
obs_spatial <- combined_presences %>%
      vect(geom = c("x", "y"), crs = crs(km_proj), keepgeom = TRUE)

forest_covariates <- rast(paste0("covariates/processed/2forest_", smoother, "_", n_knots, "k.tif")) %>% 
      terra::subset(stringr::str_detect(names(.), "spline"))

topo_covariates <- rast(paste0("covariates/processed/2topo_", smoother, "_", n_knots, "k.tif"))  %>% 
      terra::subset(stringr::str_detect(names(.), "spline")) %>% 
      terra::project(forest_covariates, method = 'max', threads = TRUE)

clim_covariates <- rast(paste0("covariates/processed/4climate_", smoother, "_", n_knots, "k.tif"))  %>% 
      terra::subset(stringr::str_detect(names(.), "spline")) %>% 
      terra::project(forest_covariates, method = 'max', threads = TRUE)

all_preds <- stack(stack(forest_covariates), stack(topo_covariates), stack(clim_covariates))

fold <- kfold(combined_presences, k=5)
occtest <- sporadic_filtered[fold == 1, ]
occtrain <- sporadic_filtered[fold != 1, ]

# fit model, biome is a categorical variable
me <- maxent(all_preds, occtrain)

# predict to entire dataset
suitability_raster <- predict(me, all_preds, progress='text') %>% rast()
plot(suitability_raster)
writeRaster(suitability_raster, filename = "model_output/rasters/maxent_rufa_all250m_sporadic.tif", overwrite = T)

pdf("figures/maxent_rufa_all250m_sporadic.pdf", width = 14, height = 9)
par(mfrow=c(1, 2))
plot(me)
plot(suitability_raster, main = paste0("Maxent suitability: ", species_choice))
dev.off()

