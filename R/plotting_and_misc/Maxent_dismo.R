rm(list = ls())
library(dismo)
library(terra)
library(tidyterra)
library(raster)
source('source/misc_functions.R')

# SET PARAMETERS ------------------------------------------
smoother = 'cr'
n_knots = 3

species_choice <- "Formica lugubris"

covars_selection <- c("clim_topo_PC1_spline1", "clim_topo_PC1_spline2",
                      "clim_topo_PC2_spline1", "clim_topo_PC2_spline2",
                      "clim_topo_PC3_spline1", "clim_topo_PC3_spline2",
                      "clim_topo_PC4_spline1", "clim_topo_PC4_spline2",
                      "clim_topo_PC5_spline1", "clim_topo_PC5_spline2",
                      "clim_topo_PC6_spline1", "clim_topo_PC6_spline2",
                      "forest_PC2_spline1", "forest_PC2_spline2",
                      "forest_PC3_spline1", "forest_PC3_spline2",
                      "forest_PC4_spline1", "forest_PC4_spline2",
                      "lat_raster",
                      "forest_mask_buff")

# DATA FILES ------------------------------------------
ROI <- vect('spatial_other/ROI_outline_27700.shp') %>% 
      terra::project(crs(km_proj))
sporadic <- read.csv('species_data/processed_csv/sporadic_combined.csv') %>% 
      filter(species == species_choice) %>% 
      dplyr::select(x, y)
exhaustive <- read.csv('species_data/processed_csv/exhaustive_combined.csv') %>% 
      filter(species == species_choice) %>% 
      dplyr::select(x, y)

clim_topo_covariates <- rast(paste0("covariates/processed/6clim_topo_300m_", smoother, "_", n_knots, "k.tif")) 
forest_covariates <- rast(paste0("covariates/processed/4forest_300m_", smoother, "_", n_knots, "k.tif")) %>% 
      resample(clim_topo_covariates)
covariates <- stack(stack(forest_covariates), stack(clim_topo_covariates))
predictors <- raster::subset(covariates, subset = covars_selection)      


presences <- rbind(sporadic, exhaustive)
res(predictors)  # Each cell is 333m wide
thinned_presences <- gridSample(presences, predictors, n=1)

# PARTITION ------------------------------------------
fold <- kfold(thinned_presences, k = 5)
occtest <- thinned_presences[fold == 1, ]
occtrain <- thinned_presences[fold != 1, ]

# TRAINING ------------------------------------------
## FIT MAXENT MODEL & PREDICT ------------------------------------------
me <- maxent(x = predictors, p = thinned_presences, nbg = 1e04)

suitability_raster <- predict(me, predictors, progress='text') %>% 
      rast()
plot(suitability_raster)

# TESTING ------------------------------------------
bg <- randomPoints(predictors, 10000)

e1 <- evaluate(me, p = occtest, a = bg, x = predictors)

pvtest <- data.frame(extract(predictors, occtest))
avtest <- data.frame(extract(predictors, bg))

e2 <- evaluate(me, p = pvtest, a = avtest)

testp <- predict(me, pvtest) 
testa <- predict(me, avtest) 

e3 <- evaluate(p = testp, a = testa)
e3
threshold(e3)

plot(e3, 'ROC')

## EXPORT ------------------------------------------
writeRaster(suitability_raster, filename = paste0("model_out/", gsub(" ", "_", species_choice), "/maxent/maxent_", gsub(" ", "_", species_choice), "_all300m_thinned.tif"), overwrite = T)

pdf(paste0("figures/", gsub(" ", "_", species_choice), "/maxent/maxent_", gsub(" ", "_", species_choice), "_all300m_thinned.pdf"), width = 14, height = 9)
par(mfrow=c(1, 2))
plot(me)
plot(suitability_raster, main = paste0("Maxent suitability: ", species_choice))
dev.off()
par(mfrow=c(1, 1))





