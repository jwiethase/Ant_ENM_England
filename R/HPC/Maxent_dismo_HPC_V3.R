rm(list = ls())
options(java.parameters = "-Xmx350g") 
library(rJava)
library(dismo)
library(terra)
library(tidyterra)
library(raster)
library(sf)
library(tidyverse)
library(ggpmisc)
library(ecospat)
library(usdm)
library(data.table)

setwd('/users/jhw538/scratch/ant_modelling')
source('source/misc_functions.R')

args <- commandArgs(trailingOnly=TRUE)
job <- as.integer(args[1])

# SET PARAMETERS ------------------------------------------
species_choices <- c("Formica rufa", "Formica lugubris")
thin_list = c(0, 100, 250) 

mult_combs <- crossing(species_choices, thin_list) # 6
comb_values <- mult_combs[job, ]
print(comb_values)

species_choice = comb_values$species_choices
thin_dist = comb_values$thin_list 

dir.create(paste0("model_out/", gsub(" ", "_", species_choice)), showWarnings = F)
dir.create(paste0("model_out/", gsub(" ", "_", species_choice), "/maxent_V3"), showWarnings = F)
dir.create(paste0("figures/", gsub(" ", "_", species_choice)), showWarnings = F)
dir.create(paste0("figures/", gsub(" ", "_", species_choice), "/maxent_V3"), showWarnings = F)

# DATA FILES ------------------------------------------
bg_bias <- read.csv('data/maxent_bias_points_1000.csv') %>% 
      dplyr::select(x, y)

forest_covariates <- rast("data/forest_PCA_30m.tif")

forest_mask_buff <- rast("data/forest_mask_buff_30m.tif") %>%
      terra::resample(forest_covariates) %>% 
      subst(0, NA)

climate_covariates <- rast("data/climate_stack_1km.tif") %>% 
      terra::resample(forest_covariates) 

topo_covariates <- rast("data/topo_stack_30m.tif") %>%
      resample(forest_covariates)

distance_forest <- rast("data/distance_forest_30m.tif") %>%
      terra::resample(forest_covariates) %>% 
      scale()
names(distance_forest) <- "distance_forest"

covariates_stack <- c(forest_covariates, climate_covariates, topo_covariates, distance_forest) %>% 
      terra::mask(forest_mask_buff)

# df <- climate_covariates %>%
#       as.data.table(xy=T) %>%
#       select(-lat_raster, -x, -y) %>%
#       drop_na()
# 
# while(TRUE) {
#       # Calculate VIF
#       vif_result <- vif(df)
# 
#       # Find the maximum VIF value and the corresponding variable
#       max_vif <- max(vif_result$VIF)
#       max_vif_var <- vif_result$Variables[which.max(vif_result$VIF)]
# 
#       # Check if the maximum VIF is below the threshold (10)
#       if (max_vif < 10) {
#             break
#       }
# 
#       # Remove the variable with the highest VIF
#       df <- df[ , !(colnames(df) %in% max_vif_var)]
# }
# 
# # The resulting dataframe (df) contains only the variables with VIF < 10
# keep <- c(names(df), names(topo_covariates), "distance_forest", "forest_PC1", "forest_PC2")

keep <- c('SBIO3_Isothermality',
          'SBIO4_Temperature_Seasonality',
          'SBIO6_Min_Temperature_of_Coldest_Month',
          'median_total_rain_coldest',
          'median_total_rain_hottest',
          'median_temp_coldest',
          'median_temp_hottest',
          'dry_duration_09perc',
          'northness',
          'eastness',
          'hillshade',
          'slope',
          'distance_forest',
          'forest_PC1',
          'forest_PC2')

predictors <- raster::subset(raster::stack(covariates_stack), 
                             subset = keep)  

print("predictors done")

sporadic <- read.csv('data/sporadic_combined.csv') %>% 
      filter(species == species_choice) %>% 
      dplyr::select(x, y)

thinned_presences <- sporadic %>% 
      vect(geom = c('x', 'y'), crs = crs(km_proj)) %>% 
      thin_spatial(., dist_meters = thin_dist, seed = 42) %>% 
      as.data.frame(geom = "XY")

unseen_records_vect <- vect('data/exhaustive_combined.shp') %>%
      filter(species == species_choice) 

unseen_records <- unseen_records_vect %>%
      as.data.frame(geom = "XY") %>%
      dplyr::select(x, y) 

bg_bias <- read.csv('data/maxent_bias_points_1000.csv') %>% 
      dplyr::select(x, y)

gc()
# FIT MAXENT MODEL ------------------------------------------
# Fit 5 times 
me_list <- maxent(x = predictors, p = thinned_presences, 
                  a = bg_bias, args = c("hinge=false", "threads=8", "replicates=5", "replicatetype=crossvalidate"))

# TESTING ----------------------------------------------
# Check how each replicate model performs against unseen data, 
# then choose the best performing model
bg <- randomPoints(predictors, 10000)

# Get AUC for each model, choose the one with the highest AUC
test_evaluation_results <- lapply(me_list@models, function(model) {
      evaluate(p = unseen_records, a = bg, model = model, x = predictors)@auc
})

# Access the best model
best_model <- me_list@models[[ which.max(test_evaluation_results)]]

# Get response curves
pdf(paste0("figures/", gsub(" ", "_", species_choice), "/maxent_V3/V3_effectPlots_", 
           gsub(" ", "_", species_choice), "_all30m_thin", thin_dist, "m.pdf"), width = 14, height = 14)
par(mfrow=c(1, 1))
response(best_model, expand = 0)
dev.off()

# Predict using the best model
suitability_preds <- predict(best_model, predictors, progress = 'text')

suitability_raster <- rast(suitability_preds) 
plot(suitability_raster)

writeRaster(suitability_raster, 
            filename = paste0("model_out/", gsub(" ", "_", species_choice), "/maxent_V3/V3_sporadic_", 
                              gsub(" ", "_", species_choice), "_all30m_thin", thin_dist, "m.tif"), overwrite = T)

pdf(paste0("figures/", gsub(" ", "_", species_choice), "/maxent_V3/V3_sporadic_",
           gsub(" ", "_", species_choice), "_all30m_thin", thin_dist, "m.pdf"), width = 14, height = 9)
par(mfrow=c(1, 2))
plot(best_model)
plot(suitability_raster, main = paste0("Maxent suitability: ", species_choice))
dev.off()
par(mfrow=c(1, 1))

unseen_records_points <- SpatialPoints(cbind(unseen_records$x, unseen_records$y),
                                       proj4string=CRS(proj4string(predictors)))

unseen_presence_env <- raster::extract(predictors, unseen_records_points, ID = F)

unseen_pred <- predict(best_model, unseen_presence_env) 

obs <- terra::extract(suitability_raster, unseen_records_vect, ID=F) %>% 
      drop_na() %>% 
      unique()

avtest <- data.frame(raster::extract(predictors, bg))
random_absence_p <- predict(best_model, avtest) 

e_unseen <- evaluate(p = unseen_pred, a = random_absence_p)

pdf(paste0("figures/", gsub(" ", "_", species_choice), "/maxent_V3/V3_sporadic_Unseen_", 
           gsub(" ", "_", species_choice), "_all30m_thin", thin_dist, "m.pdf"), width = 9, height = 7)
par(mfrow=c(2, 2))
hist(unseen_pred, main = "Raw extracted suitability", xlab = "Extracted suitability at unseen presences")
boxplot(unseen_pred, main = paste0("Raw extracted median: ", round(median(unseen_pred, na.rm = T), digits = 2)))
plot(e_unseen, 'ROC')
boyce_test <- ecospat::ecospat.boyce(fit = suitability_raster, obs = obs$layer)
title(paste0("Boyce test cor: ", boyce_test$cor))
dev.off()







