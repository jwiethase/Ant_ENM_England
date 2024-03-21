rm(list = ls())
options(java.parameters = "-Xmx250g") 
library(rJava)
library(dismo)
library(terra)
library(tidyterra)
library(raster)
library(sf)
library(tidyverse)
library(ggpmisc)

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

n_knots = 3
smoother = "tp"

covars_selection <- c("clim_topo_PC1_spline1", "clim_topo_PC1_spline2",
                      "clim_topo_PC2_spline1", "clim_topo_PC2_spline2",
                      "clim_topo_PC3_spline1", "clim_topo_PC3_spline2",
                      "clim_topo_PC4_spline1", "clim_topo_PC4_spline2",
                      "clim_topo_PC5_spline1", "clim_topo_PC5_spline2",
                      "clim_topo_PC6_spline1", "clim_topo_PC6_spline2",
                      "forest_PC1_spline1", "forest_PC1_spline2",
                      "forest_PC2_spline1", "forest_PC2_spline2")

dir.create(paste0("model_out/", gsub(" ", "_", species_choice)), showWarnings = F)
dir.create(paste0("model_out/", gsub(" ", "_", species_choice), "/maxent"), showWarnings = F)
dir.create(paste0("figures/", gsub(" ", "_", species_choice)), showWarnings = F)
dir.create(paste0("figures/", gsub(" ", "_", species_choice), "/maxent"), showWarnings = F)

# DATA FILES ------------------------------------------
bg_bias <- read.csv('data/maxent_bias_points_1000.csv') %>% 
      dplyr::select(x, y)
forest_stack <- rast("data/forest_stack_30m.tif")

forest_covariates <- rast(paste0("data/2forest_30m_", smoother, "_", n_knots, "k.tif"))

forest_mask_buff <- rast("data/forest_mask_buff_30m.tif") %>%
      terra::resample(forest_covariates) %>% 
      subst(0, NA)

clim_topo_covariates <- rast(paste0("data/6clim_topo_", smoother, "_", n_knots, "k.tif")) %>% 
      terra::resample(forest_covariates) 

print("clim_topo_covariates done")

# distance_ancient <- forest_stack %>% 
#       tidyterra::select(distance_ancient) %>% 
#       terra::resample(forest_covariates) 
# 
# distance_ancient <- distance_ancient %>% 
#       terra::scale()

# print("distance_ancient done")

covariates_stack <- c(forest_covariates, clim_topo_covariates) %>% 
      terra::mask(forest_mask_buff)

predictors <- raster::subset(stack(covariates_stack), subset = covars_selection)     

print("predictors done")

sporadic <- read.csv('data/sporadic_combined.csv') %>% 
      filter(species == species_choice) %>% 
      dplyr::select(x, y)

thinned_presences <- sporadic %>% 
      vect(geom = c('x', 'y'), crs = crs(km_proj)) %>% 
      thin_spatial(., dist_meters = thin_dist, seed = 42) %>% 
      as.data.frame(geom = "XY")

gc()

# PARTITION ------------------------------------------
fold <- kfold(thinned_presences, k = 5)
occtest <- thinned_presences[fold == 1, ]
occtrain <- thinned_presences[fold != 1, ]

# TRAINING ------------------------------------------
## FIT MAXENT MODEL & PREDICT ------------------------------------------
print("starting maxent")
me <- maxent(x = predictors, p = thinned_presences, 
             a = bg_bias)

print("Maxent done")

gc()

suitability_preds <- predict(me, predictors, progress = 'text')

suitability_raster <- rast(suitability_preds) 

writeRaster(suitability_raster, 
            filename = paste0("model_out/", gsub(" ", "_", species_choice), "/maxent/sporadic_", 
                              gsub(" ", "_", species_choice), "_", n_knots, "k_", smoother, "_all30m_thin", thin_dist, "m.tif"), overwrite = T)

print("suitability_raster done")

pdf(paste0("figures/", gsub(" ", "_", species_choice), "/maxent/sporadic_", 
           gsub(" ", "_", species_choice), "_", n_knots, "k_", smoother, "_all30m_thin", thin_dist, "m.pdf"), width = 14, height = 9)
par(mfrow=c(1, 2))
plot(me)
plot(suitability_raster, main = paste0("Maxent suitability: ", species_choice))
dev.off()
par(mfrow=c(1, 1))

# TESTING ------------------------------------------
bg <- randomPoints(predictors, 10000)

pvtest <- data.frame(raster::extract(predictors, occtest))
avtest <- data.frame(raster::extract(predictors, bg))

testp <- predict(me, pvtest) 
testa <- predict(me, avtest) 

e <- evaluate(p = testp, a = testa)

pdf(paste0("figures/", gsub(" ", "_", species_choice), "/maxent/sporadic_Eval_", 
           gsub(" ", "_", species_choice), "_", n_knots, "k_", smoother, "_all30m_thin", thin_dist/1000, "km.pdf"), width = 14, height = 9)
par(mfrow=c(2, 2))
plot(e, 'ROC')
plot(e, 'TPR')
boxplot(e)
density(e)
dev.off()







