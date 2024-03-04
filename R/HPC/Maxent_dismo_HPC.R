rm(list = ls())
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
thin_list = c(100, 250, 500, 1000)

mult_combs <- crossing(species_choices, thin_list) # 6
comb_values <- mult_combs[job, ]
print(comb_values)

species_choice = comb_values$species_choices
thin_dist = comb_values$thin_list 

n_knots = 3
smoother = "cr"

covars_selection <- c("clim_topo_PC1_spline1", "clim_topo_PC1_spline2",
                      "clim_topo_PC2_spline1", "clim_topo_PC2_spline2",
                      "clim_topo_PC3_spline1", "clim_topo_PC3_spline2",
                      "clim_topo_PC4_spline1", "clim_topo_PC4_spline2",
                      "clim_topo_PC5_spline1", "clim_topo_PC5_spline2",
                      "clim_topo_PC6_spline1", "clim_topo_PC6_spline2",
                      "forest_PC2_spline1", "forest_PC2_spline2",
                      "forest_PC3_spline1", "forest_PC3_spline2",
                      "lat_raster",
                      "forest_mask_buff")

dir.create(paste0("model_out/", gsub(" ", "_", species_choice)), showWarnings = F)
dir.create(paste0("model_out/", gsub(" ", "_", species_choice), "/maxent"), showWarnings = F)
dir.create(paste0("figures/", gsub(" ", "_", species_choice)), showWarnings = F)
dir.create(paste0("figures/", gsub(" ", "_", species_choice), "/maxent"), showWarnings = F)


# DATA FILES ------------------------------------------
ROI <- vect('data/ROI_outline_27700.shp') %>% 
      terra::project(crs(km_proj))
sporadic <- read.csv('data/sporadic_combined.csv') %>% 
      filter(species == species_choice) %>% 
      dplyr::select(x, y)
exhaustive <- read.csv('data/exhaustive_combined.csv') %>% 
      filter(species == species_choice) %>% 
      dplyr::select(x, y)

clim_topo_covariates <- rast(paste0("data/6clim_topo_", smoother, "_", n_knots, "k.tif"))
forest_covariates <- rast(paste0("data/2forest_", smoother, "_", n_knots, "k.tif")) %>% 
      resample(clim_topo_covariates)
covariates <- stack(stack(forest_covariates), stack(clim_topo_covariates))
predictors <- raster::subset(covariates, subset = covars_selection)      

effort_rast_10km <- raster('data/effort_rast_10km.tif') %>% 
      disaggregate(fact = res(.)/res(clim_topo_covariates)) %>% 
      raster::resample(., raster(clim_topo_covariates), method = "bilinear")

thinned_presences <- rbind(sporadic, exhaustive) %>% 
      vect(geom = c('x', 'y'), crs = crs(km_proj)) %>% 
      thin_spatial(., dist_meters = thin_dist) %>% 
      as.data.frame(geom = "XY")

# PARTITION ------------------------------------------
fold <- kfold(thinned_presences, k = 5)
occtest <- thinned_presences[fold == 1, ]
occtrain <- thinned_presences[fold != 1, ]

# TRAINING ------------------------------------------
## FIT MAXENT MODEL & PREDICT ------------------------------------------
me <- maxent(x = predictors, p = thinned_presences, 
             nbg = 5e04)

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
writeRaster(suitability_raster, 
            filename = paste0("model_out/", gsub(" ", "_", species_choice), "/maxent/", 
                              gsub(" ", "_", species_choice), "_", n_knots, "k_", smoother, "_all300m_thinned.tif"), overwrite = T)

pdf(paste0("figures/", gsub(" ", "_", species_choice), "/maxent/", 
           gsub(" ", "_", species_choice), "_", n_knots, "k_", smoother, "_all300m_thinned.pdf"), width = 14, height = 9)
par(mfrow=c(1, 2))
plot(me)
plot(suitability_raster, main = paste0("Maxent suitability: ", species_choice))
dev.off()
par(mfrow=c(1, 1))





