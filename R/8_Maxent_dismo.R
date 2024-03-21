# HEADER --------------------------------------------
#
# Author: Joris Wiethase
# Email: j.wiethase@gmail.com
# 
# Script Description:  
# Runs a Maxent suitability model using the 'dismo' package   

rm(list = ls())
library(dismo)
library(terra)
library(tidyterra)
library(raster)
source('source/misc_functions.R')

# SET PARAMETERS ------------------------------------------
smoother = 'tp'
n_knots = 3
thin_dist = 0

covars_selection <- c("clim_topo_PC1_spline1", "clim_topo_PC1_spline2",
                      "clim_topo_PC2_spline1", "clim_topo_PC2_spline2",
                      "clim_topo_PC3_spline1", "clim_topo_PC3_spline2",
                      "clim_topo_PC4_spline1", "clim_topo_PC4_spline2",
                      "clim_topo_PC5_spline1", "clim_topo_PC5_spline2",
                      "clim_topo_PC6_spline1", "clim_topo_PC6_spline2",
                      "forest_PC1_spline1", "forest_PC1_spline2",
                      "forest_PC2_spline1", "forest_PC2_spline2",
                      "distance_ancient")

# DATA FILES ------------------------------------------
## Covariates ------------------------------------------
ROI <- vect('spatial_other/ROI_kmproj.shp')
forest_stack <- rast("covariates/processed/forest_stack_300m.tif")

forest_covariates <- rast(paste0("covariates/processed/2forest_300m_", smoother, "_", n_knots, "k.tif"))

forest_mask_buff <- rast("covariates/processed/forest_mask_buff_30m.tif") %>%
      terra::resample(forest_covariates, method = "sum") 
forest_mask_buff <- ifel(forest_mask_buff > 10, 1, 0)

clim_topo_covariates <- rast(paste0("covariates/processed/6clim_topo_300m_", smoother, "_", n_knots, "k.tif")) %>%
      resample(forest_covariates)
clim_topo_covariates$lat_raster <- scale(clim_topo_covariates$lat_raster)

distance_ancient <- forest_stack %>% 
      tidyterra::select(distance_ancient) %>% 
      resample(forest_covariates)

covariates_stack <- c(forest_covariates, clim_topo_covariates, distance_ancient) %>% 
      terra::mask(forest_mask_buff)

predictors <- raster::subset(stack(covariates_stack), subset = covars_selection)  

for(i in c("Formica rufa", "Formica lugubris")){
      species_choice = i

      ## Observations ------------------------------------------
      sporadic <- read.csv('species_data/processed_csv/sporadic_combined.csv') %>% 
            filter(species == species_choice) %>% 
            dplyr::select(x, y)
      
      # exhaustive <- vect('species_data/processed_shp/exhaustive_combined.shp') %>% 
      #       # thin_spatial(., dist_meters = thin_dist, seed = 42) %>% 
      #       as.data.frame(geom = "XY") %>% 
      #       filter(species == species_choice) %>% 
      #       dplyr::select(x, y)
      
      thinned_presences <- sporadic %>% 
            vect(geom = c('x', 'y'), crs = crs(km_proj)) %>% 
            thin_spatial(., dist_meters = thin_dist, seed = 42) %>% 
            as.data.frame(geom = "XY")
      
      # thinned_presences <- rbind(sporadic, exhaustive)
      
      # Bias field ------------------------------------------
      #Select background points based on sampling bias raster
      # effort_rast_10km <- rast('covariates/processed/effort_rast_lgcp_10km.tif') %>%
      #       disagg(fact = 33) %>%
      #       resample(forest_mask_buff) %>%
      #       mask(forest_mask_buff) %>%
      #       normalise_raster() %>%
      #       raster()
      # 
      # bg_bias <- randomPoints(mask = effort_rast_10km, p = thinned_presences, n = 1000, prob = TRUE)
      # write.csv(bg_bias, 'covariates/processed/maxent_bias_points_1000.csv')
      bg_bias <- read.csv('covariates/processed/maxent_bias_points_1000.csv') %>% 
            dplyr::select(x, y)
      
      # PARTITION ------------------------------------------
      fold <- kfold(thinned_presences, k = 5)
      occtest <- thinned_presences[fold == 1, ]
      occtrain <- thinned_presences[fold != 1, ]
      
      # TRAINING ------------------------------------------
      ## FIT MAXENT MODEL & PREDICT ------------------------------------------
      me <- maxent(x = predictors, p = thinned_presences, 
                   a = bg_bias)
      
      suitability_preds <- predict(me, predictors, progress = 'text')
      
      suitability_raster <- rast(suitability_preds) 
      
      writeRaster(suitability_raster, 
                  filename = paste0("model_out/", gsub(" ", "_", species_choice), "/maxent/sporadic_", 
                                    gsub(" ", "_", species_choice), "_", n_knots, "k_", smoother, "_all300m_thin", thin_dist, "m.tif"), overwrite = T)
      
      pdf(paste0("figures/", gsub(" ", "_", species_choice), "/maxent/sporadic_",
                 gsub(" ", "_", species_choice), "_", n_knots, "k_", smoother, "_all300m_thin", thin_dist, "m.pdf"), width = 14, height = 9)
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
                 gsub(" ", "_", species_choice), "_", n_knots, "k_", smoother, "_all300m_thin", thin_dist/1000, "km.pdf"), width = 14, height = 9)
      par(mfrow=c(2, 2))
      plot(e, 'ROC')
      plot(e, 'TPR')
      boxplot(e)
      density(e)
      dev.off()
      
      unseen_records <- vect('species_data/processed_shp/exhaustive_combined.shp') %>%
            filter(species == species_choice) %>%
            as.data.frame(geom = "XY") %>%
            dplyr::select(x, y) 
      
      unseen_records_points <- SpatialPoints(cbind(unseen_records$x, unseen_records$y),
                                             proj4string=CRS(proj4string(predictors)))
      
      unseen_presence_env <- raster::extract(predictors, unseen_records_points)
      
      unseen_pred <- predict(me, unseen_presence_env)
      
      hist(unseen_pred)
      
}






