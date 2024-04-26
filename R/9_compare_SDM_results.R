# HEADER --------------------------------------------
#
# Author: Joris Wiethase
# Email: j.wiethase@gmail.com
# 
# Script Description:  
# Takes the Maxent and Point Process model outputs and processes them, then export and visualise.
# In the case of Point Process model outputs, these are adjusted to a more meaningful scale
# using actual nest location records

# LOAD PACKAGES -----------------------------------
rm(list = ls())
library(terra)
library(tidyverse)
library(tidyterra)
library(viridis)
library(patchwork)
source('source/misc_functions.R')

# DATA FILES ---------------------------------------
ROI <- vect('spatial_other/ROI_outline_27700.shp') %>% 
      terra::project(crs(km_proj))

FE_managed <- vect('spatial_other/Forestry_England_managed_forest.shp') %>% 
      terra::project(crs(km_proj))

# Example areas
new_forest <- FE_managed %>% 
      filter(extent == "The Open Forest") %>% 
      fillHoles

ennerdale <- FE_managed %>% 
      filter(extent == "Ennerdale") %>% 
      fillHoles

cropton <- FE_managed %>% 
      filter(extent == "Cropton") %>% 
      fillHoles

# The model outputs
# lgcp_rufa <- lapply(list.files(path = 'model_out/Formica_rufa/lgcp/', pattern = 'slice', full.names = T), rast) %>%
#       sprc() %>%
#       mosaic(fun = 'max')
# writeRaster(lgcp_rufa,
#             filename = 'model_out/Formica_rufa/lgcp/sporadic_Formica_rufa_lgcp_all30m_merged.tif',
#             overwrite=T)
lgcp_rufa <- rast('model_out/Formica_rufa/lgcp/Formica_rufa_lgcp_all30m_mosaic.tif')

# lgcp_lugubris <- lapply(list.files(path = 'model_out/Formica_lugubris/lgcp/', pattern = 'slice', full.names = T), rast) %>%
#       sprc() %>%
#       mosaic(fun = 'max')
# writeRaster(lgcp_lugubris,
#             filename = 'model_out/Formica_lugubris/lgcp/sporadic_Formica_lugubris_lgcp_all30m_merged.tif',
#             overwrite=T)
lgcp_lugubris <- rast('model_out/Formica_lugubris/lgcp/Formica_lugubris_lgcp_all30m_mosaic.tif')

maxent_rufa <- rast('model_out/Formica_rufa/maxent/V3_sporadic_Formica_rufa_all30m_thin0m.tif')
maxent_lugubris <- rast('model_out/Formica_lugubris/maxent/V3_sporadic_Formica_lugubris_all30m_thin250m.tif')

# PROCESS MODEL RESULTS ---------------------------------------
for(i in c("Formica rufa", "Formica lugubris")){
      i="Formica rufa"
      species_choice = i
      
      if(species_choice == "Formica rufa"){
            lgcp_suit_result = lgcp_rufa
            maxent_result = maxent_rufa
      } else {
            lgcp_suit_result = lgcp_lugubris
            maxent_result = maxent_lugubris
      }
      
      crs(lgcp_suit_result) <- crs(km_proj)
      
      unseen_records_vect <- vect('species_data/processed_shp/exhaustive_combined.shp') %>%
            filter(species == species_choice) 
      
      # SUITABILITY THRESHOLDS ---------------------------------------
      # Suitability is on intensity scale. Identify thresholds of suitability estimates 
      # at which the ants are actually found.
      sporadic <- read.csv('species_data/processed_csv/sporadic_combined.csv') %>% 
            filter(species == species_choice,
                   source != "dallimore", source != "nym", source != "gaitbarrows", source != "hardcastle") %>% 
            dplyr::select(x, y)
      exhaustive <- read.csv('species_data/processed_csv/exhaustive_combined.csv') %>% 
            filter(species == species_choice) %>% 
            dplyr::select(x, y)
      
      combined_presences <- rbind(sporadic, exhaustive) %>% 
            vect(geom = c('x', 'y'), crs = crs(km_proj))
      
      # Suitability estimates on intensity scale, in locations where ants are present
      suitable_present <- terra::extract(lgcp_suit_result, combined_presences, ID = F)[, "q0.5"]
      # Get suitability threshold. Since some odd records might have been found in places that
      # are not actually suitable at all, we trim the bottom end to remove these outliers
      suitability_threshold <- quantile(suitable_present, probs = c(0.05), na.rm = T)
      hist(suitable_present, main = species_choice)
      abline(v = suitability_threshold, col = "blue")
      
      # SUITABILITY MAP ---------------------------------------
      suit_lgcp_Eng <- lgcp_suit_result %>% 
            # Remove unsuitable areas
            clamp(lower = suitability_threshold, 
                  # upper = suitability_thresholds[2], 
                  values = T) %>% 
            normalise_raster() %>% 
            mask(ROI) %>% 
            tidyterra::select(q0.5) 
      
      png(paste0("figures/", gsub(" ", "_", species_choice),"/", gsub(" ", "_", species_choice), "_suitability_lgcp.png"),
          width = 7, height = 7, units = "in", res = 300)
      plot(suit_lgcp_Eng, main = paste0(species_choice, " - suitability"))
      dev.off()
      
      # EXPORT ---------------------------------------
      writeRaster(suit_lgcp_Eng, filename = paste0('model_out/', gsub(" ", "_", species_choice), '/lgcp/', gsub(" ", "_", species_choice),'_suitability_adj_30m.tif'), overwrite=TRUE)
      
      # SELECTED SUITABILITY PLOTS ---------------------------------------
      ## LGCP ---------------------------------------
      # Clamp more to highlight the suitable areas
      new_forest_suit_lgcp <- suit_lgcp_Eng %>% 
            crop(new_forest) %>% 
            mask(new_forest)
      
      ennerdale_suit_lgcp <- suit_lgcp_Eng %>% 
            crop(ennerdale) %>% 
            mask(ennerdale)
      
      cropton_suit_lgcp <- suit_lgcp_Eng %>% 
            crop(cropton) %>% 
            mask(cropton)
      
      new_forest_lgcp <- ggplot() +
            geom_spatraster(data = new_forest_suit_lgcp) +
            geom_spatvector(data = new_forest, fill = "transparent", col = "red") +
            theme_minimal() +
            ggtitle("New forest") +
            scale_fill_viridis(na.value = "transparent", name = "Suitability")
      
      cropton_lgcp <- ggplot() +
            geom_spatraster(data = cropton_suit_lgcp) +
            geom_spatvector(data = cropton, fill = "transparent", col = "red") +
            theme_minimal() +
            ggtitle("Cropton") +
            scale_fill_viridis(na.value = "transparent", name = "Suitability")
      
      ennerdale_lgcp <- ggplot() +
            geom_spatraster(data = ennerdale_suit_lgcp) +
            geom_spatvector(data = ennerdale, fill = "transparent", col = "red") +
            theme_minimal() +
            ggtitle("Ennerdale") +
            scale_fill_viridis(na.value = "transparent", name = "Suitability")
      
      combined_lgcp <- ggpubr::ggarrange(new_forest_lgcp, cropton_lgcp, ennerdale_lgcp, 
                                         common.legend = T, ncol = 3, nrow = 1, legend = "bottom")
      combined_lgcp
      
      ggsave(plot = combined_lgcp, 
             filename = paste0("figures/", gsub(" ", "_", species_choice),"/", gsub(" ", "_", species_choice), "_lgcp_suit_examples.png"),
             width = 10, height = 3.5, units = "in", dpi = 300)
      
      # Check predictions for unseen presences
      unseen_pred <- terra::extract(suit_lgcp_Eng, unseen_records_vect, ID = F)[, 'q0.5']
      
      png(paste0("figures/", gsub(" ", "_", species_choice),"/", gsub(" ", "_", species_choice), "_lgcp_unseen_preds.png"),
          width = 10, height = 5, units = "in", res = 300)
      par(mfrow=c(1, 2))
      hist(unseen_pred, main = "Raw extracted suitability", xlab = "Extracted suitability at unseen presences")
      boxplot(unseen_pred, main = paste0("Raw extracted median: ", round(median(unseen_pred, na.rm = T), digits = 2)))
      par(mfrow=c(1, 1))
      dev.off()
      
      # SUITABILITY MAP ---------------------------------------
      new_forest_suit_maxent <- maxent_result %>%
            crop(new_forest) %>%
            mask(new_forest)
      
      ennerdale_suit_maxent <- maxent_result %>%
            crop(ennerdale) %>%
            mask(ennerdale)
      
      cropton_suit_maxent <- maxent_result %>%
            crop(cropton) %>%
            mask(cropton)
      
      new_forest_maxent <- ggplot() +
            geom_spatraster(data = new_forest_suit_maxent) +
            geom_spatvector(data = new_forest, fill = "transparent", col = "red") +
            theme_minimal() +
            ggtitle("New forest") +
            scale_fill_viridis(na.value = "transparent", name = "Suitability")
      
      cropton_maxent <- ggplot() +
            geom_spatraster(data = cropton_suit_maxent) +
            geom_spatvector(data = cropton, fill = "transparent", col = "red") +
            theme_minimal() +
            ggtitle("Cropton") +
            scale_fill_viridis(na.value = "transparent", name = "Suitability")
      
      ennerdale_maxent <- ggplot() +
            geom_spatraster(data = ennerdale_suit_maxent) +
            geom_spatvector(data = ennerdale, fill = "transparent", col = "red") +
            theme_minimal() +
            ggtitle("Ennerdale") +
            scale_fill_viridis(na.value = "transparent", name = "Suitability")
      
      combined_maxent <- ggpubr::ggarrange(new_forest_maxent, cropton_maxent, ennerdale_maxent,
                                           common.legend = T, ncol = 3, nrow = 1, legend = "bottom")
      combined_maxent
      
      ggsave(plot = combined_maxent,
             filename = paste0("figures/", gsub(" ", "_", species_choice),"/", gsub(" ", "_", species_choice), "_maxent_examples.png"),
             width = 10, height = 3.5, units = "in", dpi = 300)
      
      png(paste0("figures/", gsub(" ", "_", species_choice),"/", gsub(" ", "_", species_choice), "_suitability_lgcp_maxent.png"),
          width = 12, height = 7, units = "in", res = 300)
      par(mfrow=c(1, 2))
      plot(suit_lgcp_Eng, main = paste0(species_choice, " - lgcp suitability"))
      plot(maxent_result, main = paste0(species_choice, " - Maxent suitability"))
      dev.off()
      par(mfrow=c(1, 1))
}














