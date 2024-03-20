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
lgcp_rufa <- rast('model_out/Formica_rufa/lgcp/Formica_rufa_fixedOnly300m_300mRes.tif')
lgcp_lugubris <- rast('model_out/Formica_lugubris/lgcp/Formica_lugubris_fixedOnly300m_300mRes.tif')

maxent_rufa <- rast('model_out/Formica_rufa/maxent/sporadic_Formica_rufa_3k_tp_all30m_thin100m.tif') %>% 
      mask(mask)
maxent_lugubris <- rast('model_out/Formica_lugubris/maxent/sporadic_Formica_lugubris_3k_tp_all30m_thin0m.tif')

# PROCESS MODEL RESULTS ---------------------------------------
for(i in c("Formica rufa", "Formica lugubris")){
      species_choice = i
      
      if(species_choice == "Formica rufa"){
            lgcp_suit_result = lgcp_rufa
            maxent_result = maxent_rufa
      } else {
            lgcp_suit_result = lgcp_lugubris
            maxent_result = maxent_lugubris
      }
      
      crs(lgcp_suit_result) <- crs(km_proj)

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
      suitable_present <- terra::extract(lgcp_suit_result, combined_presences)[, "q0.5"]
      # Get suitability threshold. Since some odd records might have been found in places that
      # are not actually suitable at all, we trim the bottom end to remove these outliers
      suitability_threshold <- quantile(suitable_present, probs = 0.05, na.rm = T)
      hist(suitable_present, main = species_choice)
      abline(v = suitability_threshold, col = "blue")
      
      # SUITABILITY MAP ---------------------------------------
      suit_lgcp_Eng <- lgcp_suit_result %>% 
            # Remove unsuitable areas
            clamp(lower = suitability_threshold, values = F) %>% 
            # Change scale to go from 0 to 1
            normalise_raster() %>% 
            mask(ROI) %>% 
            tidyterra::select(q0.5) 
      
      png(paste0("figures/", gsub(" ", "_", species_choice),"/", gsub(" ", "_", species_choice), "_suitability_lgcp.png"),
           width = 7, height = 7, units = "in", res = 300)
      plot(suit_lgcp_Eng, main = paste0(species_choice, " - suitability"))
      dev.off()

      # EXPORT ---------------------------------------
      writeRaster(suit_lgcp_Eng, filename = paste0('model_out/', gsub(" ", "_", species_choice), '/lgcp/', gsub(" ", "_", species_choice),'_suitability_adj_300m.tif'), overwrite=TRUE)

      # SELECTED SUITABILITY PLOTS ---------------------------------------
      ## LGCP ---------------------------------------
      # Clamp more to highlight the suitable areas
      suit_lgcp_Eng_plots <- suit_lgcp_Eng %>% 
            clamp(lower = 0.25, values = F) 
      
      new_forest_suit_lgcp <- suit_lgcp_Eng_plots %>% 
            crop(new_forest) %>% 
            mask(new_forest)
      
      ennerdale_suit_lgcp <- suit_lgcp_Eng_plots %>% 
            crop(ennerdale) %>% 
            mask(ennerdale)
      
      cropton_suit_lgcp <- suit_lgcp_Eng_plots %>% 
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
      
      ## Maxent ---------------------------------------
      maxent_result_plots <- maxent_result %>% 
            clamp(lower = 0.25, values = F) 
      
      new_forest_suit_maxent <- maxent_result_plots %>%
            crop(new_forest) %>%
            mask(new_forest)

      ennerdale_suit_maxent <- maxent_result_plots %>%
            crop(ennerdale) %>%
            mask(ennerdale)

      cropton_suit_maxent <- maxent_result_plots %>%
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
}














