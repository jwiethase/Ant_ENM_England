# Identify Forestry England-managed patches that might be colonised in the future 
rm(list = ls())
library(landscapemetrics)
library(terra)
library(tidyterra)
library(tidyverse)
source('source/misc_functions.R')

# Parameters ----------------------------------------------------
species_choice = 'Formica rufa'
model_choice = 'lgcp'
point_buffer = 0.05                  # Buffer around nest records in meter, to account for GPS inaccuracy
area_percentile_threshold = 0.01   # Percentile value to determine which forest patch areas are rarely occupied by ants 
ON_threshold = 0.1
SNO_threshold = 0.5       
max_gap_dispersal = 0.12            # Distance in meters between likely occupied and suitable that is considered likely to be dispersed across naturally
max_gap_translocation = 0      # Distance in meters beyond which we consider suitable patches very likely not occupied, and candidates for translocation

# Load and prepare data ----------------------------------------------------
FE_managed <- vect('spatial_other/Forestry_England_managed_forest.shp') %>% 
      project(crs(km_proj))
forest_stack <- rast('covariates/processed/forest_stack_300m.tif') 

if(model_choice == 'lgcp' & species_choice == 'Formica rufa'){
      suitability_map <- rast('Colin_plots/allPreds_lgcp_Eng_300m.tif')
      crs(suitability_map) <- crs(km_proj)
      suitability_map <- suitability_map 
}

if(model_choice == 'maxent' & species_choice == 'Formica rufa'){
      suitability_map <- rast('model_out/Formica_rufa/maxent/Formica_rufa_thin500_3k_cr_all30m_thinned.tif')
}

if(model_choice == 'maxent' & species_choice == 'Formica lugubris'){
      suitability_map <- rast('model_out/Formica_lugubris/maxent/Formica_lugubris_thin1000_3k_cr_all30m_thinned.tif')
}

forest_buffer_mask <- ifel(forest_stack$forest_mask_buff, 0, NA)
forest_buffer_mask <- resample(forest_buffer_mask, suitability_map)

suitability_map <- suitability_map %>% 
      mask(forest_buffer_mask)

sporadic <- read.csv('species_data/processed_csv/sporadic_combined.csv') %>% 
      filter(source != 'dallimore', source != 'nym', source != 'gaitbarrows', source != 'hardcastle') %>% 
      dplyr::select(x, y, species)
exhaustive <- read.csv('species_data/processed_csv/exhaustive_combined.csv') %>% 
      dplyr::select(x, y, species)

combined_presences <- rbind(sporadic, exhaustive) %>% 
      vect(geom = c('x', 'y'), crs = crs(km_proj)) 
      
ant_vect <- combined_presences %>% 
      filter(species == species_choice)

ant_vect_buff <- terra::buffer(ant_vect, width = point_buffer)

# 1. Create raster of relevant patches ----------------------------------------------------
# Get suitable forest patch areas only, for ON patches. Probably best to set lower suitability threshold,
# since nests might still be in areas that are no longer very suitable
suitable_forest_forON <- suitability_map %>% 
      clamp(lower = ON_threshold, values = F)
suitable_forest_forON_mask <- ifel(is.na(suitable_forest_forON), NA, 1)

# Get suitable areas only, for SNO patches. Higher suitability threshold,
# since these patches are likely not occupied yet
suitable_forest_forSNO <- suitability_map %>% 
      clamp(lower = SNO_threshold, values = F)

# Get suitable forest patch IDs
suitable_forest_forON_ID <- get_patches(suitable_forest_forON_mask, 
                                        directions = 8)[[1]][[1]]
names(suitable_forest_forON_ID) <- 'patch_ID'

# 2. Overlay nest locations with suitable forest patches to identify forest patches likely now occupied (ON) ----------------------------------------------------
ON_patch_IDs <- terra::extract(suitable_forest_forON_ID, ant_vect_buff) %>% 
      drop_na()

# Create binary raster showing forest patches where there is at least one nest
ON_patches_binary <- suitable_forest_forON_ID$patch_ID %in% ON_patch_IDs$patch_ID

ON_patches_mask <- ifel(isFALSE(ON_patches_binary), NA, 1)

# Create graduated raster showing suitability scores in these ON patches
ON_patches_gradient <- suitable_forest_forON %>% 
      mask(ON_patches_mask)

# 3. Create raster of forest patches suitable but likely not occupied now (SNO) -------------------------------------
# Patch areas that are suitable but likely not occupied. 
SNO_patches <- suitable_forest_forSNO * (1-ON_patches_binary) # Turns the TRUE/FALSE into 0/1
SNO_patches <- ifel(SNO_patches == 0, NA, SNO_patches)

# 4. SNO patches that might be colonised naturally, or might already be colonised (narrow forest gaps) -------------------------------------
# Get IDs of SNO patches
SNO_IDs <- get_patches(ifel(!is.na(SNO_patches), 1, NA), 
                                 directions = 8)[[1]][[1]]
names(SNO_IDs) <- 'patch_ID'

# Identify SNO patches that are within the distance threshold to ON patches
ON_patches_mask <- ifel(isFALSE(ON_patches_binary), NA, 1)
ON_patches_buffered_close <- terra::buffer(ON_patches_mask, max_gap_dispersal, background = 0) %>% 
      subst(0, NA) %>% 
      as.polygons()

SNO_patches_close <- terra::extract(SNO_IDs, ON_patches_buffered_close, 
                                    ID = FALSE, 
                                    touches = T) %>% 
      drop_na() %>% 
      unique()

dispersal_patches <- SNO_IDs %>% 
      filter(patch_ID %in% SNO_patches_close$patch_ID)
dispersal_patches_mask <- ifel(is.na(dispersal_patches), NA, 1)

# Create graduated raster showing suitability scores in these dispersal patches
dispersal_patches_gradient <- suitable_forest_forSNO %>% 
      mask(dispersal_patches_mask)

# 5. SNO patches that will likely not be colonised naturally (far from ON) -------------------------------------
# These patches are candidates for translocation
dispersal_patches_binary <- ifel(is.na(dispersal_patches_mask), 0, 1)
ON_or_dispersal <- ON_patches_binary + dispersal_patches_binary 
ON_or_dispersal_mask <- ifel(ON_or_dispersal > 0, 1, NA)

if(max_gap_translocation != 0){
      ON_patches_buffered_far <- terra::buffer(ON_or_dispersal_mask, max_gap_translocation, background = 0) %>% 
            subst(0, NA) %>% 
            as.polygons()
      
      SNO_patches_far <- terra::extract(SNO_IDs, ON_patches_buffered_far, 
                                        ID = FALSE, 
                                        touches = T) %>% 
            drop_na() %>% 
            unique()
      
      translocation_patches <- SNO_IDs %>% 
            filter(!patch_ID %in% SNO_patches_far$patch_ID)
} else {
      translocation_patches <- SNO_patches %>% 
            mask(ON_or_dispersal_mask, inverse = TRUE)
}

translocation_patches_mask <- ifel(is.na(translocation_patches), NA, 1)

# Create graduated raster showing suitability scores in the translocation patches
translocation_patches_gradient <- suitable_forest_forSNO %>% 
      mask(translocation_patches_mask) %>% 
      mask(FE_managed)

writeRaster(translocation_gradient_FE, "translocation_gradient_FE.tif", overwrite=T)
# writeVector(ON_patches_buffered_far, "ON_patches_buffered_far.shp", overwrite=T)

