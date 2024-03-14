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
point_buffer = 50                  # Buffer around nest records in meter, to account for GPS inaccuracy
area_percentile_threshold = 0.01   # Percentile value to determine which forest patch areas are rarely occupied by ants 
suitability_threshold = 0.25       # Low threshold value for suitability (consider suitable patches as anything above)
max_gap_dispersal = 1200            # Distance in meters between likely occupied and suitable that is considered likely to be dispersed across naturally
max_gap_translocation = 25000      # Distance in meters beyond which we consider suitable patches very likely not occupied, and candidates for translocation


# Load and prepare data ----------------------------------------------------
ROI_27700 <- vect('spatial_other/ROI_outline_27700.shp') 
FE_managed <- vect('spatial_other/Forestry_England_managed_forest.shp')
forest_stack <- rast('covariates/processed/forest_stack_300m.tif') %>% 
      terra::project(crs('epsg:27700')) 

if(model_choice == 'lgcp'){
      suitability_map <- rast('Colin_plots/allPreds_lgcp_Eng_300m.tif')
      crs(suitability_map) <- crs(km_proj)
      suitability_map <- suitability_map %>% terra::project(crs('epsg:27700')) 
}

if(model_choice == 'maxent' & species_choice == 'Formica rufa'){
      suitability_map <- rast('model_out/Formica_rufa/maxent/Formica_rufa_thin500_3k_cr_all30m_thinned.tif') %>% 
            terra::project(crs('epsg:27700')) 
}

if(model_choice == 'maxent' & species_choice == 'Formica lugubris'){
      suitability_map <- rast('model_out/Formica_lugubris/maxent/Formica_lugubris_thin1000_3k_cr_all30m_thinned.tif') %>% 
            terra::project(crs('epsg:27700')) 
}

sporadic <- read.csv('species_data/processed_csv/sporadic_combined.csv') %>% 
      filter(source != 'dallimore', source != 'nym', source != 'gaitbarrows', source != 'hardcastle') %>% 
      dplyr::select(x, y, species)
exhaustive <- read.csv('species_data/processed_csv/exhaustive_combined.csv') %>% 
      dplyr::select(x, y, species)

combined_presences <- rbind(sporadic, exhaustive) %>% 
      vect(geom = c('x', 'y'), crs = crs(km_proj)) %>% 
      terra::project(crs('epsg:27700')) 
      
ant_vect <- combined_presences %>% 
      filter(species == species_choice)

ant_vect_buff <- terra::buffer(ant_vect, width = point_buffer)

# 1. Create raster of forest patches for England ----------------------------------------------------
forest_mask <- forest_stack %>% 
      tidyterra::select(cover_VOM) 

forest_mask <- ifel(forest_mask < 0.3, NA, 1)

# Get forest area
forest_patch_area <- spatialize_lsm(forest_mask, 
                                level = 'patch',
                                metric = 'area')[[1]][[1]]

# What patch area do the ants prefer?
preferred_area <- terra::extract(forest_patch_area, ant_vect_buff) %>% 
      drop_na()

lower_threshold <- quantile(preferred_area$value, probs = area_percentile_threshold)

# Create forest patch layer above area threshold
forest_patches_sub <- ifel(forest_patch_area <= lower_threshold, NA, 1)

# Get forest patch IDs, mask for area above threshold
forest_patch_ID <- get_patches(forest_patches_sub, 
                             directions = 8)[[1]][[1]]
names(forest_patch_ID) <- 'patch_ID'

# 2. Overlay nest locations with forest patches to identify forest patches likely now occupied (ON) ----------------------------------------------------
ON_patch_IDs <- terra::extract(forest_patch_ID, ant_vect_buff) %>% 
      drop_na()

# Create binary raster showing forest patches where there is at least one nest
ON_patches <- forest_patch_ID$patch_ID %in% ON_patch_IDs$patch_ID

# 3. Create raster of forest patches suitable but likely not occupied now (SNO) -------------------------------------
# Mask suitable areas to forest patches
suitable_patches <- suitability_map %>% 
      resample(forest_patch_ID) %>% 
      mask(forest_patch_ID) %>% 
      clamp(lower = suitability_threshold, values = F)

# Patch areas that are suitable but likely not occupied. 
SNO_patches <- suitable_patches * (1-ON_patches) # Turns the TRUE/FALSE into 0/1
SNO_patches <- ifel(SNO_patches == 0, NA, SNO_patches)

# 4. SNO patches that might be colonised naturally, or might already be colonised (narrow forest gaps) -------------------------------------
# Create a mask to only keep relevant forest patches (ON or SNO)

# Get IDs of SNO patches
SNO_IDs <- get_patches(ifel(!is.na(SNO_patches), 1, NA), 
                                 directions = 8)[[1]][[1]]
names(SNO_IDs) <- 'patch_ID'

# Identify SNO patches that are within the distance threshold to ON patches
ON_patches_mask <- ifel(isFALSE(ON_patches), NA, 1)
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
dispersal_patches_binary <- ifel(is.na(dispersal_patches), 0, 1)

# 5. SNO patches that will likely not be colonised naturally (far from ON) -------------------------------------
# These patches are candidates for translocation
ON_patches_binary <- ifel(isFALSE(ON_patches), 0, 1)

ON_or_dispersal <- ON_patches_binary + dispersal_patches_binary 
ON_or_dispersal_mask <- ifel(ON_or_dispersal > 0, 1, NA)

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
translocation_patches_binary <- ifel(is.na(translocation_patches), 0, 1)

writeRaster(translocation_patches_binary, "translocation_patches_binary.tif", overwrite=T)
writeVector(ON_patches_buffered_far, "ON_patches_buffered_far.shp", overwrite=T)
