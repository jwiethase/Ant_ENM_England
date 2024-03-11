#######################################################################################
# Identify Forestry England-managed patches that might be colonised in the future 
#######################################################################################
rm(list = ls())
library(landscapemetrics)
library(terra)
library(tidyterra)
source('source/misc_functions.R')

# Parameters ----------------------------------------------------
species_choice = 'Formica rufa'
model_choice = 'lgcp'
point_buffer = 5                    # Buffer around nest records in meter, to account for accuracy
suitability_threshold = 0.2         # Threshold value for suitability
max_gap_dist = 600                  # Distance in meters that is considered likely to be dispersed across

# Load and prepare data ----------------------------------------------------
ROI_27700 <- vect('spatial_other/ROI_outline_27700.shp') 
FE_managed <- vect('spatial_other/Forestry_England_managed_forest.shp')
forest_stack <- rast('covariates/processed/forest_stack_30m.tif')  %>% 
      terra::project(crs('epsg:27700')) 

if(model_choice == 'lgcp'){
      suitability_map <- rast(paste0('model_out/', gsub(' ', '_', species_choice), '/lgcp/', gsub(' ', '_', species_choice), '_suitability_adj.tif')) %>% 
            terra::project(crs('epsg:27700')) 
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

lower_threshold <- quantile(preferred_area$value, probs = 0.005)

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

# 4. SNO patches that might be colonised (narrow forest gaps) -------------------------------------
# Create a mask to only keep relevant forest patches (ON or SNO)
SNO_patches_binary <- ifel(is.na(SNO_patches) | SNO_patches == 0, 0, 1)
ON_patches_binary <- ifel(isFALSE(ON_patches), 0, 1)
all_relevant_patches <- SNO_patches_binary + ON_patches_binary

# Get patch IDs of relevant patches
relevant_patch_ID <- get_patches(subst(all_relevant_patches, 0, NA), 
                               directions = 8)[[1]][[1]]
names(relevant_patch_ID) <- 'patch_ID'

# Identify patch IDs corresponding to ON or SNO
SNO_IDs <- relevant_patch_ID %>% 
      mask(subst(SNO_patches_binary, 0, NA))

ON_IDs <- relevant_patch_ID %>% 
      mask(subst(ON_patches_binary, 0, NA))

# Get all distances
relevant_patch_dist <- get_nearestneighbour(relevant_patch_ID, return_id = TRUE) %>% 
      mutate(dist = round(dist))

# Identify if paired patches are SNO or ON
relevant_patch_dist$patch_A[relevant_patch_dist$id %in% values(SNO_IDs$patch_ID, na.rm = T)] <- 'SNO'
relevant_patch_dist$patch_A[relevant_patch_dist$id %in% values(ON_IDs$patch_ID, na.rm = T)] <- 'ON'

relevant_patch_dist$patch_B[relevant_patch_dist$id_neighbour %in% values(SNO_IDs$patch_ID, na.rm = T)] <- 'SNO'
relevant_patch_dist$patch_B[relevant_patch_dist$id_neighbour %in% values(ON_IDs$patch_ID, na.rm = T)] <- 'ON'

# Get paired patches that were below distance threshold and SNO with ON
possible_dispersal_pairs <- relevant_patch_dist %>% 
      filter(dist <= max_gap_dist,
             patch_A != patch_B)

possible_dispersal_pairs$pair_id <- apply(possible_dispersal_pairs[, c('id', 'id_neighbour')], 1, function(x) paste(sort(x), collapse = '_'))
SNO_IDs_for_dispersal <- possible_dispersal_pairs[!duplicated(possible_dispersal_pairs$pair_id), ] %>% 
      mutate(sno_id = case_when(
            patch_A == 'SNO' ~ as.character(id),
            patch_B == 'SNO' ~ as.character(id_neighbour)
      ))

SNO_patches_for_dispersal <- SNO_IDs %>% 
      filter(patch_ID %in% SNO_IDs_for_dispersal$sno_id) %>% 
      mask(FE_managed)

# 5. Export -------------------------------------
writeRaster(SNO_patches_for_dispersal, paste0('model_out/', gsub(' ', '_', species_choice), '/', model_choice, '/', gsub(' ', '_', species_choice), '_SNO_patches_for_dispersal.tif'))
