#######################################################################################
# Identify Forestry England-managed patches that might be colonised in the future 
#######################################################################################
library(landscapemetrics)
library(terra)
library(tidyterra)
library(tidyverse)
library(data.table)
setwd('/users/jhw538/scratch/ant_modelling')
source('source/misc_functions.R')

args <- commandArgs(trailingOnly=TRUE)
job <- as.integer(args[1])

# Parameters ----------------------------------------------------
species_choices = c('Formica rufa', 'Formica lugubris')
model_choices <- c('maxent_V3', 'lgcp')
max_trans_distances <- c(0, 2, 5, 25)
ON_thresholds <- c(0.25, 0.5)
SNO_thresholds <- c(0.5, 0.75)

mult_combs <- crossing(species_choices, max_trans_distances, ON_thresholds, SNO_thresholds, model_choices) # 64
comb_values <- mult_combs[job, ]
print(comb_values)

species_choice = comb_values$species_choices
max_gap_translocation = comb_values$max_trans_distances
ON_threshold = comb_values$ON_thresholds
SNO_threshold = comb_values$SNO_thresholds
print(paste('Starting', species_choice))

model_choice = comb_values$model_choices
point_buffer = 0.05                  
max_gap_dispersal = 0.1    
translocation_threshold = 0.1

dir.create(paste0('model_out/', gsub(' ', '_', species_choice),'/', model_choice, '/maxTransDist_', max_gap_translocation, 'km'), showWarnings = F)

# Load and prepare data ----------------------------------------------------
ROI <- vect('data/ROI_outline_27700.shp') %>% 
      terra::project(crs(km_proj))
FE_managed <- vect('data/Forestry_England_managed_forest.shp') %>% 
      project(crs(km_proj))
forest_stack <- rast('data/forest_stack_30m.tif') %>%
      terra::project(crs('epsg:27700'))

if(model_choice == 'maxent_V3' & species_choice == 'Formica rufa'){
      suitability_map <- rast('model_out/Formica_rufa/maxent_V3/V3_sporadic_Formica_rufa_all30m_thin0m.tif') 
}

if(model_choice == 'maxent_V3' & species_choice == 'Formica lugubris'){
      suitability_map <- rast('model_out/Formica_lugubris/maxent_V3/V3_sporadic_Formica_lugubris_all30m_thin100m.tif')
}

if(model_choice == 'lgcp' & species_choice == 'Formica rufa'){
      suitability_map <- rast('model_out/Formica_rufa/lgcp/Formica_rufa_lgcp_all30m_mosaic.tif') 
}

if(model_choice == 'lgcp' & species_choice == 'Formica lugubris'){
      suitability_map <- rast('model_out/Formica_lugubris/lgcp/Formica_lugubris_lgcp_all30m_mosaic.tif') 
}

if(model_choice == 'lgcp'){
      forest_mask_buff <- rast("data/forest_mask_buff_30m.tif") %>%
            terra::resample(suitability_map, method = "max") %>% 
            subst(0, NA)
      
      suitability_map <- suitability_map %>% 
            mask(forest_mask_buff) %>% 
            select(q0.5)     
}

sporadic <- read.csv('data/sporadic_combined.csv') %>% 
      filter(source != 'dallimore', source != 'nym', source != 'gaitbarrows', source != 'hardcastle') %>% 
      dplyr::select(x, y, species)
exhaustive <- read.csv('data/exhaustive_combined.csv') %>% 
      dplyr::select(x, y, species)

combined_presences <- rbind(sporadic, exhaustive) %>% 
      vect(geom = c('x', 'y'), crs = crs(km_proj)) 

ant_vect <- combined_presences %>% 
      filter(species == species_choice)
ant_vect_buff <- terra::buffer(ant_vect, width = point_buffer)

if(model_choice == 'lgcp'){
      suitability_map <- suitability_map %>% 
            normalise_raster() %>% 
            mask(ROI) %>% 
            tidyterra::select(q0.5) 
}


# 1. Identify areas likely now occupied (ON) ----------------------------------------------------
suitable_forest_forON <- suitability_map %>% 
      clamp(lower = ON_threshold, values = F)
suitable_forest_forON_mask <- ifel(is.na(suitable_forest_forON), NA, 1)

suitable_forest_forON_ID <- get_patches(suitable_forest_forON_mask, 
                                        directions = 8)[[1]][[1]]
names(suitable_forest_forON_ID) <- 'patch_ID'

ON_patch_IDs <- terra::extract(suitable_forest_forON_ID, ant_vect_buff) %>% 
      drop_na()

# Create binary raster showing forest patches where there is at least one nest
ON_patches_binary <- suitable_forest_forON_ID$patch_ID %in% ON_patch_IDs$patch_ID

ON_patches_mask <- ifel(isFALSE(ON_patches_binary), NA, 1)

# Create graduated raster showing suitability scores in these ON patches
ON_patches_gradient <- suitable_forest_forON %>% 
      mask(ON_patches_mask)

writeRaster(ON_patches_mask, 
            paste0('model_out/', gsub(' ', '_', species_choice), '/', model_choice, '/', gsub(' ', '_', species_choice), 
                   '_', model_choice,
                   '_pointBuff_', point_buffer,
                   '_ON_Thresh', ON_threshold,
                   '_SNO_Thresh', SNO_threshold,
                   '_ON_patches_mask.tif'), 
            overwrite = T)

writeRaster(ON_patches_gradient, 
            paste0('model_out/', gsub(' ', '_', species_choice), '/', model_choice, '/', gsub(' ', '_', species_choice), 
                   '_', model_choice,
                   '_pointBuff_', point_buffer,
                   '_ON_Thresh', ON_threshold,
                   '_SNO_Thresh', SNO_threshold,
                   '_ON_patches_gradient.tif'), 
            overwrite = T)

print('ON patches done')

# 2. Create raster of forest patches suitable but likely not occupied now (SNO) -------------------------------------
# What forest patch area do the ants prefer?
forest_mask <- ifel(forest_stack$cover_VOM < 0.3, NA, 1)

forest_patch_area <- spatialize_lsm(forest_mask,
                                    level = 'patch',
                                    metric = 'area')[[1]][[1]]

preferred_area <- terra::extract(forest_patch_area, ant_vect_buff) %>%
      drop_na()

# Find forest patch area that 1% or less of nest records occur below
lower_threshold <- quantile(preferred_area$value, probs = 0.01)  

forest_area_mask <- ifel(forest_patch_area <= lower_threshold, 0, 1) %>% 
      terra::project(suitability_map) %>% 
      subst(0, NA)

# Get suitable areas only, for SNO patches. 
suitable_forest_forSNO <- suitability_map %>% 
      clamp(lower = SNO_threshold, values = F) %>% 
      mask(forest_area_mask)

# Patch areas that are suitable but likely not occupied. 
SNO_patches <- suitable_forest_forSNO * (1-ON_patches_binary) # Turns the TRUE/FALSE into 0/1
SNO_patches <- ifel(SNO_patches == 0, NA, SNO_patches)

print('SNO patches done')

# 3. SNO patches that might be colonised naturally, or might already be colonised (narrow forest gaps) -------------------------------------
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

dispersal_patches_gradient <- suitable_forest_forSNO %>% 
      mask(dispersal_patches_mask)

writeRaster(dispersal_patches_mask, 
            paste0('model_out/', gsub(' ', '_', species_choice), '/', model_choice, '/', gsub(' ', '_', species_choice), 
                   '_', model_choice,
                   '_pointBuff_', point_buffer,
                   '_ON_Thresh', ON_threshold,
                   '_SNO_Thresh', SNO_threshold,
                   '_disp_mask_', max_gap_dispersal*1000, 'm.tif'), 
            overwrite = T)

writeRaster(dispersal_patches_gradient, 
            paste0('model_out/', gsub(' ', '_', species_choice), '/', model_choice, '/', gsub(' ', '_', species_choice),
                   '_', model_choice,
                   '_pointBuff_', point_buffer,
                   '_ON_Thresh', ON_threshold,
                   '_SNO_Thresh', SNO_threshold,
                   '_disp_gradient_', max_gap_dispersal*1000, 'm.tif'), 
            overwrite = T)

print('Dispersal patches done')

# 4. SNO patches that will likely not be colonised naturally (far from ON) -------------------------------------
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

translocation_patches_gradient <- suitability_map %>% 
      clamp(lower = translocation_threshold, values = F) %>% 
      mask(translocation_patches_mask)

writeRaster(translocation_patches_mask, 
            paste0('model_out/', gsub(' ', '_', species_choice), '/', model_choice, '/maxTransDist_', max_gap_translocation, 'km/', gsub(' ', '_', species_choice), 
                   '_', model_choice,
                   '_pointBuff_', point_buffer,
                   '_ON_Thresh', ON_threshold,
                   '_SNO_Thresh', SNO_threshold,
                   '_trans_mask_', max_gap_translocation, 'km.tif'), 
            overwrite = T)
writeRaster(translocation_patches_gradient, 
            paste0('model_out/', gsub(' ', '_', species_choice), '/', model_choice, '/maxTransDist_', max_gap_translocation, 'km/', gsub(' ', '_', species_choice), 
                   '_', model_choice,
                   '_pointBuff_', point_buffer,
                   '_ON_Thresh', ON_threshold,
                   '_SNO_Thresh', SNO_threshold,
                   '_trans_gradient_', max_gap_translocation, 'km.tif'), 
            overwrite = T)

print('Translocation patches done')

translocation_FE <- translocation_patches %>% 
      crop(FE_managed) %>% 
      mask(FE_managed)

translocation_mask_FE <- translocation_patches_mask %>% 
      crop(FE_managed) %>% 
      mask(FE_managed)

translocation_gradient_FE <- translocation_patches_gradient %>% 
      crop(FE_managed) %>% 
      mask(FE_managed)

writeRaster(translocation_FE, 
            paste0('model_out/', gsub(' ', '_', species_choice), '/', model_choice, '/maxTransDist_', max_gap_translocation, 'km/', gsub(' ', '_', species_choice), 
                   '_', model_choice,
                   '_pointBuff_', point_buffer,
                   '_ON_Thresh', ON_threshold,
                   '_SNO_Thresh', SNO_threshold,
                   '_transFE_', max_gap_translocation, 'km.tif'), 
            overwrite = T)
writeRaster(translocation_mask_FE, 
            paste0('model_out/', gsub(' ', '_', species_choice), '/', model_choice, '/maxTransDist_', max_gap_translocation, 'km/', gsub(' ', '_', species_choice), 
                   '_', model_choice,
                   '_pointBuff_', point_buffer,
                   '_ON_Thresh', ON_threshold,
                   '_SNO_Thresh', SNO_threshold,
                   '_transMask_FE_', max_gap_translocation, 'km.tif'), 
            overwrite = T)

writeRaster(translocation_gradient_FE, 
            paste0('model_out/', gsub(' ', '_', species_choice), '/', model_choice, '/maxTransDist_', max_gap_translocation, 'km/', gsub(' ', '_', species_choice), 
                   '_', model_choice,
                   '_pointBuff_', point_buffer,
                   '_ON_Thresh', ON_threshold,
                   '_SNO_Thresh', SNO_threshold,
                   '_transGradient_FE_', max_gap_translocation, 'km.tif'), 
            overwrite = T)
print('Translocation patches FE done')


