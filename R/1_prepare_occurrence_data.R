rm(list = ls())
library(tidyverse)
library(rnrfa)
library(terra)
library(lubridate)
library(tidyterra)
library(geosphere)
source('source/misc_functions.R')

ROI <- vect('spatial_other/ROI_outline_27700.shp')

# EXHAUSTIVE SURVEY DATA -------------------------------------------------
## Procter et al. 2015 data -------------------------------------------------------------
procter_df <- read.csv('species_data/raw/Procter_nests.csv') %>% 
      vect(geom = c('x_proj', 'y_proj'), crs = crs('epsg:27700'), keepgeom = FALSE) %>% 
      terra::crop(ROI) %>%
      terra::project(crs(km_proj)) %>%
      as.data.frame(geom = 'XY') %>% 
      mutate(source = "procter")

## David and Elva 2019 data -------------------------------------------------------------
gaitbarrows_df <- read.csv('species_data/raw/Formica_rufa_Survey_Data_2019_David_Lamin_Elva_Robinson.csv') %>% 
      mutate(Grid.Reference = gsub(" ", "", Grid.Reference), 
             easting = rnrfa::osg_parse(Grid.Reference)$easting,
             northing = rnrfa::osg_parse(Grid.Reference)$northing,
             accuracy = sapply(Grid.Reference, determine_accuracy),
             species = "Formica rufa",
             source = "gaitbarrows",
             date = dmy(Date)) %>% 
      vect(geom = c('easting', 'northing'), crs = crs('epsg:27700'), keepgeom = FALSE) %>% 
      terra::crop(ROI) %>%
      terra::project(crs(km_proj)) %>%  
      as.data.frame(geom = 'XY') %>% 
      dplyr::select(species, x, y, source)

## Hardcastle data -------------------------------------------------------------
hardcastle_df <- vect('species_data/raw/2023_Anthills.gpkg') %>% 
      terra::crop(ROI) %>%
      terra::project(crs(km_proj)) %>%  
      as.data.frame(geom = 'XY') %>% 
      mutate(species = "Formica lugubris",
             source = "hardcastle") %>% 
      dplyr::select(species, x, y, source)

## NYM data -------------------------------------------------------------
nym_df <- vect('species_data/raw/NYM_27700.shp') %>% 
      terra::crop(ROI) %>%
      terra::project(crs(km_proj)) %>% 
      as.data.frame(geom = 'XY') %>% 
      mutate(species = "Formica lugubris",
             source = "nym") %>% 
      dplyr::select(species, x, y, source)

## Longshaw data -------------------------------------------------------------
longshaw_df <- read.csv('species_data/raw/FormicaLugubrisNestsLongshawEJHR2022.csv') %>% 
      vect(geom = c('lon', 'lat'), crs = crs('epsg:4326'), keepgeom = FALSE) %>% 
      terra::project(crs('epsg:27700')) %>% 
      terra::crop(ROI) %>%       
      terra::project(crs(km_proj)) %>% 
      as.data.frame(geom = 'XY') %>% 
      mutate(species = "Formica lugubris",source = 'Longshaw') %>% 
      dplyr::select(species, x, y, source)

## Combine exhaustive data -------------------------------------------------------------
exhaustive_combined <- rbind(procter_df, gaitbarrows_df, hardcastle_df, nym_df, longshaw_df)

# SPORADIC SURVEY DATA -------------------------------------------------
## BWARS data -------------------------------------------------------------
bwars_formica_df <- read.csv('species_data/raw/20240105_bwars_data.csv') %>% 
      mutate(easting = rnrfa::osg_parse(gridref)$easting,
             northing = rnrfa::osg_parse(gridref)$northing,
             accuracy = sapply(gridref, determine_accuracy),
             species = sapply(strsplit(binomial, ':'), function(x) x[1]),
             date = ymd(upper_date),
             source = "BWARS") %>% 
      rowwise() %>% 
      mutate(OS_Tile = OSGB(c(easting, northing), "10km")) %>% 
      dplyr::select(-binomial, -lower_date, -upper_date, -gridref) %>% 
      # Apply filters, only high resolution and fairly recent
      filter(accuracy <= 100, date > '2010-01-01', species %in% c('Formica rufa', 'Formica lugubris')) %>% 
      # Exclude records outside England
      vect(geom = c('easting', 'northing'), crs = crs('epsg:27700'), keepgeom = FALSE) %>% 
      terra::crop(ROI) %>% 
      terra::project(crs(km_proj)) %>% 
      as.data.frame(geom = 'XY') %>% 
      dplyr::select(species, x, y, date, OS_Tile, source)

## Hymettus survey 2011 data -------------------------------------------------------------
# These are historical sites re-sampled, so extension of BWARS data.
# Definitely previously only sampled pre-2010, so there should not be duplicates with BWARS data
hymettus_df <- read.csv('species_data/raw/HymettusSurveyEJHR.csv') %>% 
      mutate(easting = rnrfa::osg_parse(Gridref)$easting,
             northing = rnrfa::osg_parse(Gridref)$northing,
             accuracy = sapply(Gridref, determine_accuracy),
             date = dmy(Date),
             species = Taxon,
             source = "Hymettus") %>% 
      rowwise() %>% 
      mutate(OS_Tile = OSGB(c(easting, northing), "10km")) %>% 
      filter(accuracy <= 100) %>% 
      vect(geom = c('easting', 'northing'), crs = crs('epsg:27700'), keepgeom = FALSE) %>% 
      terra::crop(ROI) %>% 
      terra::project(crs(km_proj)) %>% 
      as.data.frame(geom = 'XY') %>% 
      dplyr::select(species, x, y, date, OS_Tile, source)

## Josie Monaghan survey 2018-2019 data -------------------------------------------------------------
# Also re-sampled BWARS sites. These are uncharacteristically close together sometimes though, 
# thin with 100m buffer
josie_1 <- read.csv('species_data/raw/josie_survey.csv') %>% 
      filter(Position != '', str_detect(Position, '[a-zA-Z0-9]')) %>% 
      dplyr::select(Species, Position, Date)

josie_df <- read.csv('species_data/raw/josie_survey_2.csv') %>% 
      filter(Position != '', str_detect(Position, '[a-zA-Z0-9]')) %>% 
      dplyr::select(Species, Position, Date) %>% 
      rbind(josie_1) %>%  
      mutate(Position = gsub('[^a-zA-Z0-9]', '', Position),
             easting = rnrfa::osg_parse(Position)$easting,
             northing = rnrfa::osg_parse(Position)$northing,
             Date = ifelse(is.na(Date) | Date == "", '20.07.18', Date),  # Doesn't have to be accurate, probably all from same day
             date = dmy(Date),
             species = case_when(Species == 'ru' ~ 'Formica rufa',
                                 Species == 'rufa' ~ 'Formica rufa',
                                 Species == 'lu' ~ 'Formica lugubris'),
             source = "Josie") %>% 
      rowwise() %>% 
      mutate(OS_Tile = OSGB(c(easting, northing), "10km")) %>% 
      filter(!is.na(species)) %>% 
      vect(geom = c('easting', 'northing'), crs = crs('epsg:27700'), keepgeom = FALSE) %>% 
      terra::crop(ROI) %>% 
      terra::project(crs(km_proj)) %>% 
      as.data.frame(geom = 'XY') %>% 
      dplyr::select(species, x, y, date, OS_Tile, source)

## BWARS F. nitidulus data -------------------------------------------------------------
# F. nitidulus always occurs in wood ant nest. Can tell from latitude which wood ant, and derive their presence
# Find limits of F. lugubris and F. rufa presence records, add a somewhat arbitrary 25 km buffer, just in case
F_lugubris_lowest <- min(rbind(procter_df$y, 
                               longshaw_df$y,
                               bwars_formica_df$y[bwars_formica_df$species == "Formica lugubris"],
                               hymettus_df$y[hymettus_df$species == "Formica lugubris"],
                               josie_df$y[josie_df$species == "Formica lugubris"])) - 25 

F_rufa_highest <- max(rbind(bwars_formica_df$y[bwars_formica_df$species == "Formica rufa"],
                            hymettus_df$y[hymettus_df$species == "Formica rufa"],
                            josie_df$y[josie_df$species == "Formica rufa"])) + 25 

bwars_nitidulus_df <- read.csv('species_data/raw/Formicoxenus nitidulus records BWARS to 2021.csv') %>% 
      mutate(easting = rnrfa::osg_parse(gridref)$easting,
             northing = rnrfa::osg_parse(gridref)$northing,
             accuracy = sapply(gridref, determine_accuracy),
             species = sapply(strsplit(binomial, ':'), function(x) x[1]),
             date = ymd(upper_date),
             source = "Nitidulus") %>% 
      rowwise() %>% 
      mutate(OS_Tile = OSGB(c(easting, northing), "10km")) %>% 
      dplyr::select(-binomial, -lower_date, -upper_date, -gridref) %>% 
      filter(accuracy <= 100, date > '2010-01-01') %>% 
      vect(geom = c('easting', 'northing'), crs = crs('epsg:27700'), keepgeom = FALSE) %>% 
      terra::crop(ROI) %>% 
      terra::project(crs(km_proj)) %>% 
      as.data.frame(geom = 'XY') %>% 
      dplyr::select(species, x, y, date, OS_Tile, source) %>% 
      mutate(species = case_when(y > F_rufa_highest ~ "Formica lugubris",
                                 y < F_lugubris_lowest ~ "Formica rufa")) %>% 
      filter(!is.na(species))
  
## Create sampling effort -------------------------------------------------------------
# Count number of individual days a site was sampled. For this, also add the separate samples
bwars_all_eff <- read.csv('species_data/raw/20230404 BWARS public data.csv') %>% 
      mutate(date = ymd(upper_date)) %>% 
      filter(date > '2009-01-01') %>% 
      rename(OS_Tile = osgr_to_gridref) %>% 
      dplyr::select(OS_Tile, date)

all_eff <- rbind(bwars_all_eff, 
                 hymettus_df %>% dplyr::select(OS_Tile, date), 
                 josie_df %>% dplyr::select(OS_Tile, date)) %>% 
      group_by(OS_Tile) %>% 
      summarise(days_sampled = length(unique(date)))

OS_grid <- vect('spatial_other/OSGB_Grid_10km.shp') %>% 
      terra::crop(ROI) %>%       
      terra::project(crs(km_proj))

effort_shapefile <- merge(OS_grid, all_eff, by.x = 'TILE_NAME', by.y = 'OS_Tile', all.x = TRUE) %>% 
      dplyr::select(days_sampled, TILE_NAME)

writeVector(effort_shapefile, 'spatial_other/effort_10km.shp', overwrite=TRUE)

## Combine sporadic data -------------------------------------------------------------
sporadic_combined <- rbind(bwars_formica_df, hymettus_df, josie_df, bwars_nitidulus_df) %>% 
      merge(all_eff, all.x = T)

# Export -------------------------------------------------------------
writeVector(exhaustive_combined %>% vect(geom = c('x', 'y'), crs = crs(km_proj)),
            'species_data/processed_shp/exhaustive_combined.shp', overwrite=T)
writeVector(sporadic_combined %>% vect(geom = c('x', 'y'), crs = crs(km_proj)),
            'species_data/processed_shp/sporadic_combined.shp', overwrite=T)

write.csv(exhaustive_combined, 'species_data/processed_csv/exhaustive_combined.csv')
write.csv(sporadic_combined, 'species_data/processed_csv/sporadic_combined.csv')



