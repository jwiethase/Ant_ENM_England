# HEADER --------------------------------------------
#
# Author: Joris Wiethase
# Email: j.wiethase@gmail.com
# 
# Script Description:  
# Collate and process available species records for the study species. Extract sampling effort.

rm(list = ls())
library(tidyverse)
library(rnrfa)
library(terra)
library(lubridate)
library(tidyterra)
library(geosphere)
source('source/misc_functions.R')

ROI_27700 <- vect('spatial_other/ROI_outline_27700.shp')

# EXHAUSTIVE SURVEY DATA -------------------------------------------------
# Surveys where every single nest (or close to every nest) was recorded within a study area
## Procter et al. 2015 data -------------------------------------------------------------
procter_df <- read.csv('species_data/raw/Procter_nests.csv') %>% 
      vect(geom = c('x_proj', 'y_proj'), crs = crs('epsg:27700'), keepgeom = FALSE) %>% 
      terra::crop(ROI_27700) %>%
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
      terra::crop(ROI_27700) %>%
      terra::project(crs(km_proj)) %>%  
      as.data.frame(geom = 'XY') %>% 
      dplyr::select(species, x, y, source)

## Hardcastle data -------------------------------------------------------------
hardcastle_df <- vect('species_data/raw/2023_Anthills.gpkg') %>% 
      terra::crop(ROI_27700) %>%
      terra::project(crs(km_proj)) %>%  
      as.data.frame(geom = 'XY') %>% 
      mutate(species = "Formica lugubris",
             source = "hardcastle") %>% 
      dplyr::select(species, x, y, source)

## NYM data -------------------------------------------------------------
nym_1 <- vect('species_data/raw/NYM/All 2018 field season data_all sites.shp') %>% 
      project(crs("epsg:27700")) %>% 
      crop(ROI_27700) %>% 
      as.data.frame(geom = 'XY') %>% 
      mutate(date = date(ymd_hms(time))) %>% 
      dplyr::select(x, y, date)

nym_2 <- vect('species_data/raw/NYM/Broxa 2019.shp') %>% 
      project(crs("epsg:27700")) %>% 
      crop(ROI_27700) %>% 
      as.data.frame(geom = 'XY') %>% 
      mutate(date = date(ymd_hms(time))) %>% 
      dplyr::select(x, y, date)

nym_3 <- vect('species_data/raw/NYM/Cropton margins 2019_1.shp') %>% 
      project(crs("epsg:27700")) %>% 
      crop(ROI_27700) %>% 
      as.data.frame(geom = 'XY') %>% 
      mutate(date = date(ymd_hms(time))) %>% 
      dplyr::select(x, y, date)

nym_4 <- vect('species_data/raw/NYM/Cropton margins 2019_2.shp') %>% 
      project(crs("epsg:27700")) %>% 
      crop(ROI_27700) %>% 
      as.data.frame(geom = 'XY') %>% 
      mutate(date = date(ymd_hms(time))) %>% 
      dplyr::select(x, y, date)

nym_5 <- read.csv('species_data/raw/NYM/cropton_margins_nest_data_2020.csv') %>% 
      dplyr::select(Lng, Lat, Date.visited) %>% 
      drop_na() %>% 
      vect(geom = c('Lat', 'Lng'), crs = crs('epsg:4326')) %>% 
      terra::project(crs("epsg:27700")) %>% 
      crop(ROI_27700) %>% 
      as.data.frame(geom = 'XY') %>% 
      mutate(date = date(ymd_hms(Date.visited))) %>% 
      dplyr::select(x, y, date)

nym_6 <- read.csv('species_data/raw/NYM/X000 to C000 with lat and long.csv') %>% 
      vect(geom = c('Longitude', 'Latitude'), crs = crs('epsg:4326')) %>% 
      terra::project(crs("epsg:27700")) %>% 
      crop(ROI_27700) %>% 
      as.data.frame(geom = 'XY') %>% 
      mutate(date = NA) %>% 
      dplyr::select(x, y, date)

nym_7 <- vect('species_data/raw/NYM/nests_remapped.shp') %>% 
      project(crs("epsg:27700")) %>% 
      crop(ROI_27700) %>% 
      as.data.frame(geom = 'XY') %>% 
      mutate(date = NA) %>% 
      dplyr::select(x, y, date)

nym_8 <- vect('species_data/raw/NYM/Cropton nests 2020.gpx') %>% 
      project(crs("epsg:27700")) %>% 
      crop(ROI_27700) %>% 
      as.data.frame(geom = 'XY') %>% 
      mutate(date = NA) %>% 
      dplyr::select(x, y, date)

nym_9 <- vect('species_data/raw/NYM/Broxa margins 2020.gpx') %>% 
      project(crs("epsg:27700")) %>% 
      crop(ROI_27700) %>% 
      as.data.frame(geom = 'XY') %>% 
      mutate(date = NA) %>% 
      dplyr::select(x, y, date)

nym_10 <- vect('species_data/raw/NYM/NYM_nest_2022.gpx') %>% 
      project(crs("epsg:27700")) %>% 
      crop(ROI_27700) %>% 
      as.data.frame(geom = 'XY') %>% 
      mutate(date = NA) %>% 
      dplyr::select(x, y, date)

nym_df <- rbind(nym_1, nym_2, nym_3, nym_4, nym_5, nym_6, nym_7, nym_8, nym_9, nym_10) %>% 
      mutate(species = "Formica lugubris",
             source = "nym") %>% 
      vect(geom = c('x', 'y'), crs = crs('epsg:27700'), keepgeom = FALSE) %>% 
      terra::project(crs(km_proj)) %>% 
      as.data.frame(geom = 'XY') %>% 
      dplyr::select(species, x, y, source)

## Longshaw data -------------------------------------------------------------
## Small forest patches resampled frequently. Not technically exhaustive for the wider area,
## but for these small patches. Most suitable in exhaustive category due to dense sampling
longshaw_df <- read.csv('species_data/raw/FormicaLugubrisNestsLongshawEJHR2022.csv') %>% 
      vect(geom = c('lon', 'lat'), crs = crs('epsg:4326'), keepgeom = FALSE) %>% 
      terra::project(crs('epsg:27700')) %>% 
      terra::crop(ROI_27700) %>%       
      terra::project(crs(km_proj)) %>% 
      as.data.frame(geom = 'XY') %>% 
      mutate(species = "Formica lugubris", 
             source = 'Longshaw') %>% 
      dplyr::select(species, x, y, source)

## Tom Dallimore data -------------------------------------------------------------
## Quite dense. We have detailed effort information (days visited), but effort differs from sporadic data
## (sampled numerous nests in one forest in one day)
dallimore_df <- do.call(rbind, 
                        lapply(list.files("species_data/raw/dallimore_data/", pattern = ".csv", full.names = T), 
                               read_and_select, type = "csv", cols = c('Longitude', 'Latitude', 'Nest.ID')
                        )) %>% 
      vect(geom = c('Longitude', 'Latitude'), crs = crs('epsg:4326'), keepgeom = FALSE) %>% 
      terra::project(crs('epsg:27700')) %>% 
      terra::crop(ROI_27700) %>%       
      terra::project(crs(km_proj)) %>% 
      as.data.frame(geom = 'XY') %>% 
      mutate(species = "Formica rufa",
             source = 'dallimore') %>% 
      dplyr::select(species, x, y, source)

## Combine exhaustive data -------------------------------------------------------------
exhaustive_combined <- rbind(procter_df, gaitbarrows_df, hardcastle_df, nym_df, longshaw_df, dallimore_df)

# SPORADIC SURVEY DATA -------------------------------------------------
# Surveys where volunteers mostly recorded one nest and then moved on
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
      terra::crop(ROI_27700) %>% 
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
      terra::crop(ROI_27700) %>% 
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
      terra::crop(ROI_27700) %>% 
      terra::project(crs(km_proj)) %>% 
      as.data.frame(geom = 'XY') %>% 
      dplyr::select(species, x, y, date, OS_Tile, source)

# Apply thinning to date groups, ensuring that effort information is retained 
josie_thinned <- group_thin_spatial(josie_df, group_var = "date", dist_meters = 100, crs_spec = km_proj) %>% 
      dplyr::select(species, x, y, date, OS_Tile, source) 

## BWARS F. nitidulus data -------------------------------------------------------------
# F. nitidulus always occurs in wood ant nests. Can tell from latitude which wood ant, and derive their presence
# Find limits of F. lugubris and F. rufa presence records, add a somewhat arbitrary 25 km buffer, just in case
F_lugubris_lowest <- min(rbind(procter_df$y, 
                               longshaw_df$y,
                               bwars_formica_df$y[bwars_formica_df$species == "Formica lugubris"],
                               hymettus_df$y[hymettus_df$species == "Formica lugubris"],
                               josie_thinned$y[josie_thinned$species == "Formica lugubris"],
                               nym_df$y)) - 25 

F_rufa_highest <- max(rbind(bwars_formica_df$y[bwars_formica_df$species == "Formica rufa"],
                            hymettus_df$y[hymettus_df$species == "Formica rufa"],
                            josie_thinned$y[josie_thinned$species == "Formica rufa"],
                            dallimore_df$y), 
                            gaitbarrows_df$y) + 25 

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
      terra::crop(ROI_27700) %>% 
      terra::project(crs(km_proj)) %>% 
      thin_spatial(., 50) %>% 
      as.data.frame(geom = 'XY') %>% 
      dplyr::select(species, x, y, date, OS_Tile, source) %>% 
      mutate(species = case_when(y > F_rufa_highest ~ "Formica lugubris",
                                 y < F_lugubris_lowest ~ "Formica rufa")) %>% 
      filter(!is.na(species))

## Dallimore sporadic -------------------------------------------------------------
## Since we have detailed effort information, we can filter this to match the BWARS effort.
## This is to supplement the lgcp models based only on sporadic data
dallimore_effort <- do.call(rbind,
                            lapply(list.files("species_data/raw/dallimore_data/", pattern = ".xls", full.names = T),
                                   read_and_select, type = "xls", cols = c('Date visited', 'Nest ID')
                            )) %>%
      unique() %>%
      rename("Nest.ID" = "Nest ID")

dallimore_sporadic_df <- do.call(rbind, 
                        lapply(list.files("species_data/raw/dallimore_data/", pattern = ".csv", full.names = T), 
                               read_and_select, type = "csv", cols = c('Longitude', 'Latitude', 'Nest.ID')
                        )) %>% 
      merge(., dallimore_effort, by = 'Nest.ID') %>% 
      mutate(date = ymd(`Date visited`)) %>% 
      vect(geom = c('Longitude', 'Latitude'), crs = crs('epsg:4326'), keepgeom = FALSE) %>% 
      terra::project(crs('epsg:27700')) %>% 
      terra::crop(ROI_27700) %>%   
      as.data.frame(geom = 'XY') %>% 
      rowwise() %>% 
      mutate(OS_Tile = OSGB(c(x, y), "10km")) %>%  
      vect(geom = c('x', 'y'), crs = crs('epsg:27700'), keepgeom = FALSE) %>% 
      terra::project(crs(km_proj)) %>% 
      as.data.frame(geom = 'XY') %>% 
      mutate(species = "Formica rufa",
             source = 'dallimore') %>% 
      dplyr::select(species, x, y, date, OS_Tile, source) 

# Thin spatially within date groups. This ensures that sampling effort is properly captured.
dallimore_sporadic_thinned <- group_thin_spatial(dallimore_sporadic_df, group_var = "date", dist_meters = 500, crs_spec = km_proj) %>% 
      dplyr::select(species, x, y, date, OS_Tile, source) 

## Gaitbarrows sporadic -------------------------------------------------------------
## Same as dallimore, we have detailed effort information, so we can filter this to match the BWARS effort.
gaitbarrows_sporadic <- read.csv('species_data/raw/Formica_rufa_Survey_Data_2019_David_Lamin_Elva_Robinson.csv') %>% 
      mutate(Grid.Reference = gsub(" ", "", Grid.Reference), 
             easting = rnrfa::osg_parse(Grid.Reference)$easting,
             northing = rnrfa::osg_parse(Grid.Reference)$northing,
             accuracy = sapply(Grid.Reference, determine_accuracy),
             species = "Formica rufa",
             source = "gaitbarrows",
             date = dmy(Date)) %>% 
      rowwise() %>% 
      mutate(OS_Tile = OSGB(c(easting, northing), "10km")) %>%  
      vect(geom = c('easting', 'northing'), crs = crs('epsg:27700'), keepgeom = FALSE) %>% 
      terra::crop(ROI_27700) %>%
      terra::project(crs(km_proj)) %>%  
      as.data.frame(geom = 'XY') %>% 
      dplyr::select(species, x, y, date, OS_Tile, source)

gaitbarrows_sporadic_thinned <- group_thin_spatial(gaitbarrows_sporadic, group_var = "date", dist_meters = 250, crs_spec = km_proj) %>% 
      dplyr::select(species, x, y, date, OS_Tile, source) 

## Hardcastle sporadic -------------------------------------------------------------
hardcastle_sporadic <- vect('species_data/raw/2023_Anthills.gpkg') %>% 
      as.data.frame(geom = 'XY') %>% 
      mutate(species = "Formica lugubris",
             source = "hardcastle",
             date = date(ymd_hms(`Date and Time`))) %>% 
      rowwise() %>% 
      mutate(OS_Tile = OSGB(c(x, y), "10km")) %>% 
      # Some missing dates, fill with previous. Doesn't have to be exact
      fill(date, .direction = "down") %>%  
      vect(geom = c('x', 'y'), crs = crs('epsg:27700'), keepgeom = FALSE) %>% 
      terra::crop(ROI_27700) %>%
      terra::project(crs(km_proj)) %>%  
      as.data.frame(geom = 'XY') %>% 
      dplyr::select(species, x, y, date, OS_Tile, source)

hardcastle_sporadic_thinned <- group_thin_spatial(hardcastle_sporadic, group_var = "date", dist_meters = 250, crs_spec = km_proj) %>% 
      dplyr::select(species, x, y, date, OS_Tile, source) 

## NYM sporadic -------------------------------------------------------------
## Multiple separate data files that have nests and date information, convert to sporadic
nym_sporadic <- rbind(nym_1, nym_2, nym_3, nym_4, nym_5) %>% 
      as.data.frame(geom = 'XY') %>% 
      drop_na() %>% 
      rowwise() %>% 
      mutate(OS_Tile = OSGB(c(x, y), "10km"),
             species = "Formica lugubris",
             source = "NYM") %>% 
      vect(geom = c('x', 'y'), crs = crs('epsg:27700'), keepgeom = FALSE) %>% 
      terra::crop(ROI_27700) %>%
      terra::project(crs(km_proj)) %>%  
      as.data.frame(geom = 'XY') %>% 
      dplyr::select(species, x, y, date, OS_Tile, source) 

nym_sporadic_thinned <- group_thin_spatial(nym_sporadic, group_var = "date", dist_meters = 250, crs_spec = km_proj) %>% 
      vect(geom = c('x', 'y'), crs = crs(km_proj)) %>% 
      # Remove potential re-maps
      thin_spatial(., dist_meters = 20) %>%  
      as.data.frame(geom = 'XY') %>% 
      dplyr::select(species, x, y, date, OS_Tile, source) 

## Create sampling effort -------------------------------------------------------------
# Count number of individual days a site was sampled. For this, also add the separate samples
OS_grid <- vect('spatial_other/OSGB_Grid_10km.shp') %>% 
      terra::project(crs(km_proj))

bwars_all_eff <- read.csv('species_data/raw/20230404 BWARS public data.csv') %>% 
      mutate(date = ymd(upper_date)) %>% 
      filter(date > '2009-01-01') %>% 
      rename(OS_Tile = osgr_to_gridref) %>% 
      dplyr::select(OS_Tile, date)

# The dataset for lgcp models contains exhaustive survey data that was thinned to 
# resemble BWARS data. The sampling effort from this is only suitable for lgcp models
sporadic_combined <- rbind(bwars_formica_df, hymettus_df, josie_df, bwars_nitidulus_df, 
                           dallimore_sporadic_thinned, gaitbarrows_sporadic_thinned, 
                           hardcastle_sporadic_thinned, nym_sporadic_thinned) %>%
      vect(geom = c('x', 'y'), crs = crs(km_proj)) %>%  
      thin_spatial(., 20) %>% # Thin closely clustered points, potential duplicates
      as.data.frame(geom = 'XY') 
      
all_eff_lgcp <- rbind(sporadic_combined %>% 
                            filter(source != "Nitidulus") %>% 
                            dplyr::select(OS_Tile, date),
                      bwars_all_eff %>% 
                            dplyr::select(OS_Tile, date)) %>% 
      group_by(OS_Tile) %>% 
      summarise(days_sampled = length(unique(date)))

effort_shapefile_lgcp <- merge(OS_grid, all_eff_lgcp, by.x = 'TILE_NAME', by.y = 'OS_Tile', all.x = TRUE) %>% 
      dplyr::select(days_sampled, TILE_NAME)
effort_shapefile_lgcp$days_sampled[is.na(effort_shapefile_lgcp$days_sampled)] <- 0

writeVector(effort_shapefile_lgcp, 'spatial_other/effort_lgcp_10km.shp', overwrite=TRUE)

all_eff_integrated <- rbind(bwars_all_eff, 
                      hymettus_df %>% dplyr::select(OS_Tile, date), 
                      josie_df %>% dplyr::select(OS_Tile, date)) %>% 
      group_by(OS_Tile) %>% 
      summarise(days_sampled = length(unique(date)))

effort_shapefile_integrated <- merge(OS_grid, all_eff_integrated, by.x = 'TILE_NAME', by.y = 'OS_Tile', all.x = TRUE) %>% 
      dplyr::select(days_sampled, TILE_NAME)
effort_shapefile_integrated$days_sampled[is.na(effort_shapefile_integrated$days_sampled)] <- 0

writeVector(effort_shapefile_integrated, 'spatial_other/effort_integrated_10km.shp', overwrite=TRUE)

sporadic_effort_added <- sporadic_combined %>% 
      merge(all_eff_lgcp, all.x = T) 
      
# Export -------------------------------------------------------------
writeVector(exhaustive_combined %>% vect(geom = c('x', 'y'), crs = crs(km_proj)),
            'species_data/processed_shp/exhaustive_combined.shp', overwrite=T)
writeVector(sporadic_effort_added %>% vect(geom = c('x', 'y'), crs = crs(km_proj)),
            'species_data/processed_shp/sporadic_combined.shp', overwrite=T)

write.csv(exhaustive_combined, 'species_data/processed_csv/exhaustive_combined.csv')
write.csv(sporadic_effort_added, 'species_data/processed_csv/sporadic_combined.csv')

writeVector(nym_df %>% vect(geom = c('x', 'y'), crs = crs(km_proj)),
            'nym_df.shp', overwrite=T)
