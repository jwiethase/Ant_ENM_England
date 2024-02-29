library(tidyverse)
library(tidyterra)
library(terra)
library(viridis)
library(patchwork)
library(ggtext)
source('source/misc_functions.R')

ROI <- vect('spatial_other/ROI_kmproj.shp')
small_ROI <- vect('spatial_other/ROI_small_WGS84.shp') %>% 
      terra::project(crs(km_proj)) 

FE <- vect('spatial_other/Forestry_England_managed_forest.shp')
effort_rast_10km <- rast('covariates/raw/effort_raster.tif') %>% 
      clamp(lower = 0, values = T) %>% 
      terra::project(crs(km_proj)) %>% 
      crop(ROI) %>% 
      mask(ROI) %>% 
      +1 %>% 
      log()

climate_stack <- rast("covariates/processed/climate_stack_1km.tif")
stack_30m <- rast("covariates/processed/stack_30m.tif") %>% 
      crop(small_ROI)

# Observations raw ------
sporadic <- read.csv('species_data/processed_csv/sporadic_combined.csv') %>% 
      mutate(structure = "sporadic") %>% 
      dplyr::select(x, y, species, source, structure)

exhaustive <- read.csv('species_data/processed_csv/exhaustive_combined.csv') %>% 
      mutate(structure = "exhaustive") %>% 
      dplyr::select(x, y, species, source, structure)

all_obs <- rbind(sporadic, exhaustive) %>% 
      filter(species != 'Formicoxenus nitidulus') %>% 
      vect(geom = c("x", "y"), crs = crs(km_proj), keepgeom = TRUE)

all_obs_sub <- rbind(sporadic, exhaustive) %>% 
      filter(species != 'Formicoxenus nitidulus') %>% 
      vect(geom = c("x", "y"), crs = crs(km_proj), keepgeom = TRUE) %>% 
      crop(small_ROI)

obs_overview <- ggplot() +
      geom_spatvector(data = ROI, fill = "transparent") +
      # geom_spatvector(data = FE, fill = "transparent", alpha = 0.2) +
      geom_spatvector(data = all_obs, aes(pch = structure, col = source)) +
      theme_bw() +
      facet_grid(~ species) +
      scale_shape_manual(values = c(0, 3), name = "Sampling type") +
      scale_color_discrete(name = "Data source")

obs_overview_sub <- ggplot() +
      geom_spatvector(data = all_obs_sub, aes(pch = structure, col = source)) +
      theme_bw() +
      facet_grid(~ species) +
      scale_shape_manual(values = c(0, 3), name = "Sampling type") +
      scale_color_discrete(name = "Data source")

effort_plot <- ggplot() +
      geom_spatraster(data = effort_rast_10km) +
      theme_bw() +
      scale_fill_viridis(na.value = "transparent", name = "log(Visits)") +
      ggtitle("<span style='font-size: 10pt;'>Sampling effort: Number of volunteer visits</font>") +
      theme(plot.title = element_markdown())
      
ggsave(plot = obs_overview, filename = "figures/presentation/obs_overview.png",
       dpi = 300, width = 20, height = 15, units = "cm")
ggsave(plot = obs_overview_sub, filename = "figures/presentation/obs_overview_sub.png",
       dpi = 300, width = 15, height = 15, units = "cm")
ggsave(plot = effort_plot, filename = "figures/presentation/effort_plot.png",
       dpi = 300, width = 10, height = 15, units = "cm")

# Climate raw ------
soil_temp <- ggplot() +
      geom_spatraster(data = climate_stack %>% subset(1)) +
      theme_bw() +
      scale_fill_viridis(na.value = "transparent", name = "°C", option = "A") +
      ggtitle("<span style='font-size: 10pt;'>Soil temperature maximum of hottest month</font>") +
      theme(plot.title = element_markdown())

soil_temp_season <- ggplot() +
      geom_spatraster(data = climate_stack %>% subset(2)) +
      theme_bw() +
      scale_fill_viridis(na.value = "transparent", name = "sd(°C)", option = "A") +
      ggtitle("<span style='font-size: 10pt;'>Soil temperature seasonality</font>") +
      theme(plot.title = element_markdown())

rain_coldest <- ggplot() +
      geom_spatraster(data = climate_stack %>% subset(3)) +
      theme_bw() +
      scale_fill_viridis(na.value = "transparent", name = "mm", option = "G", direction = -1) +
      ggtitle("<span style='font-size: 10pt;'>Total rain coldest quarter</font>") +
      theme(plot.title = element_markdown())

dryspell_hottest <- ggplot() +
      geom_spatraster(data = climate_stack %>% subset(5)) +
      theme_bw() +
      scale_fill_viridis(na.value = "transparent", name = "days", option = "E") +
      ggtitle("<span style='font-size: 10pt;'>Dryspell duration hottest quarter</font>") +
      theme(plot.title = element_markdown())

climate_plot <- (soil_temp | soil_temp_season) / (rain_coldest | dryspell_hottest)

ggsave(plot = climate_plot, filename = "figures/presentation/climate_covars.png",
       dpi = 300, width = 20, height = 15, units = "cm")

# Forest raw ------
names(stack_30m)

cover <- ggplot() +
      geom_spatraster(data = stack_30m %>% subset(2)) +
      theme_bw() +
      ggtitle("<span style='font-size: 10pt;'>Forest cover</font>") +
      theme(plot.title = element_markdown()) +
      scale_fill_gradient(low = "lightgreen", high = "darkgreen", name = "Percentage") 

height <- ggplot() +
      geom_spatraster(data = stack_30m %>% subset(3)) +
      theme_bw() +
      ggtitle("<span style='font-size: 10pt;'>Vegetation height median</font>") +
      theme(plot.title = element_markdown()) +
      scale_fill_gradient(low = "lightgreen", high = "darkgreen", name = "Meter") 

sd_height <- ggplot() +
      geom_spatraster(data = stack_30m %>% subset(4)) +
      theme_bw() +
      ggtitle("<span style='font-size: 10pt;'>Vegetation height stdev</font>") +
      theme(plot.title = element_markdown()) +
      scale_fill_gradient(low = "lightgreen", high = "darkgreen", name = "sd") 

ancient <- ggplot() +
      geom_spatraster(data = stack_30m %>% subset(1)) +
      theme_bw() +
      ggtitle("<span style='font-size: 10pt;'>Distance to ancient woodland</font>") +
      theme(plot.title = element_markdown()) +
      scale_fill_viridis(name = "Meter") 

forest_comb <- (height | sd_height) / (cover | ancient)
ggsave(plot = forest_comb, filename = "figures/presentation/forest_comb.png",
       dpi = 300, width = 25, height = 15, units = "cm")

# Topography ------
aspect <- ggplot() +
      geom_spatraster(data = stack_30m %>% subset(5)) +
      theme_bw() +
      ggtitle("<span style='font-size: 12pt;'>Terrain aspect</font>") +
      theme(plot.title = element_markdown()) +
      scale_fill_gradientn(colours = terrain.colors(100), name = "Degrees")

slope <- ggplot() +
      geom_spatraster(data = stack_30m %>% subset(7)) +
      theme_bw() +
      ggtitle("<span style='font-size: 12pt;'>Slope</font>") +
      theme(plot.title = element_markdown()) +
      scale_fill_gradientn(colours = terrain.colors(100), name = "Percentage")

topo_comb <- (aspect | slope)
ggsave(plot = topo_comb, filename = "figures/presentation/topo_comb.png",
       dpi = 300, width = 25, height = 10, units = "cm")


