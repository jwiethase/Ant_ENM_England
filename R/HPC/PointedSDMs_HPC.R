rm(list = ls())
library(sf)
library(terra)
library(INLA)
library(PointedSDMs)
library(tidyverse)
library(viridis)
library(tidyterra)
library(ggpmisc)
library(patchwork)
bru_options_set(control.compute = list(cpo = T)) 
setwd('/users/jhw538/scratch/ant_modelling')
source('source/misc_functions.R')

args <- commandArgs(trailingOnly=TRUE)
job <- as.integer(args[1])

# SET PARAMETERS ---------------------------------------
predictions_resolution <- 5

species_choices <- c("Formica rufa", "Formica lugubris")
max_edges <- c(7, 9)
smoother_list = c('tp')
n_knots_list = c(3)
shared_ranges = c(100, 150, 200, 250)
separate_ranges = c(25, 50, 100)
probs = c(0.1, 0.5)

mult_combs <- crossing(max_edges, smoother_list, n_knots_list, species_choices, shared_ranges, separate_ranges, probs) # 96
comb_values <- mult_combs[job, ]
print(comb_values)

max.edge = comb_values$max_edges
n_knots = comb_values$n_knots_list
smoother = comb_values$smoother_list
species_choice = comb_values$species_choices
shared_range = comb_values$shared_ranges
separat_range = comb_values$separate_ranges
prior_prob = comb_values$probs 

prior_sigma = c(0.1, 0.01)
prior_range_shared = c(shared_range, prior_prob)
prior_range_separate = c(separat_range, prior_prob)

if(species_choice == "Formica rufa"){
      thin_dist = 100
}

if(species_choice == "Formica lugubris"){
      thin_dist = 50
}

dir.create(paste0("model_out/", gsub(" ", "_", species_choice)), showWarnings = F)
dir.create(paste0("model_out/", gsub(" ", "_", species_choice), "/integrated"), showWarnings = F)
dir.create(paste0("figures/", gsub(" ", "_", species_choice)), showWarnings = F)
dir.create(paste0("figures/", gsub(" ", "_", species_choice), "/integrated"), showWarnings = F)

km_proj = '+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=km +no_defs'

covars_selection <- c("clim_topo_PC1_spline1", "clim_topo_PC1_spline2",
                      "clim_topo_PC2_spline1", "clim_topo_PC2_spline2",
                      "clim_topo_PC3_spline1", "clim_topo_PC3_spline2",
                      "clim_topo_PC4_spline1", "clim_topo_PC4_spline2",
                      "clim_topo_PC5_spline1", "clim_topo_PC5_spline2",
                      "clim_topo_PC6_spline1", "clim_topo_PC6_spline2",
                      "forest_PC1_spline1", "forest_PC1_spline2",
                      "forest_PC2_spline1", "forest_PC2_spline2",
                      "lat_raster",
                      "distance_ancient",
                      "forest_mask_buff"
)

# Data files ------------------------------------------
ROI <- vect('data/ROI_kmproj.shp')

clim_topo_covariates <- rast(paste0("data/6clim_topo_300m_", smoother, "_", n_knots, "k.tif"))
clim_topo_covariates$lat_raster <- scale(clim_topo_covariates$lat_raster)

forest_covariates <- rast(paste0("data/2forest_300m_", smoother, "_", n_knots, "k.tif")) %>% 
      resample(clim_topo_covariates)

distance_ancient <- rast("data/forest_stack_300m.tif") %>% 
      tidyterra::select(distance_ancient) %>% 
      terra::resample(clim_topo_covariates) %>% 
      terra::scale()

sporadic_sf <- read.csv('data/sporadic_combined.csv') %>% 
      filter(species == species_choice,
             source != "dallimore", source != "nym", source != "gaitbarrows", source != "hardcastle") %>% 
      dplyr::select(x, y, days_sampled) %>% 
      mutate(days_sampled_log = log(days_sampled + 1)) %>% 
      vect(geom = c("x", "y"), crs = crs(km_proj), keepgeom = TRUE) %>% 
      thin_spatial(., dist_meters = 100) %>% 
      st_as_sf() %>% 
      cbind(terra::extract(forest_covariates$forest_mask_buff, ., ID = FALSE)) %>%
      filter(forest_mask_buff == max(values(forest_covariates$forest_mask_buff, na.rm = T)))

if(species_choice == "Formica lugubris"){
      sporadic_sf <- sporadic_sf %>% 
            filter(days_sampled < 79)
}

exhaustive_sf <- read.csv('data/exhaustive_combined.csv') %>% 
      filter(species == species_choice) %>% 
      vect(geom = c("x", "y"), crs = crs(km_proj), keepgeom = TRUE) %>% 
      thin_spatial(., 50) %>%
      st_as_sf() %>% 
      cbind(terra::extract(forest_covariates$forest_mask_buff, ., ID = FALSE)) %>%
      filter(forest_mask_buff == max(values(forest_covariates$forest_mask_buff, na.rm = T)))

covariates <- c(forest_covariates, clim_topo_covariates, distance_ancient) %>% 
      tidyterra::select(all_of(covars_selection))

# Mesh -----------------------------------
all_data_spatial <- rbind(st_coordinates(sporadic_sf), 
                          st_coordinates(exhaustive_sf))
boundary <- st_as_sf(ROI)
mesh <- inla.mesh.2d(boundary = boundary,
                     loc = all_data_spatial,
                     max.edge = c(1, 5) * max.edge,
                     offset = c(1, 2) * max.edge, 
                     cutoff = max.edge/3,
                     crs = fm_CRS(km_proj),
                     min.angle = 26)

# Model setup -----------------------------------
model_setup <- intModel(sporadic_sf, exhaustive_sf,
                        Mesh = mesh, 
                        Projection = crs(km_proj),
                        spatialCovariates = covariates,
                        pointCovariates = 'days_sampled_log',
                        # responseCounts = 'pseudo_count',
                        pointsSpatial = 'shared',
                        Coordinates = c('x', 'y'),
                        Boundary = boundary)

model_setup$specifySpatial(sharedSpatial = TRUE, 
                           prior.range = prior_range_shared,
                           prior.sigma = prior_sigma)

model_setup$addBias('exhaustive_sf')
model_setup$specifySpatial(Bias = 'exhaustive_sf',
                           prior.range = prior_range_separate,
                           prior.sigma = prior_sigma)

# Run model -----------------------------------
print("Starting model...")
model <- fitISDM(model_setup, options = list(control.inla = list(strategy = "laplace",
                                                                 int.strategy = "auto")))  

if(model$summary.hyperpar[1, 1] < model$summary.hyperpar[2, 1] |
      sum(abs(model$summary.fixed$`0.5quant`)) < 1){
      stop("Random field estimation failed, reconsider priors/mesh specification.")
}

# Check output-----------------------------------
summary(model)
model$summary.hyperpar

logCPO_vect = log(model$cpo$cpo[model$cpo$cpo != 0])
logCPO_vect = logCPO_vect[is.finite(logCPO_vect)]
logCPO = round(-sum(logCPO_vect, na.rm = T), digits = 2)

boundary = st_transform(boundary, mesh$crs)
required_nx <- round((max(mesh$loc[,1]) - min(mesh$loc[,1])) / predictions_resolution)
required_ny <- round((max(mesh$loc[,2]) - min(mesh$loc[,2])) / predictions_resolution)

fixed_effects <- paste(model[["names.fixed"]], collapse = " + ")
# fixed_effects_effort <- gsub("days_sampled_log", "quantile(days_sampled_log, probs = 0.95, na.rm=T)", fixed_effects)

proj_grid <- fm_pixels(mesh, 
                       dims = c(required_nx, required_ny), 
                       format = 'sf')
proj_grid_effort <- fm_cprod(proj_grid, data.frame(days_sampled_log = quantile(sporadic_sf$days_sampled_log, probs = 0.95)[[1]]))

all_preds <- predict(object = model,
                     data = proj_grid_effort,
                     mask = boundary,
                     spatial = TRUE,
                     fun = 'linear',
                     formula = as.formula(paste0("~ shared_spatial + ", fixed_effects)))

suit_preds <- predict(object = model,
                      data = proj_grid_effort,
                      mask = boundary,
                      spatial = TRUE,
                      fun = 'linear',
                      formula = as.formula(paste0("~ ", fixed_effects)))

shared_random <- predict(object = model, 
                         data = proj_grid_effort,
                         mask = boundary,  
                         spatial = TRUE,
                         fun = 'linear', 
                         formula = ~ shared_spatial)

exhaustive_random <- predict(object = model,
                             data = proj_grid_effort,
                             mask = boundary,
                             spatial = TRUE,
                             fun = 'linear',
                             formula = ~ exhaustive_sf_biasField)

all_preds_df <- as.data.frame(all_preds$predictions)
suit_preds_df <- as.data.frame(suit_preds$predictions)
shared_random_df <- as.data.frame(shared_random$predictions)
exhaustive_random_df <- as.data.frame(exhaustive_random$predictions)

all_preds_rast <- do.call(rbind, st_geometry(all_preds_df$geometry)) %>%
      as_tibble() %>% setNames(c("x","y")) %>%
      cbind(all_preds_df) %>%
      dplyr::select(x, y, q0.025, q0.5, q0.975, sd) %>%
      rast(type = "xyz") %>%
      mask(ROI)

suit_preds_rast <- do.call(rbind, st_geometry(suit_preds_df$geometry)) %>%
      as_tibble() %>% setNames(c("x","y")) %>%
      cbind(all_preds_df) %>%
      dplyr::select(x, y, q0.025, q0.5, q0.975, sd) %>%
      rast(type = "xyz") %>%
      mask(ROI)

shared_random_rast <- do.call(rbind, st_geometry(shared_random_df$geometry)) %>% 
      as_tibble() %>% setNames(c("x","y")) %>% 
      cbind(shared_random_df) %>% 
      dplyr::select(x, y, q0.025, q0.5, q0.975, sd) %>% 
      rast(type = "xyz") %>% 
      mask(ROI)

exhaustive_random_rast <- do.call(rbind, st_geometry(exhaustive_random_df$geometry)) %>%
      as_tibble() %>% setNames(c("x","y")) %>%
      cbind(exhaustive_random_df) %>%
      dplyr::select(x, y, q0.025, q0.5, q0.975, sd) %>%
      rast(type = "xyz") %>%
      mask(ROI)

median_plot <- ggplot() +
      geom_spatraster(data = (all_preds_rast$q0.5)) +
      geom_spatvector(data = ROI, fill = NA, lwd = 1) +
      # geom_spatvector(data = obs_spatial_sporadic, cex = 0.4, col = 'red', alpha = 0.2) +
      # geom_spatvector(data = obs_spatial_exhaustive, cex = 0.4, col = 'purple', alpha = 0.2) +
      theme_minimal() +
      ggtitle("Median estimate", subtitle = "(Response Scale)") +
      scale_fill_viridis(na.value = "transparent", name = "Median", option = "D");

sd_plot <- ggplot() +
      geom_spatraster(data = (all_preds_rast$sd)) +
      geom_spatvector(data = ROI, fill = NA, lwd = 1) +
      # geom_spatvector(data = obs_spatial_sporadic, cex = 0.4, col = 'red', alpha = 0.2) +
      # geom_spatvector(data = obs_spatial_exhaustive, cex = 0.4, col = 'purple', alpha = 0.2) +
      theme_minimal() +
      ggtitle("Standard deviation", subtitle = "(Response Scale)") +
      scale_fill_viridis(na.value = "transparent", name = "SD");

suitability_plot <- ggplot() +
      geom_spatraster(data = (suit_preds_rast$q0.5)) +
      geom_spatvector(data = ROI, fill = NA, lwd = 1) +
      # geom_spatvector(data = obs_spatial_sporadic, cex = 0.4, col = 'red', alpha = 0.2) +
      # geom_spatvector(data = obs_spatial_exhaustive, cex = 0.4, col = 'purple', alpha = 0.2) +
      theme_minimal() +
      ggtitle("Suitability estimate", subtitle = "(Response Scale)") +
      scale_fill_viridis(na.value = "transparent", name = "Median", option = "D");

shared_plot <- ggplot() +
      geom_spatraster(data = (shared_random_rast$q0.5)) +
      geom_spatvector(data = ROI, fill = NA, lwd = 1) +
      geom_spatvector(data = sporadic_vect, 
                      cex = 0.01, col = 'red', alpha = 0.2) +
      geom_spatvector(data = exhaustive_vect, 
                      cex = 0.01, col = 'purple', alpha = 0.2) +
      theme_minimal() +
      ggtitle("Shared random field", subtitle = "(Response Scale)") +
      scale_fill_viridis(na.value = "transparent", name = "Median"); 

exhaustive_plot <- ggplot() +
      geom_spatraster(data = (exhaustive_random_rast$q0.5)) +
      geom_spatvector(data = ROI, fill = NA, lwd = 1) +
      geom_spatvector(data = exhaustive_vect, 
                      cex = 0.01, col = 'purple', alpha = 0.2) +
      theme_minimal() +
      ggtitle("Separate field", subtitle = "(Response Scale)") +
      scale_fill_viridis(na.value = "transparent", name = "Median"); 

hyper_df <- model$summary.hyperpar %>%
      round(digits = 2) %>%
      mutate(var = rownames(.)) %>%
      dplyr::select(var, everything()) %>%
      dplyr::select(var, mean, sd)

fixed_df <- model$summary.fixed %>%
      round(digits = 2) %>%
      mutate(var = rownames(.)) %>%
      dplyr::select(var, everything()) %>%
      dplyr::select(-mode, -kld, -mean)

hyper_plot <- ggplot() +
      annotate(geom='table',
               x=4,
               y=0,
               label=list(hyper_df)) +
      theme_void()

fixed_plot <- ggplot() +
      annotate(geom='table',
               x=4,
               y=0,
               label=list(fixed_df)) +
      theme_void()


diagnostic_plot <- ((median_plot | suitability_plot) / (shared_plot | exhaustive_plot)) | (hyper_plot / fixed_plot)

pdf(paste0("figures/", gsub(" ", "_", species_choice), "/integrated/", 
           logCPO, "_all300m_", n_knots, "k_", smoother, "_", 
           sub(" ", "_", species_choice), "_E", max.edge, "_thin", thin_dist,
          "_ShR", prior_range_shared[1], 
          "_SepR", prior_range_separate[1],
          ".pdf"),
    width = 12, height = 7)
diagnostic_plot
dev.off()
