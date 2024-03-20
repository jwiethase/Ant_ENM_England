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
source('source/misc_functions.R')

# Set parameters ---------------------------------------
predictions_resolution <- 5
species_choice <- "Formica rufa"

max.edge = 8

if(species_choice == "Formica rufa"){
      prior_range_shared = c(150, 0.1)
      prior_range_separate = c(150, 0.1)
}

if(species_choice == "Formica lugubris"){
      prior_range_shared = c(200, 0.1)
      prior_range_separate = c(20, 0.1)
      thin_dist = 100
}

prior_sigma = c(0.1, 0.01)
prior_sigma_separate = c(0.1, 0.01)
smoother = 'tp'
n_knots = 3

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
ROI <- vect('spatial_other/ROI_kmproj.shp')

clim_topo_covariates <- rast(paste0("covariates/processed/6clim_topo_300m_", smoother, "_", n_knots, "k.tif")) 
clim_topo_covariates$lat_raster <- scale(clim_topo_covariates$lat_raster)

forest_covariates <- rast(paste0("covariates/processed/2forest_300m_", smoother, "_", n_knots, "k.tif")) %>%
      terra::resample(clim_topo_covariates)

distance_ancient <- rast("covariates/processed/forest_stack_300m.tif") %>% 
      tidyterra::select(distance_ancient) %>% 
      terra::resample(clim_topo_covariates) %>% 
      terra::scale()

sporadic_sf <- read.csv('species_data/processed_csv/sporadic_combined.csv') %>% 
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
      # A small subset of records for F. lugubris (11) fell into areas of very high sampling effort. 
      # These outliers change the otherwise linear relationship between density and effort. Remove these.
      sporadic_sf <- sporadic_sf %>% 
            filter(days_sampled < 79)
}

exhaustive_sf <- read.csv('species_data/processed_csv/exhaustive_combined.csv') %>% 
      filter(species == species_choice) %>% 
      vect(geom = c("x", "y"), crs = crs(km_proj), keepgeom = TRUE) %>% 
      # Some mild spatial thinning to help fitting the model
      thin_spatial(., 50) %>%
      st_as_sf() %>% 
      cbind(terra::extract(forest_covariates$forest_mask_buff, ., ID = FALSE)) %>%
      filter(forest_mask_buff == max(values(forest_covariates$forest_mask_buff, na.rm = T)))

covariates <- c(forest_covariates, clim_topo_covariates, distance_ancient) %>% 
      tidyterra::select(all_of(covars_selection))

# coords <- as.matrix(geom(exhaustive_vect))[, c("x", "y")]
# # Make shapefiles for clusters, this is roughly where exhaustive surveys took place.
# dbscan_result <- dbscan::dbscan(coords, eps = 1, minPts = 5)
# dbscan_result
# exhaustive_vect$cluster_id <- dbscan_result$cluster
# 
# small_regions <- exhaustive_vect %>%
#       st_as_sf() %>% 
#       group_by(cluster_id) %>%
#       summarize(geometry = st_convex_hull(st_combine(geometry))) %>% 
#       vect() %>% 
#       buffer(width = 0.1) %>% 
#       st_as_sf()

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
                           prior.sigma = prior_sigma_separate)

# model_setup$addSamplers(datasetName = 'sporadic_sf', Samplers = boundary)
# model_setup$addSamplers(datasetName = 'exhaustive_sf', Samplers = small_regions)

# Run model -----------------------------------
model <- fitISDM(model_setup, options = list(control.inla = list(strategy = "laplace",
                                                                 int.strategy = "auto")))  

# if(model$summary.hyperpar[1, 2] > model$summary.hyperpar[1, 1] |
#    model$summary.hyperpar[1, 1] < model$summary.hyperpar[2, 1] |
#    sum(abs(model$summary.fixed$`0.5quant`)) < 1){
#       stop("Random field estimation failed, reconsider priors/mesh specification.")
# }

# Check output-----------------------------------
summary(model)

logCPO_vect = log(model$cpo$cpo[model$cpo$cpo != 0])
logCPO_vect = logCPO_vect[is.finite(logCPO_vect)]
logCPO = round(-sum(logCPO_vect, na.rm = T), digits = 2)

# Plot results -----------------------------------
boundary = st_transform(boundary, mesh$crs)
required_nx <- round((max(mesh$loc[,1]) - min(mesh$loc[,1])) / predictions_resolution)
required_ny <- round((max(mesh$loc[,2]) - min(mesh$loc[,2])) / predictions_resolution)

fixed_effects <- paste(model[["names.fixed"]], collapse = " + ")
fixed_effects_effort <- gsub("days_sampl_sc", "max(days_sampl_sc)", fixed_effects)

proj_grid <- fm_pixels(mesh, 
                     dims = c(required_nx, required_ny), 
                     format = 'sf')

all_preds <- predict(object = model,
                 data = proj_grid_effort,
                 mask = boundary,
                 spatial = TRUE,
                 fun = 'linear',
                 formula = as.formula(paste0("~ shared_spatial + ", fixed_effects)))

suit_preds <- predict(object = model,
                     data = proj_grid,
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

median_plot | shared_plot | exhaustive_plot
model$summary.hyperpar

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

# # Effects plot
# original_values <- data.frame(orig_values = unlist(all.seq)) %>%
#       mutate(covariate = sub("*.seq\\d+", "", rownames(.)),
#              sequence = as.numeric(gsub("\\D", "", rownames(.))))
# 
# effect_combs <- data.frame(covariate = gsub('[[:digit:]]+', '', sub("*_lc\\d+", "", rownames(model$summary.lincomb.derived))),
#                            sequence = as.numeric(gsub("\\D", "", rownames(model$summary.lincomb.derived))),
#                            quant_05 = model$summary.lincomb.derived$`0.5quant`,
#                            quant_0025 = model$summary.lincomb.derived$`0.025quant`,
#                            quant_0975 = model$summary.lincomb.derived$`0.975quant`)
# effect_combs$covariate <- gsub('_sortBase', '', effect_combs$covariate)
# 
# effect_combs_m <- merge(original_values, effect_combs) %>%
#       mutate(orig_values = if_else(str_detect(covariate, "log"), exp(orig_values) - 1, orig_values),
#              species = species_choice)
# 
# combined_data <- data.frame(values = c(effect_combs_m$quant_05,
#                                        effect_combs_m$quant_0025,
#                                        effect_combs_m$quant_0975))
# 
# # combined_data$rescaled_value <- range01(combined_data$values)
# # n <- nrow(effect_combs_m)
# # effect_combs_m$median_rescaled <- combined_data$rescaled_value[1:n]
# # effect_combs_m$quant_0025_rescaled <- combined_data$rescaled_value[(n + 1):(2 * n)]
# # effect_combs_m$quant_0975_rescaled <- combined_data$rescaled_value[(2 * n + 1):(3 * n)]
# 
# effects_plot <- ggplot(effect_combs_m) +
#       geom_line(aes(x = orig_values, y = quant_05)) +
#       geom_line(aes(x = orig_values, y = quant_0025), lty = 2, alpha = .5) +
#       geom_line(aes(x = orig_values, y = quant_0975), lty = 2, alpha = .5) +
#       theme_bw() +
#       facet_wrap( ~ covariate, scale = 'free_x') +
#       xlab("Covariate value") +
#       ylab("Intensity")

# diagnostic_plot <- (shared_plot | exhaustive_plot) | (hyper_plot / fixed_plot)
# diagnostic_plot
# # 
# diagnostic_plot <- ((median_plot | suitability_plot) / (shared_plot | exhaustive_plot)) | (hyper_plot / fixed_plot)
# diagnostic_plot

# pdf(paste0("figures/pointed_integrated_output_", sub(" ", "_", species_choice), "_E", max.edge, "_r",
#            prior_range_shared[1], "_", prior_range_shared[2], "_s",
#            prior_sigma_shared[1], "_", prior_sigma_shared[2], ".pdf"),
#     width = 12, height = 7)
# diagnostic_plot
# dev.off()
