# HEADER --------------------------------------------
#
# Author: Joris Wiethase
# Email: j.wiethase@gmail.com
# 
# Script Description:  
# Runs a Point Process model using the 'inlabru' package   

# LOAD PACKAGES -----------------------------------
rm(list = ls())
library(INLA) 
library(inlabru)
library(sf)
library(terra)
library(raster)
library(tidyverse)
library(tidyterra)
library(viridis)
library(patchwork)
library(ggpmisc)
library(sp)
library(data.table)
bru_safe_sp(force = TRUE)
bru_options_set(control.compute = list(cpo = T)) 
source('source/misc_functions.R')

# SET PARAMETERS ---------------------------------------
species_choice <- "Formica rufa"
predictions_resolution <- 2
predictions_resolution_fine <- 0.3
max.edge = 8 # Not much higher than 30

# range_multiplier = 0.1
prior_range = c(150, 0.1)
prior_sigma = c(0.1, 0.5)

smoother = 'tp'
n_knots = 3

# DATA FILES ------------------------------------------
ROI <- vect('spatial_other/ROI_kmproj.shp')
FE_managed <- vect('spatial_other/Forestry_England_managed_forest.shp') %>% 
      terra::project(crs(km_proj))
forest_stack <- rast("covariates/processed/forest_stack_300m.tif")

clim_topo_covariates <- rast(paste0("covariates/processed/6clim_topo_300m_", smoother, "_", n_knots, "k.tif"))
clim_topo_covariates$lat_raster <- terra::scale(clim_topo_covariates$lat_raster)

forest_covariates <- rast(paste0("covariates/processed/2forest_300m_", smoother, "_", n_knots, "k.tif")) 

distance_forest <- rast("covariates/processed/distance_forest_30m.tif") %>%
      terra::resample(forest_covariates) %>% 
      scale()
names(distance_forest) <- "distance_forest"

# forestdummy <- rast("covariates/processed/forest_mask_buff_300m.tif") %>% 
#       mask(ROI)
      
effort_rast_10km <- rast('covariates/processed/effort_rast_lgcp_10km.tif')  %>% 
      terra::scale()

sporadic_vect <- read.csv('species_data/processed_csv/sporadic_combined.csv') %>% 
      filter(species == species_choice) %>% 
      dplyr::select(x, y, days_sampled) %>% 
      vect(geom = c("x", "y"), crs = crs(km_proj), keepgeom = TRUE) 

sporadic_sf <- sporadic_vect %>% 
      thin_spatial(., dist_meters = 0, seed = 123) %>% 
      st_as_sf()

# test <- terra::extract(effort_rast_10km, sporadic_sf)
# hist(test$days_sampl, breaks = 20)

if(species_choice == "Formica lugubris"){
      # A small subset of records for F. lugubris (11) fell into areas of very high sampling effort. 
      # These outliers change the otherwise linear relationship between density and effort. Remove these.
      sporadic_sf <- sporadic_sf %>% 
            filter(days_sampled < 79)
}

# CREATE MESH -------------------------------------------------------
boundary <- st_as_sf(ROI) %>% as("Spatial") # sp format
sporadic_spdf <- sporadic_sf %>% as("Spatial")
mesh <- inla.mesh.2d(boundary = boundary,
                     loc = unique(geom(sporadic_spdf))[, c("x", "y")],
                     max.edge = c(1 ,5) * max.edge,
                     offset = c(1, 2) * max.edge, 
                     cutoff = max.edge/3,
                     crs = fm_crs(km_proj),
                     min.angle = 26)

# FIT MODEL ---------------------------------
# Define spatial SPDE priors
# 
# prior_range = c(range, 0.1)
# prior_sigma = c(1, 0.01)
matern <- inla.spde2.pcmatern(
      mesh,
      prior.range = prior_range,
      prior.sigma = prior_sigma,
      constr = F)

# Construct model formula
# Make the forest dummy variable a factor, where non-forest is "1", and therefore the reference level
base_terms <- c("coordinates ~ Intercept(1)",
                "distance_forest(main = distance_forest, model = 'linear')",
                "lat_raster(main = clim_topo_covariates$lat_raster, model = 'linear')",
                "effort(main = effort_rast_10km$days_sampl, model = 'linear')",
                "mySPDE(main = coordinates, model = matern)")

if(n_knots == 2){n_splines = 2} else {n_splines = n_knots-1}

clim_topo_spline_terms <- make_formula_terms(clim_topo_covariates, "clim_topo_PC", n_splines)
forest_spline_terms <- make_formula_terms(forest_covariates, "forest_PC", n_splines)

form <- as.formula(paste(c(base_terms, forest_spline_terms, clim_topo_spline_terms), collapse = " + "))

# Fit model
model <- lgcp(form,
              data = sporadic_sf %>% as("Spatial"),  
              samplers = boundary,
              domain = list(coordinates = mesh),
              options = list(control.inla = list(strategy = "laplace",
                                                 int.strategy = "auto"),
                             verbose = FALSE))

int.plot <- plot(model, "Intercept")
spde.range <- spde.posterior(model, "mySPDE", what = "range")
spde.logvar <- spde.posterior(model, "mySPDE", what = "log.variance")
range.plot <- plot(spde.range)
var.plot <- plot(spde.logvar)

multiplot(range.plot, var.plot, int.plot)

summary(model)

if(model$summary.hyperpar[1, 1] < model$summary.hyperpar[2, 1] |
   sum(abs(model$summary.fixed$`0.5quant`)) < 1){
      stop("Random field estimation failed, reconsider priors/mesh specification.")
}

logCPO_vect = log(model$cpo$cpo[model$cpo$cpo != 0])
logCPO_vect = logCPO_vect[is.finite(logCPO_vect)]
logCPO = round(-sum(logCPO_vect, na.rm = T), digits = 2)

# PREDICT -----------------------------------------
# Make prediction plots of modeled intensity. In a point process model, 
# intensity describes the point density

# Create grid prediction pixels, use full area, including non-forest patches, to make 
# random field better visible
grid_points <- rast(res = predictions_resolution, ext = ext(ROI), crs = crs(ROI), vals = 1) %>% 
      mask(ROI) %>% 
      as.data.frame(xy = T) %>% 
      SpatialPoints(proj4string = CRS(km_proj))

fixed_effects <- paste(model[["names.fixed"]], collapse = " + ")
fixed_effects_effort <- gsub("effort", "quantile(effort, probs = 0.95, na.rm=T)", fixed_effects)

# Add high constant effort for predictions
all_preds <- predict(object = model, 
                     newdata = grid_points,
                     formula = as.formula(paste0("~ mySPDE + ", fixed_effects_effort)),
                     n.samples = 200,
                     num.threads = 7,
                     seed = 42)

random_field <- predict(object = model, 
                     newdata = grid_points,
                     formula = ~ mySPDE,
                     n.samples = 200,
                     num.threads = 7,
                     seed = 42)

fixed_only <- predict(object = model, 
                       newdata = grid_points,
                       formula = as.formula(paste0("~ ", fixed_effects_effort)),
                      n.samples = 200,
                      num.threads = 7,
                      seed = 42)

all_preds_raster <- as.data.frame(all_preds) %>% 
      dplyr::select(x, y, q0.025, q0.5, q0.975, sd) %>% 
      rast(type = "xyz")

random_field_raster <- as.data.frame(random_field) %>% 
      dplyr::select(x, y, q0.025, q0.5, q0.975, sd) %>% 
      rast(type = "xyz") 

fixed_only_raster <- as.data.frame(fixed_only) %>% 
      dplyr::select(x, y, q0.025, q0.5, q0.975, sd) %>% 
      rast(type = "xyz", crs = crs(km_proj))

median_plot <- ggplot() +
      geom_spatraster(data = (all_preds_raster$q0.5)) +
      geom_spatvector(data = ROI, fill = NA, lwd = 1) +
      theme_minimal() +
      ggtitle("Median estimated point density", subtitle = "(Linear Scale)") +
      scale_fill_viridis(na.value = "transparent", name = "Median", option = 'D'); 

sd_plot <- ggplot() +
      geom_spatraster(data = (all_preds_raster$sd)) +
      geom_spatvector(data = ROI, fill = NA, lwd = 1) +
      theme_minimal() +
      ggtitle("Standard deviation", subtitle = "(Linear Scale)") +
      scale_fill_viridis(na.value = "transparent", name = "SD"); 

field_plot <- ggplot() +
      geom_spatraster(data = (random_field_raster$q0.5)) +
      geom_spatvector(data = ROI, fill = NA, lwd = 1) +
      geom_spatvector(data = vect(sporadic_sf), cex = 0.4, col = 'red', alpha = 0.2) +
      theme_minimal() +
      ggtitle("Spatial random field", subtitle = "(Linear Scale)") +
      scale_fill_viridis(na.value = "transparent", name = "Median")

fixed_only_plot <- ggplot() +
      geom_spatraster(data = (fixed_only_raster$q0.5)) +
      geom_spatvector(data = ROI, fill = NA, lwd = 1) +
      theme_minimal() +
      ggtitle("Fixed effects", subtitle = "(Linear Scale)") +
      scale_fill_viridis(na.value = "transparent", name = "Median", option = 'D'); 

hyper_df <- model$summary.hyperpar %>% 
      round(digits = 3) %>% 
      mutate(var = rownames(.)) %>% 
      dplyr::select(var, mean, sd)

fixed_df <- model$summary.fixed %>% 
      round(digits = 3) %>% 
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

diagnostic_plot <- ((median_plot | sd_plot) / (field_plot | fixed_only_plot)) | (hyper_plot / fixed_plot)
diagnostic_plot

# pdf(paste0("figures/effort_choice/", logCPO, "_fixedEffortGridMaxMaxFixed_4clim_", n_knots, "k_", smoother, "_", sub(" ", "_", species_choice), "_E", max.edge, "_r",
#            prior_range[1], "_", prior_range[2], "_s",
#            prior_sigma[1], "_", prior_sigma[2], ".pdf"),
#     width = 12, height = 7)
# diagnostic_plot
# dev.off()

# Selected areas and model fit check
new_forest <- FE_managed %>% 
      filter(extent == "The Open Forest") %>% 
      fillHoles

cropton <- FE_managed %>% 
      filter(extent == "Cropton") %>% 
      fillHoles

predictions_resolution_fine <- 0.1

grid_points_new_forest <- rast(res = predictions_resolution_fine, ext = ext(new_forest), crs = crs(new_forest), vals = 1) %>% 
      mask(new_forest) %>% 
      as.data.table(xy = T) %>% 
      SpatialPoints(proj4string = CRS(km_proj))

grid_points_cropton <- rast(res = predictions_resolution_fine, ext = ext(cropton), crs = crs(cropton), vals = 1) %>% 
      mask(cropton) %>% 
      as.data.table(xy = T) %>% 
      SpatialPoints(proj4string = CRS(km_proj))

suitability_new_forest <- predict(object = model, 
                                  newdata = grid_points_new_forest,
                                  formula = as.formula(paste0("~ mySPDE + ", fixed_effects_effort)),
                                  n.samples = 100)
suitability_cropton <- predict(object = model, 
                               newdata = grid_points_cropton,
                               formula = as.formula(paste0("~ mySPDE + ", fixed_effects_effort)),
                               n.samples = 100)

suitability_raster_new_forest <- as.data.frame(suitability_new_forest) %>% 
      dplyr::select(x, y, q0.025, q0.5, q0.975, sd) %>% 
      rast(type = "xyz") 

suitability_raster_cropton <- as.data.frame(suitability_cropton) %>% 
      dplyr::select(x, y, q0.025, q0.5, q0.975, sd) %>% 
      rast(type = "xyz") 

suitability_plot_new_forest <- ggplot() +
      geom_spatraster(data = (suitability_raster_new_forest$q0.5)) +
      geom_spatvector(data = ROI, fill = NA, lwd = 1) +
      theme_minimal() +
      ggtitle("Suitability", subtitle = "(Linear Scale)") +
      scale_fill_viridis(na.value = "transparent", name = "Median")

suitability_plot_cropton <- ggplot() +
      geom_spatraster(data = (suitability_raster_cropton$q0.5)) +
      geom_spatvector(data = ROI, fill = NA, lwd = 1) +
      theme_minimal() +
      ggtitle("Suitability", subtitle = "(Linear Scale)") +
      scale_fill_viridis(na.value = "transparent", name = "Median")

# Predict to presence locations
sporadic <- read.csv('species_data/processed_csv/sporadic_combined.csv') %>% 
      filter(species == species_choice,
             source != "dallimore", source != "nym", source != "gaitbarrows", source != "hardcastle") %>% 
      dplyr::select(x, y)
exhaustive <- read.csv('species_data/processed_csv/exhaustive_combined.csv') %>% 
      filter(species == species_choice) %>% 
      dplyr::select(x, y)

combined_presences <- rbind(sporadic, exhaustive) %>% 
      vect(geom = c('x', 'y'), crs = crs(km_proj))

combined_presences_df <- as.data.frame(combined_presences, geom='XY')

preds_env_df <- combined_presences_df %>% 
      cbind(terra::extract(clim_topo_covariates, combined_presences, ID = F)) %>% 
      cbind(terra::extract(forest_covariates, combined_presences, ID = F)) %>% 
      cbind(terra::extract(distance_forest, combined_presences, ID = F)) %>% 
      cbind(terra::extract(effort_rast_10km, combined_presences, ID = F)) %>% 
      drop_na() %>% 
      SpatialPoints(proj4string = CRS(km_proj))

suitable_present <- predict(object = model, 
                            newdata = preds_env_df,
                            formula = as.formula(paste0("~ + ", fixed_effects_effort)),
                            n.samples = 100) %>% 
      as.data.frame()

extracted_hist <- ggplot() +
      geom_rect(aes(xmin = min(values(suitability_raster$q0.5), na.rm = T),
                    xmax = max(values(suitability_raster$q0.5), na.rm = T),
                    ymin = -Inf, 
                    ymax = Inf), fill = "blue", alpha = 0.25) +
      geom_histogram(data = suitable_present, aes(q0.5)) +
      theme_bw() +
      xlim(range(values(suitability_raster$q0.5), na.rm = T))

extracted_box <- ggplot() +
      geom_rect(aes(ymin = min(values(suitability_raster$q0.5), na.rm = T),
                    ymax = max(values(suitability_raster$q0.5), na.rm = T),
                    xmin = -Inf, 
                    xmax = Inf), fill = "blue", alpha = 0.25) +
      geom_boxplot(data = suitable_present, aes(y = q0.5)) +
      ggtitle(paste0("Raw predicted median: ", round(median(suitable_present$q0.5, na.rm = T), digits = 2))) +
      theme_bw() +
      ylim(range(values(suitability_raster$q0.5), na.rm = T))

boyce_test <- ecospat::ecospat.boyce(fit = suitability_raster$q0.5, obs = suitable_present$q0.5)

par(mfrow=c(2, 3))
suitability_plot_new_forest
suitability_plot_cropton

diagnostic_plot_2 <- (suitability_plot_new_forest | suitability_plot_cropton) / (extracted_hist | extracted_box)

pdf(paste0("model_out/Formica_rufa/lgcp/diagnostic_", logCPO, "_all30m_", 
           n_knots, "k_", smoother, "_", sub(" ", "_", species_choice), 
           "_E", max.edge, "_thin", thin_dist,
           "_r", prior_range[1], "_", prior_range[2], 
           "_s", prior_sigma[1], "_", prior_sigma[2], ".pdf"),
    width = 12, height = 9)
diagnostic_plot_2
dev.off()
















# Variable importance -----------------------------------------
forest_buffer_mask <- ifel(forestdummy > 0, 1, NA)

# Use smaller area here, only for forest plus edges. We don't need to predict to non-forest areas.
grid_points_fine <- rast(res = predictions_resolution_fine, ext = ext(forest_buffer_mask), crs = crs(forest_buffer_mask), vals = 1) %>% 
      mask(resample(forest_buffer_mask, .)) %>% 
      as.data.table(xy = T) %>% 
      SpatialPoints(proj4string = CRS(km_proj))

all_preds <- predict(object = model,
                     newdata = grid_points_fine,
                     formula = as.formula(paste0('~ ', fixed_effects_effort)),
                     n.samples = 200,
                     num.threads = 7,
                     seed = 42)

all_preds_df <- as.data.table(all_preds) %>%
      dplyr::select(x, y, q0.025, q0.5, q0.975, sd) 

writeRaster(rast(all_preds_df, type = 'xyz'), "model_out/Formica_lugubris/lgcp/Formica_lugubris_fixedOnly300m_300mRes.tif", overwrite=TRUE)

no_forest <- predict(object = model,
                     newdata = grid_points_fine,
                     formula = ~ lat_raster + forestdummy + distance_ancient + 
                           quantile(effort, probs = 0.95, na.rm=T) +
                           clim_topo_PC1_spline1 + clim_topo_PC1_spline2 +
                           clim_topo_PC2_spline1 + clim_topo_PC2_spline2 +
                           clim_topo_PC3_spline1 + clim_topo_PC3_spline2 +
                           clim_topo_PC4_spline1 + clim_topo_PC4_spline2 +
                           clim_topo_PC5_spline1 + clim_topo_PC5_spline2 +
                           clim_topo_PC6_spline1 + clim_topo_PC6_spline2,
                     n.samples = 200,
                     num.threads = 7,
                     seed = 42)

no_forest_df <- as.data.table(no_forest) %>%
      dplyr::select(x, y, q0.025, q0.5, q0.975, sd) 

no_clim_topo <- predict(object = model,
                        newdata = grid_points_fine,
                        formula = ~ lat_raster + forestdummy + distance_ancient + 
                              quantile(effort, probs = 0.95, na.rm=T) +
                              forest_PC1_spline1 + forest_PC1_spline2 +
                              forest_PC2_spline1 + forest_PC2_spline2,
                        n.samples = 200,
                        num.threads = 7,
                        seed = 42)

no_clim_topo_df <- as.data.table(no_clim_topo) %>%
      dplyr::select(x, y, q0.025, q0.5, q0.975, sd) 

forest_structure_contribution <- 1 - rsq(all_preds_df$q0.5, no_forest_df$q0.5)
clim_topo_contribution <- 1 - rsq(all_preds_df$q0.5, no_clim_topo_df$q0.5)

print(forest_structure_contribution)
print(clim_topo_contribution)

print(paste("Forest variable contribution: ", forest_structure_contribution))
print(paste("Climate and topography variable contribution: ", clim_topo_contribution))












# writeRaster(suitability_raster, "model_out/Formica_lugubris/lgcp/test_suitability.tif", overwrite = T)

# predictions_resolution_small <- 0.1
# FE_managed <- vect('spatial_other/Forestry_England_managed_forest.shp') %>% 
#       terra::project(crs(km_proj))
# cropton <- FE_managed %>% 
#       filter(extent == "Cropton") %>% 
#       fillHoles
# 
# pred_area <- cropton

# PREDICT FINE -----------------------------------------
# Get the suitable area, based on coarser predictions
# sporadic <- read.csv('species_data/processed_csv/sporadic_combined.csv') %>% 
#       filter(species == species_choice,
#              source != "dallimore", source != "nym", source != "gaitbarrows", source != "hardcastle") %>% 
#       dplyr::select(x, y)
# exhaustive <- read.csv('species_data/processed_csv/exhaustive_combined.csv') %>% 
#       filter(species == species_choice) %>% 
#       dplyr::select(x, y)
# 
# combined_presences <- rbind(sporadic, exhaustive) %>% 
#       vect(geom = c('x', 'y'), crs = crs(km_proj))
# 
# suitable_present <- terra::extract(fixed_only_raster, combined_presences)[, "q0.5"]
# suitability_threshold <- quantile(suitable_present, probs = 0.05, na.rm = T)
# hist(suitable_present, main = species_choice)
# abline(v = suitability_threshold, col = "blue")
# 
# suit_lgcp_Eng <- fixed_only_raster$q0.5 %>% 
#       # Remove unsuitable areas
#       clamp(lower = suitability_threshold, values = F) %>% 
#       # Change scale to go from 0 to 1
#       normalise_raster() %>% 
#       mask(ROI)
# 
# suit_mask <- ifel(suit_lgcp_Eng == 0, NA, 1) %>% as.polygons()
# 
# grid_points_final <- rast(res = predictions_resolution_fine, ext = ext(suit_mask), crs = crs(suit_mask), vals = 1) %>% 
#       as.points() %>% 
#       mask(suit_mask) %>% 
#       as.data.table(xy = T) %>% 
#       SpatialPoints(proj4string = CRS(km_proj))
# 
# preds_final <- predict(object = model, 
#                       newdata = grid_points,
#                       formula = as.formula(paste0("~ ", fixed_effects_effort)),
#                       n.samples = 1000,
#                       num.threads = 7,
#                       seed = 42)
# 


