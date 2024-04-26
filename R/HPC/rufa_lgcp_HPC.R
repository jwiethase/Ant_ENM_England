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
bru_options_set(control.compute = list(cpo = T)) 
setwd('/users/jhw538/scratch/ant_modelling')
source('source/misc_functions.R')

args <- commandArgs(trailingOnly=TRUE)
job <- as.integer(args[1])

# SET PARAMETERS ---------------------------------------
predictions_resolution <- 2
species_choice = "Formica rufa"
print(species_choice)

max_edges <- c(5, 7, 9, 11)
ranges <- c(100, 200, 250)
range_p_list <- c(0.1, 0.5)
smoother_list = c('tp')
n_knots_list = c(3)
sigma_list = c(0.1)
sigma_p_list = c(0.01, 0.5)
thin_dist_list = c(0, 100)

mult_combs <- crossing(max_edges, ranges, range_p_list, smoother_list, n_knots_list,  
                       sigma_list, sigma_p_list, thin_dist_list) # 24
comb_values <- mult_combs[job, ]
print(comb_values)

max.edge = comb_values$max_edges
n_knots = comb_values$n_knots_list
smoother = comb_values$smoother_list
range_p = comb_values$range_p_list
range = comb_values$ranges
sigma = comb_values$sigma_list
sigma_p = comb_values$sigma_p_list
thin_dist = comb_values$thin_dist_list

prior_sigma = c(sigma, sigma_p)
prior_range = c(range, range_p)

# DATA FILES ------------------------------------------
ROI <- vect('data/ROI_kmproj.shp')
forest_stack <- rast("data/forest_stack_30m.tif")
distance_forest <- rast("data/distance_forest_30m.tif")
names(distance_forest) <- "distance_forest"
FE_managed <- vect('data/Forestry_England_managed_forest.shp') %>% 
      terra::project(crs(km_proj))

clim_topo_covariates <- rast(paste0("data/6clim_topo_", smoother, "_", n_knots, "k.tif"))
clim_topo_covariates$lat_raster <- terra::scale(clim_topo_covariates$lat_raster)

forest_covariates <- rast(paste0("data/2forest_30m_", smoother, "_3k.tif")) 

# forestdummy <- rast("data/forest_mask_buff_30m.tif") %>% 
#       mask(ROI)
# names(forestdummy) <- "forest_mask_buff"

effort_rast_10km <- rast('data/effort_rast_lgcp_10km.tif')  %>% 
      terra::scale()

# FILTER AND STANDARDISE DATA FRAMES --------------------------
sporadic_vect <- read.csv('data/sporadic_combined.csv') %>% 
      filter(species == species_choice) %>% 
      dplyr::select(x, y, days_sampled) %>% 
      vect(geom = c("x", "y"), crs = crs(km_proj), keepgeom = TRUE) 

sporadic_sf <- sporadic_vect %>% 
      thin_spatial(., dist_meters = thin_dist, seed = 123) %>% 
      st_as_sf()# %>% 
      # cbind(terra::extract(forestdummy, ., ID = FALSE)) %>%
      # filter(forest_mask_buff == max(values(forestdummy, na.rm = T)))

# CREATE MESH -------------------------------------------------------
boundary <- st_as_sf(ROI) %>% as("Spatial") 
sporadic_spdf <- sporadic_sf %>% as("Spatial")
mesh <- inla.mesh.2d(boundary = boundary,
                     loc = unique(geom(sporadic_spdf))[, c("x", "y")],
                     max.edge = c(1 ,5) * max.edge,
                     offset = c(1, 2) * max.edge, 
                     cutoff = max.edge/3,
                     crs = fm_CRS(km_proj),
                     min.angle = 26)
print("Mesh done")

# FIT MODEL ---------------------------------
matern <- inla.spde2.pcmatern(
      mesh,
      prior.range = prior_range,
      prior.sigma = prior_sigma,
      constr = F)

# Construct model formula
base_terms <- c("coordinates ~ Intercept(1)",
                "distance_forest(main = distance_forest, model = 'linear')",
                "lat_raster(main = clim_topo_covariates$lat_raster, model = 'linear')",
                "effort(main = effort_rast_10km$days_sampl, model = 'linear')",
                "mySPDE(main = coordinates, model = matern)")

if(n_knots == 2){n_splines = 2} else {n_splines = n_knots-1}

clim_topo_spline_terms <- make_formula_terms(clim_topo_covariates, "clim_topo_PC", n_splines)
forest_spline_terms <- make_formula_terms(forest_covariates, "forest_PC", n_splines)

form <- as.formula(paste(c(base_terms, clim_topo_spline_terms, forest_spline_terms), collapse = " + "))
print("Setup done")

# Fit model
model <- lgcp(form,
              data = sporadic_spdf,  
              samplers = boundary,
              domain = list(coordinates = mesh),
              options = list(control.inla = list(strategy = "laplace",
                                                 int.strategy = "auto"),
                             verbose = FALSE))
summary(model)
if(sum(abs(model$summary.fixed$`0.5quant`)) < 1 | 
   model$summary.hyperpar[1, 1] < model$summary.hyperpar[2, 1]){
      stop("Random field estimation failed, reconsider priors/mesh specification.")
}

logCPO_vect = log(model$cpo$cpo[model$cpo$cpo != 0])
logCPO_vect = logCPO_vect[is.finite(logCPO_vect)]
logCPO = round(-sum(logCPO_vect, na.rm = T), digits = 2)

# PREDICT -----------------------------------------
# Create grid prediction pixels
grid_points <- rast(res = predictions_resolution, ext = ext(ROI), crs = crs(km_proj), vals = 1) %>% 
      mask(ROI) %>% 
      as.data.frame(xy = T) %>% 
      SpatialPoints(proj4string = CRS(km_proj))

fixed_effects <- paste(model[["names.fixed"]], collapse = " + ")
fixed_effects_effort <- gsub("effort", "quantile(effort, probs = 0.95, na.rm=T)", fixed_effects)

save(model, grid_points, mesh, fixed_effects_effort,
     file = paste0("model_out/", sub(" ", "_", species_choice), "/lgcp/", logCPO, "_all30m_", n_knots, "k_", smoother, "_", 
                   sub(" ", "_", species_choice), "_E", max.edge, "_thin", thin_dist,
                   "_r", prior_range[1], "_", prior_range[2], 
                   "_s", prior_sigma[1], "_", prior_sigma[2], ".RData"))

print("Model saved")

# Diagnostic plots
suitability <- predict(object = model, 
                       newdata = grid_points,
                       formula = as.formula(paste0("~ + ", fixed_effects_effort)),
                       n.samples = 100)

random_field <- predict(object = model, 
                        newdata = grid_points,
                        formula = ~ mySPDE,
                        n.samples = 100)

suitability_raster <- as.data.frame(suitability) %>% 
      dplyr::select(x, y, q0.025, q0.5, q0.975, sd) %>% 
      rast(type = "xyz") 

random_field_raster <- as.data.frame(random_field) %>% 
      dplyr::select(x, y, q0.025, q0.5, q0.975, sd) %>% 
      rast(type = "xyz") 

suitability_plot <- ggplot() +
      geom_spatraster(data = (suitability_raster$q0.5)) +
      geom_spatvector(data = ROI, fill = NA, lwd = 1) +
      theme_minimal() +
      ggtitle("Suitability", subtitle = "(Linear Scale)") +
      scale_fill_viridis(na.value = "transparent", name = "Median")

sd_plot <- ggplot() +
      geom_spatraster(data = (suitability_raster$sd)) +
      geom_spatvector(data = ROI, fill = NA, lwd = 1) +
      theme_minimal() +
      ggtitle("Standard deviation", subtitle = "(Linear Scale)") +
      scale_fill_viridis(na.value = "transparent", name = "SD")

field_plot <- ggplot() +
      geom_spatraster(data = (random_field_raster$q0.5)) +
      geom_spatvector(data = ROI, fill = NA, lwd = 1) +
      geom_spatvector(data = vect(sporadic_spdf), cex = 0.4, col = 'red', alpha = 0.2) +
      theme_minimal() +
      ggtitle("Spatial random field", subtitle = "(Linear Scale)") +
      scale_fill_viridis(na.value = "transparent", name = "Median")

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

fixed_est_plot <- ggplot() + 
      annotate(geom='table',
               x=4,
               y=0,
               label=list(fixed_df)) +
      theme_void()

diagnostic_plot <- ((suitability_plot | sd_plot) / (field_plot | hyper_plot)) | fixed_est_plot

pdf(paste0("figures/", sub(" ", "_", species_choice), "/lgcp/", logCPO, "_all30m_", 
           n_knots, "k_", smoother, "_", sub(" ", "_", species_choice), 
           "_E", max.edge, "_thin", thin_dist,
           "_r", prior_range[1], "_", prior_range[2], 
           "_s", prior_sigma[1], "_", prior_sigma[2], ".pdf"),
    width = 12, height = 9)
diagnostic_plot
dev.off()

# Selected areas and model fit check
# Predict to presence locations
sporadic <- read.csv('data/sporadic_combined.csv') %>% 
      filter(species == species_choice,
             source != "dallimore", source != "nym", source != "gaitbarrows", source != "hardcastle") %>% 
      dplyr::select(x, y)
exhaustive <- read.csv('data/exhaustive_combined.csv') %>% 
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

new_forest <- FE_managed %>% 
      filter(extent == "The Open Forest") %>% 
      fillHoles

cropton <- FE_managed %>% 
      filter(extent == "Cropton") %>% 
      fillHoles

ennerdale <- FE_managed %>% 
      filter(extent == "Ennerdale") %>% 
      fillHoles

predictions_resolution_fine <- 0.03

grid_points_new_forest <- rast(res = predictions_resolution_fine, ext = ext(new_forest), crs = crs(new_forest), vals = 1) %>% 
      mask(new_forest) %>% 
      as.data.table(xy = T) %>% 
      SpatialPoints(proj4string = CRS(km_proj))

grid_points_cropton <- rast(res = predictions_resolution_fine, ext = ext(cropton), crs = crs(cropton), vals = 1) %>% 
      mask(cropton) %>% 
      as.data.table(xy = T) %>% 
      SpatialPoints(proj4string = CRS(km_proj))

grid_points_ennerdale <- rast(res = predictions_resolution_fine, ext = ext(ennerdale), crs = crs(ennerdale), vals = 1) %>% 
      mask(ennerdale) %>% 
      as.data.table(xy = T) %>% 
      SpatialPoints(proj4string = CRS(km_proj))

suitability_new_forest <- predict(object = model, 
                                  newdata = grid_points_new_forest,
                                  formula = as.formula(paste0("~ + ", fixed_effects_effort)),
                                  n.samples = 100)
suitability_cropton <- predict(object = model, 
                               newdata = grid_points_cropton,
                               formula = as.formula(paste0("~ + ", fixed_effects_effort)),
                               n.samples = 100)
suitability_ennerdale <- predict(object = model, 
                                 newdata = grid_points_ennerdale,
                                 formula = as.formula(paste0("~ + ", fixed_effects_effort)),
                                 n.samples = 100)

suitability_raster_new_forest <- as.data.frame(suitability_new_forest) %>% 
      dplyr::select(x, y, q0.025, q0.5, q0.975, sd) %>% 
      rast(type = "xyz") 

suitability_raster_cropton <- as.data.frame(suitability_cropton) %>% 
      dplyr::select(x, y, q0.025, q0.5, q0.975, sd) %>% 
      rast(type = "xyz") 

suitability_raster_ennerdale <- as.data.frame(suitability_ennerdale) %>% 
      dplyr::select(x, y, q0.025, q0.5, q0.975, sd) %>% 
      rast(type = "xyz") 

pdf(paste0("figures/", sub(" ", "_", species_choice), "/lgcp/diagnostic_", logCPO, "_all30m_", 
           n_knots, "k_", smoother, "_", sub(" ", "_", species_choice), 
           "_E", max.edge, "_thin", thin_dist,
           "_r", prior_range[1], "_", prior_range[2], 
           "_s", prior_sigma[1], "_", prior_sigma[2], ".pdf"),
    width = 17, height = 10)
par(mfrow=c(2, 3))
plot(suitability_raster_new_forest$q0.5)
plot(suitability_raster_cropton$q0.5)
plot(suitability_raster_ennerdale$q0.5)
hist(suitable_present$q0.5, xlim = range(values(suitability_raster$q0.5), na.rm = T))
boxplot(suitable_present$q0.5, ylim = range(values(suitability_raster$q0.5), na.rm = T))
boyce_test <- ecospat::ecospat.boyce(fit = suitability_raster$q0.5, obs = suitable_present$q0.5, method = 'pearson')
title(paste0("Boyce test cor: ", boyce_test$cor))
par(mfrow=c(1, 1))
dev.off()

