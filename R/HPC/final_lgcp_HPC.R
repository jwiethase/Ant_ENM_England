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
bru_options_set(control.compute = list(cpo = F)) 
setwd('/users/jhw538/scratch/ant_modelling')
source('source/misc_functions.R')

predictions_resolution <- 1
species_choice <- "Formica rufa"

# Get parameters of best fitting model ---------------------------------------
best_model_name <- list.files(paste0('figures/', species_choice, '/lgcp'))[1]

max.edge = sub(".*E([0-9]+\\.[0-9]+)_.*", "\\1", best_model_name)
range_multiplier = sub(".*mult([0-9]+\\.[0-9]+)_.*", "\\1", best_model_name)
n_knots = sub(".*_([0-9]+)k_.*", "\\1", best_model_name)
smoother = sub(".*k_([a-z]{2}).*", "\\1", best_model_name)
species_choice = paste("Formica", sub(".*_Formica_([^_]+)_.*", "\\1", best_model_name))
prior_sigma = c(as.numeric(sub(".*_s([0-9]+\\.[0-9]+)_([0-9]+\\.[0-9]+)_.*", "\\1", best_model_name)),
                as.numeric(sub(".*_s([0-9]+\\.[0-9]+)_([0-9]+\\.[0-9]+)_.*", "\\2", best_model_name)))

dir.create(paste0("model_out/", gsub(" ", "_", species_choice)), showWarnings = F)
dir.create(paste0("model_out/", gsub(" ", "_", species_choice), "/lgcp"), showWarnings = F)
dir.create(paste0("figures/", gsub(" ", "_", species_choice)), showWarnings = F)
dir.create(paste0("figures/", gsub(" ", "_", species_choice), "/lgcp"), showWarnings = F)

# DATA FILES ------------------------------------------
ROI <- vect('data/ROI_kmproj.shp')
sporadic <- read.csv('data/sporadic_combined.csv')

clim_topo_covariates <- rast(paste0("data/6clim_topo_", smoother, "_", n_knots, "k.tif"))
forest_covariates <- rast(paste0("data/2forest_", smoother, "_", n_knots, "k.tif"))

# clim_topo_covariates <- rast(paste0("data/6clim_topo_300m_", smoother, "_", n_knots, "k.tif"))
# forest_covariates <- rast(paste0("data/4forest_300m_", smoother, "_", n_knots, "k.tif"))

effort_rast_10km <- rast('data/effort_rast_10km.tif') 

# FILTER AND STANDARDISE DATA FRAMES --------------------------
sporadic_filtered <- sporadic %>% 
      filter(species == species_choice) %>% 
      dplyr::select(x, y)

obs_spatial <- sporadic_filtered %>% 
      vect(geom = c("x", "y"), crs = crs(km_proj), keepgeom = TRUE) %>% 
      crop(subst(forest_covariates$forest_mask_buff, 0, NA))

obs_spdf <- st_as_sf(obs_spatial) %>% 
      as("Spatial")

# CREATE MESH -------------------------------------------------------
boundary <- st_as_sf(ROI) %>% as("Spatial") 
mesh <- inla.mesh.2d(boundary = boundary,
                     loc = unique(geom(obs_spatial))[, c("x", "y")],
                     max.edge = c(1 ,5) * max.edge,
                     offset = c(1, 2) * max.edge, 
                     cutoff = max.edge/3,
                     crs = fm_CRS(km_proj),
                     min.angle = 26)

# FIT MODEL ---------------------------------
# Define spatial SPDE priors
spatial_extent <- max(sporadic_filtered$y) - min(sporadic_filtered$y)
range <- as.numeric(round(spatial_extent*range_multiplier))
prior_range = c(range, 0.1)

matern <- inla.spde2.pcmatern(
      mesh,
      prior.range = prior_range,
      prior.sigma = prior_sigma,
      constr = F)

# Construct model formula
base_terms <- c("coordinates ~ Intercept(1)",
                "forestdummy(main = forest_covariates$forest_mask_buff, model = 'linear')",
                "lat_raster(main = clim_topo_covariates$lat_raster, model = 'linear')",
                "effort_rast_10km(main = effort_rast_10km, model = 'linear')",
                "mySPDE(main = coordinates, model = matern)")

if(n_knots == 2){n_splines = 2} else {n_splines = n_knots-1}

clim_topo_spline_terms <- make_formula_terms(clim_topo_covariates, "clim_topo_PC", n_splines)
forest_spline_terms <- make_formula_terms(forest_covariates, "forest_PC", n_splines)

form <- as.formula(paste(c(base_terms, clim_topo_spline_terms, forest_spline_terms), collapse = " + "))

# Fit model
model <- lgcp(form,
              data = obs_spdf,  
              samplers = boundary,
              domain = list(coordinates = mesh),
              options = list(control.inla = list(strategy = "laplace",
                                                 int.strategy = "auto"),
                             verbose = FALSE))

if(model$summary.hyperpar[1, 2] < model$summary.hyperpar[1, 1]*0.01 | 
   model$summary.hyperpar[1, 2] > model$summary.hyperpar[1, 1] |
   model$summary.hyperpar[1, 1] < model$summary.hyperpar[2, 1]){
      stop("Random field estimation failed, reconsider priors/mesh specification.")
}

# PREDICT -----------------------------------------
# Create grid prediction pixels
required_nx <- round((max(mesh$loc[,1]) - min(mesh$loc[,1])) / predictions_resolution)
required_ny <- round((max(mesh$loc[,2]) - min(mesh$loc[,2])) / predictions_resolution)

maskborder <- st_as_sf(ROI)
maskborder = st_transform(maskborder, mesh$crs)

proj_grid <- fm_pixels(mesh, 
                       dims = c(required_nx, required_ny), 
                       format = 'sp')

fixed_effects <- paste(model[["names.fixed"]], collapse = " + ")

# Add high constant effort for predictions
fixed_effects_effort <- gsub("effort_rast_10km", "max(effort_rast_10km, na.rm = T)", fixed_effects)
all_preds <- predict(object = model, 
                     newdata = proj_grid,
                     mask = maskborder,
                     formula = as.formula(paste0("~ mySPDE + ", fixed_effects_effort)))

random_field <- predict(object = model, 
                        newdata = proj_grid,
                        mask = maskborder,
                        formula = ~ Intercept + mySPDE)

suitability <- predict(object = model, 
                       newdata = proj_grid,
                       mask = proj_grid_effort,
                       formula = as.formula(paste0("~ ", fixed_effects_effort)))

all_preds_raster <- as.data.frame(all_preds) %>% 
      rename(x = coords.x1, y = coords.x2) %>% 
      dplyr::select(x, y, q0.025, q0.5, q0.975, sd) %>% 
      rast(type = "xyz") %>% 
      mask(ROI)

random_field_raster <- as.data.frame(random_field) %>% 
      rename(x = coords.x1, y = coords.x2) %>% 
      dplyr::select(x, y, q0.025, q0.5, q0.975, sd) %>% 
      rast(type = "xyz") %>% 
      mask(ROI)

suitability_raster <- as.data.frame(suitability) %>% 
      rename(x = coords.x1, y = coords.x2) %>% 
      dplyr::select(x, y, q0.025, q0.5, q0.975, sd) %>% 
      rast(type = "xyz") %>% 
      mask(ROI)

median_plot <- ggplot() +
      geom_spatraster(data = (all_preds_raster$q0.5)) +
      geom_spatvector(data = ROI, fill = NA, lwd = 1) +
      theme_minimal() +
      ggtitle("Median estimate", subtitle = "(Response Scale)") +
      scale_fill_viridis(na.value = "transparent", name = "Median")

sd_plot <- ggplot() +
      geom_spatraster(data = (all_preds_raster$sd)) +
      geom_spatvector(data = ROI, fill = NA, lwd = 1) +
      theme_minimal() +
      ggtitle("Standard deviation", subtitle = "(Response Scale)") +
      scale_fill_viridis(na.value = "transparent", name = "SD")

field_plot <- ggplot() +
      geom_spatraster(data = (random_field_raster$q0.5)) +
      geom_spatvector(data = ROI, fill = NA, lwd = 1) +
      geom_spatvector(data = obs_spatial, cex = 0.4, col = 'red', alpha = 0.2) +
      theme_minimal() +
      ggtitle("Spatial random field", subtitle = "(Response Scale)") +
      scale_fill_viridis(na.value = "transparent", name = "Median")

suitability_plot <- ggplot() +
      geom_spatraster(data = (suitability_raster$q0.5)) +
      geom_spatvector(data = ROI, fill = NA, lwd = 1) +
      theme_minimal() +
      ggtitle("Suitability", subtitle = "(Response Scale)") +
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

fixed_plot <- ggplot() + 
      annotate(geom='table',
               x=4,
               y=0,
               label=list(fixed_df)) +
      theme_void()

diagnostic_plot <- ((median_plot | sd_plot) / (field_plot | suitability_plot)) | (hyper_plot / fixed_plot)

save(model, all_preds, suitability, random_field, comb_values, mesh, fixed_effects_effort, proj_grid, 
     file = paste0("model_out/", gsub(" ", "_", species_choice),"/lgcp/final_all30m_", n_knots, "k_", smoother, "_", 
                   sub(" ", "_", species_choice), "_E", max.edge, "_mult", range_multiplier,
                   "_r", prior_range[1], "_", prior_range[2], 
                   "_s", prior_sigma[1], "_", prior_sigma[2], ".RData"))

pdf(paste0("figures/", gsub(" ", "_", species_choice), "/lgcp/final_all30m_", n_knots, "k_", smoother, "_", 
           sub(" ", "_", species_choice), "_E", max.edge, "_mult", range_multiplier,"_r",
           prior_range[1], "_", prior_range[2], "_s",
           prior_sigma[1], "_", prior_sigma[2], "_", ".pdf"),
    width = 12, height = 7)
diagnostic_plot
dev.off()




