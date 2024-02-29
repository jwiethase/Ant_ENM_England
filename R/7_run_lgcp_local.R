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
bru_safe_sp(force = TRUE)
bru_options_set(control.compute = list(cpo = T)) 
source('source/misc_functions.R')

# SET PARAMETERS ---------------------------------------
species_choice <- "Formica rufa"
predictions_resolution <- 5
max.edge = 10 # Not much higher than 30

range_multiplier = 0.2
prior_sigma = c(0.1, 0.01)

smoother = 'cr'
n_knots = 3

# DATA FILES ------------------------------------------
ROI <- vect('spatial_other/ROI_outline_27700.shp') %>% 
      terra::project(crs(km_proj))

sporadic <- read.csv('species_data/processed_csv/sporadic_combined.csv')

clim_topo_covariates <- rast(paste0("covariates/processed/6clim_topo_300m_", smoother, "_", n_knots, "k.tif"))
# %>% 
#       terra::subset(!stringr::str_detect(names(.), "PC3"))
forest_covariates <- rast(paste0("covariates/processed/4forest_300m_", smoother, "_", n_knots, "k.tif")) %>%
      terra::subset(!stringr::str_detect(names(.), "PC1|PC4"))
effort_rast_10km <- rast('covariates/processed/effort_rast_10km.tif') 

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
boundary <- st_as_sf(ROI) %>% as("Spatial")  # sp format
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
beepr::beep()
logCPO_vect = log(model$cpo$cpo[model$cpo$cpo != 0])
logCPO_vect = logCPO_vect[is.finite(logCPO_vect)]
logCPO = round(-sum(logCPO_vect, na.rm = T), digits = 2)

int.plot <- plot(model, "Intercept")
spde.range <- spde.posterior(model, "mySPDE", what = "range")
spde.logvar <- spde.posterior(model, "mySPDE", what = "log.variance")
range.plot <- plot(spde.range)
var.plot <- plot(spde.logvar)

multiplot(range.plot, var.plot, int.plot)

# PREDICT -----------------------------------------
# Make prediction plots of modeled intensity. In a point process model, 
# intensity describes the point density
# Create grid prediction pixels
required_nx <- round((max(mesh$loc[,1]) - min(mesh$loc[,1])) / predictions_resolution)
required_ny <- round((max(mesh$loc[,2]) - min(mesh$loc[,2])) / predictions_resolution)

maskborder <- st_as_sf(ROI)
maskborder = st_transform(maskborder, mesh$crs)

proj_grid <- fm_pixels(mesh, 
                       dims = c(required_nx, required_ny), 
                       format = 'sp')

# proj_grid_effort <- fm_cprod(proj_grid,
#                              data.frame(effort_rast_10km = quantile(values(effort_rast_10km), probs = 0.95, na.rm=T)[[1]]))

fixed_effects <- paste(model[["names.fixed"]], collapse = " + ")

# Add high constant effort for predictions
# values(effort_rast_10km) <- quantile(values(effort_rast_10km), probs = 0.95, na.rm = TRUE)
fixed_effects_effort <- gsub("effort_rast_10km", "max(effort_rast_10km, na.rm=T)", fixed_effects)
all_preds <- predict(object = model, 
                     newdata = proj_grid,
                     mask = maskborder,
                     formula = as.formula(paste0("~ (mySPDE + ", fixed_effects_effort, ")")))

random_field <- predict(object = model, 
                     newdata = proj_grid,
                     mask = maskborder,
                     formula = ~ Intercept + mySPDE)

suitability <- predict(object = model, 
                       newdata = proj_grid,
                       mask = maskborder,
                       formula = as.formula(paste0("~ (", fixed_effects_effort, ")")))

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
      ggtitle("Median estimated point density", subtitle = "(Linear Scale)") +
      scale_fill_viridis(na.value = "transparent", name = "log(Median)", option = 'D'); median_plot

sd_plot <- ggplot() +
      geom_spatraster(data = (all_preds_raster$sd)) +
      geom_spatvector(data = ROI, fill = NA, lwd = 1) +
      theme_minimal() +
      ggtitle("Standard deviation", subtitle = "(Linear Scale)") +
      scale_fill_viridis(na.value = "transparent", name = "SD"); sd_plot

field_plot <- ggplot() +
      geom_spatraster(data = (random_field_raster$q0.5)) +
      geom_spatvector(data = ROI, fill = NA, lwd = 1) +
      geom_spatvector(data = obs_spatial, cex = 0.4, col = 'red', alpha = 0.2) +
      theme_minimal() +
      ggtitle("Spatial random field", subtitle = "(Linear Scale)") +
      scale_fill_viridis(na.value = "transparent", name = "Median")

suitability_plot <- ggplot() +
      geom_spatraster(data = (suitability_raster$q0.5)) +
      geom_spatvector(data = ROI, fill = NA, lwd = 1) +
      theme_minimal() +
      ggtitle("Suitability", subtitle = "(Response Scale)") +
      scale_fill_viridis(na.value = "transparent", name = "Median", option = 'D'); suitability_plot

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
diagnostic_plot
# pdf(paste0("figures/effort_choice/", logCPO, "_fixedEffortGridMaxMaxFixed_4clim_", n_knots, "k_", smoother, "_", sub(" ", "_", species_choice), "_E", max.edge, "_r",
#            prior_range[1], "_", prior_range[2], "_s",
#            prior_sigma[1], "_", prior_sigma[2], ".pdf"),
#     width = 12, height = 7)
# diagnostic_plot
# dev.off()




