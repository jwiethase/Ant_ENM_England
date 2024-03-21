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
bru_options_set(control.compute = list(cpo = T)) 
setwd('/users/jhw538/scratch/ant_modelling')
source('source/misc_functions.R')

args <- commandArgs(trailingOnly=TRUE)
job <- as.integer(args[1])

# SET PARAMETERS ---------------------------------------
predictions_resolution <- 5
species_choice = "Formica lugubris"
print(species_choice)

max_edges <- c(9, 11)
ranges <- c(50, 150, 200)
range_p_list <- c(0.1)
smoother_list = c('tp')
n_knots_list = c(3)
sigma_list = c(0.1, 1)
sigma_p_list = c(0.01, 0.5)
thin_dist_list = c(0, 100)

mult_combs <- crossing(max_edges, ranges, range_p_list, smoother_list, n_knots_list,  
                       sigma_list, sigma_p_list, thin_dist_list) # 48
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

clim_topo_covariates <- rast(paste0("data/6clim_topo_", smoother, "_", n_knots, "k.tif"))
clim_topo_covariates$lat_raster <- terra::scale(clim_topo_covariates$lat_raster)

forest_covariates <- rast(paste0("data/2forest_30m_cr_3k.tif")) 

forestdummy <- rast("data/forest_mask_buff_30m.tif") %>% 
      terra::scale() %>% 
      mask(ROI)
names(forestdummy) <- "forest_mask_buff"

effort_rast_10km <- rast('data/effort_rast_lgcp_10km.tif')  %>% 
      terra::scale()

# FILTER AND STANDARDISE DATA FRAMES --------------------------
sporadic_sf <- read.csv('data/sporadic_combined.csv') %>% 
      filter(species == species_choice) %>% 
      dplyr::select(x, y, days_sampled) %>% 
      vect(geom = c("x", "y"), crs = crs(km_proj), keepgeom = TRUE) %>% 
      thin_spatial(., dist_meters = thin_dist, seed = 123) %>% 
      st_as_sf() %>% 
      cbind(terra::extract(forestdummy, ., ID = FALSE)) %>%
      filter(forest_mask_buff == max(values(forestdummy, na.rm = T)))

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
                "forestdummy(main = forestdummy$forest_mask_buff, model = 'linear')",
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
     file = paste0("model_out/Formica_lugubris/lgcp/", logCPO, "_all30m_", n_knots, "k_", smoother, "_", 
                   sub(" ", "_", species_choice), "_E", max.edge, "_thin", thin_dist,
                   "_r", prior_range[1], "_", prior_range[2], 
                   "_s", prior_sigma[1], "_", prior_sigma[2], ".RData"))

print("Model saved")

# Diagnostic plots
suitability <- predict(object = model, 
                       newdata = grid_points,
                       formula = as.formula(paste0("~ mySPDE + ", fixed_effects_effort)),
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

pdf(paste0("model_out/Formica_lugubris/lgcp/", logCPO, "_all30m_", 
           n_knots, "k_", smoother, "_", sub(" ", "_", species_choice), 
           "_E", max.edge, "_thin", thin_dist,
           "_r", prior_range[1], "_", prior_range[2], 
           "_s", prior_sigma[1], "_", prior_sigma[2], ".pdf"),
    width = 12, height = 9)
diagnostic_plot
dev.off()