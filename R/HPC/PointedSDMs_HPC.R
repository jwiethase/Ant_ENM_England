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

# SET PARAMETERS ---------------------------------------
predictions_resolution <- 2

species_choices <- c("Formica rufa", "Formica lugubris")
max_edges <- c(5, 7, 10)
range_multipliers_sporadic <- c(0.1, 0.5, 0.7)
range_dividers_exhaustive <- c(1, 2, 3)
sigma_list <- c(0.1)
smoother_list = c('tp', 'ps', 'cr')
n_knots_list = c(3, 4)

mult_combs <- crossing(max_edges, range_multipliers_sporadic, range_dividers_exhaustive, smoother_list, n_knots_list, species_choices, sigma_list) # 432
comb_values <- mult_combs[job, ]
print(comb_values)

max.edge = comb_values$max_edges
range_multiplier_sporadic = comb_values$range_multipliers_sporadic
range_divider_exhaustive = comb_values$range_dividers_exhaustive
n_knots = comb_values$n_knots_list
smoother = comb_values$smoother_list
species_choice = comb_values$species_choices
sigma_p = comb_values$sigma_list

prior_sigma = c(sigma_p, 0.01)

dir.create(paste0("model_out/", gsub(" ", "_", species_choice)), showWarnings = F)
dir.create(paste0("model_out/", gsub(" ", "_", species_choice), "/integrated"), showWarnings = F)
dir.create(paste0("figures/", gsub(" ", "_", species_choice)), showWarnings = F)
dir.create(paste0("figures/", gsub(" ", "_", species_choice), "/integrated"), showWarnings = F)

km_proj = '+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=km +no_defs'
covars_selection <- c("clim_topo_PC1_spline1", "clim_topo_PC1_spline2",
                      "clim_topo_PC2_spline1", "clim_topo_PC2_spline2",
                      # "clim_topo_PC3_spline1", "clim_topo_PC3_spline2",
                      # "clim_topo_PC4_spline1", "clim_topo_PC4_spline2",
                      # "clim_topo_PC5_spline1", "clim_topo_PC5_spline2",
                      # "clim_topo_PC6_spline1", "clim_topo_PC6_spline2",
                      "forest_PC2_spline1", "forest_PC2_spline2",
                      "forest_PC3_spline1", "forest_PC3_spline2",
                      "lat_raster",
                      "forest_mask_buff",
                      "days_sampled_sc")

# Data import -----------------------------------
ROI <- vect('data/ROI_kmproj.shp')

sporadic <- read.csv('data/sporadic_combined.csv')  %>% 
      filter(!is.na(days_sampled), 
             species == species_choice) %>% 
      mutate(days_sampled_sc = scale(days_sampled))

exhaustive <- read.csv('data/exhaustive_combined.csv') %>% 
      filter(species == species_choice)

clim_topo_covariates <- rast(paste0("data/6clim_topo_", smoother, "_", n_knots, "k.tif")) 

forest_covariates <- rast(paste0("data/4forest_", smoother, "_", n_knots, "k.tif")) %>% 
      resample(clim_topo_covariates)

effort_rast_10km <- rast('data/effort_rast_10km.tif') %>% 
      disagg(fact = res(.)/res(clim_topo_covariates)) %>% 
      terra::resample(clim_topo_covariates, method = "max", threads = TRUE)
names(effort_rast_10km) <- "days_sampled_sc"

covariates <- c(forest_covariates, clim_topo_covariates, effort_rast_10km) %>% 
      tidyterra::select(all_of(covars_selection))
      
obs_spatial_sporadic <- sporadic %>%
      vect(geom = c("x", "y"), crs = crs(km_proj), keepgeom = TRUE) %>% 
      crop(subst(forest_covariates$forest_mask_buff, 0, NA))

obs_spatial_exhaustive <- exhaustive %>%
      vect(geom = c("x", "y"), crs = crs(km_proj), keepgeom = TRUE) %>% 
      crop(subst(forest_covariates$forest_mask_buff, 0, NA))

sporadic_sf <- st_as_sf(obs_spatial_sporadic)
exhaustive_sf <-  st_as_sf(obs_spatial_exhaustive)

small_regions <- st_as_sf(obs_spatial_exhaustive) %>%
      group_by(source) %>%
      summarize(geometry = st_convex_hull(st_combine(geometry))) %>%
      st_as_sf()

# Mesh -----------------------------------
all_data_spatial <- rbind(sporadic_sf[, c('x', 'y')], exhaustive_sf[, c('x', 'y')])
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
                        Projection = CRS(km_proj),
                        spatialCovariates = covariates,
                        pointsSpatial = 'shared',
                        Coordinates = c('x', 'y'),
                        Boundary = boundary)

spatial_extent_sporadic <- max(sporadic_sf$y) - min(sporadic_sf$y)
range_sporadic <- as.numeric(round(spatial_extent_sporadic*range_multiplier_sporadic))
prior_range_shared = c(range_sporadic, 0.1)

range_exhaustive <- as.numeric(round(range_sporadic/range_divider_exhaustive, digits = 2))
prior_range_separate = c(range_exhaustive, 0.1)

model_setup$specifySpatial(sharedSpatial = TRUE, 
                           prior.range = prior_range_shared,
                           prior.sigma = prior_sigma)

model_setup$addBias('exhaustive_sf')
model_setup$specifySpatial(Bias = 'exhaustive_sf',
                           prior.range = prior_range_separate,
                           prior.sigma = prior_sigma)
 
# model_setup$addSamplers(datasetName = 'sporadic_sf', Samplers = st_as_sf(ROI))
# model_setup$addSamplers(datasetName = 'exhaustive_sf', Samplers = small_regions)

# Run model -----------------------------------
model <- fitISDM(model_setup, options = list(control.inla = list(strategy = "laplace",
                                                                 int.strategy = "auto")))  
if(model$summary.hyperpar[1, 2] < model$summary.hyperpar[1, 1]*0.01 | 
   model$summary.hyperpar[1, 2] > model$summary.hyperpar[1, 1] |
   model$summary.hyperpar[1, 1] < model$summary.hyperpar[2, 1]){
      stop("Random field estimation failed, reconsider priors/mesh specification.")
}

logCPO_vect = log(model$cpo$cpo[model$cpo$cpo != 0])
logCPO_vect = logCPO_vect[is.finite(logCPO_vect)]
logCPO = round(-sum(logCPO_vect, na.rm = T), digits = 2)

# Check output-----------------------------------
summary(model)  
model$summary.hyperpar

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
                 data = proj_grid,
                 mask = boundary,  
                 spatial = TRUE,
                 fun = 'linear', 
                 formula = as.formula(paste0("~ shared_spatial + ", fixed_effects_effort)))

suit_preds <- predict(object = model, 
                     data = proj_grid,
                     mask = boundary,  
                     spatial = TRUE,
                     fun = 'linear', 
                     formula = as.formula(paste0("~ ", fixed_effects_effort)))

shared_random <- predict(object = model, 
                     data = proj_grid,
                     mask = boundary,  
                     spatial = TRUE,
                     fun = 'linear', 
                     formula = ~ shared_spatial)

exhaustive_random <- predict(object = model,
                         data = proj_grid,
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
      scale_fill_viridis(na.value = "transparent", name = "Median", option = "H"); 

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
      scale_fill_viridis(na.value = "transparent", name = "Median", option = "H"); 

shared_plot <- ggplot() +
      geom_spatraster(data = (shared_random_rast$q0.5)) +
      geom_spatvector(data = ROI, fill = NA, lwd = 1) +
      geom_spatvector(data = obs_spatial_sporadic, cex = 0.01, col = 'red', alpha = 0.2) +
      geom_spatvector(data = obs_spatial_exhaustive, cex = 0.01, col = 'purple', alpha = 0.2) +
      theme_minimal() +
      ggtitle("Shared random field", subtitle = "(Response Scale)") +
      scale_fill_viridis(na.value = "transparent", name = "Median"); 

exhaustive_plot <- ggplot() +
      geom_spatraster(data = (exhaustive_random_rast$q0.5)) +
      geom_spatvector(data = ROI, fill = NA, lwd = 1) +
      geom_spatvector(data = obs_spatial_exhaustive, cex = 0.4, col = 'purple', alpha = 0.2) +
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

diagnostic_plot <- ((median_plot | suitability_plot) / (shared_plot | exhaustive_plot)) | (hyper_plot / fixed_plot)
diagnostic_plot

pdf(paste0("figures/", gsub(" ", "_", species_choice), "/integrated/", 
           logCPO, "_all30m_", n_knots, "k_", smoother, "_", 
           sub(" ", "_", species_choice), "_E", max.edge, 
           "_ShMult", range_multiplier_sporadic,"_ShR", range_sporadic, 
           "_SepMult", range_multiplier_exhaustive,"_SepR", range_exhaustive, 
           "_s", sigma_p, ".pdf"),
    width = 12, height = 7)
diagnostic_plot
dev.off()
