# install.packages('INLA',repos=c(getOption('repos'),INLA='https://inla.r-inla-download.org/R/testing'), dep=TRUE)
# remotes::install_github('inlabru-org/inlabru', ref = 'stable')
# devtools::install_github('PhilipMostert/PointedSDMs')


#Define Hazard model
obs.f = function(effort, lsig){
      distance <- 1 / effort
      1-exp(-(distance/(exp(lsig)))^-1)}
# Change effort 0 values to 0.0000001

### a basic model with recoding effort, using point samplers with sqrt(1/pi) radius offset:
cmp.eff  <-  ~ mySmooth(map = coordinates, model = matern) +
      beta1.clim(map = f.clim.1(x,y) , model = "linear") +
      beta1.sq.clim(map = ( f.clim.1(x,y) )^2 , model = "linear") +
      beta2.clim(map = f.clim.2(x,y) , model = "linear") +
      beta2.sq.clim(map = ( f.clim.2(x,y) )^2 , model = "linear") +
      beta3.clim(map = f.clim.3(x,y) , model = "linear") +
      beta3.sq.clim(map = ( f.clim.3(x,y) )^2 , model = "linear") +
      beta4.clim(map = f.clim.4(x,y) , model = "linear") +
      beta4.sq.clim(map = ( f.clim.4(x,y) )^2 , model = "linear") +
      df.lsigma + Intercept

for.eff <- coordinates + eff.scaled ~
      mySmooth + Intercept +
      beta1.clim + beta1.sq.clim +
      beta2.clim + beta2.sq.clim +
      beta3.clim + beta3.sq.clim +
      beta4.clim + beta4.sq.clim +
      log(obs.f(eff.scaled, df.lsigma))

fit <- lgcp(components = cmp.eff,
            formula = for.eff,
            data = sp.spdf,
            samplers = GB.buffer.points,
            options = list( max.iter = 2,
                            offset = log( 1 / sp.spdf$eff.scaled))) ### !!! ###

# Effects plot
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

# combined_data$rescaled_value <- range01(combined_data$values)
# n <- nrow(effect_combs_m)
# effect_combs_m$median_rescaled <- combined_data$rescaled_value[1:n]
# effect_combs_m$quant_0025_rescaled <- combined_data$rescaled_value[(n + 1):(2 * n)]
# effect_combs_m$quant_0975_rescaled <- combined_data$rescaled_value[(2 * n + 1):(3 * n)]

# effects_plot <- ggplot(effect_combs_m) +
#       geom_line(aes(x = orig_values, y = quant_05)) +
#       geom_line(aes(x = orig_values, y = quant_0025), lty = 2, alpha = .5) +
#       geom_line(aes(x = orig_values, y = quant_0975), lty = 2, alpha = .5) +
#       theme_bw() +
#       facet_wrap( ~ covariate, scale = 'free_x')

r1 <- r2 <- r3 <- rast(temp_stack_qtr, nlyr=1)

for (row in 1:nrow(temp_stack_qtr)) {
      for (col in 1:ncol(temp_stack_qtr)) {
            values <- sapply(1:nlyr(temp_stack_qtr), function(i) temp_stack_qtr[[i]][row, col])
            
            if(all(is.na(values))) {
                  r1[row, col] <- NA
                  r2[row, col] <- NA
                  r3[row, col] <- NA
            } else {
                  sorted_values <- base::sort(t(data.frame(values))[, 1], decreasing = TRUE, na.last = TRUE)
                  sorted_df <- data.frame(sorted_values) %>% 
                        mutate(qtr = as.integer(gsub("X", "", row.names(.))))
                  
                  r1[row, col] <- sorted_df$qtr[1]
                  r2[row, col] <- sorted_df$qtr[2]
                  r3[row, col] <- sorted_df$qtr[3]
            }
      }
}

hottest_qtr_stack <- c(r1, r2, r3)
names(new_stack) <- c("Highest", "Second Highest", "Third Highest")





library(inlabru)
library(terra)
library(tidyterra)
library(sf)
library(sp)
library(data.table)
library(landscapemetrics)
setwd('/users/jhw538/scratch/ant_modelling')
source('source/misc_functions.R')
rsq <- function(x, y) summary(lm(y~x))$r.squared

args <- commandArgs(trailingOnly=TRUE)
job <- as.integer(args[1])

# SET PARAMETERS ---------------------------------------
predictions_resolution <- 30 # In meters

point_buffer = 50
area_percentile_threshold = 0.01

species_choices = c('Formica rufa', 'Formica lugubris')
species_choice = 'Formica lugubris'
#species_choice = species_choices[job]

print(paste("Starting", species_choice))

smoother = 'tp'
n_knots = 3

# DATA FILES ------------------------------------------
if(species_choice == 'Formica lugubris'){
      load('model_out/Formica_lugubris/lgcp/30082.7_all30m_3k_tp_Formica_lugubris_E9_thin100_r150_0.1_s1_0.01.RData')
}

if(species_choice == 'Formica rufa'){
      load('model_out/Formica_rufa/lgcp/4199081.02_all30m_3k_tp_Formica_rufa_E7_r50_0.1_s0.1_0.01.RData')
}

ROI <- vect('data/ROI_kmproj.shp')
forest_stack <- rast("data/forest_stack_30m.tif")

clim_topo_covariates <- rast(paste0("data/6clim_topo_", smoother, "_", n_knots, "k.tif")) %>% 
      tidyterra::select(contains(c("spline", "lat_raster")))
clim_topo_covariates$lat_raster <- terra::scale(clim_topo_covariates$lat_raster)

forest_covariates <- rast(paste0("data/2forest_", smoother, "_", n_knots, "k.tif"))

forest_stack$distance_ancient <- terra::scale(forest_stack$distance_ancient)

distance_ancient <- rast("data/forest_stack_30m.tif") %>% 
      tidyterra::select(distance_ancient) %>% 
      terra::scale()

forest_covariates$forest_mask_buff <- ifel(forest_covariates$forest_mask_buff == 0, 1, 2)

effort_rast_10km <- rast('data/effort_rast_lgcp_10km.tif')  %>% 
      terra::scale()

print("Covariate prep done")

sporadic <- read.csv('data/sporadic_combined.csv') %>% 
      filter(source != 'dallimore', source != 'nym', source != 'gaitbarrows', source != 'hardcastle') %>% 
      dplyr::select(x, y, species)
exhaustive <- read.csv('data/exhaustive_combined.csv') %>% 
      dplyr::select(x, y, species)

combined_presences <- rbind(sporadic, exhaustive) %>% 
      vect(geom = c('x', 'y'), crs = crs(km_proj)) %>% 
      terra::project(crs('epsg:27700')) 

ant_vect <- combined_presences %>% 
      filter(species == species_choice)

ant_vect_buff <- terra::buffer(ant_vect, width = point_buffer)

# PREDICT -----------------------------------------
forest_mask <- ifel(cover_VOM < 0.3, NA, 1)

# Get forest area
forest_patch_area <- spatialize_lsm(forest_mask, 
                                    level = 'patch',
                                    metric = 'area')[[1]][[1]]

# What patch area do the ants prefer?
preferred_area <- terra::extract(forest_patch_area, ant_vect_buff) %>% 
      drop_na()

lower_threshold <- quantile(preferred_area$value, probs = area_percentile_threshold)

# Create forest patch layer above area threshold
forest_patches_sub <- ifel(forest_patch_area <= lower_threshold, NA, 1) %>% 
      terra::project(crs(km_proj)) 

print("forest_patches_sub done")

grid_points <- rast(res = predictions_resolution, ext = ext(forest_patches_sub), crs = crs(forest_patches_sub), vals = 1) %>% 
      resample(forest_patches_sub) %>% 
      mask(forest_patches_sub) 

presence_raster <- rast(res = res(grid_points), ext=ext(forest_patches_sub), crs=crs(forest_patches_sub))
presence_raster <- rasterize(combined_presences, presence_raster, vals=1, background=NA)

combined_points <- cover(grid_points, presence_raster) %>% 
      terra::project(crs(km_proj)) %>%
      as.data.table(xy = T) %>%
      dplyr::select(x, y) %>%
      SpatialPoints(proj4string = CRS(km_proj))

all_preds <- predict(object = model, 
                     newdata = combined_points,
                     formula = as.formula(paste0("~ mySPDE + ", fixed_effects_effort)),
                     n.samples = 1000,
                     num.threads = 60,
                     seed = 42)

all_preds_raster <- as.data.frame(all_preds) %>% 
      dplyr::select(x, y, q0.025, q0.5, q0.975, sd) %>% 
      rast(type = "xyz", crs = crs(km_proj))

writeRaster(all_preds_raster, 
            filename = paste0("model_out/", gsub(" ", "_", species_choice),"/lgcp/final_allPreds_all", predictions_resolution*1000, "m_", n_knots, "k_", smoother, "_", 
                              sub(" ", "_", species_choice), "_E", max.edge,  "_thin", thin_dist,
                              "_r", prior_range[1], "_", prior_range[2], 
                              "_s", prior_sigma[1], "_", prior_sigma[2], ".tif"), overwrite=TRUE)

print("all_preds_raster done")

fixed_only <- predict(object = model, 
                      newdata = combined_points,
                      formula = as.formula(paste0("~ ", fixed_effects_effort)),
                      n.samples = 1000,
                      num.threads = 60,
                      seed = 42)

fixed_only_raster <- as.data.frame(fixed_only) %>% 
      dplyr::select(x, y, q0.025, q0.5, q0.975, sd) %>% 
      rast(type = "xyz", crs = crs(km_proj))

writeRaster(fixed_only_raster, 
            filename = paste0("model_out/", gsub(" ", "_", species_choice),"/lgcp/final_fixedOnly_all", predictions_resolution*1000, "m_", n_knots, "k_", smoother, "_", 
                              sub(" ", "_", species_choice), "_E", max.edge,  "_thin", thin_dist,
                              "_r", prior_range[1], "_", prior_range[2], 
                              "_s", prior_sigma[1], "_", prior_sigma[2], ".tif"), overwrite=TRUE)

print("fixed_only_raster done")

# VARIABLE IMPORTANCE -----------------------------------------
# Don't need quite as high resolution for this
grid_points <- rast(res = 300, ext = ext(ROI), crs = crs(km_proj), vals = 1) %>% 
      mask(ROI) %>% 
      as.data.frame(xy = T) %>% 
      SpatialPoints(proj4string = CRS(km_proj))

no_forest <- predict(object = model, 
                     newdata = grid_points,
                     formula = ~ mySPDE + lat_raster + forestdummy + quantile(effort, probs = 0.95, na.rm=T) + 
                           clim_topo_PC1_spline1 + clim_topo_PC1_spline2 + 
                           clim_topo_PC2_spline1 + clim_topo_PC2_spline2 + 
                           clim_topo_PC3_spline1 + clim_topo_PC3_spline2 + 
                           clim_topo_PC4_spline1 + clim_topo_PC4_spline2 + 
                           clim_topo_PC5_spline1 + clim_topo_PC5_spline2 + 
                           clim_topo_PC6_spline1 + clim_topo_PC6_spline2,
                     n.samples = 1000,
                     num.threads = 60,
                     seed = 42)

no_clim_topo <- predict(object = model, 
                        newdata = grid_points,
                        formula = ~ mySPDE + lat_raster + forestdummy + quantile(effort, probs = 0.95, na.rm=T) + 
                              forest_PC2_spline1 + forest_PC2_spline2 + 
                              forest_PC3_spline1 + forest_PC3_spline2,
                        n.samples = 1000,
                        num.threads = 60,
                        seed = 42)

forest_contribution <- 1 - rsq(all_preds$q0.5, no_forest$q0.5)
clim_topo_contribution <- 1 - rsq(all_preds$q0.5, no_clim_topo$q0.5)

print(forest_contribution)
print(clim_topo_contribution)

print(paste("Forest variable contribution: ", forest_contribution))
print(paste("Climate and topography variable contribution: ", clim_topo_contribution))



# PREDICT -----------------------------------------
# Create grid prediction pixels
point_buffer = 50
area_percentile_threshold = 0.01

sporadic <- read.csv('species_data/processed_csv/sporadic_combined.csv') %>%
      filter(species == species_choice,
             source != "dallimore", source != "nym", source != "gaitbarrows", source != "hardcastle") %>%
      dplyr::select(x, y)
exhaustive <- read.csv('species_data/processed_csv/exhaustive_combined.csv') %>%
      filter(species == species_choice) %>%
      dplyr::select(x, y)

combined_presences <- rbind(sporadic, exhaustive) %>%
      vect(geom = c('x', 'y'), crs = crs(km_proj)) %>%
      terra::project(crs('epsg:27700'))

ant_vect_buff <- terra::buffer(combined_presences, width = point_buffer)

forest_stack <- rast('covariates/processed/forest_stack_300m.tif') %>%
      terra::project(crs('epsg:27700'))

forest_mask <- forest_stack %>%
      tidyterra::select(cover_VOM)

forest_mask <- ifel(forest_mask < 0.3, NA, 1)

# Get forest area
forest_patch_area <- spatialize_lsm(forest_mask,
                                    level = 'patch',
                                    metric = 'area')[[1]][[1]]

preferred_area <- terra::extract(forest_patch_area, ant_vect_buff) %>%
      drop_na()

lower_threshold <- quantile(preferred_area$value, probs = area_percentile_threshold)

# Create forest patch layer above area threshold
forest_patches_sub <- ifel(forest_patch_area <= lower_threshold, NA, 1)

# Make raster of regular pixels across forest patches
grid_points <- rast(res = 300, ext = ext(forest_patches_sub), crs = crs(forest_patches_sub), vals = 1) %>%
      resample(forest_patches_sub) %>%
      mask(forest_patches_sub)

# Also make a raster of ant nest records, we need to have predictions for all of these locations
# later to adjust the outputs
presence_raster <- rast(res = res(grid_points), ext=ext(forest_patches_sub), crs=crs(forest_patches_sub))
presence_raster <- rasterize(combined_presences, presence_raster, vals=1, background=NA)

# Now combine both and make spatial points for prediction
combined_points <- cover(grid_points, presence_raster) %>%
      terra::project(crs(km_proj)) %>%
      as.data.table(xy = T) %>%
      dplyr::select(x, y) %>%
      SpatialPoints(proj4string = CRS(km_proj))

# writeRaster(grid_points, "grid_points.tif")
all_preds <- predict(object = model,
                     newdata = combined_points,
                     formula = as.formula(paste0("~ mySPDE + ", fixed_effects_effort)),
                     n.samples = 100,
                     seed = 42)

all_preds_raster <- as.data.frame(all_preds) %>%
      dplyr::select(x, y, q0.025, q0.5, q0.975, sd) %>%
      rast(type = "xyz", crs = crs(km_proj))
plot(all_preds_raster$q0.5)
writeRaster(all_preds_raster$q0.5, "all_preds_raster.tif", overwrite = T)





