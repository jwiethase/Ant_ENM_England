library(terra)
library(tidyterra)
setwd('/users/jhw538/scratch/ant_modelling')
source('source/misc_functions.R')

args <- commandArgs(trailingOnly=TRUE)
i <- as.integer(args[1])

VOM_tile_list <- list.files("VOM_chunks_smaller", full.names = TRUE)
print(paste0("Starting file: ", basename(VOM_tile_list[i])))
VOM_rast <- trim(rast(VOM_tile_list[i]))
print("VOM trim done")

# Crop the raster to the VOM chunk extent
new_ext <- ext(c(ext(VOM_rast)$xmin-100, ext(VOM_rast)$xmax+100, ext(VOM_rast)$ymin-100, ext(VOM_rast)$ymax+100))

# 300m grid on the raster, to align to
ref_rast_grid_300m <- vect('data/ref_grid_300m_ROI.shp')

# Only keep those 300m grids that are in the VOM chunk extent. Add a 100m buffer around the extent,
# so we don't get gaps later
indices <- terra::relate(ref_rast_grid_300m, new_ext, "intersects")
selected_grids <- ref_rast_grid_300m[indices, ]

intersected_grid_rast_30 <- rast(ext = ext(selected_grids), res = 30, crs = crs(selected_grids), vals = 1)
intersected_grid_VOM <- as.polygons(intersected_grid_rast_30, aggregate = F) 
print("Ref grid VOM done")

# Now aggregate to the grid polygons
#mean_result <- terra::extract(VOM_rast, intersected_grid_VOM, fun=mean_fun, bind = T)
#names(mean_result) <- c('lyr', 'mean_VOM')
#mean_result_rast <- rasterize(mean_result, intersected_grid_rast_30, 
#                              field = "mean_VOM")

#writeRaster(mean_result_rast, 
#            paste0("VOM_processed/", 
#                   gsub("VOM_Mosaic_Phase", "mean", basename(VOM_tile_list[i]))), 
#            overwrite = TRUE)
#print("mean done")

#sd_result <- terra::extract(VOM_rast, intersected_grid_VOM, fun=sd_fun, bind = T)
#names(sd_result) <- c("lyr", "sd_VOM")
#sd_result_rast <- rasterize(sd_result, intersected_grid_rast_30, 
#                              field = "sd_VOM")

#writeRaster(sd_result_rast, 
#            paste0("VOM_processed/", 
#                   gsub("VOM_Mosaic_Phase", "sd", basename(VOM_tile_list[i]))), 
#            overwrite = TRUE)
#print("sd done")

#perc09_result <- terra::extract(VOM_rast, intersected_grid_VOM, fun=calc_quantiles, bind = T)
#names(perc09_result) <- c("lyr", "perc09_VOM")
#perc09_result_rast <- rasterize(perc09_result, intersected_grid_rast_30, 
#                            field = "perc09_VOM")

#writeRaster(perc09_result_rast, 
#            paste0("VOM_processed/", 
#                   gsub("VOM_Mosaic_Phase", "perc09", basename(VOM_tile_list[i]))), 
#            overwrite = TRUE)
#print("perc09 done")

cover_result <- terra::extract(VOM_rast, intersected_grid_VOM, fun=calculate_cover, bind = T)
names(cover_result) <- c("lyr", "cover_VOM")
cover_result_rast <- rasterize(cover_result, intersected_grid_rast_30, 
                            field = "cover_VOM")
 
writeRaster(cover_result_rast, 
            paste0("VOM_processed/", 
                   gsub("VOM_Mosaic_Phase", "cover", basename(VOM_tile_list[i]))), 
           overwrite = TRUE)
print("cover done")
     

     
