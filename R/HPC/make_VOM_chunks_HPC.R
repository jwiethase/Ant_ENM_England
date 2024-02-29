library(terra)
library(tidyterra)
library(parallel)
setwd('/users/jhw538/scratch/ant_modelling')
source('source/misc_functions.R')

args <- commandArgs(trailingOnly=TRUE)
i <- as.integer(args[1])

VOM_tile_list <- list.files("VOM_tiles", full.names = TRUE)

# Function to subdivide and process a raster
make_raster_chunks <- function(raster_path, chunks = 10) {
      raster_path <- VOM_tile_list[i]
      VOM_tile <- rast(raster_path)
      total_extent <- ext((VOM_tile))
      chunk_size <- (total_extent$ymax - total_extent$ymin) / chunks
      
      for (chunk_index in 1:chunks) {
            print(paste0("Processing chunk ", chunk_index))
            ymin_chunk <- total_extent$ymin + (chunk_index - 1) * chunk_size
            ymax_chunk <- ymin_chunk + chunk_size
            chunk_extent <- c(total_extent$xmin, total_extent$xmax, ymin_chunk, ymax_chunk)
            chunk <- VOM_tile %>% terra::crop(chunk_extent)
            writeRaster(chunk, paste0("VOM_chunks_smaller/", gsub("VOM_Mosaic_Phase", paste0("VOM_Mosaic_Phase_chunk", chunk_index), basename(raster_path))), overwrite = TRUE)
      }
}

print(paste0("Starting ", basename(VOM_tile_list[i])))
make_raster_chunks(VOM_tile_list[i])
print(paste0("Finished ", basename(VOM_tile_list[i])))
