rm(list = ls())
library(terra)
library(tidyterra)
source('source/misc_functions.R')

# metrics <- c("sd", "mean", "perc09", "cover")
metrics <- c("cover")
for(i in metrics){
      rast_list <- lapply(list.files('/Users/joriswiethase/Downloads/VOM_processed', pattern = i, full.names = T), rast) %>% 
            sprc()
      rast_mosaic <- terra::mosaic(rast_list, fun = "max")
      writeRaster(rast_mosaic,
                  filename = paste0("covariates/raw/", i, "_VOM_30m_27700.tif"),
                  overwrite=T)
      print(paste0("Finished metric: ", i))
}



