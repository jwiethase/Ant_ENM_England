# HEADER --------------------------------------------
#
# Author: Joris Wiethase
# Email: j.wiethase@gmail.com
# 
# Script Description:  
# Apply further processing to the pre-processed climate layers
# - Rainfall and temperature summary of the coldest and hottest quarters of the years 2010-2013
# - Longest dry spell duration

rm(list = ls())
library(tidyverse)
library(tidyterra)
library(terra)
library(lubridate)
source('source/misc_functions.R')

# Climate downloaded from: https://data.ceda.ac.uk/badc/ukmo-hadobs/data/insitu/MOHC/HadOBS/HadUK-Grid/v1.2.0.ceda/1km

# Import daily climate and set parameters ----------------------------------------------------
temp_files <- list.files(path = "covariates/raw/clim_daily/", pattern = "tas_hadukgrid_uk_1km_mon_.*\\.nc\\.tif$", full.names = TRUE)
rain_files_all <- list.files(path = "covariates/raw/clim_daily", pattern = "rainfall_hadukgrid_uk_1km_day_.*\\.nc\\.tif$", full.names = TRUE)
temp_stack <- rast(temp_files)
names(temp_stack) <- sprintf("%d%02d", rep(2010:2013, each=12), 1:12)
start_year <- as.integer(min(regmatches(rain_files_all, regexpr("\\d{4}", rain_files_all))))
rainfall_threshold <- 1 # Rain threshold in mm, used to calculate dry spell duration
n_quarters = 4  # The number of individual coldest and hottest quarters to look for, 4 since we have data for 4 years

# Identify the hottest and coldest quarters ----------------------------------------------------
quarters <- ceiling(1:nlyr(temp_stack) / 3)
temp_stack_qtr <- tapp(temp_stack, quarters, fun = median)
avg_qtr_temp <- global(temp_stack_qtr, stat = 'median', na.rm = T)

## Coldest quarters ----------------------------------------------------
coldest_qtrs_indices <- order(avg_qtr_temp$mean)[1:n_quarters]
coldest_qtrs_dates <- lapply(coldest_qtrs_indices, function(q) {
      year <- start_year + ((q - 1) %/% 4)
      qtr <- (q - 1) %% 4 + 1
      return(paste(year, qtr, sep = "-Q"))
})

## Hottest quarters ----------------------------------------------------
hottest_qtrs_indices <- order(avg_qtr_temp$mean, decreasing = TRUE)[1:n_quarters]
hottest_qtrs_dates <- lapply(hottest_qtrs_indices, function(q) {
      year <- start_year + ((q - 1) %/% 4)
      qtr <- (q - 1) %% 4 + 1
      return(paste(year, qtr, sep = "-Q"))
})

# Get total rainfall of coldest and hottest quarters ----------------------------------------------------
get_quarterly_rainfall <- function(qtrs_dates){
      quarterly_rainfall_rasters <- list()
      for (i in seq_along(qtrs_dates)) {
            q = qtrs_dates[[i]]
            year_qtr <- strsplit(q, "-Q")[[1]]
            year <- as.numeric(year_qtr[1])
            qtr <- as.numeric(year_qtr[2])
            
            q_start <- (qtr - 1) * 3 + 1
            q_end <- q_start + 2

            q_rainfall_sum <- NULL
            
            for (m in q_start:q_end) {
                  month_pattern <- sprintf("%04d%02d", year, m)
                  month_rain_files <- grep(month_pattern, rain_files_all, value = TRUE)
                  if (length(month_rain_files) > 0) {
                        month_rain_stack <- rast(month_rain_files)
                        month_rain_sum <- sum(month_rain_stack, na.rm = TRUE)
                        if (is.null(q_rainfall_sum)) {
                              q_rainfall_sum <- month_rain_sum
                        } else {
                              q_rainfall_sum <- q_rainfall_sum + month_rain_sum
                        }
                  }
            }
            quarterly_rainfall_rasters[[i]] <- q_rainfall_sum
      }

            total_rainfall_09perc <- app(rast(quarterly_rainfall_rasters), fun = calc_quantiles)
      return(total_rainfall_09perc)
}

get_quarterly_temperature <- function(qtrs_dates){
      quarterly_temp_rasters <- list()
      for (i in seq_along(qtrs_dates)) {
            q = qtrs_dates[[i]]
            year_qtr <- strsplit(q, "-Q")[[1]]
            year <- as.numeric(year_qtr[1])
            qtr <- as.numeric(year_qtr[2])
            
            q_start <- (qtr - 1) * 3 + 1
            q_end <- q_start + 2
            
            month_temp_rasters <- list()
            
            for (k in 1:length(q_start:q_end)) {
                  month <- (q_start:q_end)[k]
                  month_pattern <- sprintf("%04d%02d", year, month)
                  month_temp_sub <- temp_stack %>% subset(month_pattern)
                  month_temp_rasters[[k]] <- month_temp_sub
            }
            quarterly_temp_rasters[[i]] <- rast(month_temp_rasters) %>% median()
      }

      temp_qtr_09perc <- app(rast(month_temp_rasters), fun=calc_quantiles) 
      return(temp_qtr_09perc)
}

perc09_total_rain_coldest <- get_quarterly_rainfall(coldest_qtrs_dates)
perc09_total_rain_hottest <- get_quarterly_rainfall(hottest_qtrs_dates)

perc09_temp_coldest <- get_quarterly_temperature(coldest_qtrs_dates)
perc09_temp_hottest <- get_quarterly_temperature(hottest_qtrs_dates)

# Get dryspell duration ----------------------------------------------------
year_files <- list()
dryspell_rasters <- list()

# Organize files by year
for (file in rain_files_all) {
      year <- substr(basename(file),  31, 34) 
      if (!year %in% names(year_files)) {
            year_files[[year]] <- c(file)
      } else {
            year_files[[year]] <- c(year_files[[year]], file)
      }
}

calculate_longest_dryspell <- function(year_stack, threshold) {
      dry_day_count <- terra::init(year_stack %>% subset(1), 0) %>% 
            mask(year_stack %>% subset(1))

      for (i in 1:nlyr(year_stack)) {
            daily_rain <- terra::subset(year_stack, i)
            is_dry <- daily_rain < threshold
            # Increment dry day count where it's dry, reset to 0 where it's not
            dry_day_count <- ifel(is_dry, dry_day_count + 1, 0)

            if (i == 1) {
                  max_dry_spell <- dry_day_count
            } else {
                  max_dry_spell <- max(max_dry_spell, dry_day_count)
            }
      }
      return(max_dry_spell)
}

# Process each year
for (year in names(year_files)) {
      year_stack <- rast(year_files[[year]])
      names(year_stack) <- paste0("day", 1:nlyr(year_stack))
      hottest_qtr <- hottest_qtrs_dates[grep(year, hottest_qtrs_dates)]
      qtr <- as.numeric(substr(hottest_qtr, 7, 7))
      start_month <- (qtr - 1) * 3 + 1
      end_month <- start_month + 2
      start_date <- as.Date(paste(year, start_month, "01", sep="-"))
      end_date <- as.Date(paste(year, end_month + 1, "01", sep="-")) - 1
      days_in_hottest_qtr <- yday(seq.Date(start_date, end_date, by="day"))
      year_stack_hottest <- year_stack %>% subset(days_in_hottest_qtr)
      longest_dryspell <- calculate_longest_dryspell(year_stack_hottest, rainfall_threshold)
      dryspell_rasters[[year]] <- longest_dryspell
}

dry_duration_09perc <- app(rast(dryspell_rasters), fun=calc_quantiles)

CEDA_clim_stack <- c(perc09_total_rain_coldest, perc09_total_rain_hottest,
                     perc09_temp_coldest, perc09_temp_hottest,
                     dry_duration_09perc)
names(CEDA_clim_stack) <- c('perc09_total_rain_coldest', 'perc09_total_rain_hottest',
                            'perc09_temp_coldest', 'perc09_temp_hottest',
                            'dry_duration_09perc')

# Export ----------------------------------------------------
writeRaster(CEDA_clim_stack, "covariates/raw/CEDA_clim_stack.tif", overwrite=TRUE)
