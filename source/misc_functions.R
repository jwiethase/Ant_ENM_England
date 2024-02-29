km_proj <- "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=km +no_defs"

determine_accuracy <- function(grid_ref) {
      accuracy_map <- c('2' = 10000, '4' = 1000, '6' = 100, '8' = 10, '10' = 1)
      num_part <- gsub("[^0-9]", "", grid_ref)
      num_length <- as.character(nchar(num_part))
      return(accuracy_map[num_length])
}

rsq <- function (x, y) cor(x, y) ^ 2

normalise_raster <- function(spatRaster){
      nx <- minmax(spatRaster)
      rn <- (spatRaster - nx[1,]) / (nx[2,] - nx[1,])
      return(rn)
} 

prepareMGCVsplines <- function(df, colname, n_knots = 3, smooth = 'tp'){
      #' Function based on discussions in: https://groups.google.com/g/r-inla-discussion-group/c/CiA4l9zhCMw
      #' See ?smooth.terms for a list of available spline functions
      require(mgcv)
            expr_str <- paste0("smoothCon(s(", colname, ", bs = '", smooth, "', k = ", n_knots, ", fx = TRUE), data = df, absorb.cons = TRUE)[[1]]$X")
            base_functions <- eval(parse(text = expr_str))
            
            num_columns <- if (n_knots > 2) { n_knots - 1 } else { n_knots }
            
            for(i in 1:num_columns){
                  df[, paste0(colname, "_spline", i)] <- base_functions[, i]
            }
      return(df)
}

range01 <- function(x){(x - min(x, na.rm = T))/(max(x, na.rm = T) - min(x, na.rm = T))}

calc_quantiles <- function(x, probs = 0.9) {
      return(stats::quantile(x, probs = probs, na.rm = TRUE))
}

calculate_cover <- function(x, threshold = 2) {
      valid_pixels <- !is.na(x)
      cover <- sum(x[valid_pixels] > threshold, na.rm = TRUE) / sum(valid_pixels)
      return(cover)
}

mean_fun <- function(x) mean(x, na.rm = TRUE)
sd_fun <- function(x) sd(x, na.rm = TRUE)

make_formula_terms <- function(rasterstack, name_leading, n_splines){
      unique_vars <- grep(paste0("^", name_leading, "[0-9]$"), names(rasterstack), value = TRUE)
      unique_vars_num <- as.numeric(gsub("\\D", "", unique_vars))
      spline_terms <- paste0(name_leading, 
                             rep(unique_vars_num, each = n_splines), "_spline", rep(1:n_splines, times = length(unique_vars)), 
                             "(main = ", deparse(substitute(rasterstack)), "$", name_leading, "", rep(unique_vars_num, each = n_splines), "_spline", rep(1:n_splines, times = length(unique_vars)), 
                             ", model = 'linear')")
      return(spline_terms)
}

disagg_to_target <- function(coarse, fine, resample = FALSE, method_disagg = "bilinear"){
      factor <- round(res(coarse)[1]/res(fine)[1])
      print(paste0("Disaggregating by factor of ", factor, "..."))
      out_rast <- terra::disagg(coarse, fact = factor, method = method_disagg)
      if(resample){
            out_rast <- out_rast %>% terra::resample(fine, threads = T, method = "bilinear")  
      }
      return(out_rast)
}

custom_scientific <- function(x) {
      formatted <- formatC(x, format = "e", digits = 1) %>%
            gsub("e\\+0*", "e", .) %>%
            gsub("e\\-", "e-", .)
      formatted <- gsub("\\.0", "", formatted)
      return(formatted)
}

rast_yeo_johnson <- function(rasterlayer, plot = FALSE){
      require(bestNormalize)
      yj_obj <- bestNormalize::yeojohnson(values(rasterlayer), standardize = T)
      transformed_lyr <- rasterlayer
      values(transformed_lyr) <- predict(yj_obj)
      if(plot){
            par(mfrow = c(2, 1))
            plot(density(values(rasterlayer, na.rm = T)), main = "Original density distr.")
            plot(density(values(transformed_lyr, na.rm = T)), main = "Transformed density distr.")
            par(mfrow = c(1, 1))
      }
      return(transformed_lyr)
}
