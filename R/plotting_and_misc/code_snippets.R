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