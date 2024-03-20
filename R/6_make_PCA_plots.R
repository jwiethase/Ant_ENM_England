rm(list = ls())
library(tidyverse)
library(viridis)
library(data.table)
library(ggpubr)
source('source/misc_functions.R')

# Climate & topography ----------------------------------------------------
clim_topo_cumulative_df <- read.csv("PCA_output/csv/HPC_clim_topo_cumulative_df.csv")
clim_topo_scores_df <- fread("PCA_output/csv/HPC_clim_topo_scores_df.csv")
clim_topo_loadings_df <- read.csv("PCA_output/csv/HPC_clim_topo_loadings_df.csv") %>% dplyr::select(-X)

# nclim_topo_included <- sum(clim_topo_cumulative_df$cum_var < 0.9)
nclim_topo_included <- 6

for (i in seq(1, nclim_topo_included - 1, by = 2)) {
      j <- i + 1
      if (j <= nclim_topo_included) {
            p <- ggplot(data = clim_topo_scores_df, aes_string(x = colnames(clim_topo_scores_df)[i], 
                                                               y = colnames(clim_topo_scores_df)[j])) +
                  geom_hex(bins = 40) +  
                  theme_minimal() +
                  scale_fill_viridis(name = "Number of pixels", option = 'H', labels = custom_scientific) +
                  xlab(paste0(colnames(clim_topo_scores_df)[i],
                              " (",
                              round(clim_topo_cumulative_df$prop_var[clim_topo_cumulative_df$PCA_comp == colnames(clim_topo_scores_df)[i]]*100),
                              "%)")) +
                  ylab(paste0(colnames(clim_topo_scores_df)[j],
                              " (",
                              round(clim_topo_cumulative_df$prop_var[clim_topo_cumulative_df$PCA_comp == colnames(clim_topo_scores_df)[j]]*100),
                              "%)")) +
                  theme(plot.title = element_text(hjust = 0.5),
                        legend.position = "none",
                        legend.text = element_text(size = 12), 
                        legend.title = element_text(size = 14), 
                        legend.key.size = unit(1.3, 'cm'))
            print(p)
            plot_name <- paste("climTopo_PC", i, "_vs_PC", j, sep = "")
            ggsave(paste0("PCA_output/figures/", plot_name, ".pdf"), p, width = 5, height = 5)
            assign(plot_name, p, envir = .GlobalEnv)
      }
}

climtopo_PC_loadings <- clim_topo_loadings_df %>% 
      dplyr::select(1:(nclim_topo_included+1)) %>% 
      # Highlight the three most important loadings for each PCA axis
      mutate_if(is.numeric, function(x){ifelse(abs(x) < sort(abs(x), decreasing=T)[4], NA, x)})

# Combined plot to get common legend for Adobe Illustrator
climTopo_comb <- ggpubr::ggarrange(climTopo_PC1_vs_PC2, climTopo_PC3_vs_PC4, climTopo_PC5_vs_PC6,
                                   nrow = 2, ncol = 2,
                                   common.legend = T, legend = "bottom")
ggsave(paste0("PCA_output/figures/climTopo_comb.pdf"), climTopo_comb, width = 15, height = 15)


# Forest characteristics ----------------------------------------------------
forest_cumulative_df <- read.csv("PCA_output/csv/forest_cumulative_df_30m.csv")
forest_scores_df <- fread("PCA_output/csv/forest_scores_df_30m.csv") 
forest_loadings_df <- read.csv("PCA_output/csv/forest_loadings_df_30m.csv")

forest_PC1_vs_PC2 <- ggplot(data = forest_scores_df, aes_string(x = "PC1", 
                                           y = "PC2")) +
      geom_hex(bins = 40) +  
      theme_minimal() +
      scale_fill_viridis(name = "Number of pixels", option = 'H', labels = custom_scientific) +
      xlab(paste0("PC1", 
                  " (", 
                  round(forest_cumulative_df$prop_var[forest_cumulative_df$PCA_comp == "PC1"]*100),
                  "%)")) +
      ylab(paste0("PC2", 
                  " (", 
                  round(forest_cumulative_df$prop_var[forest_cumulative_df$PCA_comp == "PC2"]*100),
                  "%)")) +
      theme(plot.title = element_text(hjust = 0.5),
            legend.position = "bottom")

ggsave(paste0("PCA_output/figures/forest_PC1_vs_PC2.pdf"), forest_PC1_vs_PC2, width = 5, height = 5)

                         