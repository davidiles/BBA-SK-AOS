# * NEED TO CONFIRM THAT MASSIVE ERROR IN PC-ONLY MODEL IS NOT DUE TO POOR CONVERGENCE
#    ( INCREASE NUMBER OF ITERATIONS FOR A FEW SPECIES WITH CRAZY ERROR)

# - examine spatial prediction surfaces for species with crazy error (e.g., CLSW fold 1, vs folds 2-5)



# ------------------------------------------------
# Load/install packages and set graphical themes / working directory
# ------------------------------------------------

my_packs <- c(
  
  # Data management
  'tidyverse',
  
  # Spatial functions
  'sf','raster','ggspatial','rgeos','fasterize','exactextractr','circular',
  'lubridate',
  
  # For PCA
  'factoextra',
  
  # Spatial analysis
  'inlabru','ebirdst',
  
  # For plotting
  'viridis','scales','ggpubr','ggtext','ggrepel',
  
  # Cross-validation
  'pROC')

if (any(!my_packs %in% installed.packages()[, 'Package'])) {install.packages(my_packs[which(!my_packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)}
lapply(my_packs, require, character.only = TRUE)

library(INLA) # install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep = TRUE)
library(napops) # For detectability offsets  # devtools::install_github("na-pops/napops") 
# napops::fetch_data()  # - only needs to be run once

rm(list=ls())

# Import rasters, data, and covariates from "standard analysis"
setwd("D:/Working_Files/1_Projects/Landbirds/SK_BBA_analysis/Standard_Analysis/")
`%!in%` <- Negate(`%in%`)

load("../AOS_precision/output/xval_df_integrated_20km.RData")

xval_df <- xval_df %>%
  group_by(Species) %>%
  summarize_all(mean)

xval_df$n_obs_CL <- xval_df$n_obs_SC + xval_df$n_obs_LT + xval_df$n_obs_BBA
xval_df$delta_cor <- xval_df$cor_integrated - xval_df$cor_PConly
xval_df$delta_MSE <- xval_df$MSE_integrated - xval_df$MSE_PConly
xval_df$delta_AUC <- xval_df$AUC_integrated - xval_df$AUC_PConly

mean(xval_df$delta_cor>0)
mean(xval_df$delta_AUC>0)
mean(xval_df$delta_MSE<0)

ggplot()+
  geom_text(data = xval_df, aes(x = n_obs_PC, 
                                y = cor_integrated + sign(delta_cor)*0.01, label = Species,col = delta_cor > 0))+
  geom_segment(data = xval_df, aes(x = n_obs_PC,xend = n_obs_PC, 
                                   y = (cor_PConly ),yend = cor_integrated, col = delta_cor > 0),
               size = 2,
               arrow = arrow(length = unit(0.05, "inches")))+
  scale_color_manual(values=c("red","blue"), name = "Change in cor",
                     labels = c("Decrease","Increase"), guide = "none")+
  scale_x_continuous(trans="log10", name = "Number of Detections in Point Counts")+
  ylab("Correlation\n(Crossvalidation)")+
  theme_bw()

ggplot()+
  geom_text(data = xval_df, aes(x = n_obs_PC, 
                                y = AUC_integrated + sign(delta_AUC)*0.01, label = Species,col = delta_AUC > 0))+
  geom_segment(data = xval_df, aes(x = n_obs_PC,xend = n_obs_PC, 
                                   y = (AUC_PConly ),yend = AUC_integrated, col = delta_AUC > 0),
               size = 2,
               arrow = arrow(length = unit(0.05, "inches")))+
  scale_color_manual(values=c("red","blue"), name = "Change in AUC",
                     labels = c("Decrease","Increase"))+
  scale_x_continuous(trans="log10", name = "Number of Detections in Point Counts")+
  ylab("AUC\n(Crossvalidation)")+
  theme_bw()

ggplot()+
  geom_text(data = xval_df, aes(x = n_obs_PC, 
                                y = MSE_integrated + sign(delta_MSE)*MSE_integrated*0.2, label = Species,col = delta_MSE < 0))+
  geom_segment(data = xval_df, aes(x = n_obs_PC,xend = n_obs_PC, 
                                   y = MSE_PConly,yend = MSE_integrated, col = delta_MSE < 0),
               size = 2,
               arrow = arrow(length = unit(0.05, "inches")))+
  scale_color_manual(values=c("red","blue"), name = "Change in MSE",
                     labels = c("Increase","Decrease"))+
  scale_x_continuous(trans="log10", name = "Number of Detections in Point Counts")+
  scale_y_continuous(trans="log10", name = "MSE\n(Crossvalidation)")+
  theme_bw()

# -----------------------------------------------------------
# Examine results across each cross-validation fold
# -----------------------------------------------------------

load("../AOS_precision/output/xval_df_integrated_20km.RData")

xval_df$n_obs_CL <- xval_df$n_obs_SC + xval_df$n_obs_LT + xval_df$n_obs_BBA
xval_df$delta_cor <- xval_df$cor_integrated - xval_df$cor_PConly
xval_df$delta_MSE <- xval_df$MSE_integrated - xval_df$MSE_PConly
xval_df$delta_AUC <- xval_df$AUC_integrated - xval_df$AUC_PConly

ggplot()+
 geom_segment(data = xval_df, aes(x = xval_fold,xend = xval_fold, 
                                   y = (cor_PConly ),yend = cor_integrated, col = delta_cor > 0),
               size = 2,
               arrow = arrow(length = unit(0.05, "inches")))+
  scale_color_manual(values=c("red","blue"), name = "Change in cor",
                     labels = c("Decrease","Increase"), guide = "none")+
  #scale_x_continuous(trans="log10", name = "Number of Detections in Point Counts")+
  ylab("Correlation\n(Crossvalidation)")+
  theme_bw()+
  facet_wrap(Species~.)

ggplot()+
  geom_segment(data = xval_df, aes(x = xval_fold,xend = xval_fold, 
                                   y = (AUC_PConly ),yend = AUC_integrated, col = delta_AUC > 0),
               size = 2,
               arrow = arrow(length = unit(0.05, "inches")))+
  scale_color_manual(values=c("red","blue"), name = "Change in AUC",
                     labels = c("Decrease","Increase"), guide = "none")+
  #scale_x_continuous(trans="log10", name = "Number of Detections in Point Counts")+
  ylab("AUC\n(Crossvalidation)")+
  theme_bw()+
  facet_wrap(Species~.)

ggplot()+
  geom_segment(data = xval_df, aes(x = xval_fold,xend = xval_fold, 
                                   y = (MSE_PConly ),yend = MSE_integrated, col = delta_MSE < 0),
               size = 2,
               arrow = arrow(length = unit(0.05, "inches")))+
  scale_color_manual(values=c("red","blue"), name = "Change in MSE",
                     labels = c("Decrease","Increase"), guide = "none")+
  #scale_x_continuous(trans="log10", name = "Number of Detections in Point Counts")+
  ylab("MSE\n(Crossvalidation)")+
  theme_bw()+
  facet_wrap(Species~., scales = "free_y")

# -----------------------------------------------------------
# Compare surfaces generated from point count only or integrated models
# -----------------------------------------------------------

load(file = "../AOS_precision/output/surface_comparison_100km.RData")

surface_comparison

#ggplot(surface_comparison, aes())+
  


# # -----------------------------------------------------------
# # Compare "Two Step" analysis with fully integrated analysis
# # -----------------------------------------------------------
# 
# rm(list=ls())
# load("../AOS_precision/output/xval_df_integrated_100km_noRE.RData")
# xval_df_integrated <- xval_df
# 
# load("../AOS_precision/output/xval_df_integrated_100km_TwoStep.RData")
# xval_df_TwoStep <- xval_df
# 
# xval_comparison <- full_join(xval_df_integrated[,c("Species","xval_fold","AUC_integrated","cor_integrated","MSE_integrated")],xval_df_TwoStep)
# xval_comparison$delta_cor <- xval_comparison$cor_integrated - xval_comparison$cor_TwoStep
# xval_comparison$delta_MSE <- xval_comparison$MSE_integrated - xval_comparison$MSE_TwoStep
# xval_comparison$delta_AUC <- xval_comparison$AUC_integrated - xval_comparison$AUC_TwoStep
# 
# xval_comparison[,c("Species","xval_fold","cor_integrated","cor_TwoStep","AUC_integrated","AUC_TwoStep","MSE_integrated","MSE_TwoStep")] %>% na.omit()
# 
# ggplot()+
#   geom_segment(data = xval_comparison, aes(x = xval_fold,xend = xval_fold, 
#                                    y = cor_TwoStep,yend = cor_integrated, col = delta_cor > 0),
#                size = 2,
#                arrow = arrow(length = unit(0.05, "inches")))+
#   scale_color_manual(values=c("red","blue"), name = "Change in cor",
#                      labels = c("Decrease","Increase"), guide = "none")+
#   scale_x_continuous(name = "Crossval fold")+
#   ylab("Correlation\n(Crossvalidation)")+
#   theme_bw()+
#   facet_wrap(Species~.)
# 
# ggplot()+
#   geom_segment(data = xval_comparison, aes(x = xval_fold,xend = xval_fold, 
#                                            y = AUC_TwoStep,yend = AUC_integrated, col = delta_AUC > 0),
#                size = 2,
#                arrow = arrow(length = unit(0.05, "inches")))+
#   scale_color_manual(values=c("red","blue"), name = "Change in AUC",
#                      labels = c("Decrease","Increase"), guide = "none")+
#   scale_x_continuous(name = "Crossval fold")+
#   ylab("AUC\n(Crossvalidation)")+
#   theme_bw()+
#   facet_wrap(Species~.)
# 
# ggplot()+
#   geom_segment(data = xval_comparison, aes(x = xval_fold,xend = xval_fold, 
#                                            y = MSE_TwoStep,yend = MSE_integrated, col = delta_MSE > 0),
#                size = 2,
#                arrow = arrow(length = unit(0.05, "inches")))+
#   scale_color_manual(values=c("blue","red"), name = "Change in MSE",
#                      labels = c("Decrease","Increase"), guide = "none")+
#   scale_x_continuous(name = "Crossval fold")+
#   ylab("MSE\n(Crossvalidation)")+
#   theme_bw()+
#   facet_wrap(Species~.)

