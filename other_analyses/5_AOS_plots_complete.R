# ------------------------------------------------
# Load/install packages and set graphical themes / working directory
# ------------------------------------------------

my_packs <- c(
  
  # Data management
  'tidyverse',
  
  # For plotting
  'viridis','scales','ggpubr','ggtext','ggrepel','ggthemes')

if (any(!my_packs %in% installed.packages()[, 'Package'])) {install.packages(my_packs[which(!my_packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)}
lapply(my_packs, require, character.only = TRUE)

rm(list=ls())

# Import rasters, data, and covariates from "standard analysis"
setwd("D:/Working_Files/1_Projects/Landbirds/SK_BBA_analysis/AOS_precision/script")
`%!in%` <- Negate(`%in%`)

# **********************************************************
# **********************************************************
# LOAD AND SUMMARIZE RESULTS OF ANALYSIS
# **********************************************************
# **********************************************************

load("../output/xval_PC_CL_complete.RData")

# Number of species for each fold
nspec <- xval_PC_CL_complete %>%
  group_by(xval_fold) %>%
  summarize(nspec = n())

nfold <- xval_PC_CL_complete %>%
  group_by(Species) %>%
  summarize(nfold = n()) %>%
  subset(nfold >=2)

# Arrange species by relative abundance
xval_PC_CL_complete <- xval_PC_CL_complete %>%
  subset(xval_fold %in% c(1,2)) %>%
  group_by(Species) %>%
  summarize_all(mean) %>% 
  arrange(n_det_PC) %>%
  subset(Species %in% nfold$Species)

# -----------------------------------
# Summary of AUC scores
# -----------------------------------
sum(xval_PC_CL_complete$AUC_integrated_50_50 > xval_PC_CL_complete$AUC_PC_50)
sum(xval_PC_CL_complete$AUC_PC_100 > xval_PC_CL_complete$AUC_PC_50)
sum(xval_PC_CL_complete$AUC_integrated_50_100 > xval_PC_CL_complete$AUC_PC_50)

# # -----------------------------------
# # Summary of lppd scores (CALCULATION IS MESSED UP)
# # -----------------------------------
# sum(xval_PC_CL_complete$lppd_integrated_50_50 > xval_PC_CL_complete$lppd_PC_50)
# sum(xval_PC_CL_complete$lppd_PC_100 > xval_PC_CL_complete$lppd_PC_50)
# sum(xval_PC_CL_complete$lppd_integrated_50_100 > xval_PC_CL_complete$lppd_PC_50)

# -----------------------------------
# Summary of correlation scores
# -----------------------------------
sum(xval_PC_CL_complete$cor_integrated_50_50 > xval_PC_CL_complete$cor_PC_50)
sum(xval_PC_CL_complete$cor_PC_100 > xval_PC_CL_complete$cor_PC_50)
sum(xval_PC_CL_complete$cor_integrated_50_100 > xval_PC_CL_complete$cor_PC_50)

# -----------------------------------
# Summary of MSE scores
# -----------------------------------
sum(xval_PC_CL_complete$MSE_integrated_50_50 < xval_PC_CL_complete$MSE_PC_50)
sum(xval_PC_CL_complete$MSE_PC_100 < xval_PC_CL_complete$MSE_PC_50)
sum(xval_PC_CL_complete$MSE_integrated_50_100 < xval_PC_CL_complete$MSE_PC_50)

# -----------------------------------
# Is doubling checklist coverage better than doubling point counts?
# -----------------------------------

mean(xval_PC_CL_complete$cor_integrated_50_100 > xval_PC_CL_complete$cor_PC_100)
mean(xval_PC_CL_complete$AUC_integrated_50_100 > xval_PC_CL_complete$AUC_PC_100)
mean(xval_PC_CL_complete$MSE_integrated_50_100 < xval_PC_CL_complete$MSE_PC_100)
mean(xval_PC_CL_complete$lppd_integrated_50_100 > xval_PC_CL_complete$lppd_PC_100)

# -----------------------------------
# Is integrated model better than point count only model?
# -----------------------------------

mean(xval_PC_CL_complete$cor_integrated_100_100 > xval_PC_CL_complete$cor_PC_100)
mean(xval_PC_CL_complete$AUC_integrated_100_100 > xval_PC_CL_complete$AUC_PC_100)
mean(xval_PC_CL_complete$MSE_integrated_100_100 < xval_PC_CL_complete$MSE_PC_100)
mean(xval_PC_CL_complete$lppd_integrated_100_100 > xval_PC_CL_complete$lppd_PC_100)

# **********************************************************
# **********************************************************
# GENERATE SOME PLOTS TO ILLUSTRATE DIFFERENCES AMONG SPECIES
# **********************************************************
# **********************************************************

# -----------------------------------------------------------------
# AUC
# -----------------------------------------------------------------

xval_PC_CL_complete$Species_label <- paste0(xval_PC_CL_complete$Species," (",round(xval_PC_CL_complete$n_det_PC),")")

xval_PC_CL_complete$Species <- factor(xval_PC_CL_complete$Species, levels = xval_PC_CL_complete$Species)
xval_PC_CL_complete$Species_label <- factor(xval_PC_CL_complete$Species_label, levels = xval_PC_CL_complete$Species_label)


lim <- c(0.5,max(xval_PC_CL_complete[,c("AUC_PC_50","AUC_PC_100","AUC_integrated_50_50","AUC_integrated_50_100","AUC_integrated_100_100")]))

{ labels1 <- data.frame(Species_label = xval_PC_CL_complete$Species_label,
                        Improvement = (xval_PC_CL_complete$AUC_integrated_50_50 > xval_PC_CL_complete$AUC_PC_50))
  labels1$sign <- "+"
  labels1$sign[labels1$Improvement == FALSE] <- ""
  
  labels2 <- data.frame(Species_label = xval_PC_CL_complete$Species_label,
                        Improvement = (xval_PC_CL_complete$AUC_PC_100 > xval_PC_CL_complete$AUC_PC_50))
  labels2$sign <- "+"
  labels2$sign[labels2$Improvement == FALSE] <- ""
  
  labels3 <- data.frame(Species_label = xval_PC_CL_complete$Species_label,
                        Improvement = (xval_PC_CL_complete$AUC_integrated_50_100 > xval_PC_CL_complete$AUC_PC_50))
  labels3$sign <- "+"
  labels3$sign[labels3$Improvement == FALSE] <- ""
  
}


plot_AUC_0 <- ggplot()+
  
  # AUC between predictions and validation data
  geom_point(data = xval_PC_CL_complete, aes(y = Species_label,
                                             x = AUC_PC_50),
             size = 2)+
  
  ylab("Species")+
  xlab("AUC between predictions and validation counts")+
  scale_x_continuous(limits = lim)+
  theme_few()+
  ggtitle("PC vs Integrated\n\nAUC with validation data")

plot_AUC_1 <- ggplot()+
  
  geom_text(data = labels1, aes(x = min(lim),
                                y = Species_label,
                                label = sign),
            fontface = "bold")+
  # Change in AUC with integrated model
  geom_segment(data = xval_PC_CL_complete, aes(y = Species_label,
                                               yend = Species_label, 
                                               x = AUC_PC_50,
                                               xend = AUC_integrated_50_50),
               size = 1.5,
               lineend = "butt",
               linejoin = "mitre",
               col = "black",
               alpha = 1,
               arrow = arrow(length = unit(0.1, "inches")))+
  
  # AUC between predictions and validation data
  geom_point(data = xval_PC_CL_complete, aes(y = Species_label,
                                             x = AUC_PC_50),
             size = 2)+
  
  ylab("Species")+
  xlab("AUC between predictions and validation counts")+
  scale_x_continuous(limits = lim)+
  theme_few()+
  ggtitle("PC vs Integrated\n\nAUC with validation data")

plot_AUC_2 <- ggplot()+
  
  geom_text(data = labels1, aes(x = min(lim),
                                y = Species_label,
                                label = sign),
            fontface = "bold")+
  
  geom_text(data = labels2, aes(x = min(lim)+0.01,
                                y = Species_label,
                                label = sign),
            fontface = "bold",
            col = "deepskyblue")+
  
  # Change in AUC when doubling point count coverage
  geom_segment(data = xval_PC_CL_complete, aes(y = Species_label,
                                               yend = Species_label,
                                               x = AUC_PC_50,
                                               xend = AUC_PC_100),
               size = 3,
               lineend = "butt",
               linejoin = "mitre",
               col = "deepskyblue",
               arrow = arrow(length = unit(0.1, "inches")))+
  
  # Change in AUC with integrated model
  geom_segment(data = xval_PC_CL_complete, aes(y = Species_label,
                                               yend = Species_label, 
                                               x = AUC_PC_50,
                                               xend = AUC_integrated_50_50),
               size = 1.5,
               lineend = "butt",
               linejoin = "mitre",
               col = "black",
               alpha = 1,
               arrow = arrow(length = unit(0.1, "inches")))+
  
  # AUC between predictions and validation data
  geom_point(data = xval_PC_CL_complete, aes(y = Species_label,
                                             x = AUC_PC_50),
             size = 2)+
  
  ylab("Species")+
  xlab("AUC between predictions and validation counts")+
  scale_x_continuous(limits = lim)+
  theme_few()+
  ggtitle("PC vs Integrated\n\nAUC with validation data")

plot_AUC_3 <- ggplot()+
  geom_text(data = labels1, aes(x = min(lim),
                                y = Species_label,
                                label = sign),
            fontface = "bold")+
  
  geom_text(data = labels2, aes(x = min(lim)+0.01,
                                y = Species_label,
                                label = sign),
            fontface = "bold",
            col = "deepskyblue")+
  
  geom_text(data = labels3, aes(x = min(lim)+0.02,
                                y = Species_label,
                                label = sign),
            fontface = "bold",
            col = "forestgreen")+
  
  # Change in AUC when doubling checklist coverage
  geom_segment(data = xval_PC_CL_complete, aes(y = Species_label,
                                               yend = Species_label,
                                               x = AUC_PC_50,
                                               xend = AUC_integrated_50_100),
               size = 6,
               lineend = "butt",
               linejoin = "mitre",
               col = "forestgreen",
               arrow = arrow(length = unit(0.1, "inches")))+
  
  # Change in AUC when doubling point count coverage
  geom_segment(data = xval_PC_CL_complete, aes(y = Species_label,
                                               yend = Species_label,
                                               x = AUC_PC_50,
                                               xend = AUC_PC_100),
               size = 3,
               lineend = "butt",
               linejoin = "mitre",
               col = "deepskyblue",
               arrow = arrow(length = unit(0.1, "inches")))+
  
  # Change in AUC with integrated model
  geom_segment(data = xval_PC_CL_complete, aes(y = Species_label,
                                               yend = Species_label, 
                                               x = AUC_PC_50,
                                               xend = AUC_integrated_50_50),
               size = 1.5,
               lineend = "butt",
               linejoin = "mitre",
               col = "black",
               alpha = 1,
               arrow = arrow(length = unit(0.1, "inches")))+
  
  # AUC between predictions and validation data
  geom_point(data = xval_PC_CL_complete, aes(y = Species_label,
                                             x = AUC_PC_50),
             size = 2)+
  
  ylab("Species")+
  xlab("AUC between predictions and validation counts")+
  scale_x_continuous(limits = lim)+
  theme_few()+
  ggtitle("PC vs Integrated\n\nAUC with validation data")
plot_AUC_3

plot_AUC_3 <- ggplot()+
  geom_text(data = labels1, aes(x = min(lim),
                                y = Species_label,
                                label = sign),
            fontface = "bold")+
  
  geom_text(data = labels2, aes(x = min(lim)+0.01,
                                y = Species_label,
                                label = sign),
            fontface = "bold",
            col = "deepskyblue")+
  
  geom_text(data = labels3, aes(x = min(lim)+0.02,
                                y = Species_label,
                                label = sign),
            fontface = "bold",
            col = "forestgreen")+
  
  # Change in AUC when doubling checklist coverage
  geom_segment(data = xval_PC_CL_complete, aes(y = Species_label,
                                               yend = Species_label,
                                               x = AUC_PC_50,
                                               xend = AUC_integrated_50_100),
               size = 6,
               lineend = "butt",
               linejoin = "mitre",
               col = "forestgreen",
               arrow = arrow(length = unit(0.1, "inches")))+
  
  # Change in AUC when doubling point count coverage
  geom_segment(data = xval_PC_CL_complete, aes(y = Species_label,
                                               yend = Species_label,
                                               x = AUC_PC_50,
                                               xend = AUC_PC_100),
               size = 3,
               lineend = "butt",
               linejoin = "mitre",
               col = "deepskyblue",
               arrow = arrow(length = unit(0.1, "inches")))+
  
  # Change in AUC with integrated model
  geom_segment(data = xval_PC_CL_complete, aes(y = Species_label,
                                               yend = Species_label, 
                                               x = AUC_PC_50,
                                               xend = AUC_integrated_50_50),
               size = 1.5,
               lineend = "butt",
               linejoin = "mitre",
               col = "black",
               alpha = 1,
               arrow = arrow(length = unit(0.1, "inches")))+
  
  # AUC between predictions and validation data
  geom_point(data = xval_PC_CL_complete, aes(y = Species_label,
                                             x = AUC_PC_50),
             size = 2)+
  
  ylab("Species")+
  xlab("AUC between predictions and validation counts")+
  scale_x_continuous(limits = lim)+
  theme_few()+
  ggtitle("PC vs Integrated\n\nAUC with validation data")
plot_AUC_3

# # Output figures
# png("../output/figures/xval_plot_AUC_0.png", width=12, height=6, units="in", res=300, type="cairo")
# print(plot_AUC_0)
# dev.off()
# 
# png("../output/figures/xval_plot_AUC_1.png", width=12, height=6, units="in", res=300, type="cairo")
# print(plot_AUC_1)
# dev.off()
# 
# png("../output/figures/xval_plot_AUC_2.png", width=12, height=6, units="in", res=300, type="cairo")
# print(plot_AUC_2)
# dev.off()
# 
# png("../output/figures/xval_plot_AUC_3.png", width=12, height=6, units="in", res=300, type="cairo")
# print(plot_AUC_3)
# dev.off()

sum(xval_PC_CL_complete$AUC_integrated_50_50 > xval_PC_CL_complete$AUC_PC_50)
sum(xval_PC_CL_complete$AUC_PC_100 > xval_PC_CL_complete$AUC_PC_50)
sum(xval_PC_CL_complete$AUC_integrated_50_100 > xval_PC_CL_complete$AUC_PC_50)
nrow(xval_PC_CL_complete)

# -----------------------------------------------------------------
# Correlation
# -----------------------------------------------------------------

xval_PC_CL_complete$Species = factor(xval_PC_CL_complete$Species,levels = xval_PC_CL_complete$Species)

lim <- range(xval_PC_CL_complete[,c("cor_PC_50","cor_PC_100","cor_integrated_50_50","cor_integrated_50_100")])

{ labels1 <- data.frame(Species = xval_PC_CL_complete$Species,
                        Improvement = (xval_PC_CL_complete$cor_integrated_50_50 > xval_PC_CL_complete$cor_PC_50))
  labels1$sign <- "+"
  labels1$sign[labels1$Improvement == FALSE] <- ""
  
  labels2 <- data.frame(Species = xval_PC_CL_complete$Species,
                        Improvement = (xval_PC_CL_complete$cor_PC_100 > xval_PC_CL_complete$cor_PC_50))
  labels2$sign <- "+"
  labels2$sign[labels2$Improvement == FALSE] <- ""
  
  labels3 <- data.frame(Species = xval_PC_CL_complete$Species,
                        Improvement = (xval_PC_CL_complete$cor_integrated_50_100 > xval_PC_CL_complete$cor_PC_50))
  labels3$sign <- "+"
  labels3$sign[labels3$Improvement == FALSE] <- ""
  
}

plot_cor_1 <- ggplot()+
  
  # Change in correlation with integrated model
  geom_segment(data = xval_PC_CL_complete, aes(y = Species,
                                               yend = Species, 
                                               x = cor_PC_50,
                                               xend = cor_integrated_50_50),
               size = 1.5,
               lineend = "butt",
               linejoin = "mitre",
               col = "black",
               alpha = 1,
               arrow = arrow(length = unit(0.1, "inches")))+
  
  # Correlation between predictions and validation data
  geom_point(data = xval_PC_CL_complete, aes(y = Species,
                                             x = cor_PC_50),
             size = 2)+
  
  ylab("Species")+
  xlab("Correlation between predictions and validation counts")+
  scale_x_continuous(limits = lim)+
  theme_few()+
  ggtitle("PC vs Integrated\n\nCorrelation with validation data")

plot_cor_2 <- ggplot()+
  
  # Change in correlation when doubling point count coverage
  geom_segment(data = xval_PC_CL_complete, aes(y = Species,
                                               yend = Species,
                                               x = cor_PC_50,
                                               xend = cor_PC_100),
               size = 3,
               lineend = "butt",
               linejoin = "mitre",
               col = "deepskyblue",
               arrow = arrow(length = unit(0.1, "inches")))+
  
  # Change in correlation with integrated model
  geom_segment(data = xval_PC_CL_complete, aes(y = Species,
                                               yend = Species, 
                                               x = cor_PC_50,
                                               xend = cor_integrated_50_50),
               size = 1.5,
               lineend = "butt",
               linejoin = "mitre",
               col = "black",
               alpha = 1,
               arrow = arrow(length = unit(0.1, "inches")))+
  
  # Correlation between predictions and validation data
  geom_point(data = xval_PC_CL_complete, aes(y = Species,
                                             x = cor_PC_50),
             size = 2)+
  
  ylab("Species")+
  xlab("Correlation between predictions and validation counts")+
  scale_x_continuous(limits = lim)+
  theme_few()+
  ggtitle("PC vs Integrated\n\nCorrelation with validation data")

plot_cor_3 <- ggplot()+
  
  # Change in correlation when doubling checklist coverage
  geom_segment(data = xval_PC_CL_complete, aes(y = Species,
                                               yend = Species,
                                               x = cor_PC_50,
                                               xend = cor_integrated_50_100),
               size = 6,
               lineend = "butt",
               linejoin = "mitre",
               col = "forestgreen",
               arrow = arrow(length = unit(0.1, "inches")))+
  
  # Change in correlation when doubling point count coverage
  geom_segment(data = xval_PC_CL_complete, aes(y = Species,
                                               yend = Species,
                                               x = cor_PC_50,
                                               xend = cor_PC_100),
               size = 3,
               lineend = "butt",
               linejoin = "mitre",
               col = "deepskyblue",
               arrow = arrow(length = unit(0.1, "inches")))+
  
  # Change in correlation with integrated model
  geom_segment(data = xval_PC_CL_complete, aes(y = Species,
                                               yend = Species, 
                                               x = cor_PC_50,
                                               xend = cor_integrated_50_50),
               size = 1.5,
               lineend = "butt",
               linejoin = "mitre",
               col = "black",
               alpha = 1,
               arrow = arrow(length = unit(0.1, "inches")))+
  
  # Correlation between predictions and validation data
  geom_point(data = xval_PC_CL_complete, aes(y = Species,
                                             x = cor_PC_50),
             size = 2)+
  
  ylab("Species")+
  xlab("Correlation between predictions and validation counts")+
  scale_x_continuous(limits = lim)+
  theme_few()+
  ggtitle("PC vs Integrated\n\nCorrelation with validation data")

# plot_cor_1
# plot_cor_2
#plot_cor_3


# # -----------------------------------------------------------------
# # RMSE
# # -----------------------------------------------------------------
# 
# xval_PC_CL_complete <- xval_PC_CL_complete %>%
#   mutate(RMSE_PC_50 = sqrt(MSE_PC_50),
#          RMSE_PC_100 = sqrt(MSE_PC_100),
#          RMSE_integrated_50_50 = sqrt(MSE_integrated_50_50),
#          RMSE_integrated_50_100 = sqrt(MSE_integrated_50_100))
# lim <- range(xval_PC_CL_complete[,c("RMSE_PC_50","RMSE_PC_100","RMSE_integrated_50_50","RMSE_integrated_50_100")])
# 
# plot_RMSE_1 <- ggplot()+
#   
#   # Change in RMSE with integrated model
#   geom_segment(data = xval_PC_CL_complete, aes(y = Species,
#                                              yend = Species, 
#                                              x = RMSE_PC_50,
#                                              xend = RMSE_integrated_50_50),
#                size = 1.5,
#                lineend = "butt",
#                linejoin = "mitre",
#                col = "black",
#                alpha = 1,
#                arrow = arrow(length = unit(0.1, "inches")))+
#   
#   # RMSE between predictions and validation data
#   geom_point(data = xval_PC_CL_complete, aes(y = Species,
#                                            x = RMSE_PC_50),
#              size = 2)+
#   
#   ylab("Species")+
#   xlab("RMSE between predictions and validation counts")+
#   scale_x_continuous(limits = lim)+
#   theme_few()+
#   ggtitle("PC vs Integrated\n\nRMSE with validation data")
# 
# plot_RMSE_2 <- ggplot()+
#   
#   # Change in RMSE when doubling point count coverage
#   geom_segment(data = xval_PC_CL_complete, aes(y = Species,
#                                              yend = Species,
#                                              x = RMSE_PC_50,
#                                              xend = RMSE_PC_100),
#                size = 3,
#                lineend = "butt",
#                linejoin = "mitre",
#                col = "deepskyblue",
#                arrow = arrow(length = unit(0.1, "inches")))+
#   
#   # Change in RMSE with integrated model
#   geom_segment(data = xval_PC_CL_complete, aes(y = Species,
#                                              yend = Species, 
#                                              x = RMSE_PC_50,
#                                              xend = RMSE_integrated_50_50),
#                size = 1.5,
#                lineend = "butt",
#                linejoin = "mitre",
#                col = "black",
#                alpha = 1,
#                arrow = arrow(length = unit(0.1, "inches")))+
#   
#   # RMSE between predictions and validation data
#   geom_point(data = xval_PC_CL_complete, aes(y = Species,
#                                            x = RMSE_PC_50),
#              size = 2)+
#   
#   ylab("Species")+
#   xlab("RMSE between predictions and validation counts")+
#   scale_x_continuous(limits = lim)+
#   theme_few()+
#   ggtitle("PC vs Integrated\n\nRMSE with validation data")
# 
# plot_RMSE_3 <- ggplot()+
#   
#   # Change in RMSE when doubling checklist coverage
#   geom_segment(data = xval_PC_CL_complete, aes(y = Species,
#                                              yend = Species,
#                                              x = RMSE_PC_50,
#                                              xend = RMSE_integrated_50_100),
#                size = 6,
#                lineend = "butt",
#                linejoin = "mitre",
#                col = "forestgreen",
#                arrow = arrow(length = unit(0.1, "inches")))+
#   
#   # Change in RMSE when doubling point count coverage
#   geom_segment(data = xval_PC_CL_complete, aes(y = Species,
#                                              yend = Species,
#                                              x = RMSE_PC_50,
#                                              xend = RMSE_PC_100),
#                size = 3,
#                lineend = "butt",
#                linejoin = "mitre",
#                col = "deepskyblue",
#                arrow = arrow(length = unit(0.1, "inches")))+
#   
#   # Change in RMSE with integrated model
#   geom_segment(data = xval_PC_CL_complete, aes(y = Species,
#                                              yend = Species, 
#                                              x = RMSE_PC_50,
#                                              xend = RMSE_integrated_50_50),
#                size = 1.5,
#                lineend = "butt",
#                linejoin = "mitre",
#                col = "black",
#                alpha = 1,
#                arrow = arrow(length = unit(0.1, "inches")))+
#   
#   # RMSE between predictions and validation data
#   geom_point(data = xval_PC_CL_complete, aes(y = Species,
#                                            x = RMSE_PC_50),
#              size = 2)+
#   
#   ylab("Species")+
#   xlab("RMSE between predictions and validation counts")+
#   scale_x_continuous(limits = lim,trans="log10")+
#   theme_few()+
#   ggtitle("PC vs Integrated\n\nRMSE with validation data")
# 
# plot_RMSE_1
# plot_RMSE_2
# plot_RMSE_3

# **********************************************************
# **********************************************************
# MAPS ILLUSTRATING CROSS-VALIDATION APPROACH
# **********************************************************
# **********************************************************
library(sf)
load("../output/AOS_data_package_test.RData")

for (fold in 1:5){
  set.seed(fold*100)
  
  # --------------------------------
  # Load spatial data
  # --------------------------------
  PC_sf <- PC_surveyinfo
  CL_sf <- DO_surveyinfo
  
  # Intersect with SaskSquares dataframe
  PC_sf <- st_intersection(PC_sf, SaskSquares)
  CL_sf <- st_intersection(CL_sf, SaskSquares) %>%
    subset(ProtocolType %in% c("Stationary count","Linear transect"))
  
  # --------------------------------
  # Separate different types of checklists
  # --------------------------------
  
  # Stationary counts
  SC_sf <- subset(CL_sf, ProtocolType == "Stationary count")
  
  # Linear transects
  LT_sf <- subset(CL_sf, ProtocolType == "Linear transect")
  
  # --------------------------------
  # Convert to spatial objects
  # --------------------------------
  
  PC_sp <- as(PC_sf,'Spatial')
  SC_sp <- as(SC_sf,'Spatial')
  LT_sp <- as(LT_sf,'Spatial')
  
  # --------------------------------
  # Withhold 20% of point count data for cross-validation
  # --------------------------------
  
  PC_xval <- PC_sp[PC_sp$fold == fold,]
  PC_sp <- PC_sp[PC_sp$fold != fold,]
  
  # --------------------------------
  # Withhold checklist data in validation squares as well
  # --------------------------------
  
  SC_xval <- SC_sp[SC_sp$fold == fold,]
  SC_sp <- SC_sp[SC_sp$fold != fold,]
  
  LT_xval <- LT_sp[LT_sp$fold == fold,]
  LT_sp <- LT_sp[LT_sp$fold != fold,]
  
  # --------------------------------
  # Select a random 50% of remaining squares to withhold either PC or CL data during fitting
  # --------------------------------
  
  squares_to_withhold <- SaskSquares[SaskSquares$fold != fold,]$SQUARE_ID
  squares_to_withhold <- sample(squares_to_withhold,round(length(squares_to_withhold)/2))
  
  PC_sf <- st_as_sf(PC_sp)
  SC_sf <- st_as_sf(SC_sp)
  LT_sf <- st_as_sf(LT_sp)
  
  # --------------------------------
  # PLOT DATA WITH ONLY 50% OF POINT COUNTS
  # --------------------------------
  
  xval_fold = fold
  SaskSquares_centroids <- st_centroid(SaskSquares)
  
  map_val <- ggplot() +
    geom_sf(data = SaskBoundary,colour="black",fill="#f2f2f2",lwd=0.3,show.legend = F) +
    
    geom_sf(data = SaskWater,colour=NA,fill="#59F3F3",show.legend = F)+
    geom_sf(data = SaskBoundary,colour="black",fill=NA,lwd=0.3,show.legend = F)+
    
    # Plot validation locations
    geom_sf(data = st_as_sf(PC_xval),pch=19, size=0.1,show.legend = F, col = "blue")+
    
    coord_sf(clip = "off",xlim = c(min(SaskPoints$x), max(SaskPoints$x)))+
    theme(panel.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(5,10,5,-20),legend.title.align=0.5,
          legend.title = element_markdown(lineheight=.9,hjust = "left"))
  
  map_PC_50 <- ggplot() +
    geom_sf(data = SaskBoundary,colour="black",fill="#f2f2f2",lwd=0.3,show.legend = F) +
    
    geom_sf(data = SaskWater,colour=NA,fill="#59F3F3",show.legend = F)+
    geom_sf(data = SaskBoundary,colour="black",fill=NA,lwd=0.3,show.legend = F)+
    
    # Plot validation locations
    geom_sf(data = st_as_sf(PC_xval),pch=19, size=0.1,show.legend = F, col = "blue")+
    
    # Plot 50% of point count locations
    geom_sf(data = st_as_sf(PC_sp[PC_sp$sq_id %!in% squares_to_withhold,]),pch=19, size=0.2,show.legend = F, col = "black")+
    
    coord_sf(clip = "off",xlim = c(min(SaskPoints$x), max(SaskPoints$x)))+
    theme(panel.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(5,10,5,-20),legend.title.align=0.5,
          legend.title = element_markdown(lineheight=.9,hjust = "left"))
  
  map_integrated_50_50 <- ggplot() +
    geom_sf(data = SaskBoundary,colour="black",fill="#f2f2f2",lwd=0.3,show.legend = F) +
    
    geom_sf(data = SaskWater,colour=NA,fill="#59F3F3",show.legend = F)+
    geom_sf(data = SaskBoundary,colour="black",fill=NA,lwd=0.3,show.legend = F)+
    
    # Plot validation locations
    geom_sf(data = st_as_sf(PC_xval),pch=19, size=0.1,show.legend = F, col = "blue")+
    
    
    # Plot 50% of SC locations
    geom_sf(data = st_as_sf(SC_sp[SC_sp$sq_id %!in% squares_to_withhold,]),pch=5, size=1,show.legend = F, col = "red")+
    
    # Plot 50% of LT locations
    geom_sf(data = st_as_sf(LT_sp[LT_sp$sq_id %!in% squares_to_withhold,]),pch=5, size=1,show.legend = F, col = "red")+
    
    # Plot 50% of point count locations
    geom_sf(data = st_as_sf(PC_sp[PC_sp$sq_id %!in% squares_to_withhold,]),pch=19, size=0.2,show.legend = F, col = "black")+
    
    coord_sf(clip = "off",xlim = c(min(SaskPoints$x), max(SaskPoints$x)))+
    theme(panel.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(5,10,5,-20),legend.title.align=0.5,
          legend.title = element_markdown(lineheight=.9,hjust = "left"))
  
  map_PC_100 <- ggplot() +
    geom_sf(data = SaskBoundary,colour="black",fill="#f2f2f2",lwd=0.3,show.legend = F) +
    
    geom_sf(data = SaskWater,colour=NA,fill="#59F3F3",show.legend = F)+
    geom_sf(data = SaskBoundary,colour="black",fill=NA,lwd=0.3,show.legend = F)+
    
    # Plot validation locations
    geom_sf(data = st_as_sf(PC_xval),pch=19, size=0.1,show.legend = F, col = "blue")+
    
    # Plot 100% of point count locations
    geom_sf(data = st_as_sf(PC_sp),pch=19, size=0.2,show.legend = F, col = "black")+
    
    coord_sf(clip = "off",xlim = c(min(SaskPoints$x), max(SaskPoints$x)))+
    theme(panel.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(5,10,5,-20),legend.title.align=0.5,
          legend.title = element_markdown(lineheight=.9,hjust = "left"))
  
  map_integrated_50_100 <- ggplot() +
    geom_sf(data = SaskBoundary,colour="black",fill="#f2f2f2",lwd=0.3,show.legend = F) +
    
    geom_sf(data = SaskWater,colour=NA,fill="#59F3F3",show.legend = F)+
    geom_sf(data = SaskBoundary,colour="black",fill=NA,lwd=0.3,show.legend = F)+
    
    # Plot validation locations
    geom_sf(data = st_as_sf(PC_xval),pch=19, size=0.1,show.legend = F, col = "blue")+
    
    # Plot 100% of SC locations
    geom_sf(data = st_as_sf(SC_sp),pch=5, size=1,show.legend = F, col = "red")+
    
    # Plot 100% of LT locations
    geom_sf(data = st_as_sf(SC_sp),pch=5, size=1,show.legend = F, col = "red")+
    
    # Plot 50% of point count locations
    geom_sf(data = st_as_sf(PC_sp[PC_sp$sq_id %!in% squares_to_withhold,]),pch=19, size=0.2,show.legend = F, col = "black")+
    
    coord_sf(clip = "off",xlim = c(min(SaskPoints$x), max(SaskPoints$x)))+
    theme(panel.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(5,10,5,-20),legend.title.align=0.5,
          legend.title = element_markdown(lineheight=.9,hjust = "left"))
  
  # Output maps
  png(paste0("../output/figures/xval_maps/fold",xval_fold,"_val.png"), width=6.5, height=8, units="in", res=300, type="cairo")
  print(map_val)
  dev.off()
  
  png(paste0("../output/figures/xval_maps/fold",xval_fold,"_PC_50.png"), width=6.5, height=8, units="in", res=300, type="cairo")
  print(map_PC_50)
  dev.off()
  
  png(paste0("../output/figures/xval_maps/fold",xval_fold,"_PC_100.png"), width=6.5, height=8, units="in", res=300, type="cairo")
  print(map_PC_100)
  dev.off()
  
  png(paste0("../output/figures/xval_maps/fold",xval_fold,"_integrated_50_50.png"), width=6.5, height=8, units="in", res=300, type="cairo")
  print(map_integrated_50_50)
  dev.off()
  
  png(paste0("../output/figures/xval_maps/fold",xval_fold,"_integrated_50_100.png"), width=6.5, height=8, units="in", res=300, type="cairo")
  print(map_integrated_50_100)
  dev.off()
  
}

# ****************************************************
# ****************************************************
# Calculate/plot mean survey effort per square
# ****************************************************
# ****************************************************

# --------------------------------
# Load spatial data
# --------------------------------
PC_sf <- PC_surveyinfo
CL_sf <- DO_surveyinfo

# Intersect with SaskSquares dataframe
PC_sf <- st_intersection(PC_sf, SaskSquares)
CL_sf <- st_intersection(CL_sf, SaskSquares) %>%
  subset(ProtocolType %in% c("Stationary count","Linear transect"))

# --------------------------------
# Calculate mean survey effort per square
# --------------------------------
PC_effort <- PC_sf %>%
  as.data.frame() %>%
  group_by(sq_id) %>%
  summarize(PC_n = n(),
            PC_hours = sum(DurationInMinutes)/60)

CL_effort <- CL_sf %>%
  as.data.frame() %>%
  group_by(sq_id,ProtocolType) %>%
  summarize(CL_n = n(),
            CL_hours = sum(DurationInHours),
            CL_km = sum(TravelDistance_m)/1000)


mean(PC_effort$PC_n)


# ------------------------------------------------------------------------------------------
# Remaining questions to potentially explore:
# ------------------------------------------------------------------------------------------

# - is it better to double point count effort within each square, compared to adding overlapping checklists
#       - within a square, can checklists substitute for point counts?

# - what is the relative effort of point counts versus checklists within each square?
#       - multiple ways to assess 'checklist effort', not clear how to directly compare
#       - time spent surveying (total survey hours) is not directly comparable
#            - instead, compare average hours/survey, and also number of surveys

# - include offset for "survey duration", in addition to survey duration covariate (or does flexible survey duration effect already soak up all the variation)?
#      - see if this improves cross-validation scores?
#      