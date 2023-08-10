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

load("../output/xval_PC_CL_subset.RData")

# Only run species where all 5 folds have been completed
nfolds <- xval_PC_CL_subset %>%
  group_by(Species) %>%
  summarize(n = n()) %>%
  subset(n >=2)

# Arrange species by relative abundance
xval_PC_CL_subset <- xval_PC_CL_subset %>%
  subset(Species %in% nfolds$Species) %>%
  group_by(Species) %>%
  summarize_all(mean) %>% 
  arrange(n_det_PC)

xval_PC_CL_subset$Species <- factor(xval_PC_CL_subset$Species, levels = xval_PC_CL_subset$Species)

# -----------------------------------
# Summary of AUC scores
# -----------------------------------
sum(xval_PC_CL_subset$AUC_integrated_50_50 > xval_PC_CL_subset$AUC_PConly_50)
sum(xval_PC_CL_subset$AUC_PConly_100 > xval_PC_CL_subset$AUC_PConly_50)
sum(xval_PC_CL_subset$AUC_integrated_50_100 > xval_PC_CL_subset$AUC_PConly_50)

# -----------------------------------
# Summary of lppd scores
# -----------------------------------
sum(xval_PC_CL_subset$lppd_integrated_50_50 > xval_PC_CL_subset$lppd_PConly_50)
sum(xval_PC_CL_subset$lppd_PConly_100 > xval_PC_CL_subset$lppd_PConly_50)
sum(xval_PC_CL_subset$lppd_integrated_50_100 > xval_PC_CL_subset$lppd_PConly_50)

# -----------------------------------
# Summary of correlation scores
# -----------------------------------
sum(xval_PC_CL_subset$cor_integrated_50_50 > xval_PC_CL_subset$cor_PConly_50)
sum(xval_PC_CL_subset$cor_PConly_100 > xval_PC_CL_subset$cor_PConly_50)
sum(xval_PC_CL_subset$cor_integrated_50_100 > xval_PC_CL_subset$cor_PConly_50)

# -----------------------------------
# Summary of MSE scores
# -----------------------------------
sum(xval_PC_CL_subset$MSE_integrated_50_50 < xval_PC_CL_subset$MSE_PConly_50)
sum(xval_PC_CL_subset$MSE_PConly_100 < xval_PC_CL_subset$MSE_PConly_50)
sum(xval_PC_CL_subset$MSE_integrated_50_100 < xval_PC_CL_subset$MSE_PConly_50)

# -----------------------------------
# Is doubling checklist coverage better than doubling point counts?
# -----------------------------------

mean(xval_PC_CL_subset$cor_integrated_50_100 > xval_PC_CL_subset$cor_PConly_100)
mean(xval_PC_CL_subset$AUC_integrated_50_100 > xval_PC_CL_subset$AUC_PConly_100)
mean(xval_PC_CL_subset$MSE_integrated_50_100 < xval_PC_CL_subset$MSE_PConly_100)
mean(xval_PC_CL_subset$lppd_integrated_50_100 > xval_PC_CL_subset$lppd_PConly_100)

# **********************************************************
# **********************************************************
# GENERATE SOME PLOTS TO ILLUSTRATE DIFFERENCES AMONG SPECIES
# **********************************************************
# **********************************************************

# -----------------------------------------------------------------
# AUC
# -----------------------------------------------------------------

lim <- c(0.5,max(xval_PC_CL_subset[,c("AUC_PConly_50","AUC_PConly_100","AUC_integrated_50_50","AUC_integrated_50_100")]))

{ labels1 <- data.frame(Species = xval_PC_CL_subset$Species,
                        Improvement = (xval_PC_CL_subset$AUC_integrated_50_50 > xval_PC_CL_subset$AUC_PConly_50))
  labels1$sign <- "+"
  labels1$sign[labels1$Improvement == FALSE] <- ""
  
  labels2 <- data.frame(Species = xval_PC_CL_subset$Species,
                        Improvement = (xval_PC_CL_subset$AUC_PConly_100 > xval_PC_CL_subset$AUC_PConly_50))
  labels2$sign <- "+"
  labels2$sign[labels2$Improvement == FALSE] <- ""
  
  labels3 <- data.frame(Species = xval_PC_CL_subset$Species,
                        Improvement = (xval_PC_CL_subset$AUC_integrated_50_100 > xval_PC_CL_subset$AUC_PConly_50))
  labels3$sign <- "+"
  labels3$sign[labels3$Improvement == FALSE] <- ""
  
}


plot_AUC_0 <- ggplot()+
  
  # AUC between predictions and validation data
  geom_point(data = xval_PC_CL_subset, aes(y = Species,
                                           x = AUC_PConly_50),
             size = 2)+
  
  ylab("Species")+
  xlab("AUC between predictions and validation counts")+
  scale_x_continuous(limits = lim)+
  theme_few()+
  ggtitle("PC vs Integrated\n\nAUC with validation data")

plot_AUC_1 <- ggplot()+
  
  geom_text(data = labels1, aes(x = min(lim),
                               y = Species,
                               label = sign),
            fontface = "bold")+
  # Change in AUC with integrated model
  geom_segment(data = xval_PC_CL_subset, aes(y = Species,
                                             yend = Species, 
                                             x = AUC_PConly_50,
                                             xend = AUC_integrated_50_50),
               size = 1.5,
               lineend = "butt",
               linejoin = "mitre",
               col = "black",
               alpha = 1,
               arrow = arrow(length = unit(0.1, "inches")))+
  
  # AUC between predictions and validation data
  geom_point(data = xval_PC_CL_subset, aes(y = Species,
                                           x = AUC_PConly_50),
             size = 2)+
  
  ylab("Species")+
  xlab("AUC between predictions and validation counts")+
  scale_x_continuous(limits = lim)+
  theme_few()+
  ggtitle("PC vs Integrated\n\nAUC with validation data")

plot_AUC_2 <- ggplot()+
  
  geom_text(data = labels1, aes(x = min(lim),
                                y = Species,
                                label = sign),
            fontface = "bold")+
  
  geom_text(data = labels2, aes(x = min(lim)+0.01,
                                y = Species,
                                label = sign),
            fontface = "bold",
            col = "deepskyblue")+
  
  # Change in AUC when doubling point count coverage
  geom_segment(data = xval_PC_CL_subset, aes(y = Species,
                                             yend = Species,
                                             x = AUC_PConly_50,
                                             xend = AUC_PConly_100),
               size = 3,
               lineend = "butt",
               linejoin = "mitre",
               col = "deepskyblue",
               arrow = arrow(length = unit(0.1, "inches")))+
  
  # Change in AUC with integrated model
  geom_segment(data = xval_PC_CL_subset, aes(y = Species,
                                             yend = Species, 
                                             x = AUC_PConly_50,
                                             xend = AUC_integrated_50_50),
               size = 1.5,
               lineend = "butt",
               linejoin = "mitre",
               col = "black",
               alpha = 1,
               arrow = arrow(length = unit(0.1, "inches")))+
  
  # AUC between predictions and validation data
  geom_point(data = xval_PC_CL_subset, aes(y = Species,
                                           x = AUC_PConly_50),
             size = 2)+
  
  ylab("Species")+
  xlab("AUC between predictions and validation counts")+
  scale_x_continuous(limits = lim)+
  theme_few()+
  ggtitle("PC vs Integrated\n\nAUC with validation data")

plot_AUC_3 <- ggplot()+
  geom_text(data = labels1, aes(x = min(lim),
                                y = Species,
                                label = sign),
            fontface = "bold")+
  
  geom_text(data = labels2, aes(x = min(lim)+0.01,
                                y = Species,
                                label = sign),
            fontface = "bold",
            col = "deepskyblue")+
  
  geom_text(data = labels3, aes(x = min(lim)+0.02,
                                y = Species,
                                label = sign),
            fontface = "bold",
            col = "forestgreen")+
  
  # Change in AUC when doubling checklist coverage
  geom_segment(data = xval_PC_CL_subset, aes(y = Species,
                                             yend = Species,
                                             x = AUC_PConly_50,
                                             xend = AUC_integrated_50_100),
               size = 6,
               lineend = "butt",
               linejoin = "mitre",
               col = "forestgreen",
               arrow = arrow(length = unit(0.1, "inches")))+
  
  # Change in AUC when doubling point count coverage
  geom_segment(data = xval_PC_CL_subset, aes(y = Species,
                                             yend = Species,
                                             x = AUC_PConly_50,
                                             xend = AUC_PConly_100),
               size = 3,
               lineend = "butt",
               linejoin = "mitre",
               col = "deepskyblue",
               arrow = arrow(length = unit(0.1, "inches")))+
  
  # Change in AUC with integrated model
  geom_segment(data = xval_PC_CL_subset, aes(y = Species,
                                             yend = Species, 
                                             x = AUC_PConly_50,
                                             xend = AUC_integrated_50_50),
               size = 1.5,
               lineend = "butt",
               linejoin = "mitre",
               col = "black",
               alpha = 1,
               arrow = arrow(length = unit(0.1, "inches")))+
  
  # AUC between predictions and validation data
  geom_point(data = xval_PC_CL_subset, aes(y = Species,
                                           x = AUC_PConly_50),
             size = 2)+
  
  ylab("Species")+
  xlab("AUC between predictions and validation counts")+
  scale_x_continuous(limits = lim)+
  theme_few()+
  ggtitle("PC vs Integrated\n\nAUC with validation data")
plot_AUC_3

# Output figures
png("../output/figures/xval_plot_AUC_0.png", width=12, height=6, units="in", res=300, type="cairo")
print(plot_AUC_0)
dev.off()

png("../output/figures/xval_plot_AUC_1.png", width=12, height=6, units="in", res=300, type="cairo")
print(plot_AUC_1)
dev.off()

png("../output/figures/xval_plot_AUC_2.png", width=12, height=6, units="in", res=300, type="cairo")
print(plot_AUC_2)
dev.off()

png("../output/figures/xval_plot_AUC_3.png", width=12, height=6, units="in", res=300, type="cairo")
print(plot_AUC_3)
dev.off()

sum(xval_PC_CL_subset$AUC_integrated_50_50 > xval_PC_CL_subset$AUC_PConly_50)
sum(xval_PC_CL_subset$AUC_PConly_100 > xval_PC_CL_subset$AUC_PConly_50)
sum(xval_PC_CL_subset$AUC_integrated_50_100 > xval_PC_CL_subset$AUC_PConly_50)
nrow(xval_PC_CL_subset)

# -----------------------------------------------------------------
# Correlation
# -----------------------------------------------------------------

xval_PC_CL_subset$Species = factor(xval_PC_CL_subset$Species,levels = xval_PC_CL_subset$Species)

lim <- range(xval_PC_CL_subset[,c("cor_PConly_50","cor_PConly_100","cor_integrated_50_50","cor_integrated_50_100")])

{ labels1 <- data.frame(Species = xval_PC_CL_subset$Species,
                        Improvement = (xval_PC_CL_subset$cor_integrated_50_50 > xval_PC_CL_subset$cor_PConly_50))
  labels1$sign <- "+"
  labels1$sign[labels1$Improvement == FALSE] <- ""
  
  labels2 <- data.frame(Species = xval_PC_CL_subset$Species,
                        Improvement = (xval_PC_CL_subset$cor_PConly_100 > xval_PC_CL_subset$cor_PConly_50))
  labels2$sign <- "+"
  labels2$sign[labels2$Improvement == FALSE] <- ""
  
  labels3 <- data.frame(Species = xval_PC_CL_subset$Species,
                        Improvement = (xval_PC_CL_subset$cor_integrated_50_100 > xval_PC_CL_subset$cor_PConly_50))
  labels3$sign <- "+"
  labels3$sign[labels3$Improvement == FALSE] <- ""
  
}

plot_cor_1 <- ggplot()+
  
  # Change in correlation with integrated model
  geom_segment(data = xval_PC_CL_subset, aes(y = Species,
                                             yend = Species, 
                                             x = cor_PConly_50,
                                             xend = cor_integrated_50_50),
               size = 1.5,
               lineend = "butt",
               linejoin = "mitre",
               col = "black",
               alpha = 1,
               arrow = arrow(length = unit(0.1, "inches")))+
  
  # Correlation between predictions and validation data
  geom_point(data = xval_PC_CL_subset, aes(y = Species,
                                           x = cor_PConly_50),
             size = 2)+
  
  ylab("Species")+
  xlab("Correlation between predictions and validation counts")+
  scale_x_continuous(limits = lim)+
  theme_few()+
  ggtitle("PC vs Integrated\n\nCorrelation with validation data")

plot_cor_2 <- ggplot()+
  
  # Change in correlation when doubling point count coverage
  geom_segment(data = xval_PC_CL_subset, aes(y = Species,
                                             yend = Species,
                                             x = cor_PConly_50,
                                             xend = cor_PConly_100),
               size = 3,
               lineend = "butt",
               linejoin = "mitre",
               col = "deepskyblue",
               arrow = arrow(length = unit(0.1, "inches")))+
  
  # Change in correlation with integrated model
  geom_segment(data = xval_PC_CL_subset, aes(y = Species,
                                             yend = Species, 
                                             x = cor_PConly_50,
                                             xend = cor_integrated_50_50),
               size = 1.5,
               lineend = "butt",
               linejoin = "mitre",
               col = "black",
               alpha = 1,
               arrow = arrow(length = unit(0.1, "inches")))+
  
  # Correlation between predictions and validation data
  geom_point(data = xval_PC_CL_subset, aes(y = Species,
                                           x = cor_PConly_50),
             size = 2)+
  
  ylab("Species")+
  xlab("Correlation between predictions and validation counts")+
  scale_x_continuous(limits = lim)+
  theme_few()+
  ggtitle("PC vs Integrated\n\nCorrelation with validation data")

plot_cor_3 <- ggplot()+
  
  # Change in correlation when doubling checklist coverage
  geom_segment(data = xval_PC_CL_subset, aes(y = Species,
                                             yend = Species,
                                             x = cor_PConly_50,
                                             xend = cor_integrated_50_100),
               size = 6,
               lineend = "butt",
               linejoin = "mitre",
               col = "forestgreen",
               arrow = arrow(length = unit(0.1, "inches")))+
  
  # Change in correlation when doubling point count coverage
  geom_segment(data = xval_PC_CL_subset, aes(y = Species,
                                             yend = Species,
                                             x = cor_PConly_50,
                                             xend = cor_PConly_100),
               size = 3,
               lineend = "butt",
               linejoin = "mitre",
               col = "deepskyblue",
               arrow = arrow(length = unit(0.1, "inches")))+
  
  # Change in correlation with integrated model
  geom_segment(data = xval_PC_CL_subset, aes(y = Species,
                                             yend = Species, 
                                             x = cor_PConly_50,
                                             xend = cor_integrated_50_50),
               size = 1.5,
               lineend = "butt",
               linejoin = "mitre",
               col = "black",
               alpha = 1,
               arrow = arrow(length = unit(0.1, "inches")))+
  
  # Correlation between predictions and validation data
  geom_point(data = xval_PC_CL_subset, aes(y = Species,
                                           x = cor_PConly_50),
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
# xval_PC_CL_subset <- xval_PC_CL_subset %>%
#   mutate(RMSE_PConly_50 = sqrt(MSE_PConly_50),
#          RMSE_PConly_100 = sqrt(MSE_PConly_100),
#          RMSE_integrated_50_50 = sqrt(MSE_integrated_50_50),
#          RMSE_integrated_50_100 = sqrt(MSE_integrated_50_100))
# lim <- range(xval_PC_CL_subset[,c("RMSE_PConly_50","RMSE_PConly_100","RMSE_integrated_50_50","RMSE_integrated_50_100")])
# 
# plot_RMSE_1 <- ggplot()+
#   
#   # Change in RMSE with integrated model
#   geom_segment(data = xval_PC_CL_subset, aes(y = Species,
#                                              yend = Species, 
#                                              x = RMSE_PConly_50,
#                                              xend = RMSE_integrated_50_50),
#                size = 1.5,
#                lineend = "butt",
#                linejoin = "mitre",
#                col = "black",
#                alpha = 1,
#                arrow = arrow(length = unit(0.1, "inches")))+
#   
#   # RMSE between predictions and validation data
#   geom_point(data = xval_PC_CL_subset, aes(y = Species,
#                                            x = RMSE_PConly_50),
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
#   geom_segment(data = xval_PC_CL_subset, aes(y = Species,
#                                              yend = Species,
#                                              x = RMSE_PConly_50,
#                                              xend = RMSE_PConly_100),
#                size = 3,
#                lineend = "butt",
#                linejoin = "mitre",
#                col = "deepskyblue",
#                arrow = arrow(length = unit(0.1, "inches")))+
#   
#   # Change in RMSE with integrated model
#   geom_segment(data = xval_PC_CL_subset, aes(y = Species,
#                                              yend = Species, 
#                                              x = RMSE_PConly_50,
#                                              xend = RMSE_integrated_50_50),
#                size = 1.5,
#                lineend = "butt",
#                linejoin = "mitre",
#                col = "black",
#                alpha = 1,
#                arrow = arrow(length = unit(0.1, "inches")))+
#   
#   # RMSE between predictions and validation data
#   geom_point(data = xval_PC_CL_subset, aes(y = Species,
#                                            x = RMSE_PConly_50),
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
#   geom_segment(data = xval_PC_CL_subset, aes(y = Species,
#                                              yend = Species,
#                                              x = RMSE_PConly_50,
#                                              xend = RMSE_integrated_50_100),
#                size = 6,
#                lineend = "butt",
#                linejoin = "mitre",
#                col = "forestgreen",
#                arrow = arrow(length = unit(0.1, "inches")))+
#   
#   # Change in RMSE when doubling point count coverage
#   geom_segment(data = xval_PC_CL_subset, aes(y = Species,
#                                              yend = Species,
#                                              x = RMSE_PConly_50,
#                                              xend = RMSE_PConly_100),
#                size = 3,
#                lineend = "butt",
#                linejoin = "mitre",
#                col = "deepskyblue",
#                arrow = arrow(length = unit(0.1, "inches")))+
#   
#   # Change in RMSE with integrated model
#   geom_segment(data = xval_PC_CL_subset, aes(y = Species,
#                                              yend = Species, 
#                                              x = RMSE_PConly_50,
#                                              xend = RMSE_integrated_50_50),
#                size = 1.5,
#                lineend = "butt",
#                linejoin = "mitre",
#                col = "black",
#                alpha = 1,
#                arrow = arrow(length = unit(0.1, "inches")))+
#   
#   # RMSE between predictions and validation data
#   geom_point(data = xval_PC_CL_subset, aes(y = Species,
#                                            x = RMSE_PConly_50),
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

# # **********************************************************
# # **********************************************************
# # MAPS ILLUSTRATING CROSS-VALIDATION APPROACH
# # **********************************************************
# # **********************************************************
# library(sf)
# load("../output/AOS_data_package.RData")
# 
# for (fold in 1:5){
#   
#   # --------------------------------
#   # Load spatial data
#   # --------------------------------
#   PC_sf <- PC_surveyinfo 
#   CL_sf <- DO_surveyinfo 
#   
#   # Intersect with SaskSquares dataframe
#   PC_sf <- st_intersection(PC_sf, SaskSquares)
#   CL_sf <- st_intersection(CL_sf, SaskSquares)
#   
#   
#   # --------------------------------
#   # Separate different types of checklists
#   # --------------------------------
#   
#   # Stationary counts
#   SC_sf <- subset(CL_sf, ProtocolType == "Stationary count")
#   
#   # Linear transects
#   LT_sf <- subset(CL_sf, ProtocolType == "Linear transect")
#   
#   # --------------------------------
#   # Convert to spatial objects
#   # --------------------------------
#   
#   PC_sp <- as(PC_sf,'Spatial')
#   SC_sp <- as(SC_sf,'Spatial')
#   LT_sp <- as(LT_sf,'Spatial')
#   
#   # --------------------------------
#   # Withhold 20% of point count data for cross-validation
#   # --------------------------------
#   
#   PC_xval <- PC_sp[PC_sp$fold == fold,]
#   PC_sp <- PC_sp[PC_sp$fold != fold,]
#   
#   # --------------------------------
#   # Withhold checklist data in validation squares as well
#   # --------------------------------
#   
#   SC_xval <- SC_sp[SC_sp$fold == fold,]
#   SC_sp <- SC_sp[SC_sp$fold != fold,]
#   
#   LT_xval <- LT_sp[LT_sp$fold == fold,]
#   LT_sp <- LT_sp[LT_sp$fold != fold,]
#   
#   # --------------------------------
#   # Select a random 50% of remaining squares to withhold either PC or CL data during fitting
#   # --------------------------------
#   
#   set.seed(fold*100)
#   
#   squares_to_withhold <- SaskSquares[SaskSquares$fold != fold,]$SQUARE_ID
#   squares_to_withhold <- sample(squares_to_withhold,round(length(squares_to_withhold)/2))
#   
#   # --------------------------------
#   # Calculate number of point counts and checklists per non-validation square
#   # --------------------------------
#   
#   n_samples_PC <- PC_sp %>% as.data.frame() %>% group_by(sq_id) %>% summarize(n_PC = n())
#   n_samples_SC <- SC_sp %>% as.data.frame() %>% group_by(sq_id) %>% summarize(n_SC = n())
#   n_samples_LT <- LT_sp %>% as.data.frame() %>% group_by(sq_id) %>% summarize(n_LT = n())
#   
#   n_samples <- full_join(n_samples_PC,n_samples_SC) %>%
#     rowwise() %>%
#     mutate(n_min = min(c(n_PC,n_SC))) %>%
#     na.omit()
#   
#   
#   PC_sp <- PC_sp %>% 
#     st_as_sf() %>%
#     subset(sq_id %in% n_samples$sq_id) %>%
#     full_join(n_samples) %>%
#     group_by(sq_id) %>%
#     mutate(samp = sample(n())) %>%
#     filter(samp <= n_min) %>%
#     as('Spatial')
#   
#   SC_sp <- SC_sp %>% 
#     st_as_sf() %>%
#     subset(sq_id %in% n_samples$sq_id) %>%
#     full_join(n_samples) %>%
#     group_by(sq_id) %>%
#     mutate(samp = sample(n())) %>%
#     filter(samp <= n_min) %>%
#     as('Spatial')
#   
#   PC_sf <- st_as_sf(PC_sp)
#   SC_sf <- st_as_sf(SC_sp)
#   
#   # --------------------------------
#   # PLOT DATA WITH ONLY 50% OF POINT COUNTS
#   # --------------------------------
#   xval_fold = fold
#   SaskSquares_centroids <- st_centroid(SaskSquares)
#   
#   map_val <- ggplot() +
#     geom_sf(data = SaskBoundary,colour="black",fill="#f2f2f2",lwd=0.3,show.legend = F) +
#     
#     geom_sf(data = SaskWater,colour=NA,fill="#59F3F3",show.legend = F)+
#     geom_sf(data = SaskBoundary,colour="black",fill=NA,lwd=0.3,show.legend = F)+
#     
#     # Plot validation locations
#     geom_sf(data = st_as_sf(PC_xval),pch=19, size=0.1,show.legend = F, col = "blue")+
#     
#     coord_sf(clip = "off",xlim = c(min(SaskPoints$x), max(SaskPoints$x)))+
#     theme(panel.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#     theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
#     theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
#     theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
#     theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(5,10,5,-20),legend.title.align=0.5,
#           legend.title = element_markdown(lineheight=.9,hjust = "left"))
#   
#   map_PConly_50 <- ggplot() +
#     geom_sf(data = SaskBoundary,colour="black",fill="#f2f2f2",lwd=0.3,show.legend = F) +
#     
#     geom_sf(data = SaskWater,colour=NA,fill="#59F3F3",show.legend = F)+
#     geom_sf(data = SaskBoundary,colour="black",fill=NA,lwd=0.3,show.legend = F)+
#     
#     # Plot validation locations
#     geom_sf(data = st_as_sf(PC_xval),pch=19, size=0.1,show.legend = F, col = "blue")+
#     
#     # Plot 50% of point count locations
#     geom_sf(data = st_as_sf(PC_sp[PC_sp$sq_id %!in% squares_to_withhold,]),pch=19, size=0.2,show.legend = F, col = "black")+
#     
#     coord_sf(clip = "off",xlim = c(min(SaskPoints$x), max(SaskPoints$x)))+
#     theme(panel.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#     theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
#     theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
#     theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
#     theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(5,10,5,-20),legend.title.align=0.5,
#           legend.title = element_markdown(lineheight=.9,hjust = "left"))
#   
#   map_integrated_50_50 <- ggplot() +
#     geom_sf(data = SaskBoundary,colour="black",fill="#f2f2f2",lwd=0.3,show.legend = F) +
#     
#     geom_sf(data = SaskWater,colour=NA,fill="#59F3F3",show.legend = F)+
#     geom_sf(data = SaskBoundary,colour="black",fill=NA,lwd=0.3,show.legend = F)+
#     
#     # Plot validation locations
#     geom_sf(data = st_as_sf(PC_xval),pch=19, size=0.1,show.legend = F, col = "blue")+
#     
#     
#     # Plot 50% of SC locations
#     geom_sf(data = st_as_sf(SC_sp[SC_sp$sq_id %!in% squares_to_withhold,]),pch=5, size=1,show.legend = F, col = "red")+
#     
#     # Plot 50% of point count locations
#     geom_sf(data = st_as_sf(PC_sp[PC_sp$sq_id %!in% squares_to_withhold,]),pch=19, size=0.2,show.legend = F, col = "black")+
#     
#     coord_sf(clip = "off",xlim = c(min(SaskPoints$x), max(SaskPoints$x)))+
#     theme(panel.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#     theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
#     theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
#     theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
#     theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(5,10,5,-20),legend.title.align=0.5,
#           legend.title = element_markdown(lineheight=.9,hjust = "left"))
#   
#   map_PConly_100 <- ggplot() +
#     geom_sf(data = SaskBoundary,colour="black",fill="#f2f2f2",lwd=0.3,show.legend = F) +
#     
#     geom_sf(data = SaskWater,colour=NA,fill="#59F3F3",show.legend = F)+
#     geom_sf(data = SaskBoundary,colour="black",fill=NA,lwd=0.3,show.legend = F)+
#     
#     # Plot validation locations
#     geom_sf(data = st_as_sf(PC_xval),pch=19, size=0.1,show.legend = F, col = "blue")+
#     
#     # Plot 100% of point count locations
#     geom_sf(data = st_as_sf(PC_sp),pch=19, size=0.2,show.legend = F, col = "black")+
#     
#     coord_sf(clip = "off",xlim = c(min(SaskPoints$x), max(SaskPoints$x)))+
#     theme(panel.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#     theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
#     theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
#     theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
#     theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(5,10,5,-20),legend.title.align=0.5,
#           legend.title = element_markdown(lineheight=.9,hjust = "left"))
#   
#   map_integrated_50_100 <- ggplot() +
#     geom_sf(data = SaskBoundary,colour="black",fill="#f2f2f2",lwd=0.3,show.legend = F) +
#     
#     geom_sf(data = SaskWater,colour=NA,fill="#59F3F3",show.legend = F)+
#     geom_sf(data = SaskBoundary,colour="black",fill=NA,lwd=0.3,show.legend = F)+
#     
#     # Plot validation locations
#     geom_sf(data = st_as_sf(PC_xval),pch=19, size=0.1,show.legend = F, col = "blue")+
#     
#     
#     # Plot 100% of SC locations
#     geom_sf(data = st_as_sf(SC_sp),pch=5, size=1,show.legend = F, col = "red")+
#     
#     # Plot 50% of point count locations
#     geom_sf(data = st_as_sf(PC_sp[PC_sp$sq_id %!in% squares_to_withhold,]),pch=19, size=0.2,show.legend = F, col = "black")+
#     
#     coord_sf(clip = "off",xlim = c(min(SaskPoints$x), max(SaskPoints$x)))+
#     theme(panel.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#     theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
#     theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
#     theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
#     theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(5,10,5,-20),legend.title.align=0.5,
#           legend.title = element_markdown(lineheight=.9,hjust = "left"))
#   
#   # Output maps
#   png(paste0("../output/figures/xval_maps/fold",xval_fold,"_val.png"), width=6.5, height=8, units="in", res=300, type="cairo")
#   print(map_val)
#   dev.off()
#   
#   png(paste0("../output/figures/xval_maps/fold",xval_fold,"_PConly_50.png"), width=6.5, height=8, units="in", res=300, type="cairo")
#   print(map_PConly_50)
#   dev.off()
#   
#   png(paste0("../output/figures/xval_maps/fold",xval_fold,"_PConly_100.png"), width=6.5, height=8, units="in", res=300, type="cairo")
#   print(map_PConly_100)
#   dev.off()
#   
#   png(paste0("../output/figures/xval_maps/fold",xval_fold,"_integrated_50_50.png"), width=6.5, height=8, units="in", res=300, type="cairo")
#   print(map_integrated_50_50)
#   dev.off()
#   
#   png(paste0("../output/figures/xval_maps/fold",xval_fold,"_integrated_50_100.png"), width=6.5, height=8, units="in", res=300, type="cairo")
#   print(map_integrated_50_100)
#   dev.off()
#   
# }


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