# ------------------------------------------------
# Load/install packages and set graphical themes / working directory
# ------------------------------------------------

my_packs <- c(
  
  # Data management
  'tidyverse',
  
  # Spatial functions
  #'sf','raster','ggspatial','rgeos','fasterize','exactextractr','circular',
  #'lubridate',
  
  # For PCA
  #'factoextra',
  
  # Spatial analysis
  #'inlabru','ebirdst',
  
  # For plotting
  'viridis','scales','ggpubr','ggtext','ggrepel')

if (any(!my_packs %in% installed.packages()[, 'Package'])) {install.packages(my_packs[which(!my_packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)}
lapply(my_packs, require, character.only = TRUE)

#library(INLA) # install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep = TRUE)
# library(napops) # For detectability offsets  # devtools::install_github("na-pops/napops") 
# napops::fetch_data()  # - only needs to be run once

rm(list=ls())

# Import rasters, data, and covariates from "standard analysis"
setwd("D:/Working_Files/1_Projects/Landbirds/SK_BBA_analysis/AOS_precision/script")
`%!in%` <- Negate(`%in%`)

# **********************************************************
# **********************************************************
# PART 1: ANALYSIS WHEN WITHHOLDING BOTH CHECKLISTS AND POINT COUNTS
# **********************************************************
# **********************************************************

load("../output/AOS_data_package.RData")
load("../output/xval_PC_CL.RData")
rm(mean)

nfolds <- xval_PC_CL %>%
  group_by(Species) %>%
  summarize(nfolds = n())

xval_PC_CL <- xval_PC_CL %>%
  group_by(Species) %>%
  summarise_all(.funs = mean,na.rm=T) %>%
  left_join(.,species_distribution_summary, by = c("Species"="sp_code")) %>%
  left_join(nfolds) %>%
  subset(nfolds > 1)

# -----------------------------------
# Cross-validation of POINT COUNT data
# -----------------------------------

xval_PC_CL$delta_cor_PC <- xval_PC_CL$cor_PC_integrated - xval_PC_CL$cor_PC_PConly

lim <- range(xval_PC_CL$n_PC,xval_PC_CL$n_sq)
plot1_a <- ggplot()+
  geom_point(data = xval_PC_CL, aes(x = n_PC, 
                                      y = cor_PC_PConly),
             size = 2)+
  geom_text(data = xval_PC_CL, aes(x = n_PC,
                                   y = cor_PC_PConly - sign(delta_cor_PC)*0.02,
                                   label = Species,
                                   hjust = 0.5),
            size = 4)+
  ylab("Correlation")+
  xlab("Number of atlas squares with detections")+
  scale_x_continuous(limits = lim)+
  scale_y_continuous(limits = c(0,0.5))+
  theme_bw()+
  ggtitle("PC vs Integrated\n\nCorrelation with validation data")

plot1_b <- ggplot()+
  geom_point(data = xval_PC_CL, aes(x = n_PC, 
                                    y = cor_PC_PConly),
             size = 2)+
  geom_segment(data = xval_PC_CL, aes(x = n_PC,xend = n_sq, 
                                                  y = cor_PC_PConly,
                                      yend = cor_PC_integrated, 
                                      col = delta_cor_PC > 0),
               size = 2,
               arrow = arrow(length = unit(0.05, "inches")))+
  geom_text(data = xval_PC_CL, aes(x = n_PC,
                                         y = cor_PC_PConly - sign(delta_cor_PC)*0.02, 
                                         col = delta_cor_PC > 0,
                                         label = Species,
                                         hjust = 0.5),
                  size = 4)+
  scale_color_manual(values=c("red","blue"), name = "Change in cor",
                     labels = c("Decrease","Increase"), guide = "none")+
  
  ylab("Correlation")+
  xlab("Number of atlas squares with detections")+
  scale_x_continuous(limits = lim)+
  scale_y_continuous(limits = c(0,0.5))+
  theme_bw()+
  ggtitle("PC vs Integrated\n\nCorrelation with validation data")

plot1_a
plot1_b

lim <- range(xval_PC_CL$n_PC,xval_PC_CL$n_sq)
xval_PC_CL$delta_AUC_PC <- xval_PC_CL$AUC_PC_integrated - xval_PC_CL$AUC_PC_PConly
plot2_a <- ggplot()+
  geom_point(data = xval_PC_CL, aes(x = n_PC, 
                                    y = AUC_PC_PConly),
             size = 2)+
  geom_text(data = xval_PC_CL, aes(x = n_PC,
                                   y = AUC_PC_PConly - sign(delta_AUC_PC)*0.02,
                                   label = Species,
                                   hjust = 0.5),
            size = 4)+
  ylab("AUC")+
  xlab("Number of atlas squares with detections")+
  scale_x_continuous(limits = lim)+
  scale_y_continuous(limits = c(0.5,1))+
  theme_bw()+
  ggtitle("PC vs Integrated\n\nAUC with validation data")

plot2_b <- ggplot()+
  geom_point(data = xval_PC_CL, aes(x = n_PC, 
                                    y = AUC_PC_PConly),
             size = 2)+
  geom_segment(data = xval_PC_CL, aes(x = n_PC,xend = n_sq, 
                                      y = AUC_PC_PConly,
                                      yend = AUC_PC_integrated, 
                                      col = delta_AUC_PC > 0),
               size = 2,
               arrow = arrow(length = unit(0.05, "inches")))+
  geom_text(data = xval_PC_CL, aes(x = n_PC,
                                   y = AUC_PC_PConly - sign(delta_AUC_PC)*0.02, 
                                   col = delta_AUC_PC > 0,
                                   label = Species,
                                   hjust = 0.5),
            size = 4)+
  scale_color_manual(values=c("red","blue"), name = "Change in AUC",
                     labels = c("Decrease","Increase"), guide = "none")+
  
  ylab("AUC")+
  xlab("Number of atlas squares with detections")+
  scale_x_continuous(limits = lim)+
  scale_y_continuous(limits = c(0.5,1))+
  theme_bw()+
  ggtitle("PC vs Integrated\n\nAUC with validation data")

#plot2_a
plot2_b


# 
# xval_PC_CL$delta_MSE_PC <- xval_PC_CL$MSE_PC_integrated - xval_PC_CL$MSE_PC_PConly
# lim <- range(xval_PC_CL[,c("MSE_PC_integrated","MSE_PC_PConly")])
# plot3_a <- ggplot()+
#   geom_point(data = xval_PC_CL, aes(x = n_PC, 
#                                     y = MSE_PC_PConly),
#              size = 2)+
#   geom_text(data = xval_PC_CL, aes(x = n_PC,
#                                    y = MSE_PC_PConly - sign(delta_MSE_PC)*0.02,
#                                    label = Species,
#                                    hjust = 0.5),
#             size = 4)+
#   ylab("MSE")+
#   xlab("Number of atlas squares with detections")+
#   scale_x_continuous(trans = "log10", limits = c(20,2000))+
#   scale_y_continuous(limits = lim, trans = "log10")+
#   theme_bw()+
#   ggtitle("PC vs Integrated\n\nMSE with validation data")
# 
# 
# plot3_b <- ggplot()+
#   geom_point(data = xval_PC_CL, aes(x = n_PC, 
#                                     y = MSE_PC_PConly),
#              size = 2)+
#   geom_segment(data = xval_PC_CL, aes(x = n_PC,xend = n_sq, 
#                                       y = MSE_PC_PConly,
#                                       yend = MSE_PC_integrated, 
#                                       col = delta_MSE_PC > 0),
#                size = 2,
#                arrow = arrow(length = unit(0.05, "inches")))+
#   geom_text(data = xval_PC_CL, aes(x = n_PC,
#                                    y = MSE_PC_PConly - sign(delta_MSE_PC)*0.02, 
#                                    col = delta_MSE_PC > 0,
#                                    label = Species,
#                                    hjust = 0.5),
#             size = 4)+
#   scale_color_manual(values=c("red","blue"), name = "Change in MSE",
#                      labels = c("Decrease","Increase"), guide = "none")+
#   
#   ylab("MSE")+
#   xlab("Number of atlas squares with detections")+
#   scale_x_continuous(trans = "log10", limits = c(20,2000))+
#   scale_y_continuous(limits = lim, trans = "log10")+
#   theme_bw()+
#   ggtitle("PC vs Integrated\n\nMSE with validation data")
# 
# plot3_a
# plot3_b
# 

# **********************************************************
# **********************************************************
# PART 2: ANALYSIS WHEN WITHHOLDING ONLY POINT COUNTS
# **********************************************************
# **********************************************************
rm(list=ls())

`%!in%` <- Negate(`%in%`)
load("../output/AOS_data_package.RData")
load("../output/xval_PC_CL.RData")
load("../output/xval_PC.RData")
rm(mean)

xval_PC <- xval_PC %>%
  dplyr::rename(cor_PC_integrated2 = cor_PC_integrated,
                AUC_PC_integrated2 = AUC_PC_integrated,
                MSE_PC_integrated2 = MSE_PC_integrated) %>%
  dplyr::select(Species,xval_fold,cor_PC_integrated2,AUC_PC_integrated2,MSE_PC_integrated2)

# Join
xval_PC_CL2 = full_join(xval_PC,xval_PC_CL)

nfolds <- xval_PC_CL2 %>%
  group_by(Species) %>%
  summarize(nfolds = sum(!is.na(cor_PC_integrated2)))

xval_PC_CL2 <- xval_PC_CL2 %>%
  group_by(Species) %>%
  summarise_all(.funs = mean,na.rm=T) %>%
  left_join(.,species_distribution_summary, by = c("Species"="sp_code")) %>%
  left_join(nfolds) %>%
  subset(nfolds >= 1)

# -----------------------------------
# Cross-validation of POINT COUNT data
# -----------------------------------

xval_PC_CL2$delta_cor_PC <- xval_PC_CL2$cor_PC_integrated2 - xval_PC_CL2$cor_PC_PConly

plot1_c <- ggplot()+
  geom_point(data = xval_PC_CL2, aes(x = n_PC, 
                                    y = cor_PC_PConly),
             size = 2)+
  geom_text(data = xval_PC_CL2, aes(x = n_PC,
                                   y = cor_PC_PConly - sign(delta_cor_PC)*0.02,
                                   label = Species,
                                   hjust = 0.5),
            size = 4)+
  ylab("Correlation")+
  xlab("Number of atlas squares with detections")+
  scale_x_continuous(trans = "log10", limits = c(20,2000))+
  scale_y_continuous(limits = c(0,0.5))+
  theme_bw()+
  ggtitle("PC vs Integrated\n\nCorrelation with validation data")

plot1_d <- ggplot()+
  geom_point(data = xval_PC_CL2, aes(x = n_PC, 
                                    y = cor_PC_PConly),
             size = 2)+
  geom_segment(data = xval_PC_CL2, aes(x = n_PC,xend = n_sq, 
                                      y = cor_PC_PConly,
                                      yend = cor_PC_integrated2, 
                                      col = delta_cor_PC > 0),
               size = 2,
               arrow = arrow(length = unit(0.05, "inches")))+
  geom_text(data = xval_PC_CL2, aes(x = n_PC,
                                   y = cor_PC_PConly - sign(delta_cor_PC)*0.02, 
                                   col = delta_cor_PC > 0,
                                   label = Species,
                                   hjust = 0.5),
            size = 4)+
  scale_color_manual(values=c("red","blue"), name = "Change in cor",
                     labels = c("Decrease","Increase"), guide = "none")+
  
  ylab("Correlation")+
  xlab("Number of atlas squares with detections")+
  scale_x_continuous(trans = "log10", limits = c(20,2000))+
  scale_y_continuous(limits = c(0,0.5))+
  theme_bw()+
  ggtitle("PC vs Integrated\n\nCorrelation with validation data")

plot1_c
plot1_d

xval_PC_CL2$delta_AUC_PC <- xval_PC_CL2$AUC_PC_integrated2 - xval_PC_CL2$AUC_PC_PConly
plot2_c <- ggplot()+
  geom_point(data = xval_PC_CL2, aes(x = n_PC, 
                                    y = AUC_PC_PConly),
             size = 2)+
  geom_text(data = xval_PC_CL2, aes(x = n_PC,
                                   y = AUC_PC_PConly - sign(delta_AUC_PC)*0.02,
                                   label = Species,
                                   hjust = 0.5),
            size = 4)+
  ylab("AUC")+
  xlab("Number of atlas squares with detections")+
  scale_x_continuous(trans = "log10", limits = c(20,2000))+
  scale_y_continuous(limits = c(0.5,1))+
  theme_bw()+
  ggtitle("PC vs Integrated\n\nAUC with validation data")

plot2_d <- ggplot()+
  geom_point(data = xval_PC_CL2, aes(x = n_PC, 
                                    y = AUC_PC_PConly),
             size = 2)+
  geom_segment(data = xval_PC_CL2, aes(x = n_PC,xend = n_sq, 
                                      y = AUC_PC_PConly,
                                      yend = AUC_PC_integrated2, 
                                      col = delta_AUC_PC > 0),
               size = 2,
               arrow = arrow(length = unit(0.05, "inches")))+
  geom_text(data = xval_PC_CL2, aes(x = n_PC,
                                   y = AUC_PC_PConly - sign(delta_AUC_PC)*0.02, 
                                   col = delta_AUC_PC > 0,
                                   label = Species,
                                   hjust = 0.5),
            size = 4)+
  scale_color_manual(values=c("red","blue"), name = "Change in AUC",
                     labels = c("Decrease","Increase"), guide = "none")+
  
  ylab("AUC")+
  xlab("Number of atlas squares with detections")+
  scale_x_continuous(trans = "log10", limits = c(20,2000))+
  scale_y_continuous(limits = c(0.5,1))+
  theme_bw()+
  ggtitle("PC vs Integrated\n\nAUC with validation data")

plot2_c
plot2_d

# -----------------------------------
# Cross-validation of POINT COUNT data
# -----------------------------------
ggplot(xval_PC_CL2)+
  geom_point(aes(x = cor_PC_integrated, y = cor_PC_integrated2))+
  geom_abline(intercept = 0, slope = 1)+
  theme_bw()

ggplot(xval_PC_CL2)+
  geom_point(aes(x = AUC_PC_integrated, y = AUC_PC_integrated2))+
  geom_abline(intercept = 0, slope = 1)+
  theme_bw()

# -----------------------------------
# Cross-validation of STATIONARY COUNT checklist data
# -----------------------------------
# 
# xval_PC_CL$delta_cor_SC <- xval_PC_CL$cor_SC_integrated - xval_PC_CL$cor_SC_CLonly
# lim <- max(abs(xval_PC_CL$delta_cor_SC))
# ggplot()+
#   geom_segment(data = xval_PC_CL, aes(x = 0,xend = delta_cor_SC, 
#                                                   y = Species,yend = Species, col = delta_cor_SC > 0),
#                size = 2,
#                arrow = arrow(length = unit(0.05, "inches")))+
#   scale_color_manual(values=c("red","blue"), name = "Change in cor",
#                      labels = c("Decrease","Increase"), guide = "none")+
#   
#   ylab("Species")+
#   xlab("Change in Correlation")+
#   scale_x_continuous(limits = c(-lim,lim))+
#   theme_bw()+
#   ggtitle("SC vs Integrated\n\nCorrelation with validation data")

xval_PC_CL$delta_AUC_SC <- xval_PC_CL$AUC_SC_integrated - xval_PC_CL$AUC_SC_CLonly
lim <- max(abs(xval_PC_CL$delta_AUC_SC))
ggplot()+
  geom_segment(data = xval_PC_CL, aes(x = 0,
                                      xend = delta_AUC_SC, 
                                                  y = Species,yend = Species, col = delta_AUC_SC > 0),
               size = 2,
               arrow = arrow(length = unit(0.05, "inches")))+
  scale_color_manual(values=c("red","blue"), name = "Change in AUC",
                     labels = c("Decrease","Increase"), guide = "none", drop = FALSE)+
  
  ylab("Species")+
  xlab("Change in AUC")+
  scale_x_continuous(limits = c(-lim,lim))+
  theme_bw()+
  ggtitle("SC vs Integrated\n\nAUC (presence/absence predictions)")+
  facet_grid()

# xval_PC_CL$delta_MSE_SC <- xval_PC_CL$MSE_SC_integrated - xval_PC_CL$MSE_SC_CLonly
# lim <- max(abs(xval_PC_CL$delta_MSE_SC))
# ggplot()+
#   geom_segment(data = xval_PC_CL, aes(x = 0,xend = delta_MSE_SC, 
#                                                   y = Species,yend = Species, col = delta_MSE_SC < 0),
#                size = 2,
#                arrow = arrow(length = unit(0.05, "inches")))+
#   scale_color_manual(values=c("red","blue"), name = "Change in MSE",
#                      labels = c("Decrease","Increase"), guide = "none", drop = FALSE)+
#   
#   ylab("Species")+
#   xlab("Change in MSE")+
#   scale_x_continuous(limits = c(-lim,lim))+
#   theme_bw()+
#   ggtitle("SC vs Integrated\n\nMSE with validation data")

# -----------------------------------
# Cross-validation of LINEAR TRANSECT checklist data
# -----------------------------------
# 
# xval_PC_CL$delta_cor_LT <- xval_PC_CL$cor_LT_integrated - xval_PC_CL$cor_LT_CLonly
# lim <- max(abs(xval_PC_CL$delta_cor_LT))
# ggplot()+
#   geom_segment(data = xval_PC_CL, aes(x = 0,xend = delta_cor_LT, 
#                                                   y = Species,yend = Species, col = delta_cor_LT > 0),
#                size = 2,
#                arrow = arrow(length = unit(0.05, "inches")))+
#   scale_color_manual(values=c("red","blue"), name = "Change in cor",
#                      labels = c("Decrease","Increase"), guide = "none")+
#   
#   ylab("Species")+
#   xlab("Change in Correlation")+
#   scale_x_continuous(limits = c(-lim,lim))+
#   theme_bw()+
#   ggtitle("LT vs Integrated\n\nCorrelation with validation data")

xval_PC_CL$delta_AUC_LT <- xval_PC_CL$AUC_LT_integrated - xval_PC_CL$AUC_LT_CLonly
lim <- max(abs(xval_PC_CL$delta_AUC_LT))
ggplot()+
  geom_segment(data = xval_PC_CL, aes(x = 0,xend = delta_AUC_LT, 
                                                  y = Species,yend = Species, col = delta_AUC_LT > 0),
               size = 2,
               arrow = arrow(length = unit(0.05, "inches")))+
  scale_color_manual(values=c("red","blue"), name = "Change in AUC",
                     labels = c("Decrease","Increase"), guide = "none", drop = FALSE)+
  
  ylab("Species")+
  xlab("Change in AUC")+
  scale_x_continuous(limits = c(-lim,lim))+
  theme_bw()+
  ggtitle("LT vs Integrated\n\nAUC (presence/absence predictions)")

# xval_PC_CL$delta_MSE_LT <- xval_PC_CL$MSE_LT_integrated - xval_PC_CL$MSE_LT_CLonly
# lim <- max(abs(xval_PC_CL$delta_MSE_LT))
# ggplot()+
#   geom_segment(data = xval_PC_CL, aes(x = 0,xend = delta_MSE_LT, 
#                                                   y = Species,yend = Species, col = delta_MSE_LT < 0),
#                size = 2,
#                arrow = arrow(length = unit(0.05, "inches")))+
#   scale_color_manual(values=c("red","blue"), name = "Change in MSE",
#                      labels = c("Decrease","Increase"), guide = "none", drop = FALSE)+
#   
#   ylab("Species")+
#   xlab("Change in MSE")+
#   scale_x_continuous(limits = c(-lim,lim))+
#   theme_bw()+
#   ggtitle("LT vs Integrated\n\nMSE with validation data")


# -----------------------------------------------------------
# Examine results across each cross-validation fold
# -----------------------------------------------------------

load("../AOS_precision/output/xval_PC_CL.RData")

xval_PC_CL$n_obs_CL <- xval_PC_CL$n_obs_SC + xval_PC_CL$n_obs_LT + xval_PC_CL$n_obs_BBA
xval_PC_CL$delta_cor <- xval_PC_CL$cor_integrated - xval_PC_CL$cor_PConly
xval_PC_CL$delta_MSE <- xval_PC_CL$MSE_integrated - xval_PC_CL$MSE_PConly
xval_PC_CL$delta_AUC <- xval_PC_CL$AUC_integrated - xval_PC_CL$AUC_PConly

ggplot()+
  geom_segment(data = xval_PC_CL, aes(x = xval_fold,xend = xval_fold,
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
  geom_segment(data = xval_PC_CL, aes(x = xval_fold,xend = xval_fold, 
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
  geom_segment(data = xval_PC_CL, aes(x = xval_fold,xend = xval_fold, 
                                   y = (MSE_PConly ),yend = MSE_integrated, col = delta_MSE < 0),
               size = 2,
               arrow = arrow(length = unit(0.05, "inches")))+
  scale_color_manual(values=c("red","blue"), name = "Change in MSE",
                     labels = c("Decrease","Increase"), guide = "none")+
  #scale_x_continuous(trans="log10", name = "Number of Detections in Point Counts")+
  ylab("MSE\n(Crossvalidation)")+
  theme_bw()+
  facet_wrap(Species~.)+
  scale_y_continuous(trans = "log10")

# *******************************************************************
# *******************************************************************
# Figures for AOS
# *******************************************************************
# *******************************************************************

library(sf)
load("../output/AOS_data_package.RData")

# ---------------------------------------------
# Summarize data on a square-by-square level
# ---------------------------------------------
PC_sf <- PC_surveyinfo %>% mutate(count = NA)
CL_sf <- DO_surveyinfo %>% mutate(count = NA)

SaskSquares <- SaskSquares %>% dplyr::rename(sq_id = SQUARE_ID)

# Summary of point count information in each square
PC_summary <- PC_sf %>%
  as.data.frame() %>%
  group_by(sq_id) %>%
  summarize(n_PC = n())

# Summary of checklist information in each square
CL_summary <- CL_sf %>%
  as.data.frame() %>%
  group_by(sq_id) %>%
  summarize(n_CL = n(),
            n_SC = sum(ProtocolType == "Stationary count"),
            n_LT = sum(ProtocolType == "Linear transect"),
            n_BBA = sum(ProtocolType == "Breeding Bird Atlas"))

SaskSquares_n <- SaskSquares %>%
  left_join(PC_summary) %>%
  left_join(CL_summary) %>%
  dplyr::select(sq_id,n_PC,n_CL,n_SC,n_LT,n_BBA,geometry)%>%
  rowwise() %>%
  mutate(PC_CL = sum(n_PC,n_CL,na.rm=TRUE))

SaskSquares_n_centroids <- st_centroid(SaskSquares_n)

map_PC <- ggplot() +
  geom_sf(data = SaskBoundary,colour="black",fill="#f2f2f2",lwd=0.3,show.legend = F) +

  geom_sf(data = SaskWater,colour=NA,fill="#59F3F3",show.legend = F)+
  geom_sf(data = SaskBoundary,colour="black",fill=NA,lwd=0.3,show.legend = F)+

  # Plot point count locations
  geom_sf(data = subset(SaskSquares_n_centroids,n_PC > 0 & !is.na(n_PC)),pch=19, size=0.1,show.legend = F, col = "black")+

  coord_sf(clip = "off",xlim = c(min(SaskPoints$x), max(SaskPoints$x)))+
  theme(panel.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
  theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(5,10,5,-20),legend.title.align=0.5,
        legend.title = element_markdown(lineheight=.9,hjust = "left"))

png(paste0("../output/figures/map_PC.png"), width=6.5, height=8, units="in", res=300, type="cairo")
print(map_PC)
dev.off()

map_PC_SC <- ggplot() +
  geom_sf(data = SaskBoundary,colour="black",fill="#f2f2f2",lwd=0.3,show.legend = F) +

  geom_sf(data = SaskWater,colour=NA,fill="#59F3F3",show.legend = F)+
  geom_sf(data = SaskBoundary,colour="black",fill=NA,lwd=0.3,show.legend = F)+

  # Plot point count locations
  geom_sf(data = subset(SaskSquares_n_centroids,n_PC > 0 & !is.na(n_PC)),pch=19, size=0.1,show.legend = F, col = "black")+

  # Plot black circles where SC checklists were collected
  geom_sf(data = subset(SaskSquares_n_centroids,n_SC > 0 & !is.na(n_SC)),pch=1, size=1,show.legend = F, col = "black")+

  coord_sf(clip = "off",xlim = c(min(SaskPoints$x), max(SaskPoints$x)))+
  theme(panel.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
  theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(5,10,5,-20),legend.title.align=0.5,
        legend.title = element_markdown(lineheight=.9,hjust = "left"))

png(paste0("../output/figures/map_PC_SC.png"), width=6.5, height=8, units="in", res=300, type="cairo")
print(map_PC_SC)
dev.off()

map_PC_SC_LT <- ggplot() +
  geom_sf(data = SaskBoundary,colour="black",fill="#f2f2f2",lwd=0.3,show.legend = F) +

  geom_sf(data = SaskWater,colour=NA,fill="#59F3F3",show.legend = F)+
  geom_sf(data = SaskBoundary,colour="black",fill=NA,lwd=0.3,show.legend = F)+

  # Plot point count locations
  geom_sf(data = subset(SaskSquares_n_centroids,n_PC > 0 & !is.na(n_PC)),pch=19, size=0.1,show.legend = F, col = "black")+

  # Plot black circles where SC checklists were collected
  geom_sf(data = subset(SaskSquares_n_centroids,n_SC > 0 & !is.na(n_SC)),pch=1, size=1,show.legend = F, col = "black")+

  # Plot red diamonds where LT checklists were collected
  geom_sf(data = subset(SaskSquares_n_centroids,n_LT > 0 & !is.na(n_LT)),pch=5, size=1,show.legend = F, col = "red")+

  coord_sf(clip = "off",xlim = c(min(SaskPoints$x), max(SaskPoints$x)))+
  theme(panel.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
  theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(5,10,5,-20),legend.title.align=0.5,
        legend.title = element_markdown(lineheight=.9,hjust = "left"))

png(paste0("../output/figures/map_PC_SC_LT.png"), width=6.5, height=8, units="in", res=300, type="cairo")
print(map_PC_SC_LT)
dev.off()


# --------------------------------
# --------------------------------
SaskSquares_n_centroids <- SaskSquares_n_centroids %>%
  left_join(.,as.data.frame(SaskSquares)[,c("sq_id","fold")])

# Plot of data locations 
map_xval_PC <- ggplot() +
  geom_sf(data = SaskBoundary,colour="black",fill="#f2f2f2",lwd=0.3,show.legend = F) +
  
  geom_sf(data = SaskWater,colour=NA,fill="#59F3F3",show.legend = F)+
  geom_sf(data = SaskBoundary,colour="black",fill=NA,lwd=0.3,show.legend = F)+
  
  # Plot point count locations
  geom_sf(data = subset(SaskSquares_n_centroids,n_PC > 0 & !is.na(n_PC) & fold != 1),pch=19, size=0.1,show.legend = F, col = "black")+
  
  geom_sf(data = subset(SaskSquares_n_centroids,n_PC > 0 & !is.na(n_PC) & fold == 1),pch=19, size=1,show.legend = F, col = "gray50")+
  
  # Plot black circles where SC checklists were collected
  #geom_sf(data = subset(SaskSquares_n_centroids,n_SC > 0 & !is.na(n_SC)),pch=1, size=1,show.legend = F, col = "black")+
  
  # Plot red diamonds where LT checklists were collected
  #geom_sf(data = subset(SaskSquares_n_centroids,n_LT > 0 & !is.na(n_LT)),pch=5, size=1,show.legend = F, col = "red")+
  
  coord_sf(clip = "off",xlim = c(min(SaskPoints$x), max(SaskPoints$x)))+
  theme(panel.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
  theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(5,10,5,-20),legend.title.align=0.5,
        legend.title = element_markdown(lineheight=.9,hjust = "left"))

png(paste0("../output/figures/map_xval_PC.png"), width=6.5, height=8, units="in", res=300, type="cairo")
print(map_xval_PC)
dev.off()

# Plot of data locations 
map_xval_PC_CL <- ggplot() +
  geom_sf(data = SaskBoundary,colour="black",fill="#f2f2f2",lwd=0.3,show.legend = F) +
  
  geom_sf(data = SaskWater,colour=NA,fill="#59F3F3",show.legend = F)+
  geom_sf(data = SaskBoundary,colour="black",fill=NA,lwd=0.3,show.legend = F)+
  
  # Plot point count locations
  geom_sf(data = subset(SaskSquares_n_centroids,n_PC > 0 & !is.na(n_PC) & fold != 1),pch=19, size=0.1,show.legend = F, col = "black")+
  
  geom_sf(data = subset(SaskSquares_n_centroids,n_PC > 0 & !is.na(n_PC) & fold == 1),pch=19, size=1,show.legend = F, col = "gray50")+
  
  # Plot black circles where SC checklists were collected
  geom_sf(data = subset(SaskSquares_n_centroids,n_SC > 0 & !is.na(n_SC) & fold != 1),pch=1, size=1,show.legend = F, col = "black")+
  
  # Plot red diamonds where LT checklists were collected
  geom_sf(data = subset(SaskSquares_n_centroids,n_LT > 0 & !is.na(n_LT) & fold != 1),pch=5, size=1,show.legend = F, col = "red")+
  
  coord_sf(clip = "off",xlim = c(min(SaskPoints$x), max(SaskPoints$x)))+
  theme(panel.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
  theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(5,10,5,-20),legend.title.align=0.5,
        legend.title = element_markdown(lineheight=.9,hjust = "left"))

png(paste0("../output/figures/map_xval_PC_CL.png"), width=6.5, height=8, units="in", res=300, type="cairo")
print(map_xval_PC_CL)
dev.off()