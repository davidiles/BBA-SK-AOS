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
setwd("F:/Working_Files/1_Projects/Landbirds/SK_BBA_analysis/Standard_Analysis/")
`%!in%` <- Negate(`%in%`)

load("../AOS_precision/output/xval_df_50km_TSS.RData")

nfolds <- xval_df_50km_TSS %>%
  group_by(Species) %>%
  summarize(nfold = n()) %>%
  subset(nfold >= 4)

xval_df_50km_TSS <- xval_df_50km_TSS %>%
  group_by(Species) %>%
  summarize_all(mean) %>%
  subset(Species %in% nfolds$Species)

xval_df_50km_TSS$n_obs_CL <- xval_df_50km_TSS$n_obs_SC + xval_df_50km_TSS$n_obs_LT
xval_df_50km_TSS$delta_cor <- xval_df_50km_TSS$cor_integrated - xval_df_50km_TSS$cor_PConly
xval_df_50km_TSS$delta_MSE <- xval_df_50km_TSS$MSE_integrated - xval_df_50km_TSS$MSE_PConly
xval_df_50km_TSS$delta_AUC <- xval_df_50km_TSS$AUC_integrated - xval_df_50km_TSS$AUC_PConly

mean(xval_df_50km_TSS$delta_cor>0)
mean(xval_df_50km_TSS$delta_AUC>0)
mean(xval_df_50km_TSS$delta_MSE<0)

ggplot()+
  geom_text(data = xval_df_50km_TSS, aes(x = n_obs_PC, 
                                y = cor_integrated + sign(delta_cor)*0.01, label = Species,col = delta_cor > 0))+
  geom_segment(data = xval_df_50km_TSS, aes(x = n_obs_PC,xend = n_obs_PC, 
                                   y = (cor_PConly ),yend = cor_integrated, col = delta_cor > 0),
               size = 2,
               arrow = arrow(length = unit(0.05, "inches")))+
  scale_color_manual(values=c("red","blue"), name = "Change in cor",
                     labels = c("Decrease","Increase"), guide = "none")+
  scale_x_continuous(trans="log10", name = "Number of Detections in Point Counts")+
  ylab("Correlation\n(Crossvalidation)")+
  theme_bw()

ggplot()+
  geom_text(data = xval_df_50km_TSS, aes(x = n_obs_PC, 
                                y = AUC_integrated + sign(delta_AUC)*0.01, label = Species,col = delta_AUC > 0))+
  geom_segment(data = xval_df_50km_TSS, aes(x = n_obs_PC,xend = n_obs_PC, 
                                   y = (AUC_PConly ),yend = AUC_integrated, col = delta_AUC > 0),
               size = 2,
               arrow = arrow(length = unit(0.05, "inches")))+
  scale_color_manual(values=c("red","blue"), name = "Change in AUC",
                     labels = c("Decrease","Increase"))+
  scale_x_continuous(trans="log10", name = "Number of Detections in Point Counts")+
  ylab("AUC\n(Crossvalidation)")+
  theme_bw()

ggplot()+
  geom_text(data = xval_df_50km_TSS, aes(x = n_obs_PC, 
                                y = MSE_integrated + sign(delta_MSE)*MSE_integrated*0.2, label = Species,col = delta_MSE < 0))+
  geom_segment(data = xval_df_50km_TSS, aes(x = n_obs_PC,xend = n_obs_PC, 
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

load("../AOS_precision/output/xval_df_50km_TSS.RData")

xval_df_50km_TSS$n_obs_CL <- xval_df_50km_TSS$n_obs_SC + xval_df_50km_TSS$n_obs_LT + xval_df_50km_TSS$n_obs_BBA
xval_df_50km_TSS$delta_cor <- xval_df_50km_TSS$cor_integrated - xval_df_50km_TSS$cor_PConly
xval_df_50km_TSS$delta_MSE <- xval_df_50km_TSS$MSE_integrated - xval_df_50km_TSS$MSE_PConly
xval_df_50km_TSS$delta_AUC <- xval_df_50km_TSS$AUC_integrated - xval_df_50km_TSS$AUC_PConly

ggplot()+
  geom_segment(data = xval_df_50km_TSS, aes(x = xval_fold,xend = xval_fold,
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
  geom_segment(data = xval_df_50km_TSS, aes(x = xval_fold,xend = xval_fold, 
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
  geom_segment(data = xval_df_50km_TSS, aes(x = xval_fold,xend = xval_fold, 
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


setwd("D:/Working_Files/1_Projects/Landbirds/SK_BBA_analysis/Standard_Analysis/")
`%!in%` <- Negate(`%in%`)

# Load squares to withhold for crossvalidation
SaskSquares <- read_sf("../AOS_precision/output/SaskSquares_xval_20km.shp") %>%
  st_transform(crs(PC_surveyinfo))
SaskSquares$sq_idx <- as.numeric(factor(SaskSquares$SQUARE_ID))
SaskSquares_centroids <- st_centroid(SaskSquares)


load("!Data_processed/SaskGrid_PCA.RData")
SaskPoints <- st_centroid(SaskGrid) %>% dplyr::select(pointid)
SaskPoints$y <- st_coordinates(SaskPoints)[,2]
SaskPoints$x <- st_coordinates(SaskPoints)[,1]

SaskBoundary <- read_sf("!Shapefiles/SaskBoundary/SaskBoundary_Project.shp")

SaskWater <- read_sf("!Shapefiles/Water/SaskWaterClip.shp") %>%
  st_transform(crs = st_crs(SaskBoundary))
SaskWater$area <- st_area(SaskWater) %>% as.numeric()
SaskWater <- SaskWater[SaskWater$area>2.580e+06 ,]

# --------------------------------
# Load Point count data
# --------------------------------

# Pointcount_dataset: locations and covariates associated with each point count
load("D:/Working_Files/1_Projects/Landbirds/SK_BBA_analysis/AOS_precision/data/Pointcount_data.RData")

# Covariates / spatial locations of each point count survey
PC_surveyinfo <- Pointcount_data$PC_surveyinfo

# Counts associated with each point count survey
PC_matrix <- Pointcount_data$PC_matrix

# ******
# Selection criteria for point counts to use in analysis
# *********
PC_to_use <- which(PC_surveyinfo$DurationInMinutes %in% c(3,5,10)) # Time of day, day of year????
PC_surveyinfo <- PC_surveyinfo[PC_to_use,]
PC_matrix <- PC_matrix[PC_to_use,]

rm(Pointcount_data)

# --------------------------------
# Load Checklist data
# --------------------------------

# Pointcount_dataset: locations and covariates associated with each point count
load("D:/Working_Files/1_Projects/Landbirds/SK_BBA_analysis/AOS_precision/data/Checklist_data.RData")

# Covariates / spatial locations of each checklist survey (daily observations)
DO_surveyinfo <- Checklist_data$DO_surveyinfo

# Counts associated with each checklist survey
DO_matrix <- Checklist_data$DO_matrix

# --------------------------------
# Prepare to plot
# --------------------------------
PC_sf <- PC_surveyinfo
CL_sf <- DO_surveyinfo

# Stationary counts
SC_sf <- subset(CL_sf, ProtocolType == "Stationary count" & DurationInHours <= 1)

# Linear transects (travel distance included)
LT_sf <- subset(CL_sf, ProtocolType == "Linear transect" & 
                  TravelDistance_m <= 5000 & 
                  DurationInHours > 0 &
                  DurationInHours <= 1)

# Breeding Bird Atlas (travel distance mostly missing)
BBA_sf <- subset(CL_sf, ProtocolType == "Breeding Bird Atlas" &
                   DurationInHours > 0 &
                   DurationInHours <= 1)

CL_sf2 <- rbind(SC_sf,LT_sf,BBA_sf)


# ---------------------------------------------
# Summarize data on a square-by-square level
# ---------------------------------------------

SaskSquares <- SaskSquares %>% dplyr::rename(sq_id = SQUARE_ID)

# Summary of point count information in each square
PC_summary <- PC_sf %>%
  as.data.frame() %>%
  group_by(sq_id) %>%
  summarize(n_PC = n())

# Summary of checklist information in each square
CL_summary <- CL_sf2 %>%
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

png(paste0("../AOS_precision/output/figures/map_PC.png"), width=6.5, height=8, units="in", res=300, type="cairo")
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

png(paste0("../AOS_precision/output/figures/map_PC_SC.png"), width=6.5, height=8, units="in", res=300, type="cairo")
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

png(paste0("../AOS_precision/output/figures/map_PC_SC_LT.png"), width=6.5, height=8, units="in", res=300, type="cairo")
print(map_PC_SC_LT)
dev.off()