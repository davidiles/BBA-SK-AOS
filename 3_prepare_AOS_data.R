# ****************************************************************
# ****************************************************************
# CONDUCT CROSS-VALIDATION ANALYSIS ON A SUBSET OF SPECIES
#  - OMITTING BOTH POINT COUNTS AND CHECKLISTS FROM 20% OF SQUARES IN EACH CROSSVALIDATION FOLD
# ****************************************************************
# ****************************************************************

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
  
  # Cross-validation
  'pROC')

if (any(!my_packs %in% installed.packages()[, 'Package'])) {install.packages(my_packs[which(!my_packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)}
lapply(my_packs, require, character.only = TRUE)

library(INLA) # install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep = TRUE)

rm(list=ls())

# Import rasters, data, and covariates from "standard analysis"
setwd("D:/Working_Files/1_Projects/Landbirds/SK_BBA_analysis/AOS_precision/script")
`%!in%` <- Negate(`%in%`)

# ************************************************************************************
# ------------------------------------------------------------------------------------
# Load/prepare data
# ------------------------------------------------------------------------------------
# ************************************************************************************
target_crs <- "+proj=lcc +lat_0=48 +lon_0=-106 +lat_1=51 +lat_2=57 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"

# --------------------------------
# Load 1 km x 1 km covariate spatial grid (same grid onto which predictions will be made)
# --------------------------------

load("../../Standard_Analysis/!Data_processed/SaskGrid_PCA.RData") # SHOULD ALREADY BE IN TARGET CRS
SaskPoints <- st_centroid(SaskGrid) %>% dplyr::select(pointid)
SaskPoints$y <- st_coordinates(SaskPoints)[,2]
SaskPoints$x <- st_coordinates(SaskPoints)[,1]

# Raster with target resolution
SaskRaster <- raster("../../Standard_Analysis/!Shapefiles/RasterGrid/SaskRaster.tif")

# Study area boundary
SaskBoundary <- read_sf("../../Standard_Analysis/!Shapefiles/SaskBoundary/SaskBoundary_Project.shp")

# ------------------------------------------
# Process SaskWater raster
# ------------------------------------------
SaskWater <- read_sf("../../Standard_Analysis/!Shapefiles/Water/SaskWaterClip.shp") %>%
  st_transform(crs(target_crs))
SaskWater$area <- st_area(SaskWater)
SaskWater$area <- as.numeric(SaskWater$area)
SaskWater <- SaskWater[SaskWater$area>2.580e+06 ,]

Water_centroids <- SaskGrid %>%
  dplyr::select(pointid) %>%
  st_centroid() %>%
  st_intersection(.,SaskWater)

# *************************************************************
# Calculate number of unique atlas squares in which each species was detected
# (in both point counts and checklists)
# *************************************************************

# Load squares to withhold for crossvalidation
SaskSquares <- read_sf("../output/SaskSquares_xval_50km.shp") %>%
  st_transform(crs(target_crs))
SaskSquares$sq_idx <- as.numeric(factor(SaskSquares$SQUARE_ID))
SaskSquares_centroids <- st_centroid(SaskSquares)

# ******************************************************************
# ******************************************************************
# Selection criteria for point counts to use in analysis
# ******************************************************************
# ******************************************************************

# Pointcount_dataset: locations and covariates associated with each point count
load("../data/Pointcount_data_2.RData")

# Covariates / spatial locations of each point count survey
PC_surveyinfo <- Pointcount_data$PC_surveyinfo

# Counts associated with each point count survey
PC_matrix <- Pointcount_data$PC_matrix

# Acceptable days of year for surveys
start_date <- lubridate::ymd("2022-05-31") %>% yday()
end_date <- lubridate::ymd("2022-07-07") %>% yday()

PC_to_use <- subset(PC_surveyinfo,
                    DurationInMinutes %in% c(3,5,10) &
                     HSS >= -0.5 &
                     HSS <= 4 &
                     yday >= start_date &
                     yday <= end_date)

dim(PC_surveyinfo) # 19709 rows
dim(PC_to_use)     # 12058 rows

PC_surveyinfo <- PC_surveyinfo[PC_to_use$obs_index,]
PC_matrix <- PC_matrix[PC_to_use$obs_index,]

rm(Pointcount_data,PC_to_use)

# ******************************************************************
# ******************************************************************
# Selection criteria for CHECKLISTS to use in analysis
# ******************************************************************
# ******************************************************************

load("../data/Checklist_data_2.RData")

# Covariates / spatial locations of each checklist survey (daily observations)
DO_surveyinfo <- Checklist_data$DO_surveyinfo

start_date <- lubridate::ymd("2022-05-15") %>% yday()
end_date <- lubridate::ymd("2022-07-15") %>% yday()

# Select stationary counts to use
SC_to_use <- subset(DO_surveyinfo,
                    ProtocolType == "Stationary count" & 
                      HSS >= -1 &
                      HSS <= 12 &
                      yday >= start_date &
                      yday <= end_date &
                      DurationInHours <= 1)

# Select linear transects to use
LT_to_use <- subset(DO_surveyinfo,
                    ProtocolType == "Linear transect" & 
                      HSS >= -1 &
                      HSS <= 12 &
                      yday >= start_date &
                      yday <= end_date &
                      DurationInHours > (10/60) &
                      DurationInHours <= 2 &
                      TravelDistance_m <= 10000)

CL_to_use <- c(SC_to_use$obs_index,LT_to_use$obs_index)
dim(DO_surveyinfo) # 24672
DO_surveyinfo <- DO_surveyinfo[CL_to_use,]
dim(DO_surveyinfo) # 11616

# Counts associated with each checklist survey
DO_matrix <- Checklist_data$DO_matrix[CL_to_use,]

# --------------------------------
# Z-standardize covariates prior to analysis (to have mean = 0 and SD = 1)
#   (helps with model convergence, prior specification, etc)
# --------------------------------

# Covariates to be used in this analysis: PC1, PC2, PC3, Water_5km

# PC1
mean <- mean(PC_surveyinfo$PC1,na.rm = TRUE) ; sd <- sd(PC_surveyinfo$PC1,na.rm = TRUE)
PC_surveyinfo$PC1 <- (PC_surveyinfo$PC1 - mean)/sd
DO_surveyinfo$PC1 <- (DO_surveyinfo$PC1 - mean)/sd
SaskGrid$PC1 <- (SaskGrid$PC1 - mean)/sd

# PC2
mean <- mean(PC_surveyinfo$PC2,na.rm = TRUE) ; sd <- sd(PC_surveyinfo$PC2,na.rm = TRUE)
PC_surveyinfo$PC2 <- (PC_surveyinfo$PC2 - mean)/sd
DO_surveyinfo$PC2 <- (DO_surveyinfo$PC2 - mean)/sd
SaskGrid$PC2 <- (SaskGrid$PC2 - mean)/sd

# PC3
mean <- mean(PC_surveyinfo$PC3,na.rm = TRUE) ; sd <- sd(PC_surveyinfo$PC3,na.rm = TRUE)
PC_surveyinfo$PC3 <- (PC_surveyinfo$PC3 - mean)/sd
DO_surveyinfo$PC3 <- (DO_surveyinfo$PC3 - mean)/sd
SaskGrid$PC3 <- (SaskGrid$PC3 - mean)/sd

# Water_5km
mean <- mean(PC_surveyinfo$Water_5km,na.rm = TRUE) ; sd <- sd(PC_surveyinfo$Water_5km,na.rm = TRUE)
PC_surveyinfo$Water_5km <- (PC_surveyinfo$Water_5km - mean)/sd
DO_surveyinfo$Water_5km <- (DO_surveyinfo$Water_5km - mean)/sd
SaskGrid$Water_5km <- (SaskGrid$Water_5km - mean)/sd

# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# SELECT SPECIES FOR ANALYSIS
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------

# *************************************************************
# Calculate relative abundace (number of squares in which each species was observed)
# *************************************************************

species_distribution_summary <- data.frame()

for (sp_code in colnames(PC_matrix)){
  
  if (sp_code %!in% colnames(DO_matrix)) next
  
  # point counts
  PC_df <- data.frame(sp_code = sp_code,
                      obs_index = PC_surveyinfo$obs_index,
                      sq_id = PC_surveyinfo$sq_id,
                      count = PC_matrix[,sp_code])
  # checklists
  CL_df <- data.frame(sp_code = sp_code,
                      obs_index = DO_surveyinfo$obs_index,
                      sq_id = DO_surveyinfo$sq_id,
                      ProtocolType = DO_surveyinfo$ProtocolType,
                      count = DO_matrix[,sp_code])
  
  # Calculate presence/absence summaries within each atlas square
  sq_summary_PC <- PC_df %>%
    group_by(sp_code,sq_id) %>%
    summarize(n_surveys_PC = n(),
              count_sum_PC = sum(count),
              count_mean_PC = mean(count),
              count_max_PC = max(count)) %>%
    mutate(presence_PC = count_sum_PC > 0)
  
  sq_summary_CL <- CL_df %>%
    group_by(sp_code,sq_id) %>%
    summarize(n_surveys_CL = n(),
              count_sum_CL = sum(count),
              count_mean_CL = mean(count),
              count_max_CL = max(count)) %>%
    mutate(presence_CL = count_sum_CL > 0)
  
  # merge
  sq_summary_sp <- full_join(sq_summary_PC,sq_summary_CL)
  sq_summary_sp[is.na(sq_summary_sp)] <- 0
  
  species_distribution_summary <- rbind(species_distribution_summary,
                                        data.frame(sp_code = sp_code,
                                                   n_sq = sum((sq_summary_sp$presence_PC + sq_summary_sp$presence_CL)>0),
                                                   n_PC = sum(sq_summary_sp$presence_PC > 0),
                                                   n_CL = sum(sq_summary_sp$presence_CL > 0),
                                                   n_PConly = sum(sq_summary_sp$presence_PC > 0 & sq_summary_sp$presence_CL == 0),
                                                   n_CLonly = sum(sq_summary_sp$presence_CL > 0 & sq_summary_sp$presence_PC == 0),
                                                   n_both = sum(sq_summary_sp$presence_PC > 0 & sq_summary_sp$presence_CL > 0)))
}

# Proportion of total number of 'presence' squares that species was detected in point counts
species_distribution_summary$prop_PC <- species_distribution_summary$n_PC/species_distribution_summary$n_sq

# Restrict to species detected in at least 25 squares
species_distribution_summary <- subset(species_distribution_summary,
                                       n_PC >= 25 & n_CL >= 25)
ggplot(species_distribution_summary, aes(x = n_PC, y = n_CL,label=sp_code, col = prop_PC))+
  geom_point()+
  geom_text_repel()+
  scale_color_gradientn(colors = viridis(10))+
  scale_x_continuous(trans = "log10")+
  scale_y_continuous(trans = "log10")+
  theme_bw()

# *************************************************************
# Select species for analysis
# *************************************************************
set.seed(999)
species_to_fit <- sample_n(species_distribution_summary,20)

ggplot(species_to_fit, aes(x = n_PC, y = n_CL,label=sp_code, col = prop_PC))+
  geom_point()+
  geom_text_repel()+
  scale_color_gradientn(colors = viridis(10),
                        trans = "log10")+
  scale_x_continuous(trans = "log10")+
  scale_y_continuous(trans = "log10")+
  theme_bw()

save.image("../output/AOS_data_package.RData")
