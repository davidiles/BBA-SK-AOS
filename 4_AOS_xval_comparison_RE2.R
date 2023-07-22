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
  'inlabru',
  
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

# **************************************************************************************
# Functions for plotting
# **************************************************************************************

# ------------------------------------------------
# Function to rasterize
# ------------------------------------------------

cut.fn <- function(df, column_name, lower_bound = NA, upper_bound = NA){
  
  max_val <- upper_bound
  max_val <- ifelse(is.na(max_val), 0, max_val)
  
  max_lev <- ifelse(max_val > 1.6, 4, 
                    ifelse(max_val > 0.8, 4, 3))
  
  cut_levs <- signif(max_val/(2^((max_lev-1):0)), 2)
  cut_levs <- unique(cut_levs)
  cut_levs <- ifelse(is.na(cut_levs), 0, cut_levs)
  
  if (lower_bound %in% cut_levs) cut_levs <- cut_levs[-which(cut_levs == lower_bound)]
  if (lower_bound > min(cut_levs)) cut_levs = cut_levs[-which(cut_levs < lower_bound)]
  
  max_lev <- length(cut_levs) ## do this because sometimes there are fewer levels
  
  cut_levs_labs <- c(paste0("0-",lower_bound),
                     paste(lower_bound, cut_levs[1], sep="-"),
                     paste(cut_levs[-max_lev], cut_levs[-1], sep="-"), 
                     paste(cut_levs[max_lev], "+"))
  
  cut_levs <- c(-1, lower_bound, cut_levs, 1000) %>% unique()
  
  # ** MULTIPLY BY 15 FOR COMPARISON TO BSC
  df$pc_cut <- cut(as.data.frame(df)[,column_name], cut_levs, labels=cut_levs_labs, ordered=TRUE)
  
  sp_abund_raster <- fasterize(df,SaskRaster,field = "pc_cut")
  sp_abund_raster <- mask(sp_abund_raster,SaskBoundary)
  
  sp_abund_raster_data <- as.data.frame(rasterToPoints(sp_abund_raster))
  sp_abund_raster_data$layer[sp_abund_raster_data$layer==1] <- cut_levs_labs[1]
  sp_abund_raster_data$layer[sp_abund_raster_data$layer==2] <- cut_levs_labs[2]
  sp_abund_raster_data$layer[sp_abund_raster_data$layer==3] <- cut_levs_labs[3]
  sp_abund_raster_data$layer[sp_abund_raster_data$layer==4] <- cut_levs_labs[4]
  sp_abund_raster_data$layer[sp_abund_raster_data$layer==5] <- cut_levs_labs[5]
  sp_abund_raster_data$layer[sp_abund_raster_data$layer==6] <- cut_levs_labs[6]
  
  sp_abund_raster_data$layer <- factor(sp_abund_raster_data$layer,
                                       levels= c(cut_levs_labs[1],cut_levs_labs[2],cut_levs_labs[3],
                                                 cut_levs_labs[4],cut_levs_labs[5],cut_levs_labs[6]),
                                       ordered = T)
  return(sp_abund_raster_data)
}

# ------------------------------------------------
# Colour scales
# ------------------------------------------------

colscale_density <- c("#FEFEFE", "#FFF4B3", "#F5D271", "#F2B647", "#EC8E00", "#CA302A")
colpal_density <- colorRampPalette(colscale_density)

colscale_relabund <-c("#FEFEFE", "#FBF7E2", "#FCF8D0", "#EEF7C2", "#CEF2B0", "#94E5A0", "#51C987", "#18A065", "#008C59", "#007F53", "#006344")
colpal_relabund <- colorRampPalette(colscale_relabund)


# ************************************************************************************
# ------------------------------------------------------------------------------------
# Load/prepare data
# ------------------------------------------------------------------------------------
# ************************************************************************************

# --------------------------------
# Load 1 km x 1 km covariate spatial grid (same grid onto which predictions will be made)
# --------------------------------

load("!Data_processed/SaskGrid_PCA.RData")
SaskPoints <- st_centroid(SaskGrid) %>% dplyr::select(pointid)
SaskPoints$y <- st_coordinates(SaskPoints)[,2]
SaskPoints$x <- st_coordinates(SaskPoints)[,1]

# Raster with target resolution
SaskRaster <- raster("!Shapefiles/RasterGrid/SaskRaster.tif")

# Study area boundary
SaskBoundary <- read_sf("!Shapefiles/SaskBoundary/SaskBoundary_Project.shp")

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


# ------------------------------------------------
# Process species names / labels (from Birds Canada)
# ------------------------------------------------

Sask_spcd <- read.csv("!Data/SaskAtlas_Species_2022-05-16.csv", stringsAsFactors = F)
countSpaces <- function(s) { sapply(gregexpr(" ", s), function(p) { sum(p>=0) } ) }

# Process Species Names so they fit
Sask_spcd$Label <- NA
Sask_spcd$CommonName[Sask_spcd$CommonName=="Rock Pigeon (Feral Pigeon)"] <- "Rock Pigeon"

for (i in 1:nrow(Sask_spcd)) {
  Name <- Sask_spcd$CommonName[i]
  if(nchar(Name) > 13){
    if(countSpaces(Name)>0){
      Sask_spcd$Label[i] <- gsub(" "," \n",Name)
    } 
    
  }
  else {
    Sask_spcd$Label[i] <- Sask_spcd$CommonName[i]
  }
  
}

# ------------------------------------------------
# List of species in the NA-POPS database (i.e., species with detectability offsets)
#    - species without detectability offsets will only have relative abundance maps (not actual density maps)
# ------------------------------------------------

napops_species <- napops::list_species()

# *************************************************************
# Calculate number of unique atlas squares in which each species was detected
# (in both point counts and checklists)
# *************************************************************

# Load squares to withhold for crossvalidation
SaskSquares <- read_sf("../AOS_precision/output/SaskSquares_xval.shp") %>%
  st_transform(crs(PC_surveyinfo))

SaskSquares$sq_idx <- as.numeric(factor(SaskSquares$SQUARE_ID))
SaskSquares_centroids <- st_centroid(SaskSquares)

# species_distribution_summary <- data.frame()
# 
# for (sp_code in colnames(PC_matrix)){
#   
#   if (sp_code %!in% colnames(DO_matrix)) next
#   
#   # point counts
#   PC_df <- data.frame(sp_code = sp_code,
#                       obs_index = PC_surveyinfo$obs_index,
#                       sq_id = PC_surveyinfo$sq_id,
#                       count = PC_matrix[,sp_code])
#   # checklists
#   CL_df <- data.frame(sp_code = sp_code,
#                       obs_index = DO_surveyinfo$obs_index,
#                       sq_id = DO_surveyinfo$sq_id,
#                       ProtocolType = DO_surveyinfo$ProtocolType,
#                       count = DO_matrix[,sp_code])
#   
#   # Calculate presence/absence summaries within each atlas square
#   sq_summary_PC <- PC_df %>%
#     group_by(sp_code,sq_id) %>%
#     summarize(n_surveys_PC = n(),
#               count_sum_PC = sum(count),
#               count_mean_PC = mean(count),
#               count_max_PC = max(count)) %>%
#     mutate(presence_PC = count_sum_PC > 0)
#   
#   sq_summary_CL <- CL_df %>%
#     group_by(sp_code,sq_id) %>%
#     summarize(n_surveys_CL = n(),
#               count_sum_CL = sum(count),
#               count_mean_CL = mean(count),
#               count_max_CL = max(count)) %>%
#     mutate(presence_CL = count_sum_CL > 0)
#   
#   # merge
#   sq_summary_sp <- full_join(sq_summary_PC,sq_summary_CL)
#   sq_summary_sp[is.na(sq_summary_sp)] <- 0
#   
#   species_distribution_summary <- rbind(species_distribution_summary,
#                                         data.frame(sp_code = sp_code,
#                                                    n_sq = sum((sq_summary_sp$presence_PC + sq_summary_sp$presence_CL)>0),
#                                                    n_PC = sum(sq_summary_sp$presence_PC > 0),
#                                                    n_CL = sum(sq_summary_sp$presence_CL > 0),
#                                                    n_PConly = sum(sq_summary_sp$presence_PC > 0 & sq_summary_sp$presence_CL == 0),
#                                                    n_CLonly = sum(sq_summary_sp$presence_CL > 0 & sq_summary_sp$presence_PC == 0),
#                                                    n_both = sum(sq_summary_sp$presence_PC > 0 & sq_summary_sp$presence_CL > 0)))
#   print(sp_code)
# }
# 
# # *************************************************************
# # Select species for analysis
# # *************************************************************
# 
# # Proportion of total that is CL only
# species_distribution_summary$prop_CLonly <- species_distribution_summary$n_CLonly/species_distribution_summary$n_sq
# 
# # Number of unique locations where species were detected
# species_relabund_PC <- colSums(PC_matrix>0,na.rm = TRUE) %>% sort(decreasing = TRUE) %>% as.data.frame() %>% rename(PC = 1)
# species_relabund_PC$Species <- rownames(species_relabund_PC)
# species_relabund_DO <- colSums(DO_matrix>0,na.rm = TRUE) %>% sort(decreasing = TRUE) %>% as.data.frame() %>% rename(DO = 1)
# species_relabund_DO$Species <- rownames(species_relabund_DO)
# species_relabund <- full_join(species_relabund_PC,species_relabund_DO)
# species_relabund <- na.omit(species_relabund)
# 
# species_relabund <- full_join(species_relabund,species_distribution_summary,
#                               by = c("Species" = "sp_code")) %>%
#   dplyr::relocate(Species)

#save(species_relabund,file="../AOS_precision/output/species_relabund.RData")
load(file="../AOS_precision/output/species_relabund.RData")

# Only consider species detected in at least 100 squares
species_to_fit <- subset(species_relabund, 
                         n_PC > 100 & n_CL > 100 & n_sq >=250 & n_sq <= 500) %>%
  arrange(desc(prop_CLonly))


ggplot(species_to_fit, aes(x = n_sq, y = prop_CLonly,label=Species, col = prop_CLonly))+
  geom_point()+
  geom_text_repel()+
  scale_color_gradientn(colors = viridis(10))+
  theme_bw()

# ************************************************************************************
# ------------------------------------------------------------------------------------
# Analysis (note some additional data wrangling is needed on a species-by-species basis)
# ------------------------------------------------------------------------------------
# ************************************************************************************

xval_df <- data.frame()
if (file.exists("../AOS_precision/output/xval_df_integrated_RE2.RData")){
  load("../AOS_precision/output/xval_df_integrated_RE2.RData")
}

