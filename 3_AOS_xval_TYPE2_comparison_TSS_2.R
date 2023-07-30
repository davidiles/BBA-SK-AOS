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
# 
# 
# # **************************************************************************************
# # Functions for plotting
# # **************************************************************************************
# 
# # ------------------------------------------------
# # Function to rasterize
# # ------------------------------------------------
# 
# cut.fn <- function(df, column_name, lower_bound = NA, upper_bound = NA){
# 
#   max_val <- upper_bound
#   max_val <- ifelse(is.na(max_val), 0, max_val)
# 
#   max_lev <- ifelse(max_val > 1.6, 4,
#                     ifelse(max_val > 0.8, 4, 3))
# 
#   cut_levs <- signif(max_val/(2^((max_lev-1):0)), 2)
#   cut_levs <- unique(cut_levs)
#   cut_levs <- ifelse(is.na(cut_levs), 0, cut_levs)
# 
#   if (lower_bound %in% cut_levs) cut_levs <- cut_levs[-which(cut_levs == lower_bound)]
#   if (lower_bound > min(cut_levs)) cut_levs = cut_levs[-which(cut_levs < lower_bound)]
# 
#   max_lev <- length(cut_levs) ## do this because sometimes there are fewer levels
# 
#   cut_levs_labs <- c(paste0("0-",lower_bound),
#                      paste(lower_bound, cut_levs[1], sep="-"),
#                      paste(cut_levs[-max_lev], cut_levs[-1], sep="-"),
#                      paste(cut_levs[max_lev], "+"))
# 
#   cut_levs <- c(-1, lower_bound, cut_levs, 1000) %>% unique()
# 
#   # ** MULTIPLY BY 15 FOR COMPARISON TO BSC
#   df$pc_cut <- cut(as.data.frame(df)[,column_name], cut_levs, labels=cut_levs_labs, ordered=TRUE)
# 
#   sp_abund_raster <- fasterize(df,SaskRaster,field = "pc_cut")
#   sp_abund_raster <- mask(sp_abund_raster,SaskBoundary)
# 
#   sp_abund_raster_data <- as.data.frame(rasterToPoints(sp_abund_raster))
#   sp_abund_raster_data$layer[sp_abund_raster_data$layer==1] <- cut_levs_labs[1]
#   sp_abund_raster_data$layer[sp_abund_raster_data$layer==2] <- cut_levs_labs[2]
#   sp_abund_raster_data$layer[sp_abund_raster_data$layer==3] <- cut_levs_labs[3]
#   sp_abund_raster_data$layer[sp_abund_raster_data$layer==4] <- cut_levs_labs[4]
#   sp_abund_raster_data$layer[sp_abund_raster_data$layer==5] <- cut_levs_labs[5]
#   sp_abund_raster_data$layer[sp_abund_raster_data$layer==6] <- cut_levs_labs[6]
# 
#   sp_abund_raster_data$layer <- factor(sp_abund_raster_data$layer,
#                                        levels= c(cut_levs_labs[1],cut_levs_labs[2],cut_levs_labs[3],
#                                                  cut_levs_labs[4],cut_levs_labs[5],cut_levs_labs[6]),
#                                        ordered = T)
#   return(sp_abund_raster_data)
# }
# 
# # Custom cut points
# cut.fn2 <- function(df, column_name, cut_levs = cut_levs){
# 
# 
#   lower_bound <- min(cut_levs)
#   upper_bound <- max(cut_levs)
# 
#   max_lev <- length(cut_levs) ## do this because sometimes there are fewer levels
# 
#   cut_levs_labs <- c(paste0("0-",cut_levs[1]),
#                      paste(cut_levs[-max_lev], cut_levs[-1], sep="-"),
#                      paste(cut_levs[max_lev], "+"))
# 
#   cut_levs <- c(-1, cut_levs, 1000) %>% unique()
# 
#   # ** MULTIPLY BY 15 FOR COMPARISON TO BSC
#   df$pc_cut <- cut(as.data.frame(df)[,column_name],
#                    cut_levs,
#                    labels=cut_levs_labs,
#                    ordered=TRUE)
# 
#   sp_abund_raster <- fasterize(df,SaskRaster,field = "pc_cut")
#   sp_abund_raster <- mask(sp_abund_raster,SaskBoundary)
# 
#   sp_abund_raster_data <- as.data.frame(rasterToPoints(sp_abund_raster))
#   for (i in unique(sp_abund_raster_data$layer)) sp_abund_raster_data$layer[sp_abund_raster_data$layer==i] <- cut_levs_labs[i]
#   # sp_abund_raster_data$layer[sp_abund_raster_data$layer==1] <- cut_levs_labs[1]
#   # sp_abund_raster_data$layer[sp_abund_raster_data$layer==2] <- cut_levs_labs[2]
#   # sp_abund_raster_data$layer[sp_abund_raster_data$layer==3] <- cut_levs_labs[3]
#   # sp_abund_raster_data$layer[sp_abund_raster_data$layer==4] <- cut_levs_labs[4]
#   # sp_abund_raster_data$layer[sp_abund_raster_data$layer==5] <- cut_levs_labs[5]
#   # sp_abund_raster_data$layer[sp_abund_raster_data$layer==6] <- cut_levs_labs[6]
#   #
#   sp_abund_raster_data$layer <- factor(sp_abund_raster_data$layer,
#                                        levels= cut_levs_labs,
#                                        ordered = T)
#   return(sp_abund_raster_data)
# }
# 
# # ------------------------------------------------
# # Colour scales
# # ------------------------------------------------
# 
# colscale_density <- c("#FEFEFE", "#FFF4B3", "#F5D271", "#F2B647", "#EC8E00", "#CA302A")
# colpal_density <- colorRampPalette(colscale_density)
# 
# colscale_relabund <-c("#FEFEFE", "#FBF7E2", "#FCF8D0", "#EEF7C2", "#CEF2B0", "#94E5A0", "#51C987", "#18A065", "#008C59", "#007F53", "#006344")
# colpal_relabund <- colorRampPalette(colscale_relabund)
# 
# 
# # ************************************************************************************
# # ------------------------------------------------------------------------------------
# # Load/prepare data
# # ------------------------------------------------------------------------------------
# # ************************************************************************************
# 
# # --------------------------------
# # Load 1 km x 1 km covariate spatial grid (same grid onto which predictions will be made)
# # --------------------------------
# 
# load("!Data_processed/SaskGrid_PCA.RData")
# SaskPoints <- st_centroid(SaskGrid) %>% dplyr::select(pointid)
# SaskPoints$y <- st_coordinates(SaskPoints)[,2]
# SaskPoints$x <- st_coordinates(SaskPoints)[,1]
# 
# # Raster with target resolution
# SaskRaster <- raster("!Shapefiles/RasterGrid/SaskRaster.tif")
# 
# # Study area boundary
# SaskBoundary <- read_sf("!Shapefiles/SaskBoundary/SaskBoundary_Project.shp")
# 
# # --------------------------------
# # Load Point count data
# # --------------------------------
# 
# # Pointcount_dataset: locations and covariates associated with each point count
# load("D:/Working_Files/1_Projects/Landbirds/SK_BBA_analysis/AOS_precision/data/Pointcount_data_2.RData")
# 
# # Covariates / spatial locations of each point count survey
# PC_surveyinfo <- Pointcount_data$PC_surveyinfo
# 
# # Counts associated with each point count survey
# PC_matrix <- Pointcount_data$PC_matrix
# 
# # ******
# # Selection criteria for point counts to use in analysis
# # *********
# 
# # Acceptable days of year for surveys
# start_date <- lubridate::ymd("2022-05-31") %>% yday()
# end_date <- lubridate::ymd("2022-07-07") %>% yday()
# 
# PC_to_use <- which(PC_surveyinfo$DurationInMinutes %in% c(3,5,10) &
#                      PC_surveyinfo$HSS >= -0.5 &
#                      PC_surveyinfo$HSS <= 4 &
#                      PC_surveyinfo$yday >= start_date & 
#                      PC_surveyinfo$yday <= end_date) # Time of day, day of year????
# PC_surveyinfo <- PC_surveyinfo[PC_to_use,]
# PC_matrix <- PC_matrix[PC_to_use,]
# 
# rm(Pointcount_data)
# 
# # --------------------------------
# # Load Checklist data
# # --------------------------------
# 
# load("D:/Working_Files/1_Projects/Landbirds/SK_BBA_analysis/AOS_precision/data/Checklist_data_2.RData")
# 
# # Covariates / spatial locations of each checklist survey (daily observations)
# DO_surveyinfo <- Checklist_data$DO_surveyinfo
# 
# # ******
# # Selection criteria for CHECKLISTS to use in analysis
# # *********
# start_date <- lubridate::ymd("2022-05-15") %>% yday()
# end_date <- lubridate::ymd("2022-07-15") %>% yday()
# 
# DO_to_use <- which(DO_surveyinfo$HSS >= -1 &
#                      DO_surveyinfo$HSS <= 12 &
#                      DO_surveyinfo$yday >= start_date &
#                      DO_surveyinfo$yday <= end_date) # Time of day, day of year????
# DO_surveyinfo <- DO_surveyinfo[DO_to_use,]
# 
# # Counts associated with each checklist survey
# DO_matrix <- Checklist_data$DO_matrix
# DO_matrix <- DO_matrix[DO_to_use,]
# 
# # --------------------------------
# # Z-standardize covariates prior to analysis (to have mean = 0 and SD = 1)
# #   (helps with model convergence, prior specification, etc)
# # --------------------------------
# 
# # Covariates to be used in this analysis: PC1, PC2, PC3, Water_5km
# 
# # PC1
# mean <- mean(PC_surveyinfo$PC1,na.rm = TRUE) ; sd <- sd(PC_surveyinfo$PC1,na.rm = TRUE)
# PC_surveyinfo$PC1 <- (PC_surveyinfo$PC1 - mean)/sd
# DO_surveyinfo$PC1 <- (DO_surveyinfo$PC1 - mean)/sd
# SaskGrid$PC1 <- (SaskGrid$PC1 - mean)/sd
# 
# # PC2
# mean <- mean(PC_surveyinfo$PC2,na.rm = TRUE) ; sd <- sd(PC_surveyinfo$PC2,na.rm = TRUE)
# PC_surveyinfo$PC2 <- (PC_surveyinfo$PC2 - mean)/sd
# DO_surveyinfo$PC2 <- (DO_surveyinfo$PC2 - mean)/sd
# SaskGrid$PC2 <- (SaskGrid$PC2 - mean)/sd
# 
# # PC3
# mean <- mean(PC_surveyinfo$PC3,na.rm = TRUE) ; sd <- sd(PC_surveyinfo$PC3,na.rm = TRUE)
# PC_surveyinfo$PC3 <- (PC_surveyinfo$PC3 - mean)/sd
# DO_surveyinfo$PC3 <- (DO_surveyinfo$PC3 - mean)/sd
# SaskGrid$PC3 <- (SaskGrid$PC3 - mean)/sd
# 
# # Water_5km
# mean <- mean(PC_surveyinfo$Water_5km,na.rm = TRUE) ; sd <- sd(PC_surveyinfo$Water_5km,na.rm = TRUE)
# PC_surveyinfo$Water_5km <- (PC_surveyinfo$Water_5km - mean)/sd
# DO_surveyinfo$Water_5km <- (DO_surveyinfo$Water_5km - mean)/sd
# SaskGrid$Water_5km <- (SaskGrid$Water_5km - mean)/sd
# 
# # ------------------------------------------------
# # Process species names / labels (from Birds Canada)
# # ------------------------------------------------
# 
# Sask_spcd <- read.csv("!Data/SaskAtlas_Species_2022-05-16.csv", stringsAsFactors = F)
# countSpaces <- function(s) { sapply(gregexpr(" ", s), function(p) { sum(p>=0) } ) }
# 
# # Process Species Names so they fit
# Sask_spcd$Label <- NA
# Sask_spcd$CommonName[Sask_spcd$CommonName=="Rock Pigeon (Feral Pigeon)"] <- "Rock Pigeon"
# 
# for (i in 1:nrow(Sask_spcd)) {
#   Name <- Sask_spcd$CommonName[i]
#   if(nchar(Name) > 13){
#     if(countSpaces(Name)>0){
#       Sask_spcd$Label[i] <- gsub(" "," \n",Name)
#     }
# 
#   }
#   else {
#     Sask_spcd$Label[i] <- Sask_spcd$CommonName[i]
#   }
# 
# }
# 
# # ------------------------------------------------
# # List of species in the NA-POPS database (i.e., species with detectability offsets)
# #    - species without detectability offsets will only have relative abundance maps (not actual density maps)
# # ------------------------------------------------
# 
# # napops_species <- napops::list_species()
# 
# # *************************************************************
# # Calculate number of unique atlas squares in which each species was detected
# # (in both point counts and checklists)
# # *************************************************************
# 
# # Load squares to withhold for crossvalidation
# SaskSquares <- read_sf("../AOS_precision/output/SaskSquares_xval_50km.shp") %>%
#   st_transform(crs(PC_surveyinfo))
# SaskSquares$sq_idx <- as.numeric(factor(SaskSquares$SQUARE_ID))
# SaskSquares_centroids <- st_centroid(SaskSquares)
# 
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
# 
# save(species_relabund,file="../AOS_precision/output/species_relabund.RData")
# load(file="../AOS_precision/output/species_relabund.RData")
# 
# # Only consider species detected in at least 100 squares
# species_to_fit <- subset(species_relabund,
#                          n_PC > 100 & n_sq >=100) %>%
#   arrange(desc(prop_CLonly))
# 
# species_to_fit <- sample_n(species_to_fit,20)
# 
# ggplot(species_to_fit, aes(x = n_PC, y = n_CL,label=Species, col = prop_CLonly))+
#   geom_point()+
#   geom_text_repel()+
#   scale_color_gradientn(colors = viridis(10))+
#   theme_bw()
# 
# # ************************************************************************************
# # ------------------------------------------------------------------------------------
# # Analysis (note some additional data wrangling is needed on a species-by-species basis)
# # ------------------------------------------------------------------------------------
# # ************************************************************************************
# SaskWater <- read_sf("!Shapefiles/Water/SaskWaterClip.shp")
# SaskWater <- st_transform(SaskWater, crs = st_crs(SaskBoundary))
# SaskWater$area <- st_area(SaskWater)
# SaskWater$area <- as.numeric(SaskWater$area)
# SaskWater <- SaskWater[SaskWater$area>2.580e+06 ,]
# 
# Water_centroids <- SaskGrid %>%
#   dplyr::select(pointid) %>%
#   st_centroid() %>%
#   st_intersection(.,SaskWater)
# 
# save.image("D:/Working_Files/1_Projects/Landbirds/SK_BBA_analysis/AOS_precision/output/wksp_50km_2.RData")

load("D:/Working_Files/1_Projects/Landbirds/SK_BBA_analysis/AOS_precision/output/wksp_50km_2.RData")

xval_TYPE2_df_50km_TSS <- data.frame()
if (file.exists("../AOS_precision/output/xval_TYPE2_df_50km_TSS_2.RData")){
  load("../AOS_precision/output/xval_TYPE2_df_50km_TSS_2.RData")
}

# covariates to include in models
covariates_to_include <- c("PC1","PC2","PC3","Water_5km")

for (fold in sort(unique(SaskSquares$fold))){
  for (sp_code in rev(species_to_fit$Species)){
    
    if (file.exists("../AOS_precision/output/xval_TYPE2_df_50km_TSS_2.RData")){
      load("../AOS_precision/output/xval_TYPE2_df_50km_TSS_2.RData")
    }
    
    # Check if this species/xval fold have already been run. If so, skip
    if (nrow(xval_TYPE2_df_50km_TSS)>0){
      if (nrow(subset(xval_TYPE2_df_50km_TSS,Species == sp_code & xval_fold == fold))>0) next
    }
    
    print(paste(sp_code," fold ",fold))
    
    # --------------------------------
    # Prepare point count data for this species
    # --------------------------------
    
    PC_sf <- PC_surveyinfo %>% mutate(count = NA)
    if (sp_code %in% colnames(PC_matrix)) PC_sf$count <- PC_matrix[,sp_code]
    PC_sf <- subset(PC_sf, !is.na(count))
    
    # --------------------------------
    # Checklists will be treated as presence/absence only
    # --------------------------------
    
    CL_sf <- DO_surveyinfo %>% mutate(count = NA) # Checklist counts
    if (sp_code %in% colnames(DO_matrix)) CL_sf$count <- DO_matrix[,sp_code]
    CL_sf$presence <- as.numeric(CL_sf$count > 0)
    
    CL_sf$surveyID <- 1:nrow(CL_sf)
    
    # Remove NAs
    CL_sf <- subset(CL_sf, !is.na(count)) 
    
    # # --------------------------------
    # # Check if QPAD offsets exist, and generate them if so
    # # --------------------------------
    # 
    # offset_exists <- FALSE
    # offset_5min_Pointcount <- 0
    # sp_cr = sp_edr = NA
    # sp_napops <- subset(napops_species,Species == sp_code)
    # if (nrow(sp_napops)>0){
    #   if (sp_napops$Removal == 1 & sp_napops$Distance == 1){
    #     offset_exists <- TRUE
    #     sp_cr <- cue_rate(species = sp_code,od = 153, TSSr = 1, model = 1)
    #     sp_edr <- edr(species = sp_code,road = FALSE, forest = 0.5,model = 1)
    #     
    #     # Calculate A and p, which jointly determine offset
    #     PC_sf$A_metres <- c(pi*sp_edr$EDR_est^2)
    #     PC_sf$p <- 1-exp(-PC_sf$DurationInMinutes*sp_cr$CR_est[1,1])
    #     PC_sf$QPAD_offset <- log(PC_sf$A_metres * PC_sf$p)
    #     
    #     # Calculate offset for a 5-minute point count (needed later in analysis)
    #     offset_5min_Pointcount <- c(log((pi*sp_edr$EDR_est^2)*(1-exp(-5*sp_cr$CR_est[1,1]))))
    #   }
    # }
    
    # Intersect with SaskSquares dataframe
    PC_sf <- st_intersection(PC_sf, SaskSquares)
    CL_sf <- st_intersection(CL_sf, SaskSquares)
    
    # --------------------------------
    # STANDARDIZE 'TIME SINCE SUNRISE' COVARIATE
    # --------------------------------
    
    hss_mean <- mean(as.numeric(c(PC_sf$HSS,CL_sf$HSS)))
    hss_sd <- sd(as.numeric(c(PC_sf$HSS,CL_sf$HSS)))
    
    PC_sf$TSS <- (as.numeric(PC_sf$HSS)-hss_mean)/hss_sd
    CL_sf$TSS <- (as.numeric(CL_sf$HSS)-hss_mean)/hss_sd
    
    # --------------------------------
    # DEFINE SQUARE-DAY COVARIATE
    # --------------------------------
    squareday_df <- rbind(as.data.frame(PC_sf)[,c("Date","sq_id")],
                          as.data.frame(CL_sf)[,c("Date","sq_id")]) %>%
      unique()
    
    squareday_df$square_day <- 1:nrow(squareday_df)
    
    
    PC_sf <- left_join(PC_sf,squareday_df)
    CL_sf <- left_join(CL_sf,squareday_df)
    
    # --------------------------------
    # DEFINE OBSERVATION-LEVEL RANDOM EFFECT
    # --------------------------------
    
    PC_sf$surveyID <- 1:nrow(PC_sf)
    
    # --------------------------------
    # Separate different types of checklists
    # --------------------------------
    
    # Stationary counts
    SC_sf <- subset(CL_sf, ProtocolType == "Stationary count" & DurationInHours <= 1) # Almost all surveys are less than 1 hour
    
    dim(subset(CL_sf, ProtocolType == "Stationary count"))
    dim(SC_sf) # Dropped 712 surveys out of 9154
    
    # Linear transects (travel distance included)
    LT_sf <- subset(CL_sf, ProtocolType == "Linear transect"  & 
                      DurationInHours > (10/60) &
                      DurationInHours <= 2 &
                      TravelDistance_m <= 10000)
    
    dim(subset(CL_sf, ProtocolType == "Linear transect"))
    dim(LT_sf) # Dropped 1932 surveys out of 5094 (about half due to distance, other half due to time)
    
    #cor(LT_sf$DurationInHours,LT_sf$TravelDistance_m) # Only a correlation of about 0.15
    # --------------------------------
    # Construct standardized effort corrections
    # --------------------------------
    
    SC_sf$SC_duration <- scale(SC_sf$DurationInHours)
    LT_sf$LT_duration <- scale(LT_sf$DurationInHours)
    LT_sf$LT_distance <- scale(LT_sf$TravelDistance_m)
    
    # --------------------------------
    # Convert to spatial objects
    # --------------------------------
    
    PC_sp <- as(PC_sf,'Spatial')
    SC_sp <- as(SC_sf,'Spatial')
    LT_sp <- as(LT_sf,'Spatial')
    
    # --------------------------------
    # Create a spatial mesh, which is used to fit the residual spatial field
    # --------------------------------
    
    all_points <- rbind(PC_sp[,"obs_index"],SC_sp[,"obs_index"],LT_sp[,"obs_index"]) %>%
      st_as_sf()
    
    # Note: mesh developed using tutorial at: https://rpubs.com/jafet089/886687
    max.edge = diff(range(st_coordinates(all_points)[,1]))/15
    bound.outer = diff(range(st_coordinates(all_points)[,1]))/3
    cutoff = max.edge/5
    bound.outer = diff(range(st_coordinates(all_points)[,1]))/3
    
    mesh_spatial <- inla.mesh.2d(loc = st_coordinates(all_points),
                                 cutoff = max.edge/2, #max.edge/5,
                                 max.edge = c(1,2)*max.edge,
                                 offset=c(max.edge, bound.outer))
    mesh_locs <- mesh_spatial$loc[,c(1,2)] %>% as.data.frame()
    
    matern_coarse <- inla.spde2.pcmatern(mesh_spatial,
                                         #prior.range = c(max.edge*5, 0.1), # 10% chance range is smaller than 500000
                                         prior.range = c(250000,0.1),
                                         prior.sigma = c(2, 0.1), # 1% chance sd is larger than 2
                                         constr = TRUE # sum to 0 constraint
    )                
    
    # --------------------------------
    # Define random effect prior
    # --------------------------------
    
    kappa_prec <- list(prior = "pcprec", param = c(1,0.1))
    
    # --------------------------------
    # Create mesh to model effect of time since sunrise (TSS)
    # --------------------------------
    TSS_range <- range(c(PC_sf$TSS,CL_sf$TSS))
    TSS_meshpoints <- seq(TSS_range[1]-0.1,TSS_range[2]+0.1,length.out = 15)
    TSS_mesh1D = inla.mesh.1d(TSS_meshpoints,boundary="free")
    TSS_spde = inla.spde2.pcmatern(TSS_mesh1D,
                                   prior.range = c(3,0.5), 
                                   prior.sigma = c(2,0.1)) # 10% chance sd is larger than 2
    
    # --------------------------------
    # Create mesh to model effect of checklist duration
    # --------------------------------
    # 
    SC_duration_meshpoints <- seq(min(SC_sf$SC_duration)-0.1,max(SC_sf$SC_duration)+0.1,length.out = 11)
    SC_duration_mesh1D = inla.mesh.1d(SC_duration_meshpoints,boundary="free")
    SC_duration_spde = inla.spde2.pcmatern(SC_duration_mesh1D,
                                         prior.range = c(1,0.5), 
                                         prior.sigma = c(2,0.1)) # 10% chance sd is larger than 2
    
    LT_duration_meshpoints <- seq(min(LT_sf$LT_duration)-0.1,max(LT_sf$LT_duration)+0.1,length.out = 11)
    LT_duration_mesh1D = inla.mesh.1d(LT_duration_meshpoints,boundary="free")
    LT_duration_spde = inla.spde2.pcmatern(LT_duration_mesh1D,
                                         prior.range = c(1,0.5), 
                                         prior.sigma = c(2,0.1)) # 10% chance sd is larger than 2
    
    LT_distance_meshpoints <- seq(min(LT_sf$LT_distance)-0.1,max(LT_sf$LT_distance)+0.1,length.out = 11)
    LT_distance_mesh1D = inla.mesh.1d(LT_distance_meshpoints,boundary="free")
    LT_distance_spde = inla.spde2.pcmatern(LT_distance_mesh1D,
                                           prior.range = c(1,0.5), 
                                           prior.sigma = c(2,0.1)) # 10% chance sd is larger than 2
    
    # --------------------------------
    # Model components
    # --------------------------------
    
    covariates_to_include <- c("PC1","PC2","PC3","Water_5km")
    
    
    model_components = as.formula(paste0('~
  Intercept_PC(1)+
  Intercept_SC(1)+
  Intercept_LT(1)+
  TSS(main = TSS,model = TSS_spde) +
  kappa_surveyID(surveyID, model = "iid", constr = TRUE, hyper = list(prec = kappa_prec))+
  kappa_squareID(sq_idx, model = "iid", constr = TRUE, hyper = list(prec = kappa_prec))+
  kappa_squareday(square_day, model = "iid", constr = TRUE, hyper = list(prec = kappa_prec))+
  spde_coarse(main = coordinates, model = matern_coarse) +
  
  SC_duration(main = SC_duration,model = SC_duration_spde) +
  LT_duration(main = LT_duration,model = LT_duration_spde) +
  LT_distance(main = LT_distance,model = LT_distance_spde) +
    
  ',
                                         
                                         paste0("Beta1_",covariates_to_include,'(1,model="linear", mean.linear = 0, prec.linear = 4)', collapse = " + "))
    )
    
    model_components
    
    # --------------------------------
    # Model formulas
    # --------------------------------
    
    model_formula_PC = as.formula(paste0('count ~
                  Intercept_PC +
                  TSS +
                  kappa_surveyID +
                  kappa_squareID +
                  kappa_squareday +
                  spde_coarse +
                   ',
                                         paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + ")))
    
    model_formula_SC = as.formula(paste0('presence ~ log(1/exp(-exp(

                  Intercept_SC +
                  TSS +
                  kappa_squareID +
                  kappa_squareday +
                  spde_coarse +
                  SC_duration +
                                       ',
                                         paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
                                         "))-1)"))
    
    
    model_formula_LT = as.formula(paste0('presence ~ log(1/exp(-exp(

                  Intercept_LT +
                  TSS +
                  kappa_squareID +
                  kappa_squareday +
                  spde_coarse +
                  LT_distance +
                  LT_duration +
                                       ',
                                         paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
                                         "))-1)"))
    
    # --------------------------------
    # Withhold 20% of point count data for cross-validation
    # --------------------------------
    
    PC_xval <- PC_sp[PC_sp$fold == fold,]
    PC_sp <- PC_sp[PC_sp$fold != fold,]
    PC_xval$presence <- as.numeric(PC_xval$count > 0)
    
    # --------------------------------
    # Withhold 20% of checklist data for cross-validation
    # --------------------------------
    
    SC_xval <- SC_sp[SC_sp$fold == fold,]
    SC_sp <- SC_sp[SC_sp$fold != fold,]
    
    LT_xval <- LT_sp[LT_sp$fold == fold,]
    LT_sp <- LT_sp[LT_sp$fold != fold,]
    
    # --------------------------------
    # Specify model likelihoods
    # --------------------------------
    
    like_PC <- like(family = "poisson",
                    formula = model_formula_PC,
                    data = PC_sp)
    like_SC <- like(family = "binomial",
                    formula = model_formula_SC,
                    data = SC_sp)
    like_LT <- like(family = "binomial",
                    formula = model_formula_LT,
                    data = LT_sp)
    
    # &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    # FIT MODELS 
    # &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    
    
    bru_options_reset() # bru_options_set(inla.mode = "experimental")
    
    inits <- c(-5,-5,-5,rep(0,length(covariates_to_include))) %>% as.list()
    names(inits) <- c("Intercept_PC","Intercept_SC","Intercept_LT",paste0("Beta1_",covariates_to_include))
    
    # -----------------------------
    # Use only point count data
    # -----------------------------
    
    start <- Sys.time()
    fit_PConly <- bru(components = model_components,
                      like_PC,
                      options = list(
                        control.inla = list(int.strategy = "eb"),
                        bru_verbose = 4,
                        bru_max_iter = 5,
                        bru_initial = inits))
    end <- Sys.time()
    runtime_PConly <- difftime( end,start, units="mins") # 16 min
    
    # -----------------------------
    # Use checklists only
    # -----------------------------
    
    start <- Sys.time()
    fit_CLonly <- bru(components = model_components,
                      like_SC,like_LT,
                      options = list(
                        control.inla = list(int.strategy = "eb"),
                        bru_verbose = 4,
                        bru_max_iter = 5,
                        bru_initial = inits))
    end <- Sys.time()
    runtime_CLonly <- difftime( end,start, units="mins") # 9 min
    
    # -----------------------------
    # Use point counts and checklists
    # -----------------------------
    
    start <- Sys.time()
    fit_integrated <- bru(components = model_components,
                          like_PC,like_SC,like_LT,
                          options = list(
                            control.inla = list(int.strategy = "eb"),
                            bru_verbose = 4,
                            bru_max_iter = 5,
                            bru_initial = inits))
    end <- Sys.time()
    runtime_integrated <- difftime( end,start, units="mins") # 26 min
    
    # &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    # Predictions on withheld data
    # &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    
    nsamp = 500
    
    pred_formula_PC = as.formula(paste0(' ~

                  Intercept_PC +
                  TSS +
                  spde_coarse +
                   ',
                                        paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + ")))
    
    pred_formula_SC = as.formula(paste0(' ~

                  Intercept_SC +
                  TSS +
                  SC_duration +
                  spde_coarse +
                   ',
                                        paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + ")))
    
    pred_formula_LT = as.formula(paste0(' ~

                  Intercept_LT +
                  TSS +
                  LT_duration +
                  LT_distance +
                  spde_coarse +
                   ',
                                        paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + ")))
    
    # -------------------------
    # Predictions on point count data
    # -------------------------
    
    pred_PC_PConly <- generate(fit_PConly, 
                               PC_xval, 
                               formula = pred_formula_PC,
                               n.samples = nsamp) %>% 
      apply(.,1,median)
    
    pred_PC_integrated <- generate(fit_integrated, 
                                   PC_xval, 
                                   formula = pred_formula_PC,
                                   n.samples = nsamp) %>% 
      apply(.,1,median)
    
    # Add lognormal variance corrections, and exponentiate to place on count scale
    pred_PC_PConly     <- exp(pred_PC_PConly     + 0.5/summary(fit_PConly)$inla$hyperpar["Precision for kappa_surveyID",4]     + 0.5/summary(fit_PConly)$inla$hyperpar["Precision for kappa_squareID",4]     + 0.5/summary(fit_PConly)$inla$hyperpar["Precision for kappa_squareday",4])
    pred_PC_integrated <- exp(pred_PC_integrated + 0.5/summary(fit_integrated)$inla$hyperpar["Precision for kappa_surveyID",4] + 0.5/summary(fit_integrated)$inla$hyperpar["Precision for kappa_squareID",4] + 0.5/summary(fit_integrated)$inla$hyperpar["Precision for kappa_squareday",4])
    
    # -------------------------
    # Predictions on stationary checklists
    # -------------------------
    
    pred_SC_CLonly <- generate(fit_CLonly, 
                               SC_xval, 
                               formula = pred_formula_SC,
                               n.samples = nsamp) %>% 
      apply(.,1,median)
    
    pred_SC_integrated <- generate(fit_integrated, 
                                   SC_xval, 
                                   formula = pred_formula_SC,
                                   n.samples = nsamp) %>% 
      apply(.,1,median)
    
    # Add lognormal variance corrections, place on count scale
    pred_SC_CLonly     <- exp(pred_SC_CLonly + 0.5/summary(fit_CLonly)$inla$hyperpar["Precision for kappa_squareID",4]     + 0.5/summary(fit_CLonly)$inla$hyperpar["Precision for kappa_squareday",4])
    pred_SC_integrated <- exp(pred_SC_integrated + 0.5/summary(fit_integrated)$inla$hyperpar["Precision for kappa_squareID",4] + 0.5/summary(fit_integrated)$inla$hyperpar["Precision for kappa_squareday",4])
    
    
    # -------------------------
    # Predictions on linear transect checklists
    # -------------------------
    
    pred_LT_CLonly <- generate(fit_CLonly, 
                               LT_xval, 
                               formula = pred_formula_LT,
                               n.samples = nsamp) %>% 
      apply(.,1,median)
    
    pred_LT_integrated <- generate(fit_integrated, 
                                   LT_xval, 
                                   formula = pred_formula_LT,
                                   n.samples = nsamp) %>% 
      apply(.,1,median)
    
    # Add lognormal variance corrections, place on count scale
    pred_LT_CLonly     <- exp(pred_LT_CLonly + 0.5/summary(fit_CLonly)$inla$hyperpar["Precision for kappa_squareID",4]     + 0.5/summary(fit_CLonly)$inla$hyperpar["Precision for kappa_squareday",4])
    pred_LT_integrated <- exp(pred_LT_integrated + 0.5/summary(fit_integrated)$inla$hyperpar["Precision for kappa_squareID",4] + 0.5/summary(fit_integrated)$inla$hyperpar["Precision for kappa_squareday",4])
    
    # &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    # Calculate crossvalidation scores
    # &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    
    # -------------------------
    # Predictions on point count data
    # -------------------------
    
    # AUC (for presence/absence predictions)
    AUC_PC_PConly <- as.numeric(auc(PC_xval$presence, 1-exp(-pred_PC_PConly)))
    AUC_PC_integrated <- as.numeric(auc(PC_xval$presence, 1-exp(-pred_PC_integrated)))
    
    # correlation (for count predictions)
    cor_PC_PConly <- as.numeric(cor(PC_xval$count, pred_PC_PConly))
    cor_PC_integrated <- as.numeric(cor(PC_xval$count, pred_PC_integrated))
    
    # MSE (for count predictions)
    MSE_PC_PConly <- mean((PC_xval$count - pred_PC_PConly)^2)
    MSE_PC_integrated <- mean((PC_xval$count - pred_PC_integrated)^2)
    
    # -------------------------
    # Predictions on stationary checklists
    # -------------------------
    
    # AUC (for presence/absence predictions)
    AUC_SC_CLonly <- as.numeric(auc(SC_xval$presence, 1-exp(-pred_SC_CLonly)))
    AUC_SC_integrated <- as.numeric(auc(SC_xval$presence, 1-exp(-pred_SC_integrated)))
    
    # correlation (for count predictions)
    cor_SC_CLonly <- as.numeric(cor(SC_xval$count, pred_SC_CLonly))
    cor_SC_integrated <- as.numeric(cor(SC_xval$count, pred_SC_integrated))
    
    # MSE (for count predictions)
    MSE_SC_CLonly <- mean((SC_xval$count - pred_SC_CLonly)^2)
    MSE_SC_integrated <- mean((SC_xval$count - pred_SC_integrated)^2)
    
    # -------------------------
    # Predictions on linear transect checklists
    # -------------------------
    
    # AUC (for presence/absence predictions)
    AUC_LT_CLonly <- as.numeric(auc(LT_xval$presence, 1-exp(-pred_LT_CLonly)))
    AUC_LT_integrated <- as.numeric(auc(LT_xval$presence, 1-exp(-pred_LT_integrated)))
    
    # correlation (for count predictions)
    cor_LT_CLonly <- as.numeric(cor(LT_xval$count, pred_LT_CLonly))
    cor_LT_integrated <- as.numeric(cor(LT_xval$count, pred_LT_integrated))
    
    # MSE (for count predictions)
    MSE_LT_CLonly <- mean((LT_xval$count - pred_LT_CLonly)^2)
    MSE_LT_integrated <- mean((LT_xval$count - pred_LT_integrated)^2)
    
    
    # -------------------------------------------------------
    # Save results
    # -------------------------------------------------------
    
    if (file.exists("../AOS_precision/output/xval_TYPE2_df_50km_TSS_2.RData")){
      load("../AOS_precision/output/xval_TYPE2_df_50km_TSS_2.RData")
    }
    
    # Compare crossvalidation metrics
    xval_TYPE2_df_50km_TSS <- rbind(xval_TYPE2_df_50km_TSS,
                                    data.frame(Species = sp_code,
                                               xval_fold = fold,
                                               
                                               # Number of times species was observed in each dataset
                                               n_obs_PC = sum(PC_sp$count>0),
                                               n_obs_SC = sum(SC_sp$count>0),
                                               n_obs_LT = sum(LT_sp$count>0),
                                               
                                               # Amount of data in crossvalidation squares
                                               n_obs_PC_xval = sum(PC_xval$count>0),
                                               n_obs_SC_xval = sum(SC_xval$count>0),
                                               n_obs_LT_xval = sum(LT_xval$count>0),
                                               
                                               # -------------------------------
                                               # Crossvalidation metrics - point counts
                                               # -------------------------------
                                               
                                               # Correlation between predictions and observed
                                               cor_PC_PConly = cor_PC_PConly,
                                               cor_PC_integrated = cor_PC_integrated,
                                               
                                               # AUC
                                               AUC_PC_PConly = AUC_PC_PConly,
                                               AUC_PC_integrated = AUC_PC_integrated,
                                               
                                               # Mean squared error
                                               MSE_PC_PConly = MSE_PC_PConly,
                                               MSE_PC_integrated = MSE_PC_integrated,
                                               
                                               # -------------------------------
                                               # Crossvalidation metrics - point counts
                                               # -------------------------------
                                               
                                               # Correlation between predictions and observed
                                               cor_SC_CLonly = cor_SC_CLonly,
                                               cor_SC_integrated = cor_SC_integrated,
                                               
                                               # AUC
                                               AUC_SC_CLonly = AUC_SC_CLonly,
                                               AUC_SC_integrated = AUC_SC_integrated,
                                               
                                               # Mean squared error
                                               MSE_SC_CLonly = MSE_SC_CLonly,
                                               MSE_SC_integrated = MSE_SC_integrated,
                                               
                                               # -------------------------------
                                               # Crossvalidation metrics - point counts
                                               # -------------------------------
                                               
                                               # Correlation between predictions and observed
                                               cor_LT_CLonly = cor_LT_CLonly,
                                               cor_LT_integrated = cor_LT_integrated,
                                               
                                               # AUC
                                               AUC_LT_CLonly = AUC_LT_CLonly,
                                               AUC_LT_integrated = AUC_LT_integrated,
                                               
                                               # Mean squared error
                                               MSE_LT_CLonly = MSE_LT_CLonly,
                                               MSE_LT_integrated = MSE_LT_integrated
                                               
                                    )
    )
    
    save(xval_TYPE2_df_50km_TSS, file = "../AOS_precision/output/xval_TYPE2_df_50km_TSS_2.RData")
    
    rm(list = c("fit_PConly","fit_integrated","fit_CLonly"))
    
  } # close species loop
} # xval fold
