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

# Custom cut points
cut.fn2 <- function(df, column_name, cut_levs = cut_levs){
  
  
  lower_bound <- min(cut_levs)
  upper_bound <- max(cut_levs)
  
  max_lev <- length(cut_levs) ## do this because sometimes there are fewer levels
  
  cut_levs_labs <- c(paste0("0-",cut_levs[1]),
                     paste(cut_levs[-max_lev], cut_levs[-1], sep="-"), 
                     paste(cut_levs[max_lev], "+"))
  
  cut_levs <- c(-1, cut_levs, 1000) %>% unique()
  
  # ** MULTIPLY BY 15 FOR COMPARISON TO BSC
  df$pc_cut <- cut(as.data.frame(df)[,column_name], 
                   cut_levs, 
                   labels=cut_levs_labs, 
                   ordered=TRUE)
  
  sp_abund_raster <- fasterize(df,SaskRaster,field = "pc_cut")
  sp_abund_raster <- mask(sp_abund_raster,SaskBoundary)
  
  sp_abund_raster_data <- as.data.frame(rasterToPoints(sp_abund_raster))
  for (i in unique(sp_abund_raster_data$layer)) sp_abund_raster_data$layer[sp_abund_raster_data$layer==i] <- cut_levs_labs[i]
  # sp_abund_raster_data$layer[sp_abund_raster_data$layer==1] <- cut_levs_labs[1]
  # sp_abund_raster_data$layer[sp_abund_raster_data$layer==2] <- cut_levs_labs[2]
  # sp_abund_raster_data$layer[sp_abund_raster_data$layer==3] <- cut_levs_labs[3]
  # sp_abund_raster_data$layer[sp_abund_raster_data$layer==4] <- cut_levs_labs[4]
  # sp_abund_raster_data$layer[sp_abund_raster_data$layer==5] <- cut_levs_labs[5]
  # sp_abund_raster_data$layer[sp_abund_raster_data$layer==6] <- cut_levs_labs[6]
  # 
  sp_abund_raster_data$layer <- factor(sp_abund_raster_data$layer,
                                       levels= cut_levs_labs,
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

# napops_species <- napops::list_species()

# *************************************************************
# Calculate number of unique atlas squares in which each species was detected
# (in both point counts and checklists)
# *************************************************************

# Load squares to withhold for crossvalidation
SaskSquares <- read_sf("../AOS_precision/output/SaskSquares_xval_20km.shp") %>%
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


ggplot(species_to_fit, aes(x = n_PC, y = n_CL,label=Species, col = prop_CLonly))+
  geom_point()+
  geom_text_repel()+
  scale_color_gradientn(colors = viridis(10))+
  theme_bw()

# ************************************************************************************
# ------------------------------------------------------------------------------------
# Analysis (note some additional data wrangling is needed on a species-by-species basis)
# ------------------------------------------------------------------------------------
# ************************************************************************************
SaskWater <- read_sf("!Shapefiles/Water/SaskWaterClip.shp")
SaskWater <- st_transform(SaskWater, crs = st_crs(SaskBoundary))
SaskWater$area <- st_area(SaskWater)
SaskWater$area <- as.numeric(SaskWater$area)
SaskWater <- SaskWater[SaskWater$area>2.580e+06 ,]

Water_centroids <- SaskGrid %>%
  dplyr::select(pointid) %>%
  st_centroid() %>%
  st_intersection(.,SaskWater)

surface_comparison <- data.frame()
if (file.exists("../AOS_precision/output/surface_comparison_RE_effort.RData")){
  load("../AOS_precision/output/surface_comparison_RE_effort.RData")
}

# covariates to include in models
covariates_to_include <- c("PC1","PC2","PC3","Water_5km")

for (sp_code in species_to_fit$Species){

  if (file.exists("../AOS_precision/output/surface_comparison_RE_effort.RData")){
    load("../AOS_precision/output/surface_comparison_RE_effort.RData")
  }
  
  # Check if this species/xval fold have already been run. If so, skip
  if (nrow(surface_comparison)>0){
    if (nrow(subset(surface_comparison,Species == sp_code))>0) next
  }
  
  print(sp_code)
  
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
  #     sp_cr <- cue_rate(species = sp_code,od = 153, tssr = 1, model = 1)
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
  
  PC_sp <- as(PC_sf,'Spatial')
  
  # --------------------------------
  # Separate different types of checklists
  # --------------------------------
  
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
  
  SC_sp <- as(SC_sf,'Spatial')
  LT_sp <- as(LT_sf,'Spatial')
  BBA_sp <- as(BBA_sf,'Spatial')
  
  # --------------------------------
  # Create a spatial mesh, which is used to fit the residual spatial field
  # --------------------------------
  
  all_points <- rbind(PC_sp[,"obs_index"],SC_sp[,"obs_index"],LT_sp[,"obs_index"],BBA_sp[,"obs_index"]) %>%
    st_as_sf()
  
  # Note: mesh developed using tutorial at: https://rpubs.com/jafet089/886687
  max.edge = diff(range(st_coordinates(all_points)[,1]))/15
  bound.outer = diff(range(st_coordinates(all_points)[,1]))/3
  cutoff = max.edge/5
  bound.outer = diff(range(st_coordinates(all_points)[,1]))/3
  
  mesh_spatial <- inla.mesh.2d(loc = st_coordinates(all_points),
                               cutoff = max.edge/5,
                               max.edge = c(1,2)*max.edge,
                               offset=c(max.edge, bound.outer))
  
  mesh_locs <- mesh_spatial$loc[,c(1,2)] %>% as.data.frame()
  dim(mesh_locs)
  #plot(mesh_spatial)
  
  matern_coarse <- inla.spde2.pcmatern(mesh_spatial,
                                       prior.range = c(max.edge*5, 0.1), # 10% chance range is smaller than 500000
                                       prior.sigma = c(1, 0.1)      # 10% chance sd is larger than 1
  )                # sum to 1 constraint
  
  # --------------------------------
  # Define square-level random effects
  # --------------------------------
  
  # Point counts
  kappa_prec <- list(prior = "pcprec", param = c(1,0.1))
  
  # --------------------------------
  # Create 'temporal mesh' to model effect of checklist duration
  # --------------------------------
  # 
  SC_duration_meshpoints <- seq(0,max(SC_sf$DurationInHours)+0.1,length.out = 11)
  SC_duration_mesh1D = inla.mesh.1d(SC_duration_meshpoints,boundary="free")
  SC_duration_spde = inla.spde2.pcmatern(SC_duration_mesh1D,
                                         prior.range = c(1,0.9), # 90% chance range is smaller than 1
                                         prior.sigma = c(2,0.1)) # 10% chance sd is larger than 2
  # 
  # LT_duration_meshpoints <- seq(0,max(LT_sf$DurationInHours)+0.1,length.out = 11)
  # LT_duration_mesh1D = inla.mesh.1d(LT_duration_meshpoints,boundary="free")
  # LT_duration_spde = inla.spde2.pcmatern(LT_duration_mesh1D,
  #                                        prior.range = c(1,0.9), # 90% chance range is smaller than 1
  #                                        prior.sigma = c(2,0.1)) # 10% chance sd is larger than 2
  # 
  # 
  # BBA_duration_meshpoints <- seq(0,max(BBA_sf$DurationInHours)+0.1,length.out = 11)
  # BBA_duration_mesh1D = inla.mesh.1d(BBA_duration_meshpoints,boundary="free")
  # BBA_duration_spde = inla.spde2.pcmatern(BBA_duration_mesh1D,
  #                                         prior.range = c(1,0.9), # 90% chance range is smaller than 1
  #                                         prior.sigma = c(2,0.1)) # 10% chance sd is larger than 2
  # 
  
  # --------------------------------
  # NOTE: I've written the code below so we can add/remove covariates easily by including them
  #       in the vector called "covariates_to_include"
  #
  #       We could simplify this potentially to make the code easier to read/understand
  # --------------------------------
  
  # --------------------------------
  # Model components
  # --------------------------------
  
  covariates_to_include <- c("PC1","PC2","PC3","Water_5km")
  
  
  # LT_effort(main = DurationInHours,model = LT_duration_spde) +
  # BBA_effort(main = DurationInHours,model = BBA_duration_spde) +
  #Intercept_LT(1)+
  #Intercept_BBA(1)+
  #kappa_BBA(sq_idx, model = "iid", constr = TRUE, hyper = list(prec = kappa_prec))+
  
  model_components = as.formula(paste0('~
  Intercept_PC(1)+
  Intercept_SC(1)+
  
  kappa_PC(sq_idx, model = "iid", constr = TRUE, hyper = list(prec = kappa_prec))+
  kappa_SC(sq_idx, model = "iid", constr = TRUE, hyper = list(prec = kappa_prec))+
  kappa_shared(sq_idx, model = "iid", constr = TRUE, hyper = list(prec = kappa_prec))+
  SC_effort(main = DurationInHours,model = SC_duration_spde) +
  spde_coarse(main = coordinates, model = matern_coarse) + 
  ',
                                       
                                       paste0("Beta1_",covariates_to_include,'(1,model="linear", mean.linear = 0, prec.linear = 4)', collapse = " + "),
                                       " + ",
                                       paste0("Beta2_",covariates_to_include,'(1,model="linear", mean.linear = 0, prec.linear = 16)', collapse = " + "))
  )
  
  model_components
  
  # --------------------------------
  # Model formulas
  # --------------------------------
  
  model_formula_PConly = as.formula(paste0('count ~
                  Intercept_PC +
                  kappa_PC +
                  spde_coarse +
                   ',
                                           paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
                                           " + ",
                                           paste0("Beta2_",covariates_to_include,'*',covariates_to_include,"^2", collapse = " + ")))
  
  model_formula_PC = as.formula(paste0('count ~
                  Intercept_PC +
                  spde_coarse +
                  kappa_PC +
                  kappa_shared +
                   ',
                                       paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
                                       " + ",
                                       paste0("Beta2_",covariates_to_include,'*',covariates_to_include,"^2", collapse = " + ")))
  
  model_formula_SC = as.formula(paste0('presence ~ log(1/exp(-exp(

                  Intercept_SC +
                  spde_coarse +
                  kappa_SC +
                  kappa_shared +
                  SC_effort +
                                       ',
                                       paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
                                       " + ",
                                       paste0("Beta2_",covariates_to_include,'*',covariates_to_include,"^2", collapse = " + "),
                                       "))-1)"))
  
  model_formula_LT = as.formula(paste0('presence ~ log(1/exp(-exp(

                  Intercept_LT +
                  spde_coarse +
                  kappa_shared +
                                       ',
                                       paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
                                       " + ",
                                       paste0("Beta2_",covariates_to_include,'*',covariates_to_include,"^2", collapse = " + "),
                                       "))-1)"))
  
  model_formula_BBA = as.formula(paste0('presence ~ log(1/exp(-exp(

                  Intercept_BBA +
                  spde_coarse +
                  kappa_shared +
                                       ',
                                        paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
                                        " + ",
                                        paste0("Beta2_",covariates_to_include,'*',covariates_to_include,"^2", collapse = " + "),
                                        "))-1)"))
  
  # --------------------------------
  # Specify model likelihoods
  # --------------------------------
  like_PConly <- like(family = "poisson",
                      formula = model_formula_PConly,
                      data = PC_sp)
  like_PC <- like(family = "poisson",
                  formula = model_formula_PC,
                  data = PC_sp)
  like_SC <- like(family = "binomial",
                  formula = model_formula_SC,
                  data = SC_sp)
  like_LT <- like(family = "binomial",
                  formula = model_formula_LT,
                  data = LT_sp)
  like_BBA <- like(family = "binomial",
                   formula = model_formula_BBA,
                   data = BBA_sp)
  
  # --------------------------------
  # Select reasonable initial values (should not affect inference, but affects model convergence)
  # --------------------------------
  
  inits <- c(-5,-5,rep(0,length(covariates_to_include)*2)) %>% as.list()
  #,"Intercept_LT","Intercept_BBA"
  names(inits) <- c("Intercept_PC","Intercept_SC",paste0("Beta1_",covariates_to_include), paste0("Beta2_",covariates_to_include))
  
  # --------------------------------
  # Fit models
  # --------------------------------
  
  #bru_options_set(inla.mode = "experimental")
  #bru_options_reset()
  
  start <- Sys.time()
  
  fit_PConly <- bru(components = model_components,
                    like_PConly,
                    options = list(
                      control.inla = list(int.strategy = "eb"),
                      bru_verbose = 4,
                      bru_max_iter = 25,
                      bru_initial = inits)) 
  
  fit_integrated <- bru(components = model_components, 
                        like_PC,like_SC,#like_LT,like_BBA,
                        options = list(#control.compute = list(waic = TRUE, cpo = TRUE, config = TRUE),
                          control.inla = list(int.strategy = "eb"),
                          bru_verbose = 4,
                          bru_max_iter = 25,
                          bru_initial = inits))
  
  # &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  # PREDICTION SURFACES
  # &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  
  pred_formula_PC = as.formula(paste0(' ~

                  Intercept_PC +
                  spde_coarse +
                   ',
                                      paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
                                      " + ",
                                      paste0("Beta2_",covariates_to_include,'*',covariates_to_include,"^2", collapse = " + ")))
  
  species_name = Sask_spcd$CommonName[which(Sask_spcd$spcd == sp_code)]
  species_label = Sask_spcd$Label[which(Sask_spcd$spcd == sp_code)]
  
  # ********
  # PREDICTIONS FROM MODEL FIT TO POINT COUNTS ONLY
  # ********
  
  nsamp = 100
  pred_surface_PConly <- generate(fit_PConly,
                                  as(SaskGrid,'Spatial'),
                                  formula =  pred_formula_PC,
                                  n.samples = nsamp) %>%
    apply(.,1,median) %>%
    exp()
  
  pred_surface_integrated <- generate(fit_integrated,
                                      as(SaskGrid,'Spatial'),
                                      formula =  pred_formula_PC,
                                      n.samples = nsamp)%>%
    apply(.,1,median) %>%
    exp()
  end = Sys.time()
  
  # Add lognormal variance correction
  pred_surface_PConly <- pred_surface_PConly * exp(0.5*1/summary(fit_PConly)$inla$hyperpar["Precision for kappa_PC",4])
  
  # Add lognormal variance correction
  pred_surface_integrated <- pred_surface_integrated * exp(0.5*1/summary(fit_integrated)$inla$hyperpar["Precision for kappa_PC",4]) * exp(0.5*1/summary(fit_integrated)$inla$hyperpar["Precision for kappa_shared",4])
  
  #rm(list = c("fit_integrated","fit_PConly"))
  
  pred_surface_PConly[Water_centroids$pointid] <- NA # Trim out pixels that have centroids in open water (Note: probably a better way to do this
  pred_surface_integrated[Water_centroids$pointid] <- NA # Trim out pixels that have centroids in open water (Note: probably a better way to do this)
  
  pred_grid_sp <- SaskGrid
  pred_grid_sp$pred_PConly_med  <- pred_surface_PConly
  pred_grid_sp$pred_integrated_med  <- pred_surface_integrated
  
  # ---------------------------------------------
  # Prepare to plot locations where species was observed
  # ---------------------------------------------
  
  SaskSquares_species <- SaskSquares %>%
    dplyr::rename(sq_id = SQUARE_ID)
  
  # Summary of point count information in each square
  PC_summary <- PC_sf %>%
    as.data.frame() %>%
    group_by(sq_id) %>%
    summarize(mean_count_PC = mean(count,na.rm = TRUE),
              observed_PC = as.numeric(mean(count,na.rm = TRUE)>0),
              n_PC = n())
  
  # Summary of checklist information in each square
  CL_summary <- CL_sf %>%
    as.data.frame() %>%
    group_by(sq_id) %>%
    summarize(mean_count_CL = mean(count,na.rm = TRUE),
              observed_CL = as.numeric(mean(count,na.rm = TRUE)>0),
              n_CL = n())
  
  SaskSquares_species <- SaskSquares_species %>% 
    left_join(PC_summary) %>%
    left_join(CL_summary) %>%
    dplyr::select(sq_id,observed_PC,n_PC,observed_CL,n_CL,geometry)%>%
    rowwise() %>% 
    mutate(observed_either = sum(observed_PC,observed_CL,na.rm=TRUE))
  SaskSquares_species$observed_either[is.na(SaskSquares_species$observed_PC) & is.na(SaskSquares_species$observed_CL)] <- NA
  
  SaskSquares_species_centroids <- st_centroid(SaskSquares_species)
  
  # ---------------------------------------------
  # PLOTS WITH HUMAN-RELEVANT COLOR SCALE
  # ---------------------------------------------
  
  raster_pred_surface_PConly <- cut.fn2(df = pred_grid_sp,
                                        column_name = "pred_PConly_med",
                                        cut_levs <- 1/c(100,50,25,10,2,1))
  
  raster_pred_surface_integrated <- cut.fn2(df = pred_grid_sp,
                                            column_name = "pred_integrated_med",
                                            cut_levs <- 1/c(100,50,25,10,2,1))
  
  pred_surface_map_PConly_q50 <- ggplot() +
    geom_sf(data = SaskBoundary,colour="black",fill="#f2f2f2",lwd=0.3,show.legend = F) +
    geom_raster(data = raster_pred_surface_PConly , aes(x = x, y = y, fill = layer)) +
    scale_fill_manual(name = "<span style='font-size:13pt'>Relative Abundance</span><br><span style='font-size:7pt'>Per point count</span><br><span style='font-size:7pt'>(Posterior Median)</span>",
                      values = colpal_relabund(length(levels(raster_pred_surface_PConly$layer))), drop=FALSE)+
    
    geom_sf(data = SaskWater,colour=NA,fill="#59F3F3",show.legend = F)+
    geom_sf(data = SaskBoundary,colour="black",fill=NA,lwd=0.3,show.legend = F)+
    
    # Place black dots in squares where species was observed
    #geom_sf(data = subset(SaskSquares_species_centroids,observed_PC > 0 & !is.na(observed_PC)),pch=19, size=0.1,show.legend = F, col = "black")+
    # Place small gray dots in squares where species was not observed
    #geom_sf(data = subset(SaskSquares_species_centroids,observed_PC == 0 & !is.na(observed_PC)),pch=1, size=0.1,show.legend = F, col = "gray70")+
    
    coord_sf(clip = "off",xlim = c(min(SaskPoints$x), max(SaskPoints$x)))+
    theme(panel.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    annotate(geom="text",x=346000,y=960000, label= paste0(species_label),lineheight = .85,hjust = 0,size=6,fontface =2) +
    annotate(geom="text",x=346000,y=550000, label= paste0("Prepared on ",Sys.Date()),size=3,lineheight = .75,hjust = 0,color="#3b3b3b")+
    theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(5,10,5,-20),legend.title.align=0.5,
          legend.title = element_markdown(lineheight=.9,hjust = "left"))
  
  #print(pred_surface_map_PConly_q50)
  
  png(paste0("../AOS_precision/output/figures/",sp_code,"_plot1_PConly_effort.png"), width=6.5, height=8, units="in", res=300, type="cairo")
  print(pred_surface_map_PConly_q50)
  dev.off()
  
  pred_surface_map_integrated_q50 <- ggplot() +
    geom_sf(data = SaskBoundary,colour="black",fill="#f2f2f2",lwd=0.3,show.legend = F) +
    geom_raster(data = raster_pred_surface_integrated , aes(x = x, y = y, fill = layer)) +
    scale_fill_manual(name = "<span style='font-size:13pt'>Relative Abundance</span><br><span style='font-size:7pt'>Per point count</span><br><span style='font-size:7pt'>(Posterior Median)</span>",
                      values = colpal_relabund(length(levels(raster_pred_surface_integrated$layer))), drop=FALSE)+
    
    geom_sf(data = SaskWater,colour=NA,fill="#59F3F3",show.legend = F)+
    geom_sf(data = SaskBoundary,colour="black",fill=NA,lwd=0.3,show.legend = F)+
    
    # Place black dots in squares where species was observed
    #geom_sf(data = subset(SaskSquares_species_centroids,observed_PC > 0 & !is.na(observed_PC)),pch=19, size=0.1,show.legend = F, col = "black")+
    # Place small gray dots in squares where species was not observed
    #geom_sf(data = subset(SaskSquares_species_centroids,observed_PC == 0 & !is.na(observed_PC)),pch=1, size=0.1,show.legend = F, col = "gray70")+
    
    # Also plot squares where species was observed in checklists
    #geom_sf(data = subset(SaskSquares_species,observed_CL > 0 & !is.na(observed_CL)),size=0.1,show.legend = F, col = "black", fill = "transparent")+
    
    
    coord_sf(clip = "off",xlim = c(min(SaskPoints$x), max(SaskPoints$x)))+
    theme(panel.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    annotate(geom="text",x=346000,y=960000, label= paste0(species_label),lineheight = .85,hjust = 0,size=6,fontface =2) +
    annotate(geom="text",x=346000,y=550000, label= paste0("Prepared on ",Sys.Date()),size=3,lineheight = .75,hjust = 0,color="#3b3b3b")+
    theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(5,10,5,-20),legend.title.align=0.5,
          legend.title = element_markdown(lineheight=.9,hjust = "left"))
  
  #print(pred_surface_map_integrated_q50)
  
  png(paste0("../AOS_precision/output/figures/",sp_code,"_plot1_integrated_effort.png"), width=6.5, height=8, units="in", res=300, type="cairo")
  print(pred_surface_map_integrated_q50)
  dev.off()
  
  # ---------------------------------------------
  # PLOTS WITH BIOLOGICALLY RELEVANT COLORSCALE
  # ---------------------------------------------
  
  lower_bound <- 0.01
  upper_bound <- quantile(pred_grid_sp$pred_integrated_med,0.99,na.rm = TRUE) %>% signif(2)
  if (lower_bound >= upper_bound){
    upper_bound <- quantile(pred_grid_sp$pred_integrated_med,0.99,na.rm = TRUE) %>% signif(2)
    lower_bound <- (upper_bound/5) %>% signif(2)
  }
  
  
  raster_pred_surface_PConly <- cut.fn(df = pred_grid_sp,
                                        column_name = "pred_PConly_med",
                                       lower_bound = lower_bound,
                                       upper_bound = upper_bound)
  
  raster_pred_surface_integrated <- cut.fn(df = pred_grid_sp,
                                           column_name = "pred_integrated_med",
                                           lower_bound = lower_bound,
                                           upper_bound = upper_bound)
  
  pred_surface_map_PConly_q50 <- ggplot() +
    geom_sf(data = SaskBoundary,colour="black",fill="#f2f2f2",lwd=0.3,show.legend = F) +
    geom_raster(data = raster_pred_surface_PConly , aes(x = x, y = y, fill = layer)) +
    scale_fill_manual(name = "<span style='font-size:13pt'>Relative Abundance</span><br><span style='font-size:7pt'>Per point count</span><br><span style='font-size:7pt'>(Posterior Median)</span>",
                      values = colpal_relabund(length(levels(raster_pred_surface_PConly$layer))), drop=FALSE)+
    
    geom_sf(data = SaskWater,colour=NA,fill="#59F3F3",show.legend = F)+
    geom_sf(data = SaskBoundary,colour="black",fill=NA,lwd=0.3,show.legend = F)+
    
    # Place black dots in squares where species was observed
    geom_sf(data = subset(SaskSquares_species_centroids,observed_PC > 0 & !is.na(observed_PC)),pch=19, size=0.1,show.legend = F, col = "black")+
    # Place small gray dots in squares where species was not observed
    geom_sf(data = subset(SaskSquares_species_centroids,observed_PC == 0 & !is.na(observed_PC)),pch=1, size=0.1,show.legend = F, col = "gray70")+
    
    coord_sf(clip = "off",xlim = c(min(SaskPoints$x), max(SaskPoints$x)))+
    theme(panel.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    annotate(geom="text",x=346000,y=960000, label= paste0(species_label),lineheight = .85,hjust = 0,size=6,fontface =2) +
    annotate(geom="text",x=346000,y=550000, label= paste0("Prepared on ",Sys.Date()),size=3,lineheight = .75,hjust = 0,color="#3b3b3b")+
    theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(5,10,5,-20),legend.title.align=0.5,
          legend.title = element_markdown(lineheight=.9,hjust = "left"))
  
  #print(pred_surface_map_PConly_q50)
  
  png(paste0("../AOS_precision/output/figures/",sp_code,"_plot2_PConly_effort.png"), width=6.5, height=8, units="in", res=300, type="cairo")
  print(pred_surface_map_PConly_q50)
  dev.off()
  
  pred_surface_map_integrated_q50 <- ggplot() +
    geom_sf(data = SaskBoundary,colour="black",fill="#f2f2f2",lwd=0.3,show.legend = F) +
    geom_raster(data = raster_pred_surface_integrated , aes(x = x, y = y, fill = layer)) +
    scale_fill_manual(name = "<span style='font-size:13pt'>Relative Abundance</span><br><span style='font-size:7pt'>Per point count</span><br><span style='font-size:7pt'>(Posterior Median)</span>",
                      values = colpal_relabund(length(levels(raster_pred_surface_integrated$layer))), drop=FALSE)+
    
    geom_sf(data = SaskWater,colour=NA,fill="#59F3F3",show.legend = F)+
    geom_sf(data = SaskBoundary,colour="black",fill=NA,lwd=0.3,show.legend = F)+
    
    
    # Place black dots in squares where species was observed
    geom_sf(data = subset(SaskSquares_species_centroids,observed_PC > 0 & !is.na(observed_PC)),pch=19, size=0.1,show.legend = F, col = "black")+
    
    # Place black diamonds in squares where species was observed in checklists
    geom_sf(data = subset(SaskSquares_species_centroids,observed_CL > 0 & !is.na(observed_CL)),pch=5, size=0.5,show.legend = F, col = "black")+
    
    # Place small gray dots in squares where species was not observed
    geom_sf(data = subset(SaskSquares_species_centroids,observed_either == 0 & !is.na(observed_either)),pch=1, size=0.1,show.legend = F, col = "gray70")+
    
    coord_sf(clip = "off",xlim = c(min(SaskPoints$x), max(SaskPoints$x)))+
    theme(panel.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    annotate(geom="text",x=346000,y=960000, label= paste0(species_label),lineheight = .85,hjust = 0,size=6,fontface =2) +
    annotate(geom="text",x=346000,y=550000, label= paste0("Prepared on ",Sys.Date()),size=3,lineheight = .75,hjust = 0,color="#3b3b3b")+
    theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(5,10,5,-20),legend.title.align=0.5,
          legend.title = element_markdown(lineheight=.9,hjust = "left"))
  
  #print(pred_surface_map_integrated_q50)
  
  png(paste0("../AOS_precision/output/figures/",sp_code,"_plot2_integrated_effort.png"), width=6.5, height=8, units="in", res=300, type="cairo")
  print(pred_surface_map_integrated_q50)
  dev.off()
  
  # -------------------------------------------------------
  # Calculate difference between surfaces
  # -------------------------------------------------------
  
  # Percent of landscape "in range"
  range_PConly <- mean(pred_grid_sp$pred_PConly_med>0.01,na.rm = TRUE)
  range_integrated <- mean(pred_grid_sp$pred_integrated_med>0.01,na.rm = TRUE)
  
  # Percent of landscape "in core of range"
  core_PConly <- mean(pred_grid_sp$pred_PConly_med>0.1,na.rm = TRUE)
  core_integrated <- mean(pred_grid_sp$pred_integrated_med>0.1,na.rm = TRUE)
  
  # Correlation between maps
  correlation_between_surfaces <- cor(pred_grid_sp$pred_PConly_med,pred_grid_sp$pred_integrated_med, use = "complete.obs")
  
  
  end <- Sys.time() 
  runtime_mins <- difftime( end,start, units="mins")
  # -------------------------------------------------------
  # Save results
  # -------------------------------------------------------
  
  if (file.exists("../AOS_precision/output/surface_comparison_RE_effort.RData")){
    load("../AOS_precision/output/surface_comparison_RE_effort.RData")
  }
  
  surface_comparison <- rbind(surface_comparison,
                              data.frame(Species = sp_code,
                                         range_PConly = range_PConly,
                                         range_integrated = range_integrated,
                                         core_PConly = core_PConly,
                                         core_integrated = core_integrated,
                                         correlation_between_surfaces = correlation_between_surfaces,
                                         runtime_mins = round(as.numeric(runtime_mins))
                              ))
  
  
  print(paste(sp_code," ... ",round(runtime_mins),"mins"))
  
  save(surface_comparison, file = "../AOS_precision/output/surface_comparison_RE_effort.RData")
  
  rm(list = c("pred_grid_sp","pred_surface_integrated","pred_surface_PConly","fit_integrated","fit_PConly","raster_pred_surface_PConly","raster_pred_surface_integrated"))
  
} # close species loop
