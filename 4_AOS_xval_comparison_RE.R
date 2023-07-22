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
if (file.exists("../AOS_precision/output/xval_df_integrated_RE.RData")){
  load("../AOS_precision/output/xval_df_integrated_RE.RData")
}

# covariates to include in models
covariates_to_include <- c("PC1","PC2","PC3","Water_5km")
#for (fold in sort(unique(SaskSquares$fold))){
fold = 1

for (sp_code in rev(species_to_fit$Species)){
   
    if (file.exists("../AOS_precision/output/xval_df_integrated_RE.RData")){
      load("../AOS_precision/output/xval_df_integrated_RE.RData")
    }
    
    # Check if this species/xval fold have already been run. If so, skip
    if (nrow(xval_df)>0){
      if (nrow(subset(xval_df,Species == sp_code & xval_fold == fold))>0) next
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
    
    # --------------------------------
    # Withhold 20% of point count data for cross-validation
    # --------------------------------
    
    PC_sp <- as(PC_sf,'Spatial')
    PC_xval <- PC_sp[PC_sp$fold == fold,]
    PC_sp <- PC_sp[PC_sp$fold != fold,]
    
    PC_xval$presence <- as.numeric(PC_xval$count > 0)
    
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
    
    # In meters
    #mesh_cutoff <- 100         
    #max_edge_inner <- 20000   
    #max_edge_outer <- 100000
    
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
                                         prior.sigma = c(1, 0.1), # 10% chance sd is larger than 1
                                         constr = TRUE # sum to 1 constraint
    )                
    
    # --------------------------------
    # Define square-level random effects
    # --------------------------------
    
    # Point counts
    kappa_prec <- list(prior = "pcprec", param = c(1,0.1))
    
    # --------------------------------
    # Create 'temporal mesh' to model effect of checklist duration
    # --------------------------------
    # 
    # SC_duration_meshpoints <- seq(0,max(SC_sf$DurationInHours)+0.1,length.out = 11)
    # SC_duration_mesh1D = inla.mesh.1d(SC_duration_meshpoints,boundary="free")
    # SC_duration_spde = inla.spde2.pcmatern(SC_duration_mesh1D,
    #                                        prior.range = c(1,0.9), # 90% chance range is smaller than 1
    #                                        prior.sigma = c(2,0.1)) # 10% chance sd is larger than 2
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
    
    # SC_effort(main = DurationInHours,model = SC_duration_spde) +
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
                        bru_max_iter = 5,
                        bru_initial = inits)) 
    
    fit_integrated <- bru(components = model_components, 
                          like_PC,like_SC,#like_LT,like_BBA,
                          options = list(#control.compute = list(waic = TRUE, cpo = TRUE, config = TRUE),
                            control.inla = list(int.strategy = "eb"),
                            bru_verbose = 4,
                            bru_max_iter = 5,
                            bru_initial = inits))
    
    # &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    # Cross-validation on withheld data
    # &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    
    pred_formula_PC = as.formula(paste0(' ~

                  Intercept_PC +
                  spde_coarse +
                   ',
                                        paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
                                        " + ",
                                        paste0("Beta2_",covariates_to_include,'*',covariates_to_include,"^2", collapse = " + ")))
    
    # -------------------------------------------------------
    # Generate predictions
    # -------------------------------------------------------

    pred_PConly <- generate(fit_PConly, 
                            PC_xval, 
                            formula = pred_formula_PC,
                            n.samples = 500)
    
    pred_integrated <- generate(fit_integrated, 
                                PC_xval, 
                                formula = pred_formula_PC,
                                n.samples = 500)

    # &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    # Compare crossvalidation accuracy between the two models
    # &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    
    # pointwise predictive density
    ppd_PConly <- ppd_integrated <-rep(NA, nrow(pred_PConly))
    for (i in 1:length(ppd_PConly)){
      ppd_PConly[i] <- mean(dpois(PC_xval$count[i],exp(pred_PConly[i,])))
      ppd_integrated[i] <- mean(dpois(PC_xval$count[i],exp(pred_integrated[i,])))
    }
    
    # AUC (for presence/absence predictions)
    AUC_PConly <- AUC_integrated <- rep(NA, ncol(pred_PConly))
    
    suppressMessages(
      for (i in 1:length(AUC_PConly)){
        AUC_PConly[i] <- as.numeric(auc(PC_xval$presence, 1-exp(-exp(pred_PConly[,i]))))
        AUC_integrated[i] <- as.numeric(auc(PC_xval$presence, 1-exp(-exp(pred_integrated[,i]))))
      }
    )
    
    # correlation (for count predictions)
    cor_PConly <- cor_integrated <- rep(NA, ncol(pred_PConly))
    for (i in 1:length(cor_PConly)){
      cor_PConly[i] <- as.numeric(cor(PC_xval$count, exp(pred_PConly[,i])))
      cor_integrated[i] <- as.numeric(cor(PC_xval$count, exp(pred_integrated[,i])))
    }
    
    # MSE (for count predictions)
    MSE_PConly <- MSE_integrated <- rep(NA, ncol(pred_PConly))
    for (i in 1:length(MSE_PConly)){
      MSE_PConly[i] <- mean((PC_xval$count - exp(pred_PConly[,i]))^2)
      MSE_integrated[i] <- mean((PC_xval$count - exp(pred_integrated[,i]))^2)
    }
    
    # &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    # CALCULATE MEAN FOR EACH PREDICTION, THEN AVERAGE (METHOD 2)
    # &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    
    # AUC (for presence/absence predictions)
    AUC_PConly2 <- as.numeric(auc(PC_xval$presence, apply(1-exp(-exp(pred_PConly)),1,mean)))
    AUC_integrated2 <- as.numeric(auc(PC_xval$presence, apply(1-exp(-exp(pred_integrated)),1,mean)))
    
    # correlation (for count predictions)
    cor_PConly2 <- as.numeric(cor(PC_xval$count, apply(exp(pred_PConly),1,mean)))
    cor_integrated2 <- as.numeric(cor(PC_xval$count, apply(exp(pred_integrated),1,mean)))
    
    # MSE (for count predictions)
    MSE_PConly2 <- mean((PC_xval$count - apply(exp(pred_PConly),1,mean))^2)
    MSE_integrated2 <- mean((PC_xval$count - apply(exp(pred_integrated),1,mean))^2)
    
    # -------------------------------------------------------
    # Save results
    # -------------------------------------------------------
    
    if (file.exists("../AOS_precision/output/xval_df_integrated_RE.RData")){
      load("../AOS_precision/output/xval_df_integrated_RE.RData")
    }
    
    end <- Sys.time() 
    runtime_mins <- difftime( end,start, units="mins")
    
    # Compare crossvalidation metrics
    xval_df <- rbind(xval_df,
                     data.frame(Species = sp_code,
                                xval_fold = fold,
                                
                                # Number of times species was observed in each dataset
                                n_obs_PC = sum(PC_sp$count>0),
                                n_obs_SC = sum(SC_sp$count>0),
                                n_obs_LT = sum(LT_sp$count>0),
                                n_obs_BBA = sum(BBA_sp$count>0),
                                
                                # Amount of data in crossvalidation squares
                                n_obs_PC_xval = sum(PC_xval$count>0),
                                n_obs_SC_xval = sum(SC_sp$count[SC_sp$fold == fold]>0),
                                n_obs_LT_xval = sum(LT_sp$count[LT_sp$fold == fold]>0),
                                n_obs_BBA_xval = sum(BBA_sp$count[BBA_sp$fold == fold]>0),
                                
                                # Log pointwise predictive density
                                lppd_PConly = sum(log(ppd_PConly)),
                                lppd_integrated = sum(log(ppd_integrated)),
                                
                                # Correlation between predictions and observed
                                cor_PConly = mean(cor_PConly),
                                cor_integrated = mean(cor_integrated),
                                
                                # AUC
                                AUC_PConly = mean(AUC_PConly),
                                AUC_integrated = mean(AUC_integrated),
                                
                                MSE_PConly = mean(MSE_PConly),
                                MSE_integrated = mean(MSE_integrated),
                                
                                # Correlation between predictions and observed
                                cor_PConly2 = cor_PConly2,
                                cor_integrated2 = cor_integrated2,
                                
                                # AUC
                                AUC_PConly2 = AUC_PConly2,
                                AUC_integrated2 = AUC_integrated2,
                                
                                # Mean squared error
                                MSE_PConly2 = MSE_PConly2,
                                MSE_integrated2 = MSE_integrated2,
                                runtime_mins = runtime_mins
                     )
    )
    
    print(paste(sp_code," fold ",fold," ... ",round(runtime_mins), "mins"))
    
    save(xval_df, file = "../AOS_precision/output/xval_df_integrated_RE.RData")
    
    rm(list = c("pred_integrated","pred_PConly","fit_PConly","fit_integrated"))
    
  } # close species loop
#} # xval fold

