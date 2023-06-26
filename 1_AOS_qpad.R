# *******************************************

# Do QPAD offsets improve cross-validation accuracy compared to models that omit them?
#  - Models are fit to point counts only, not checklists

# *******************************************

# ------------------------------------------------
# Load/install packages and set graphical themes / working directory
# ------------------------------------------------

my_packs <- c(
  
  # Data management
  'tidyverse',
  
  # Spatial functions
  'sf','raster','ggspatial','rgeos','fasterize','exactextractr','circular',
  
  # For PCA
  'factoextra',
  
  # Spatial analysis
  'inlabru','ebirdst',
  
  # For plotting
  'viridis','scales','ggpubr','ggtext',
  
  # Cross-validation
  'pROC')

if (any(!my_packs %in% installed.packages()[, 'Package'])) {install.packages(my_packs[which(!my_packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)}
lapply(my_packs, require, character.only = TRUE)

library(INLA) # install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep = TRUE)
library(napops) # For detectability offsets  # devtools::install_github("na-pops/napops") 
# napops::fetch_data()  # - only needs to be run once

rm(list=ls())

# Draw rasters, data, and covariates from "standard analysis"
setwd("D:/Working_Files/1_Projects/Landbirds/SK_BBA_analysis/Standard_Analysis/")
`%!in%` <- Negate(`%in%`)

# ************************************************************************************
# ------------------------------------------------------------------------------------
# Load/prepare data
# ------------------------------------------------------------------------------------
# ************************************************************************************

# --------------------------------
# Load 1 km x 1 km covariate spatial grid (same grid onto which predictions will be made)
# --------------------------------

load("!Data_processed/SaskGrid_PCA.RData")

# Study area boundary
SaskBoundary <- read_sf("!Shapefiles/SaskBoundary/SaskBoundary_Project.shp")

# --------------------------------
# Load Point count data
# --------------------------------

# Pointcount_dataset: locations and covariates associated with each point count
load("!Data_processed/Pointcount_data.RData")

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
load("!Data_processed/Checklist_data.RData")

# Covariates / spatial locations of each checklist survey (daily observations)
DO_surveyinfo <- Checklist_data$DO_surveyinfo

# Counts associated with each checklist survey
DO_matrix <- Checklist_data$DO_matrix

# ******
# Selection criteria for checklists to use in analysis
# *********
DO_to_use <- which(DO_surveyinfo$ProtocolType == "Linear transect" &
                     DO_surveyinfo$TravelDistance_m <= 5000 &
                     DO_surveyinfo$DurationInHours > 0 &
                     DO_surveyinfo$DurationInHours <= 5)

DO_surveyinfo <- DO_surveyinfo[DO_to_use,]
DO_matrix <- DO_matrix[DO_to_use,]

rm(Checklist_data)

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


# ------------------------------------------------
# Identify lists of species to fit models for
# ------------------------------------------------

species_relabund_PC <- colSums(PC_matrix[,-1]>0,na.rm = TRUE) %>% sort(decreasing = TRUE) %>% as.data.frame() %>% rename(PC = 1)
species_relabund_PC$Species <- rownames(species_relabund_PC)
species_relabund_DO <- colSums(DO_matrix[,-1]>0,na.rm = TRUE) %>% sort(decreasing = TRUE) %>% as.data.frame() %>% rename(DO = 1)
species_relabund_DO$Species <- rownames(species_relabund_DO)
species_relabund <- full_join(species_relabund_PC,species_relabund_DO)

# Must have been detected in at least 50 point counts
species_to_fit <- subset(species_relabund, PC >= 50)$Species

# Remove "unknown" species
species_to_fit <- species_to_fit[-which(species_to_fit %in% c("UNDU","UNGU","UNWO","UNYE"))]

length(species_to_fit) # 169 species

# ************************************************************************************
# ------------------------------------------------------------------------------------
# Analysis (note some additional data wrangling is needed on a species-by-species basis)
# ------------------------------------------------------------------------------------
# ************************************************************************************

# Load squares to withhold for crossvalidation
SaskSquares <- read_sf("../AOS_precision/output/SaskSquares_xval.shp") %>%
  st_transform(crs(PC_surveyinfo))


#species_to_fit <- c("RWBL")
waic_df <- xval_df <- data.frame()


for (sp_code in species_to_fit){
  
  if (sp_code %in% waic_df$sp_code) next
  
  print(sp_code)
  
  # If model has already been fit, skip it
  model_filename <- paste0("!Results/Model_Summaries/", sp_code, "_model_summary.RData")
  prediction_filename <- paste0("!Results/Density_Predictions/Prediction_Matrix_", sp_code, ".RData")
  
  # If results have already been generated, skip this species
  if (file.exists(model_filename) & file.exists(prediction_filename)) next
  
  # --------------------------------
  # Prepare data for this species
  # --------------------------------
  
  PC_sf <- PC_surveyinfo %>% mutate(count = NA) # Point counts
  CL_sf <- DO_surveyinfo %>% mutate(count = NA) # Checklist counts (linear transect only, currently)
  
  # Fill in counts
  if (sp_code %in% colnames(PC_matrix)) PC_sf$count <- PC_matrix[,sp_code]
  if (sp_code %in% colnames(DO_matrix)) CL_sf$count <- DO_matrix[,sp_code]
  
  PC_sf <- subset(PC_sf, !is.na(count))
  CL_sf <- subset(CL_sf, !is.na(count)) 
  
  # --------------------------------
  # Checklists will be treated as presence/absence only
  # --------------------------------
  
  CL_sf$presence <- as.numeric(CL_sf$count > 0)
  
  # --------------------------------
  # Check if QPAD offsets exist, and generate them if so
  # --------------------------------
  
  offset_exists <- FALSE
  offset_5min_Pointcount <- 0
  sp_cr = sp_edr = NA
  sp_napops <- subset(napops_species,Species == sp_code)
  if (nrow(sp_napops)>0){
    if (sp_napops$Removal == 1 & sp_napops$Distance == 1){
      offset_exists <- TRUE
      sp_cr <- cue_rate(species = sp_code,od = 153, tssr = 1, model = 1)
      sp_edr <- edr(species = sp_code,road = FALSE, forest = 0.5,model = 1)
      
      # Calculate A and p, which jointly determine offset
      PC_sf$A_metres <- c(pi*sp_edr$EDR_est^2)
      PC_sf$p <- 1-exp(-PC_sf$DurationInMinutes*sp_cr$CR_est[1,1])
      PC_sf$QPAD_offset <- log(PC_sf$A_metres * PC_sf$p)
      
      # Calculate offset for a 5-minute point count (needed later in analysis)
      offset_5min_Pointcount <- c(log((pi*sp_edr$EDR_est^2)*(1-exp(-5*sp_cr$CR_est[1,1]))))
    }
  }
  
  # --------------------------------
  # Withhold data for cross-validation
  # --------------------------------
  
  # Intersect with SaskSquares dataframe
  PC_sf <- st_intersection(PC_sf, SaskSquares)
  CL_sf <- st_intersection(CL_sf, SaskSquares)
  
  # --------------------------------
  # Convert to spatial objects (required by INLA)
  # --------------------------------
  
  PC_sp <- as(PC_sf,'Spatial')
  CL_sp <- as(CL_sf,'Spatial')
  
  PC_xval <- subset(PC_sp, fold == 1)
  CL_xval <- subset(CL_sp, fold == 1)
  
  PC_sp <- subset(PC_sp, fold != 1)
  CL_sp <- subset(CL_sp, fold != 1)
  
  # --------------------------------
  # Create a spatial mesh, which is used to fit the residual spatial field
  # --------------------------------
  
  # In meters
  mesh_cutoff <- 30000         
  max_edge_inner <- 30000   
  max_edge_outer <- 150000
  
  mesh_spatial <- inla.mesh.2d(loc.domain = PC_sp,
                               cutoff = mesh_cutoff, 
                               max.edge = c(max_edge_inner,max_edge_outer))
  
  mesh_locs <- mesh_spatial$loc[,c(1,2)] %>% as.data.frame()
  
  matern_coarse <- inla.spde2.pcmatern(mesh_spatial,
                                       prior.range = c(100000, 0.1), # 10% chance range is smaller than 100000
                                       prior.sigma = c(1, 0.1),      # 10% chance sd is larger than 1
                                       constr = TRUE)                # sum to 1 constraint
  
  # --------------------------------
  # Define square-level random effects
  # --------------------------------
  
  # Point counts
  kappa_PC_prec <- list(prior = "pcprec", param = c(1,0.1))
  PC_sp$sq_idx <- as.numeric(factor(PC_sp$sq_id))
  
  # Linear transects
  kappa_CL_prec <- list(prior = "pcprec", param = c(1,0.1))
  CL_sp$sq_idx <- as.numeric(factor(CL_sp$sq_id))
  
  # --------------------------------
  # Create 'temporal mesh' to model effect of checklist duration
  # --------------------------------
  
  CL_duration_meshpoints <- seq(0,max(CL_sp$DurationInHours)+0.1,length.out = 11)
  CL_duration_mesh1D = inla.mesh.1d(CL_duration_meshpoints,boundary="free")
  CL_duration_spde = inla.spde2.pcmatern(CL_duration_mesh1D,
                                         prior.range = c(1,0.9), # 90% chance range is smaller than 1
                                         prior.sigma = c(2,0.1)) # 10% chance sd is larger than 2
  
  
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
  
  model_components = as.formula(paste0('~
  Intercept_PC(1, model = "linear", mean.linear = -10, prec.linear = 0.01)+
  kappa_PC(sq_idx, model = "iid", constr = TRUE, hyper = list(prec = kappa_PC_prec))+
  spde_coarse(main = coordinates, model = matern_coarse) + 
  
  ',
                                       
                                       paste0("Beta1_",covariates_to_include,'(1,model="linear", mean.linear = 0, prec.linear = 4)', collapse = " + "),
                                       " + ",
                                       paste0("Beta2_",covariates_to_include,'(1,model="linear", mean.linear = 0, prec.linear = 16)', collapse = " + "))
  )
  
  model_components
  
  # --------------------------------
  # Model formulas (for Point count and Checklist likelihoods)
  # --------------------------------
  
  # If species has a QPAD offset calculated, use it.  Otherwise, fit without offsets
  if (!offset_exists) next
  
  
  model_formula_QPAD = as.formula(paste0('count ~
                  QPAD_offset +
                  Intercept_PC +
                  spde_coarse +
                  kappa_PC +
                   ',
                                         paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
                                         " + ",
                                         paste0("Beta2_",covariates_to_include,'*',covariates_to_include,"^2", collapse = " + ")))
  
  model_formula_NULL = as.formula(paste0('count ~
                  Intercept_PC +
                  spde_coarse +
                  kappa_PC +
                   ',
                                         paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
                                         " + ",
                                         paste0("Beta2_",covariates_to_include,'*',covariates_to_include,"^2", collapse = " + ")))
  
  
  # --------------------------------
  # Prediction formula (for creating density maps)
  # --------------------------------
  
  pred_formula = as.formula(paste0(' ~

                  Intercept_PC +
                  spde_coarse +
                   ',
                                   paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + "),
                                   " + ",
                                   paste0("Beta2_",covariates_to_include,'*',covariates_to_include,"^2", collapse = " + ")))
  
  # --------------------------------
  # Specify model likelihoods
  # --------------------------------
  
  like_QPAD <- like(family = "poisson",
                    formula = model_formula_QPAD,
                    data = PC_sp)
  
  like_NULL <- like(family = "poisson",
                    formula = model_formula_NULL,
                    data = PC_sp)
  
  # --------------------------------
  # Select reasonable initial values (should not affect inference, but affects model convergence)
  # --------------------------------
  
  inits <- c(-5,rep(0,length(covariates_to_include)*2)) %>% as.list()
  names(inits) <- c("Intercept_PC",paste0("Beta1_",covariates_to_include), paste0("Beta2_",covariates_to_include))
  
  # --------------------------------
  # Fit model with QPAD offsets
  # --------------------------------
  
  fit_QPAD <- bru(components = model_components, 
                  like_QPAD,
                  options = list(control.compute = list(waic = TRUE, cpo = TRUE, config = TRUE),
                                 bru_verbose = 4,
                                 bru_max_iter = 10,
                                 bru_initial = inits))
  
  as.numeric(sum(fit_QPAD$bru_timings$Time))/60
  
  # --------------------------------
  # Fit model without QPAD offsets
  # --------------------------------
  fit_NULL <- bru(components = model_components, 
                  like_NULL,
                  options = list(control.compute = list(waic = TRUE, cpo = TRUE, config = TRUE),
                                 bru_verbose = 4,
                                 bru_max_iter = 10,
                                 bru_initial = inits))
  
  as.numeric(sum(fit_NULL$bru_timings$Time))/60
  
  # -------------------------------------------------------
  # Store WAIC from each model
  # -------------------------------------------------------
  
  # Compare WAIC for each model
  waic_df <- rbind(waic_df,
                   data.frame(sp_code = sp_code,
                              waic_QPAD = fit_QPAD$waic$waic,
                              waic_NULL = fit_NULL$waic$waic))
  
  # -------------------------------------------------------
  # Store WAIC from each model
  # -------------------------------------------------------
  
  # Compare WAIC for each model
  waic_df <- rbind(waic_df,
                   data.frame(sp_code = sp_code,
                              waic_QPAD = fit_QPAD$waic$waic,
                              waic_NULL = fit_NULL$waic$waic))
  
  # &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  # Cross-validation on withheld data
  # &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

  PC_xval_df <- PC_xval %>%
    as.data.frame() %>% 
    mutate(presence = as.numeric(count>0))
  
  # -------------------------------------------------------
  # Generate spatial predictions from model with QPAD offset, and calculate xval metrics
  # -------------------------------------------------------
  
  start <- Sys.time()
  
  # Generate predictions based on model with QPAD offsets
  pred_QPAD <- generate(fit_QPAD, 
                        PC_xval, 
                        formula = ~ Intercept_PC + spde_coarse + QPAD_offset,
                        n.samples = 500)
  
  # Full distribution of predictions (500 samples) for each withheld data point
  fit_count_QPAD = exp(pred_QPAD)
  fit_presence_QPAD = 1-exp(-fit_count_QPAD)
  
  # Use mean of those distributions to evaluate cross-validation accuracy
  yhat_QPAD <- apply(fit_count_QPAD,1,mean)
  phat_QPAD <- apply(fit_presence_QPAD,1,mean)
  
  AUC_QPAD = auc(PC_xval_df$presence, phat_QPAD) %>% as.numeric()
  MSE_QPAD = mean(( PC_xval_df$count - yhat_QPAD)^2)
  MAE_QPAD = mean(abs( PC_xval_df$count - yhat_QPAD))
  cor_QPAD = cor( PC_xval_df$count, yhat_QPAD) %>% as.numeric()
  logLik_QPAD = sum(log((exp(-yhat_QPAD)*yhat_QPAD^PC_xval_df$count)/factorial(PC_xval_df$count)))
  
  # -------------------------------------------------------
  # Generate spatial predictions from NULL model and calculate xval metrics
  # -------------------------------------------------------
  
  # Generate predictions based on model with NULL offsets
  pred_NULL <- generate(fit_NULL, 
                        PC_xval, 
                        formula = ~ Intercept_PC + spde_coarse,
                        n.samples = 500)
  
  # Full distribution of predictions (500 samples) for each withheld data point
  fit_count_NULL = exp(pred_NULL)
  fit_presence_NULL = 1-exp(-fit_count_NULL)
  
  # Use mean of those distributions to evaluate cross-validation accuracy
  yhat_NULL <- apply(fit_count_NULL,1,mean)
  phat_NULL <- apply(fit_presence_NULL,1,mean)
  
  AUC_NULL = auc(PC_xval_df$presence, apply(fit_presence_NULL,1,mean)) %>% as.numeric()
  MSE_NULL = mean(( PC_xval_df$count - apply(fit_count_NULL,1,mean))^2)
  MAE_NULL = mean(abs( PC_xval_df$count - apply(fit_count_NULL,1,mean)))
  cor_NULL = cor( PC_xval_df$count, apply(fit_count_NULL,1,mean)) %>% as.numeric()
  logLik_NULL = sum(log((exp(-yhat_NULL)*yhat_NULL^PC_xval_df$count)/factorial(PC_xval_df$count)))
  
  end <- Sys.time()
  
  end-start
  
  # Compare WAIC for each model
  xval_df <- rbind(xval_df,
                   data.frame(sp_code = sp_code,
                              
                              # Model fit with QPAD offset
                              AUC_QPAD = AUC_QPAD,
                              MSE_QPAD = MSE_QPAD,
                              MAE_QPAD = MAE_QPAD,
                              cor_QPAD = cor_QPAD,
                              logLik_QPAD = logLik_QPAD,
                              
                              # Model fit with NULL offset
                              AUC_NULL = AUC_NULL,
                              MSE_NULL = MSE_NULL,
                              MAE_NULL = MAE_NULL,
                              cor_NULL = cor_NULL,
                              logLik_NULL = logLik_NULL
                              
                              )
                   )
  
  print(sp_code)
  
  save(xval_df, file = "../AOS_precision/output/xval_df_QPAD.RData")
  save(waic_df, file = "../AOS_precision/output/waic_df_QPAD.RData")
  
} # close species loop


# --------------------------------
# Compare crossvalidation accuracy
# --------------------------------

# Difference in log likelihood
xval_df$delta_LL <- xval_df$logLik_QPAD - xval_df$logLik_NULL

hist(xval_df$delta_LL)

mean(xval_df$delta_LL)
median(xval_df$delta_LL)

mean(xval_df$delta_LL>0)
