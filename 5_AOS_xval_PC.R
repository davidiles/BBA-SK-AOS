# ****************************************************************
# ****************************************************************
# CONDUCT CROSS-VALIDATION ANALYSIS ON A SUBSET OF SPECIES
#  - OMITTING ONLY POINT COUNTS FROM 20% OF SQUARES IN EACH CROSSVALIDATION FOLD
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


# ************************************************************************************
# ------------------------------------------------------------------------------------
# CONDUCT ANALYSIS
# ------------------------------------------------------------------------------------
# ************************************************************************************
setwd("C:/Users/ilesd/OneDrive - EC-EC/Iles/Projects/Landbirds/SK_BBA_analysis/AOS_precision/script")
load("../output/AOS_data_package.RData")

xval_PC <- data.frame()
if (file.exists("../output/xval_PC.RData")){
  load("../output/xval_PC.RData")
}

# covariates to include in models
covariates_to_include <- c("PC1","PC2","PC3","Water_5km")

for (sp_code in rev(species_to_fit$sp_code)){
  for (fold in sort(unique(SaskSquares$fold))){
    
    if (file.exists("../output/xval_PC.RData")){
      load("../output/xval_PC.RData")
    }
    
    # Check if this species/xval fold have already been run. If so, skip
    if (nrow(xval_PC)>0){
      if (nrow(subset(xval_PC,Species == sp_code & xval_fold == fold))>0) next
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
    SC_sf <- subset(CL_sf, ProtocolType == "Stationary count")
    
    # Linear transects
    LT_sf <- subset(CL_sf, ProtocolType == "Linear transect")
    
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
    
    #SC_xval <- SC_sp[SC_sp$fold == fold,]
    #SC_sp <- SC_sp[SC_sp$fold != fold,]
    
    #LT_xval <- LT_sp[LT_sp$fold == fold,]
    #LT_sp <- LT_sp[LT_sp$fold != fold,]
    
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
    # SPECIFY AND FIT MODELS 
    # &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    
    bru_options_reset() # bru_options_set(inla.mode = "experimental")
    
    inits <- c(-5,-5,-5,rep(0,length(covariates_to_include))) %>% as.list()
    names(inits) <- c("Intercept_PC","Intercept_SC","Intercept_LT",paste0("Beta1_",covariates_to_include))
    
    # -----------------------------
    # Use only point count data
    # -----------------------------
  #   
  #   model_components = as.formula(paste0('~
  # Intercept_PC(1)+
  # TSS(main = TSS,model = TSS_spde) +
  # kappa_surveyID(surveyID, model = "iid", constr = TRUE, hyper = list(prec = kappa_prec))+
  # kappa_squareID(sq_idx, model = "iid", constr = TRUE, hyper = list(prec = kappa_prec))+
  # kappa_squareday(square_day, model = "iid", constr = TRUE, hyper = list(prec = kappa_prec))+
  # spde_coarse(main = coordinates, model = matern_coarse) +
  #  
  # ',
  #                                        
  #                                        paste0("Beta1_",covariates_to_include,'(1,model="linear", mean.linear = 0, prec.linear = 4)', collapse = " + "))
  #   )
  #   
  #   start <- Sys.time()
  #   fit_PConly <- bru(components = model_components,
  #                     like_PC,
  #                     options = list(
  #                       control.inla = list(int.strategy = "eb"),
  #                       bru_verbose = 4,
  #                       bru_max_iter = 5,
  #                       bru_initial = inits))
  #   end <- Sys.time()
  #   runtime_PConly <- difftime( end,start, units="mins") # 16 min
  #   
    # -----------------------------
    # Use checklists only
    # -----------------------------
    
  #   model_components = as.formula(paste0('~
  # Intercept_SC(1)+
  # Intercept_LT(1)+
  # TSS(main = TSS,model = TSS_spde) +
  # kappa_squareID(sq_idx, model = "iid", constr = TRUE, hyper = list(prec = kappa_prec))+
  # kappa_squareday(square_day, model = "iid", constr = TRUE, hyper = list(prec = kappa_prec))+
  # spde_coarse(main = coordinates, model = matern_coarse) +
  # 
  # SC_duration(main = SC_duration,model = SC_duration_spde) +
  # LT_duration(main = LT_duration,model = LT_duration_spde) +
  # LT_distance(main = LT_distance,model = LT_distance_spde) +
  #   
  # ',
  #                                        
  #                                        paste0("Beta1_",covariates_to_include,'(1,model="linear", mean.linear = 0, prec.linear = 4)', collapse = " + "))
  #   )
  #   
  #   start <- Sys.time()
  #   fit_CLonly <- bru(components = model_components,
  #                     like_SC,like_LT,
  #                     options = list(
  #                       control.inla = list(int.strategy = "eb"),
  #                       bru_verbose = 4,
  #                       bru_max_iter = 5,
  #                       bru_initial = inits))
  #   end <- Sys.time()
  #   runtime_CLonly <- difftime( end,start, units="mins") # 9 min
  #   
    # -----------------------------
    # Use point counts and checklists
    # -----------------------------
    
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
    # 
    # pred_formula_SC = as.formula(paste0(' ~
    # 
    #               Intercept_SC +
    #               TSS +
    #               SC_duration +
    #               spde_coarse +
    #                ',
    #                                     paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + ")))
    # 
    # pred_formula_LT = as.formula(paste0(' ~
    # 
    #               Intercept_LT +
    #               TSS +
    #               LT_duration +
    #               LT_distance +
    #               spde_coarse +
    #                ',
    #                                     paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + ")))
    # 
    # -------------------------
    # Predictions on point count data
    # -------------------------
    # 
    # pred_PC_PConly <- generate(fit_PConly, 
    #                            PC_xval, 
    #                            formula = pred_formula_PC,
    #                            n.samples = nsamp) %>% 
    #   apply(.,1,median)
    
    pred_PC_integrated <- generate(fit_integrated, 
                                   PC_xval, 
                                   formula = pred_formula_PC,
                                   n.samples = nsamp) %>% 
      apply(.,1,median)
    
    # Add lognormal variance corrections, and exponentiate to place on count scale
    #pred_PC_PConly     <- exp(pred_PC_PConly     + 0.5/summary(fit_PConly)$inla$hyperpar["Precision for kappa_surveyID",4]     + 0.5/summary(fit_PConly)$inla$hyperpar["Precision for kappa_squareID",4]     + 0.5/summary(fit_PConly)$inla$hyperpar["Precision for kappa_squareday",4])
    pred_PC_integrated <- exp(pred_PC_integrated + 0.5/summary(fit_integrated)$inla$hyperpar["Precision for kappa_surveyID",4] + 0.5/summary(fit_integrated)$inla$hyperpar["Precision for kappa_squareID",4] + 0.5/summary(fit_integrated)$inla$hyperpar["Precision for kappa_squareday",4])
    
    # # -------------------------
    # # Predictions on stationary checklists
    # # -------------------------
    # 
    # pred_SC_CLonly <- generate(fit_CLonly, 
    #                            SC_xval, 
    #                            formula = pred_formula_SC,
    #                            n.samples = nsamp) %>% 
    #   apply(.,1,median)
    # 
    # pred_SC_integrated <- generate(fit_integrated, 
    #                                SC_xval, 
    #                                formula = pred_formula_SC,
    #                                n.samples = nsamp) %>% 
    #   apply(.,1,median)
    # 
    # # Add lognormal variance corrections, place on count scale
    # pred_SC_CLonly     <- exp(pred_SC_CLonly + 0.5/summary(fit_CLonly)$inla$hyperpar["Precision for kappa_squareID",4]     + 0.5/summary(fit_CLonly)$inla$hyperpar["Precision for kappa_squareday",4])
    # pred_SC_integrated <- exp(pred_SC_integrated + 0.5/summary(fit_integrated)$inla$hyperpar["Precision for kappa_squareID",4] + 0.5/summary(fit_integrated)$inla$hyperpar["Precision for kappa_squareday",4])
    # 
    # 
    # # -------------------------
    # # Predictions on linear transect checklists
    # # -------------------------
    # 
    # pred_LT_CLonly <- generate(fit_CLonly, 
    #                            LT_xval, 
    #                            formula = pred_formula_LT,
    #                            n.samples = nsamp) %>% 
    #   apply(.,1,median)
    # 
    # pred_LT_integrated <- generate(fit_integrated, 
    #                                LT_xval, 
    #                                formula = pred_formula_LT,
    #                                n.samples = nsamp) %>% 
    #   apply(.,1,median)
    # 
    # # Add lognormal variance corrections, place on count scale
    # pred_LT_CLonly     <- exp(pred_LT_CLonly + 0.5/summary(fit_CLonly)$inla$hyperpar["Precision for kappa_squareID",4]     + 0.5/summary(fit_CLonly)$inla$hyperpar["Precision for kappa_squareday",4])
    # pred_LT_integrated <- exp(pred_LT_integrated + 0.5/summary(fit_integrated)$inla$hyperpar["Precision for kappa_squareID",4] + 0.5/summary(fit_integrated)$inla$hyperpar["Precision for kappa_squareday",4])
    # 
    # &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    # Calculate crossvalidation scores
    # &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    
    # -------------------------
    # Predictions on point count data
    # -------------------------
    
    # AUC (for presence/absence predictions)
    #AUC_PC_PConly <- as.numeric(auc(PC_xval$presence, 1-exp(-pred_PC_PConly)))
    AUC_PC_integrated <- as.numeric(auc(PC_xval$presence, 1-exp(-pred_PC_integrated)))
    
    # correlation (for count predictions)
    #cor_PC_PConly <- as.numeric(cor(PC_xval$count, pred_PC_PConly))
    cor_PC_integrated <- as.numeric(cor(PC_xval$count, pred_PC_integrated))
    
    # MSE (for count predictions)
    #MSE_PC_PConly <- mean((PC_xval$count - pred_PC_PConly)^2)
    MSE_PC_integrated <- mean((PC_xval$count - pred_PC_integrated)^2)
    
    # # -------------------------
    # # Predictions on stationary checklists
    # # -------------------------
    # 
    # # AUC (for presence/absence predictions)
    # AUC_SC_CLonly <- as.numeric(auc(SC_xval$presence, 1-exp(-pred_SC_CLonly)))
    # AUC_SC_integrated <- as.numeric(auc(SC_xval$presence, 1-exp(-pred_SC_integrated)))
    # 
    # # correlation (for count predictions)
    # cor_SC_CLonly <- as.numeric(cor(SC_xval$count, pred_SC_CLonly))
    # cor_SC_integrated <- as.numeric(cor(SC_xval$count, pred_SC_integrated))
    # 
    # # MSE (for count predictions)
    # MSE_SC_CLonly <- mean((SC_xval$count - pred_SC_CLonly)^2)
    # MSE_SC_integrated <- mean((SC_xval$count - pred_SC_integrated)^2)
    # 
    # # -------------------------
    # # Predictions on linear transect checklists
    # # -------------------------
    # 
    # # AUC (for presence/absence predictions)
    # AUC_LT_CLonly <- as.numeric(auc(LT_xval$presence, 1-exp(-pred_LT_CLonly)))
    # AUC_LT_integrated <- as.numeric(auc(LT_xval$presence, 1-exp(-pred_LT_integrated)))
    # 
    # # correlation (for count predictions)
    # cor_LT_CLonly <- as.numeric(cor(LT_xval$count, pred_LT_CLonly))
    # cor_LT_integrated <- as.numeric(cor(LT_xval$count, pred_LT_integrated))
    # 
    # # MSE (for count predictions)
    # MSE_LT_CLonly <- mean((LT_xval$count - pred_LT_CLonly)^2)
    # MSE_LT_integrated <- mean((LT_xval$count - pred_LT_integrated)^2)
    # 
    
    # -------------------------------------------------------
    # Save results
    # -------------------------------------------------------
    
    if (file.exists("../output/xval_PC.RData")){
      load("../output/xval_PC.RData")
    }
    
    # Compare crossvalidation metrics
    xval_PC <- rbind(xval_PC,
                                    data.frame(Species = sp_code,
                                               xval_fold = fold,
                                               
                                               # Number of times species was observed in each dataset
                                               n_obs_PC = sum(PC_sp$count>0),
                                               n_obs_SC = sum(SC_sp$count>0),
                                               n_obs_LT = sum(LT_sp$count>0),
                                               
                                               # Amount of data in crossvalidation squares
                                               #n_obs_PC_xval = sum(PC_xval$count>0),
                                               #n_obs_SC_xval = sum(SC_xval$count>0),
                                               #n_obs_LT_xval = sum(LT_xval$count>0),
                                               
                                               # -------------------------------
                                               # Crossvalidation metrics - point counts
                                               # -------------------------------
                                               
                                               # Correlation between predictions and observed
                                               #cor_PC_PConly = cor_PC_PConly,
                                               cor_PC_integrated = cor_PC_integrated,
                                               
                                               # AUC
                                               #AUC_PC_PConly = AUC_PC_PConly,
                                               AUC_PC_integrated = AUC_PC_integrated,
                                               
                                               # Mean squared error
                                               #MSE_PC_PConly = MSE_PC_PConly
                                               MSE_PC_integrated = MSE_PC_integrated
                                               
                                               # # -------------------------------
                                               # # Crossvalidation metrics - point counts
                                               # # -------------------------------
                                               # 
                                               # # Correlation between predictions and observed
                                               # cor_SC_CLonly = cor_SC_CLonly,
                                               # cor_SC_integrated = cor_SC_integrated,
                                               # 
                                               # # AUC
                                               # AUC_SC_CLonly = AUC_SC_CLonly,
                                               # AUC_SC_integrated = AUC_SC_integrated,
                                               # 
                                               # # Mean squared error
                                               # MSE_SC_CLonly = MSE_SC_CLonly,
                                               # MSE_SC_integrated = MSE_SC_integrated,
                                               # 
                                               # # -------------------------------
                                               # # Crossvalidation metrics - point counts
                                               # # -------------------------------
                                               # 
                                               # # Correlation between predictions and observed
                                               # cor_LT_CLonly = cor_LT_CLonly,
                                               # cor_LT_integrated = cor_LT_integrated,
                                               # 
                                               # # AUC
                                               # AUC_LT_CLonly = AUC_LT_CLonly,
                                               # AUC_LT_integrated = AUC_LT_integrated,
                                               # 
                                               # # Mean squared error
                                               # MSE_LT_CLonly = MSE_LT_CLonly,
                                               # MSE_LT_integrated = MSE_LT_integrated
                                               # 
                                    )
    )
    
    save(xval_PC, file = "../output/xval_PC.RData")
    
    rm(list = c("fit_PConly","fit_integrated","fit_CLonly"))
    
  } # close species loop
} # xval fold


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
