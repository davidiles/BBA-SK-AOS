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


# ************************************************************************************
# ------------------------------------------------------------------------------------
# CONDUCT ANALYSIS
# ------------------------------------------------------------------------------------
# ************************************************************************************
defaultW <- getOption("warn") 
options(warn = -1) # Suppress warnings

setwd("D:/Working_Files/1_Projects/Landbirds/SK_BBA_analysis/AOS_precision/script")
load("../output/AOS_data_package.RData")


xval_PC_CL_subset <- data.frame()
if (file.exists("../output/xval_PC_CL_subset2.RData")){
  load("../output/xval_PC_CL_subset2.RData")
}

# covariates to include in models
covariates_to_include <- c("PC1","PC2","PC3","Water_5km")

set.seed(999)
species_to_fit <- subset(species_distribution_summary, n_sq >= 200 & n_PC >= 100 & n_CL >= 100) %>%
  sample_n(30) %>%
  arrange(n_PC)

 for (sp_code in species_to_fit$sp_code){
   for (fold in 1:5){
     
    if (file.exists("../output/xval_PC_CL_subset2.RData")){
      load("../output/xval_PC_CL_subset2.RData")
    }
    
    # Check if this species/xval fold have already been run. If so, skip
    if (nrow(xval_PC_CL_subset)>0){
      if (nrow(subset(xval_PC_CL_subset,Species == sp_code & xval_fold == fold))>0) next
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
    TSS_meshpoints <- seq(TSS_range[1]-0.1,TSS_range[2]+0.1,length.out = 11)
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
    # Withhold checklist data in validation squares as well
    # --------------------------------
    
    SC_xval <- SC_sp[SC_sp$fold == fold,]
    SC_sp <- SC_sp[SC_sp$fold != fold,]
    
    LT_xval <- LT_sp[LT_sp$fold == fold,]
    LT_sp <- LT_sp[LT_sp$fold != fold,]
    
    # --------------------------------
    # Calculate number of point counts and checklists per non-validation square
    # --------------------------------
    
    n_samples_PC <- PC_sp %>% as.data.frame() %>% group_by(sq_id) %>% summarize(n_PC = n())
    n_samples_SC <- SC_sp %>% as.data.frame() %>% group_by(sq_id) %>% summarize(n_SC = n())
    n_samples_LT <- LT_sp %>% as.data.frame() %>% group_by(sq_id) %>% summarize(n_LT = n())
    
    n_samples <- full_join(n_samples_PC,n_samples_SC) %>%
      rowwise() %>%
      mutate(n_min = min(c(n_PC,n_SC))) %>%
      na.omit()
    
    # Ensure equal sample sizes
    PC_sp <- PC_sp %>% 
      st_as_sf() %>%
      subset(sq_id %in% n_samples$sq_id) %>%
      full_join(n_samples) %>%
      group_by(sq_id) %>%
      mutate(samp = sample(n())) %>%
      filter(samp <= n_min) %>%
      as('Spatial')
    
    SC_sp <- SC_sp %>% 
      st_as_sf() %>%
      subset(sq_id %in% n_samples$sq_id) %>%
      full_join(n_samples) %>%
      group_by(sq_id) %>%
      mutate(samp = sample(n())) %>%
      filter(samp <= n_min) %>%
      as('Spatial')
    
    # --------------------------------
    # Select a random 50% of remaining squares to withhold either PC or CL data during fitting
    # --------------------------------
    
    set.seed(fold*100)
    
    squares_to_withhold <- SaskSquares[SaskSquares$fold != fold,]$SQUARE_ID
    squares_to_withhold <- sample(squares_to_withhold,round(length(squares_to_withhold)/2))
    
    # &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    # SPECIFY AND FIT MODELS 
    # &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    
    bru_options_reset() # bru_options_set(inla.mode = "experimental")
    
    inits <- c(-5,-5,rep(0,length(covariates_to_include))) %>% as.list()
    names(inits) <- c("Intercept_PC","Intercept_SC",paste0("Beta1_",covariates_to_include))
    
    # -----------------------------
    # FIT MODEL TO 50% OF POINT COUNT DATA
    # -----------------------------
    
    model_components = as.formula(paste0('~
  Intercept_PC(1)+
  TSS(main = TSS,model = TSS_spde) +
  kappa_surveyID(surveyID, model = "iid", constr = TRUE, hyper = list(prec = kappa_prec))+
  kappa_squareID(sq_idx, model = "iid", constr = TRUE, hyper = list(prec = kappa_prec))+
  kappa_squareday(square_day, model = "iid", constr = TRUE, hyper = list(prec = kappa_prec))+
  spde_coarse(main = coordinates, model = matern_coarse) +
   
  ',
                                         
                                         paste0("Beta1_",covariates_to_include,'(1,model="linear", mean.linear = 0, prec.linear = 4)', collapse = " + "))
    )
    
    start <- Sys.time()
    fit_PConly_50 <- bru(components = model_components,
                         
                         like(family = "poisson",
                              formula = model_formula_PC,
                              data = PC_sp[PC_sp$sq_id %!in% squares_to_withhold,]),
                         
                         options = list(
                           control.inla = list(int.strategy = "eb"),
                           bru_verbose = 4,
                           bru_max_iter = 5,
                           bru_initial = inits))
    end <- Sys.time()
    runtime_PConly_50 <- difftime( end,start, units="mins")
    
    # -----------------------------
    # FIT MODEL TO 100% OF POINT COUNT DATA
    # -----------------------------
    
    model_components = as.formula(paste0('~
  Intercept_PC(1)+
  TSS(main = TSS,model = TSS_spde) +
  kappa_surveyID(surveyID, model = "iid", constr = TRUE, hyper = list(prec = kappa_prec))+
  kappa_squareID(sq_idx, model = "iid", constr = TRUE, hyper = list(prec = kappa_prec))+
  kappa_squareday(square_day, model = "iid", constr = TRUE, hyper = list(prec = kappa_prec))+
  spde_coarse(main = coordinates, model = matern_coarse) +
   
  ',
                                         
                                         paste0("Beta1_",covariates_to_include,'(1,model="linear", mean.linear = 0, prec.linear = 4)', collapse = " + "))
    )
    
    start <- Sys.time()
    fit_PConly_100 <- bru(components = model_components,
                          
                          like(family = "poisson",
                               formula = model_formula_PC,
                               data = PC_sp),
                          
                          options = list(
                            control.inla = list(int.strategy = "eb"),
                            bru_verbose = 4,
                            bru_max_iter = 5,
                            bru_initial = inits))
    end <- Sys.time()
    runtime_PConly_100 <- difftime( end,start, units="mins")
    
    # -----------------------------
    # FIT MODEL TO 50% OF POINT COUNT DATA + 100% OF CHECKLIST DATA
    # -----------------------------
    
    model_components = as.formula(paste0('~
  Intercept_PC(1)+
  Intercept_SC(1)+
  TSS(main = TSS,model = TSS_spde) +
  kappa_surveyID(surveyID, model = "iid", constr = TRUE, hyper = list(prec = kappa_prec))+
  kappa_squareID(sq_idx, model = "iid", constr = TRUE, hyper = list(prec = kappa_prec))+
  kappa_squareday(square_day, model = "iid", constr = TRUE, hyper = list(prec = kappa_prec))+
  spde_coarse(main = coordinates, model = matern_coarse) +
  SC_duration(main = SC_duration,model = SC_duration_spde) +
    
  ',
                                         
                                         paste0("Beta1_",covariates_to_include,'(1,model="linear", mean.linear = 0, prec.linear = 4)', collapse = " + "))
    )
    
    start <- Sys.time()
    fit_integrated_50_50 <- bru(components = model_components,
                                
                                like(family = "poisson",
                                     formula = model_formula_PC,
                                     data = PC_sp[PC_sp$sq_id %!in% squares_to_withhold,]),
                                
                                like(family = "binomial",
                                     formula = model_formula_SC,
                                     data = SC_sp[SC_sp$sq_id %!in% squares_to_withhold,]),
                                
                                options = list(
                                  control.inla = list(int.strategy = "eb"),
                                  bru_verbose = 4,
                                  bru_max_iter = 5,
                                  bru_initial = inits))
    end <- Sys.time()
    runtime_integrated_50_50 <- difftime( end,start, units="mins")
    
    start <- Sys.time()
    fit_integrated_50_100 <- bru(components = model_components,
                                 
                                 like(family = "poisson",
                                      formula = model_formula_PC,
                                      data = PC_sp[PC_sp$sq_id %!in% squares_to_withhold,]),
                                 
                                 like(family = "binomial",
                                      formula = model_formula_SC,
                                      data = SC_sp),
                                 options = list(
                                   control.inla = list(int.strategy = "eb"),
                                   bru_verbose = 4,
                                   bru_max_iter = 5,
                                   bru_initial = inits))
    end <- Sys.time()
    runtime_integrated_50_100 <- difftime( end,start, units="mins")
    
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
    
    # -------------------------
    # Predictions on point count data
    # -------------------------
    
    pred_PConly_50 <- generate(fit_PConly_50, 
                               PC_xval, 
                               formula = pred_formula_PC,
                               n.samples = nsamp)
    
    pred_PConly_100 <- generate(fit_PConly_100,
                                PC_xval,
                                formula = pred_formula_PC,
                                n.samples = nsamp)
    
    pred_integrated_50_50 <- generate(fit_integrated_50_50, 
                                      PC_xval, 
                                      formula = pred_formula_PC,
                                      n.samples = nsamp)
    
    pred_integrated_50_100 <- generate(fit_integrated_50_100, 
                                       PC_xval, 
                                       formula = pred_formula_PC,
                                       n.samples = nsamp)
    
    # Add lognormal variance corrections, and exponentiate to place on count scale
    pred_PConly_50     <- exp(pred_PConly_50     + 0.5/summary(fit_PConly_50)$inla$hyperpar["Precision for kappa_surveyID",4]     + 0.5/summary(fit_PConly_50)$inla$hyperpar["Precision for kappa_squareID",4]     + 0.5/summary(fit_PConly_50)$inla$hyperpar["Precision for kappa_squareday",4])
    pred_PConly_100     <- exp(pred_PConly_100     + 0.5/summary(fit_PConly_100)$inla$hyperpar["Precision for kappa_surveyID",4]     + 0.5/summary(fit_PConly_100)$inla$hyperpar["Precision for kappa_squareID",4]     + 0.5/summary(fit_PConly_100)$inla$hyperpar["Precision for kappa_squareday",4])
    pred_integrated_50_50 <- exp(pred_integrated_50_50 + 0.5/summary(fit_integrated_50_50)$inla$hyperpar["Precision for kappa_surveyID",4] + 0.5/summary(fit_integrated_50_50)$inla$hyperpar["Precision for kappa_squareID",4] + 0.5/summary(fit_integrated_50_50)$inla$hyperpar["Precision for kappa_squareday",4])
    pred_integrated_50_100 <- exp(pred_integrated_50_100 + 0.5/summary(fit_integrated_50_100)$inla$hyperpar["Precision for kappa_surveyID",4] + 0.5/summary(fit_integrated_50_100)$inla$hyperpar["Precision for kappa_squareID",4] + 0.5/summary(fit_integrated_50_100)$inla$hyperpar["Precision for kappa_squareday",4])
    
    
    # &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    # Calculate crossvalidation scores
    # &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    
    cor_PConly_50 <- cor_PConly_100 <- cor_integrated_50_50 <- cor_integrated_50_100 <- rep(NA,nsamp)
    AUC_PConly_50 <- AUC_PConly_100 <- AUC_integrated_50_50 <- AUC_integrated_50_100 <- rep(NA,nsamp)
    MSE_PConly_50 <- MSE_PConly_100 <- MSE_integrated_50_50 <- MSE_integrated_50_100 <- rep(NA,nsamp)
    lppd_PConly_50 <- lppd_PConly_100 <- lppd_integrated_50_50 <- lppd_integrated_50_100 <- rep(NA,nsamp)
    
      for (i in 1:nsamp){
        cor_PConly_50[i] <- as.numeric(cor(PC_xval$count, pred_PConly_50[,i]))
        cor_PConly_100[i] <- as.numeric(cor(PC_xval$count, pred_PConly_100[,i]))
        cor_integrated_50_50[i] <- as.numeric(cor(PC_xval$count, pred_integrated_50_50[,i]))
        cor_integrated_50_100[i] <- as.numeric(cor(PC_xval$count, pred_integrated_50_100[,i]))
        
        if (length(unique(PC_xval$presence))>1){
          suppressMessages({
          AUC_PConly_50[i] <- as.numeric(auc(PC_xval$presence, 1-exp(-pred_PConly_50[,i])))
          AUC_PConly_100[i] <- as.numeric(auc(PC_xval$presence, 1-exp(-pred_PConly_100[,i])))
          AUC_integrated_50_50[i] <- as.numeric(auc(PC_xval$presence, 1-exp(-pred_integrated_50_50[,i])))
          AUC_integrated_50_100[i] <- as.numeric(auc(PC_xval$presence, 1-exp(-pred_integrated_50_100[,i])))
          })
        }
        
        MSE_PConly_50[i] <- mean((PC_xval$count - pred_PConly_50[,i])^2)
        MSE_PConly_100[i] <- mean((PC_xval$count - pred_PConly_100[,i])^2)
        MSE_integrated_50_50[i] <- mean((PC_xval$count - pred_integrated_50_50[,i])^2)
        MSE_integrated_50_100[i] <- mean((PC_xval$count - pred_integrated_50_100[,i])^2)
        
        # Due to numerical overflow, set minimum log(density) to -1000
        lppd_PConly_50 <- log(dpois(PC_xval$count,pred_PConly_50[,i]))
        lppd_PConly_100 <- log(dpois(PC_xval$count,pred_PConly_100[,i]))
        lppd_integrated_50_50 <- log(dpois(PC_xval$count,pred_integrated_50_50[,i]))
        lppd_integrated_50_100 <- log(dpois(PC_xval$count,pred_integrated_50_100[,i]))
        
        lppd_PConly_50[lppd_PConly_50 == -Inf] = -1000
        lppd_PConly_100[lppd_PConly_100 == -Inf] = -1000
        lppd_integrated_50_50[lppd_integrated_50_50 == -Inf] = -1000
        lppd_integrated_50_100[lppd_integrated_50_100 == -Inf] = -1000
        
        lppd_PConly_50[i] <- sum(lppd_PConly_50)
        lppd_PConly_100[i] <- sum(lppd_PConly_100)
        lppd_integrated_50_50[i] <- sum(lppd_integrated_50_50)
        lppd_integrated_50_100[i] <- sum(lppd_integrated_50_100)
        
      }
  
    # &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    # Save results
    # &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    
    if (file.exists("../output/xval_PC_CL_subset2.RData")){
      load("../output/xval_PC_CL_subset2.RData")
    }
    
    n_sq_det_PC <- PC_sp %>% 
      as.data.frame() %>%
      subset(count>0) %>%
      summarize(n_sq_det = length(unique(sq_id)))
    
    n_sq_det_SC <- SC_sp %>% 
      as.data.frame() %>%
      subset(count>0) %>%
      summarize(n_sq_det = length(unique(sq_id)))
    
    n_det_PC <- PC_sp %>% 
      as.data.frame() %>%
      summarize(n_det = sum(count>0))
    
    n_det_SC <- SC_sp %>% 
      as.data.frame() %>%
      summarize(n_det = sum(count>0))
    
    
    # Compare crossvalidation metrics
    xval_PC_CL_subset <- rbind(xval_PC_CL_subset,
                               data.frame(Species = sp_code,
                                          xval_fold = fold,
                                          
                                          # Number of squares in 'full' dataset
                                          n_sq = length(unique(PC_sp$sq_id)),
                                          
                                          # Number of squares in which species was detected in full PC dataset
                                          n_sq_det_PC = n_sq_det_PC$n_sq_det,
                                          # Number of detections in full PC dataset
                                          n_det_PC = n_det_PC$n_det,
                                          
                                          # Number of squares in which species was detected in full SC dataset
                                          n_sq_det_SC = n_sq_det_SC$n_sq_det,
                                          # Number of detections in full SC dataset
                                          n_det_SC = n_det_SC$n_det,
                                          
                                          # Number of detections in crossvalidation squares
                                          n_det_PC_xval = sum(PC_xval$count>0),
                                          n_det_SC_xval = sum(SC_xval$count>0),
                                          
                                          # -------------------------------
                                          # Crossvalidation metrics
                                          # -------------------------------
                                          
                                          # Correlation between predictions and observed
                                          cor_PConly_50 = mean(cor_PConly_50),
                                          cor_PConly_100 = mean(cor_PConly_100),
                                          cor_integrated_50_50 = mean(cor_integrated_50_50),
                                          cor_integrated_50_100 = mean(cor_integrated_50_100),
                                          
                                          # AUC
                                          AUC_PConly_50 = mean(AUC_PConly_50),
                                          AUC_PConly_100 = mean(AUC_PConly_100),
                                          AUC_integrated_50_50 = mean(AUC_integrated_50_50),
                                          AUC_integrated_50_100 = mean(AUC_integrated_50_100),
                                          
                                          # Mean squared error
                                          MSE_PConly_50 = mean(MSE_PConly_50),
                                          MSE_PConly_100 = mean(MSE_PConly_100),
                                          MSE_integrated_50_50 = mean(MSE_integrated_50_50),
                                          MSE_integrated_50_100 = mean(MSE_integrated_50_100),
                                          
                                          # LPPD
                                          lppd_PConly_50 = mean(lppd_PConly_50),
                                          lppd_PConly_100 = mean(lppd_PConly_100),
                                          lppd_integrated_50_50 = mean(lppd_integrated_50_50),
                                          lppd_integrated_50_100 = mean(lppd_integrated_50_100)
                                          
                               )
    )
    
    if (min(xval_PC_CL_subset[,c("lppd_PConly_50","lppd_PConly_100","lppd_integrated_50_50","lppd_integrated_50_100")]) == "-Inf") break
    
    save(xval_PC_CL_subset, file = "../output/xval_PC_CL_subset2.RData")
    
    rm(list = c("fit_PConly_50","fit_PConly_100","fit_integrated_50_50","fit_integrated_50_100",
                "pred_PConly_50","pred_PConly_100","pred_integrated_50_50","pred_integrated_50_100"))
    
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
