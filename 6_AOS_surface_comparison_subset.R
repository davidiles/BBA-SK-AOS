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

rm(list=ls())

defaultW <- getOption("warn") 
options(warn = -1) # Suppress warnings

setwd("D:/Working_Files/1_Projects/Landbirds/SK_BBA_analysis/AOS_precision/script")
load("../output/AOS_data_package.RData")


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


# *************************************************************************
# *************************************************************************
# CONDUCT ANALYSIS
# *************************************************************************
# *************************************************************************
sp_code = "BRTH"
fold = 3

# covariates to include in models
covariates_to_include <- c("PC1","PC2","PC3","Water_5km")

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
# Generate maps
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

# ------------------------------------------------
# Process species names / labels (from Birds Canada)
# ------------------------------------------------

Sask_spcd <- read.csv("../../Standard_Analysis/!Data/SaskAtlas_Species_2022-05-16.csv", stringsAsFactors = F)
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
species_name = Sask_spcd$CommonName[which(Sask_spcd$spcd == sp_code)]
species_label = Sask_spcd$Label[which(Sask_spcd$spcd == sp_code)]

nsamp = 500
pred_grid_sp <- SaskGrid

pred_formula_PC = as.formula(paste0(' ~

                  Intercept_PC +
                  TSS +
                  spde_coarse +
                   ',
                                    paste0("Beta1_",covariates_to_include,'*',covariates_to_include, collapse = " + ")))



pred_PConly_50 <- generate(fit_PConly_50,
                                as(SaskGrid,'Spatial'),
                                formula =  pred_formula_PC,
                                n.samples = nsamp) %>% apply(.,1,median)

pred_PConly_100 <- generate(fit_PConly_100,
                           as(SaskGrid,'Spatial'),
                           formula =  pred_formula_PC,
                           n.samples = nsamp) %>% apply(.,1,median)

pred_integrated_50_50 <- generate(fit_integrated_50_50,
                           as(SaskGrid,'Spatial'),
                           formula =  pred_formula_PC,
                           n.samples = nsamp) %>% apply(.,1,median)

pred_integrated_50_100 <- generate(fit_integrated_50_100,
                           as(SaskGrid,'Spatial'),
                           formula =  pred_formula_PC,
                           n.samples = nsamp) %>% apply(.,1,median)

# Lognormal variance corrections
pred_grid_sp$pred_PConly_50     <- exp(pred_PConly_50     + 0.5/summary(fit_PConly_50)$inla$hyperpar["Precision for kappa_surveyID",4]     + 0.5/summary(fit_PConly_50)$inla$hyperpar["Precision for kappa_squareID",4]     + 0.5/summary(fit_PConly_50)$inla$hyperpar["Precision for kappa_squareday",4])
pred_grid_sp$pred_PConly_100     <- exp(pred_PConly_100     + 0.5/summary(fit_PConly_100)$inla$hyperpar["Precision for kappa_surveyID",4]     + 0.5/summary(fit_PConly_100)$inla$hyperpar["Precision for kappa_squareID",4]     + 0.5/summary(fit_PConly_100)$inla$hyperpar["Precision for kappa_squareday",4])
pred_grid_sp$pred_integrated_50_50     <- exp(pred_integrated_50_50     + 0.5/summary(fit_PConly_50)$inla$hyperpar["Precision for kappa_surveyID",4]     + 0.5/summary(fit_PConly_50)$inla$hyperpar["Precision for kappa_squareID",4]     + 0.5/summary(fit_PConly_50)$inla$hyperpar["Precision for kappa_squareday",4])
pred_grid_sp$pred_integrated_50_100     <- exp(pred_integrated_50_100     + 0.5/summary(fit_PConly_50)$inla$hyperpar["Precision for kappa_surveyID",4]     + 0.5/summary(fit_PConly_50)$inla$hyperpar["Precision for kappa_squareID",4]     + 0.5/summary(fit_PConly_50)$inla$hyperpar["Precision for kappa_squareday",4])

lower_bound <- 0.01
upper_bound <- quantile(pred_grid_sp$pred_PConly_100,0.99,na.rm = TRUE) %>% signif(2)
if (lower_bound >= upper_bound){
  upper_bound <- quantile(pred_grid_sp$pred_PConly_100,0.99,na.rm = TRUE) %>% signif(2)
  lower_bound <- (upper_bound/5) %>% signif(2)
}


# Create rasters
raster_PConly_50 <- cut.fn(df = pred_grid_sp,
                                      column_name = "pred_PConly_50",
                            lower_bound = lower_bound,
                            upper_bound = upper_bound)

raster_PConly_100 <- cut.fn(df = pred_grid_sp,
                            column_name = "pred_PConly_100",
                            lower_bound = lower_bound,
                            upper_bound = upper_bound)

raster_integrated_50_50 <- cut.fn(df = pred_grid_sp,
                             column_name = "pred_integrated_50_50",
                             lower_bound = lower_bound,
                             upper_bound = upper_bound)

raster_integrated_50_100 <- cut.fn(df = pred_grid_sp,
                                   column_name = "pred_integrated_50_100",
                                   lower_bound = lower_bound,
                                   upper_bound = upper_bound)

# ---------------------------------------------
# Prepare to plot locations where species was observed
# ---------------------------------------------

PC_sf <- PC_surveyinfo %>% mutate(count = NA)
if (sp_code %in% colnames(PC_matrix)) PC_sf$count <- PC_matrix[,sp_code]
PC_sf <- subset(PC_sf, !is.na(count))

CL_sf <- DO_surveyinfo %>% mutate(count = NA) # Checklist counts
if (sp_code %in% colnames(DO_matrix)) CL_sf$count <- DO_matrix[,sp_code]
CL_sf$presence <- as.numeric(CL_sf$count > 0)
CL_sf$surveyID <- 1:nrow(CL_sf)

# Remove NAs
CL_sf <- subset(CL_sf, !is.na(count)) 

# Intersect with SaskSquares dataframe
PC_sf <- st_intersection(PC_sf, SaskSquares)
CL_sf <- st_intersection(CL_sf, SaskSquares)

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
# Generate plots
# ---------------------------------------------

pred_map_PConly_50 <- ggplot() +
  
  geom_sf(data = SaskBoundary,colour="black",fill="#f2f2f2",lwd=0.3,show.legend = F) +
  geom_raster(data = raster_PConly_50 , aes(x = x, y = y, fill = layer)) +
  scale_fill_manual(name = "<span style='font-size:13pt'>Relative Abundance</span><br><span style='font-size:7pt'>Per point count</span><br><span style='font-size:7pt'>(Posterior Median)</span>",
                    values = colpal_relabund(length(levels(raster_PConly_50$layer))), drop=FALSE)+
  
  geom_sf(data = SaskWater,colour=NA,fill="#59F3F3",show.legend = F)+
  geom_sf(data = SaskBoundary,colour="black",fill=NA,lwd=0.3,show.legend = F)+
  
  coord_sf(clip = "off",xlim = c(min(SaskPoints$x), max(SaskPoints$x)))+
  theme(panel.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
  annotate(geom="text",x=346000,y=960000, label= paste0(species_label),lineheight = .85,hjust = 0,size=6,fontface =2) +
  annotate(geom="text",x=346000,y=550000, label= paste0("Prepared on ",Sys.Date()),size=3,lineheight = .75,hjust = 0,color="#3b3b3b")+
  theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(5,10,5,-20),legend.title.align=0.5,
        legend.title = element_markdown(lineheight=.9,hjust = "left"))

png(paste0("../output/figures/xval_maps/xval_",sp_code,"_fold",fold,"_PConly_50.png"), width=6.5, height=8, units="in", res=300, type="cairo")
print(pred_map_PConly_50)
dev.off()

pred_map_PConly_100 <- ggplot() +
  
  geom_sf(data = SaskBoundary,colour="black",fill="#f2f2f2",lwd=0.3,show.legend = F) +
  geom_raster(data = raster_PConly_100 , aes(x = x, y = y, fill = layer)) +
  scale_fill_manual(name = "<span style='font-size:13pt'>Relative Abundance</span><br><span style='font-size:7pt'>Per point count</span><br><span style='font-size:7pt'>(Posterior Median)</span>",
                    values = colpal_relabund(length(levels(raster_PConly_100$layer))), drop=FALSE)+
  
  geom_sf(data = SaskWater,colour=NA,fill="#59F3F3",show.legend = F)+
  geom_sf(data = SaskBoundary,colour="black",fill=NA,lwd=0.3,show.legend = F)+
  
  coord_sf(clip = "off",xlim = c(min(SaskPoints$x), max(SaskPoints$x)))+
  theme(panel.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
  annotate(geom="text",x=346000,y=960000, label= paste0(species_label),lineheight = .85,hjust = 0,size=6,fontface =2) +
  annotate(geom="text",x=346000,y=550000, label= paste0("Prepared on ",Sys.Date()),size=3,lineheight = .75,hjust = 0,color="#3b3b3b")+
  theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(5,10,5,-20),legend.title.align=0.5,
        legend.title = element_markdown(lineheight=.9,hjust = "left"))

png(paste0("../output/figures/xval_maps/xval_",sp_code,"_fold",fold,"_PConly_100.png"), width=6.5, height=8, units="in", res=300, type="cairo")
print(pred_map_PConly_100)
dev.off()

pred_map_integrated_50_50 <- ggplot() +
  
  geom_sf(data = SaskBoundary,colour="black",fill="#f2f2f2",lwd=0.3,show.legend = F) +
  geom_raster(data = raster_integrated_50_50 , aes(x = x, y = y, fill = layer)) +
  scale_fill_manual(name = "<span style='font-size:13pt'>Relative Abundance</span><br><span style='font-size:7pt'>Per point count</span><br><span style='font-size:7pt'>(Posterior Median)</span>",
                    values = colpal_relabund(length(levels(raster_integrated_50_50$layer))), drop=FALSE)+
  
  geom_sf(data = SaskWater,colour=NA,fill="#59F3F3",show.legend = F)+
  geom_sf(data = SaskBoundary,colour="black",fill=NA,lwd=0.3,show.legend = F)+
  
  coord_sf(clip = "off",xlim = c(min(SaskPoints$x), max(SaskPoints$x)))+
  theme(panel.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
  annotate(geom="text",x=346000,y=960000, label= paste0(species_label),lineheight = .85,hjust = 0,size=6,fontface =2) +
  annotate(geom="text",x=346000,y=550000, label= paste0("Prepared on ",Sys.Date()),size=3,lineheight = .75,hjust = 0,color="#3b3b3b")+
  theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(5,10,5,-20),legend.title.align=0.5,
        legend.title = element_markdown(lineheight=.9,hjust = "left"))

png(paste0("../output/figures/xval_maps/xval_",sp_code,"_fold",fold,"_integrated_50_50.png"), width=6.5, height=8, units="in", res=300, type="cairo")
print(pred_map_integrated_50_50)
dev.off()

pred_map_integrated_50_100 <- ggplot() +
  
  geom_sf(data = SaskBoundary,colour="black",fill="#f2f2f2",lwd=0.3,show.legend = F) +
  geom_raster(data = raster_integrated_50_100 , aes(x = x, y = y, fill = layer)) +
  scale_fill_manual(name = "<span style='font-size:13pt'>Relative Abundance</span><br><span style='font-size:7pt'>Per point count</span><br><span style='font-size:7pt'>(Posterior Median)</span>",
                    values = colpal_relabund(length(levels(raster_integrated_50_100$layer))), drop=FALSE)+
  
  geom_sf(data = SaskWater,colour=NA,fill="#59F3F3",show.legend = F)+
  geom_sf(data = SaskBoundary,colour="black",fill=NA,lwd=0.3,show.legend = F)+
  
  coord_sf(clip = "off",xlim = c(min(SaskPoints$x), max(SaskPoints$x)))+
  theme(panel.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
  annotate(geom="text",x=346000,y=960000, label= paste0(species_label),lineheight = .85,hjust = 0,size=6,fontface =2) +
  annotate(geom="text",x=346000,y=550000, label= paste0("Prepared on ",Sys.Date()),size=3,lineheight = .75,hjust = 0,color="#3b3b3b")+
  theme(legend.margin=margin(0,0,0,0),legend.box.margin=margin(5,10,5,-20),legend.title.align=0.5,
        legend.title = element_markdown(lineheight=.9,hjust = "left"))

png(paste0("../output/figures/xval_maps/xval_",sp_code,"_fold",fold,"_integrated_50_100.png"), width=6.5, height=8, units="in", res=300, type="cairo")
print(pred_map_integrated_50_100)
dev.off()
