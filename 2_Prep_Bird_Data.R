# ************************************************************************
# ************************************************************************
# Code to format bird data from point counts and atlas checklists
# ************************************************************************
# ************************************************************************

# ------------------------------------------------------------
# Load/install packages and set graphical themes / working directory
# ------------------------------------------------------------

my_packs <- c('tidyverse','sf','raster','CircStats','ggspatial','spsurvey',
              'lubridate','photobiology')

if (any(!my_packs %in% installed.packages()[, 'Package'])) {install.packages(my_packs[which(!my_packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)}
lapply(my_packs, require, character.only = TRUE)

# load naturecounts r package (install if necessary)
library(naturecounts) # remotes::install_github("BirdsCanada/naturecounts")

rm(list=ls())

setwd("D:/Working_Files/1_Projects/Landbirds/SK_BBA_analysis/Standard_Analysis/")
`%!in%` <- Negate(`%in%`)


# Acceptable days of year for surveys
start_date <- lubridate::ymd("2022-05-31") %>% yday()
end_date <- lubridate::ymd("2022-07-07") %>% yday()

# NOTE: Surveys will only be included if they were conducted 30 min before sunrise to 4 hours after sunrise

# ------------------------------------------------------------
# Pull species names from naturecounts
# ------------------------------------------------------------

species <- search_species_code()

# Remove duplicates and remove incorrect codes
species[which(duplicated(species$species_id)),]
species <- subset(species, !(BSCDATA %in% c("CAGO","UNDD","GRPA","ASTK","MONP","COSN","RWBU","PEPL","UNSG","GUSP","RODO","ECDO","WPWI","TTWO","GRAJ","SIFL","HYBG","YWAR","NSTS","OTHR","DOWS","EWWR","UNWX","COMO","SBDU","JUNI","RIPH","BNOW","BDOW")))

species_vec <- sort(unique(species$BSCDATA))

# ************************************************************
# ************************************************************
# PART 1: PREPARE POINT COUNT DATA
#
# This script creats two objects:
#   1) PC_surveyinfo - an sf object with "survey information" for each point count
#                    - e.g., Date, Location, Observer, Point count duration, etc
#   2) PC_matrix - a matrix object with counts for each species for each point count
#  
#   Note that PC_surveyinfo and PC_matrix should have same number of rows
#
# ************************************************************
# ************************************************************

PCData <- read.csv("!Data/SaskAtlas_PCData_2022-05-16.csv")

# Add species codes to data
PCData$spcd <- NA
for (i in unique(PCData$species_id)) {
  if(i %in% species$species_id) PCData$spcd[which(PCData$species_id == i)] <- species$BSCDATA[species$species_id==i]
}

# Create a "surveyID" column (unique to each point count)
PCData <- PCData %>% mutate(surveyID = as.numeric(factor(paste0(latitude,longitude,ObservationDate,TimeCollected))))

# Duration in minutes for each point count
PCData$DurationInMinutes <- round(PCData$DurationInHours * 60)

# ---------------------------------------------------------
# Initial QA/QC
# ---------------------------------------------------------

# Remove species not in species list
PCData <- subset(PCData, spcd %in% species_vec & !is.na(spcd))

# Remove NA counts
PCData <- subset(PCData, !is.na(ObservationCount))

# Set any counts over 9 to 9
PCData$ObservationCount[PCData$ObservationCount > 9] <- 9

# ---------------------------------------------------------
# Dataframe with one row per survey
# ---------------------------------------------------------

# Construct a dataset that contains JUST survey information (not counts)
PC_surveyinfo <- PCData %>% 
  as.data.frame() %>%
  dplyr::select(surveyID,SiteCode,latitude,longitude,ObservationDate,TimeCollected,DurationInMinutes) %>%
  unique() %>%
  rename(sq_id = SiteCode) %>%
  st_as_sf(.,coords = c("longitude","latitude"), crs = CRS("+init=epsg:4326"), remove = FALSE)


# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# ONLY SELECT SURVEYS IN APPROPRIATE DATE AND TIME RANGE
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

# ------------------------------------------------
# RESTRICT DATES OF SURVEYS TO BETWEEN MAY 31 AND JULY 7
# ------------------------------------------------

PC_surveyinfo$Date <- lubridate::ymd(PC_surveyinfo$ObservationDate)
PC_surveyinfo$yday <- lubridate::ymd(PC_surveyinfo$ObservationDate) %>% yday()
PC_surveyinfo <- subset(PC_surveyinfo, yday >= start_date & yday <= end_date)

# ------------------------------------------------
# RESTRICT TIMES OF SURVEYS TO WITHIN 4 HOURS OF SUNRISE
# ------------------------------------------------

PC_surveyinfo$Hours <- floor(PC_surveyinfo$TimeCollected)
PC_surveyinfo$Mins <- round(60*(PC_surveyinfo$TimeCollected - PC_surveyinfo$Hours))
PC_surveyinfo$Time <- PC_surveyinfo$Date + hours(PC_surveyinfo$Hours) + minutes(PC_surveyinfo$Hours)
tz(PC_surveyinfo$Time) <-  "Canada/Central"
PC_surveyinfo$Sunrise <- PC_surveyinfo$Time

# Calculate local sunrise at each location
# For some reason this takes forever if done as vector calculation.  Use for loop (about 2 min)
start <- Sys.time()
for (i in 1:nrow(PC_surveyinfo)){
  PC_surveyinfo$Sunrise[i] <- sunrise_time(
    date = PC_surveyinfo$Date[i],
    tz = "Canada/Central",
    geocode = tibble::tibble(lon = PC_surveyinfo$longitude[i], 
                             lat = PC_surveyinfo$latitude[i]),
    twilight = "sunlight",
    unit.out = "datetime"
  )
}
end <- Sys.time()

# Hours since sunrise
PC_surveyinfo$HSS <- with(PC_surveyinfo, difftime(Time,Sunrise,units="hours") )
#PC_surveyinfo <- subset(PC_surveyinfo, HSS >= -0.5 & HSS <= 4)

# ---------------------------------------------------------
# Create matrix of species counts (each row corresponds to one in PC_surveyinfo)
# ---------------------------------------------------------

PC_matrix <- matrix(0,nrow = nrow(PC_surveyinfo),ncol = length(species_vec),
                    dimnames = list(PC_surveyinfo$surveyID,species_vec))

for (i in 1:nrow(PCData)){
  PC_matrix[which(rownames(PC_matrix) == PCData$surveyID[i]),which(colnames(PC_matrix) == PCData$spcd[i])] <- PCData$ObservationCount[i]
  print(i/nrow(PCData))
}

# Remove columns (species) that were never observed to save storage space
PC_matrix <- PC_matrix[,-which(colSums(PC_matrix)==0)]



# ************************************************************
# ************************************************************
# PART 2: PREPARE CHECKLIST DATA
#
# This script creates two objects:
#   1) DO_surveyinfo - an sf object with "survey information" for each checklist
#                    - e.g., Date, Location, Observer, Checklist effort, etc
#   2) DO_matrix - a matrix object with counts for each species on each checklist
#  
# ************************************************************
# ************************************************************

DOData <- read.csv("!Data/SaskAtlas_DOData_2022-05-16.csv")

# Add species codes to data
DOData$spcd <- NA
for (i in unique(DOData$species_id)) {
  if(i %in% species$species_id) DOData$spcd[which(DOData$species_id == i)] <- species$BSCDATA[species$species_id==i]
}

# ------------------------------------------------------------
# QA/QC
# ------------------------------------------------------------

# Only include complete checklists (AllSpeciesReported == "Y")
DOData$AllSpeciesReported[DOData$AllSpeciesReported == "y"] <- "Y"
table(DOData$AllSpeciesReported, useNA = "always")
DOData <- subset(DOData, AllSpeciesReported == "Y")

# There are many NAs in the checklist dataset.  Assign these a count of 1, assuming this was a detection
# NOTE THAT ANALYSES WILL ONLY USE CHECKLISTS AS PRESENCE/ABSENCE RECORDS
DOData$ObservationCount[is.na(DOData$ObservationCount)] <- 1

# ---------------------------------------------------------
# Dataframe with one row per checklist
# ---------------------------------------------------------

# Summarize remaining checklist information
DO_surveyinfo <- DOData %>% 
  group_by(SamplingEventIdentifier,
           SiteCode,
           ObservationDate,TimeCollected,
           CollectorNumber,ProtocolType,
           AllIndividualsReported,AllSpeciesReported,
           
           # Effort covariates
           DurationInHours,
           NumberOfObservers,
           
           # Travel distance
           EffortMeasurement1,
           
           # Does checklist include point counts (EffortMeasurement3 == 1)
           EffortMeasurement3) %>%
  
  summarize(n_species = length(unique(spcd)),
            n_birds = sum(ObservationCount),
            
            # There is one checklist with a minor discrepancy in lat/lon coordinates among rows, so
            # calculate the mean lat/lon.  All other checklists are consistent among rows
            latitude = mean(latitude),
            longitude = mean(longitude)) %>%
  
  rename(TravelDistance_m = EffortMeasurement1,
         IncludesPointCounts = EffortMeasurement3,
         sq_id = SiteCode) %>%
  ungroup()


# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# ONLY SELECT SURVEYS IN APPROPRIATE DATE AND TIME RANGE
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

# ------------------------------------------------
# RESTRICT DATES OF SURVEYS TO BETWEEN MAY 31 AND JULY 7
# ------------------------------------------------

DO_surveyinfo$Date <- lubridate::ymd(DO_surveyinfo$ObservationDate)
DO_surveyinfo$yday <- lubridate::ymd(DO_surveyinfo$ObservationDate) %>% yday()
DO_surveyinfo <- subset(DO_surveyinfo, yday >= start_date & yday <= end_date)

# ------------------------------------------------
# RESTRICT TIMES OF SURVEYS TO WITHIN 4 HOURS OF SUNRISE
# ------------------------------------------------

DO_surveyinfo$Hours <- floor(DO_surveyinfo$TimeCollected)
DO_surveyinfo$Mins <- round(60*(DO_surveyinfo$TimeCollected - DO_surveyinfo$Hours))
DO_surveyinfo$Time <- DO_surveyinfo$Date + hours(DO_surveyinfo$Hours) + minutes(DO_surveyinfo$Hours)
tz(DO_surveyinfo$Time) <-  "Canada/Central"
DO_surveyinfo$Sunrise <- DO_surveyinfo$Time

# Calculate local sunrise at each location
# For some reason this takes forever if done as vector calculation.  Use for loop (about 2 min)
start <- Sys.time()
for (i in 1:nrow(DO_surveyinfo)){
  DO_surveyinfo$Sunrise[i] <- sunrise_time(
    date = DO_surveyinfo$Date[i],
    tz = "Canada/Central",
    geocode = tibble::tibble(lon = DO_surveyinfo$longitude[i], 
                             lat = DO_surveyinfo$latitude[i]),
    twilight = "sunlight",
    unit.out = "datetime"
  )
}
end <- Sys.time()

# Hours since sunrise
DO_surveyinfo$HSS <- with(DO_surveyinfo, difftime(Time,Sunrise,units="hours") )

#DO_surveyinfo <- subset(DO_surveyinfo, HSS >= -0.5 & HSS <= 4)

# ---------------------------------------------------------
# Create matrix of species counts (each row corresponds to one in DO_surveyinfo)
# ---------------------------------------------------------

DO_matrix <- matrix(0,nrow = nrow(DO_surveyinfo),ncol = length(species_vec),
                    dimnames = list(DO_surveyinfo$SamplingEventIdentifier,species_vec))

for (i in 1:nrow(DOData)) DO_matrix[which(rownames(DO_matrix) == DOData$SamplingEventIdentifier[i]),which(colnames(DO_matrix) == DOData$spcd[i])] <- DOData$ObservationCount[i]

# Remove columns (species) that were never observed to save storage space
DO_matrix <- DO_matrix[,-which(colSums(DO_matrix)==0)]

# ************************************************************
# ************************************************************
# PART 3: Extract covariates associated with each survey location
#
#   - currently just pull values from SaskGrid
# ************************************************************
# ************************************************************

load(file = "!Data_processed/SaskGrid_PCA.RData")

# ---------------------------------------------------------
# Point counts
# ---------------------------------------------------------

# Unique numerical identifier for each point count, called 'obs_index'
PC_surveyinfo$obs_index <- 1:nrow(PC_surveyinfo)

# Add covariates using spatial intersection with SaskGrid
PC_surveyinfo <- PC_surveyinfo %>%
  st_transform(crs(SaskGrid)) %>%
  st_intersection(SaskGrid)

PC_surveyinfo <- PC_surveyinfo %>% arrange(obs_index) # Ensure this is ordered correctly, so that rows match up with PC_matrix

# Move geometry to last column
PC_surveyinfo <- PC_surveyinfo %>% relocate(geometry, .after = last_col())

Pointcount_data <- list(PC_surveyinfo = PC_surveyinfo,
                        PC_matrix = PC_matrix)

# ---------------------------------------------------------
# Checklists
# ---------------------------------------------------------

# Unique numerical identifier for each checklist, called 'obs_index'
DO_surveyinfo$obs_index <- 1:nrow(DO_surveyinfo)

# Add covariates using spatial intersection with SaskGrid
DO_surveyinfo <- DO_surveyinfo %>%
  st_as_sf(coords = c("longitude","latitude"), crs = CRS("+init=epsg:4326")) %>%
  st_transform(crs(SaskGrid)) %>%
  st_intersection(SaskGrid)

DO_surveyinfo <- DO_surveyinfo %>% arrange(obs_index) # Ensure this is ordered correctly, so that rows match up with DO_matrix

# Move geometry to last column
DO_surveyinfo <- DO_surveyinfo %>% relocate(geometry, .after = last_col())

Checklist_data <- list(DO_surveyinfo = DO_surveyinfo,
                       DO_matrix = DO_matrix)



# ---------------------------------------------------------
# Save
# ---------------------------------------------------------

save(Pointcount_data,file = "D:/Working_Files/1_Projects/Landbirds/SK_BBA_analysis/AOS_precision/data/Pointcount_data.RData")
save(Checklist_data,file = "D:/Working_Files/1_Projects/Landbirds/SK_BBA_analysis/AOS_precision/data/Checklist_data.RData")
