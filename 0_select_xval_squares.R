# *******************************************
# *******************************************
# This script sets aside a spatially balanced subset of atlas squares for cross-validation
# *******************************************
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
  
  # For plotting
  'viridis','scales','ggpubr','ggtext',
  
  # For GRTS
  'spsurvey')

if (any(!my_packs %in% installed.packages()[, 'Package'])) {install.packages(my_packs[which(!my_packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)}
lapply(my_packs, require, character.only = TRUE)


rm(list=ls())

setwd("D:/Working_Files/1_Projects/Landbirds/SK_BBA_analysis/Standard_Analysis/")
`%!in%` <- Negate(`%in%`)

# Load SaskSquares spatial object
SaskSquares <- read_sf("!Shapefiles/SaskSquares/SaskSquares.shp")

# Centroids of sask atlas squares
SK_centroids <- st_centroid(SaskSquares)

# Divide grid into 5 spatially balanced folds
grts_grid <- spsurvey::grts(SK_centroids, n_base = nrow(SK_centroids))$sites_base
grts_grid$grts_order <- 1:nrow(grts_grid)
grts_grid$fold <- cut(grts_grid$grts_order,seq(1,max(grts_grid$grts_order),length.out = 6),include.lowest = TRUE) %>% as.numeric()
table(grts_grid$fold)

SaskSquares_xval <- full_join(SaskSquares,as.data.frame(grts_grid)[,c("SQUARE_ID","fold","grts_order")])

ggplot(SaskSquares_xval, aes(fill = fold))+
  geom_sf(col = "transparent")+
  scale_fill_gradientn(colors=viridis(5))

st_write(SaskSquares_xval,"../AOS_precision/output/SaskSquares_xval.shp",append=FALSE)

