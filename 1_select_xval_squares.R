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

set.seed(999)

setwd("D:/Working_Files/1_Projects/Landbirds/SK_BBA_analysis/Standard_Analysis/")
`%!in%` <- Negate(`%in%`)

# Load SaskSquares spatial object
SaskSquares <- read_sf("!Shapefiles/SaskSquares/SaskSquares.shp")

# 
# # GRTS
# SaskSquares_cent <- st_centroid(SaskSquares)
# grts_grid <- spsurvey::grts(SaskSquares_cent, n_base = nrow(SaskSquares_cent))$sites_base
# grts_grid$grts_order <- 1:nrow(grts_grid)
# 
# # Separate into 5 folds
# SaskSquares_xval <- full_join(SaskSquares,as.data.frame(grts_grid)[,c("SQUARE_ID","grts_order")])
# SaskSquares_xval$fold <- cut(SaskSquares_xval$grts_order,seq(1,max(SaskSquares_xval$grts_order),length.out = 6),include.lowest = TRUE) %>% as.numeric()
# table(SaskSquares_xval$fold)
# 
# ggplot(SaskSquares_xval, aes(fill = fold))+
#   geom_sf(col = "transparent")+
#   scale_fill_gradientn(colors=viridis(5))
# 
# st_write(SaskSquares_xval,"../AOS_precision/output/SaskSquares_xval.shp",append=FALSE)

# ------------------------------------------------------
# Optional code to create larger blocks (e.g., 50 km blocks)
# ------------------------------------------------------

# Divide squares into spatially balanced folds
SaskSquares_centroids <- st_centroid(SaskSquares)
grts_grid <- spsurvey::grts(SaskSquares, n_base = nrow(SaskSquares))$sites_base
grts_grid$grts_order <- 1:nrow(grts_grid)
# Create 100 km grid across study region
ext <- extent(SaskSquares)
ext[1] <- ext[1] - 100000
ext[2] <- ext[2] + 100000
ext[3] <- ext[3] - 100000
ext[4] <- ext[4] + 100000

grid <- raster(ext, resolution = 50000, crs = crs(SaskSquares))
grid <- as(grid, "SpatialPolygons") %>% st_as_sf()
grid$block_id <- 1:nrow(grid)
grid_cent <- grid %>% st_centroid()

# Divide grid into 5 spatially balanced folds
grts_grid <- spsurvey::grts(grid_cent, n_base = nrow(grid_cent))$sites_base
grts_grid$grts_order <- 1:nrow(grts_grid)
grts_grid <- grts_grid %>% arrange(block_id)
grid$grts_order <- grts_grid$grts_order
grid$fold <- cut(grid$grts_order,seq(1,max(grid$grts_order),length.out = 6),include.lowest = TRUE) %>% as.numeric()
table(grid$fold)

#plot(as(SaskSquares,'Spatial'))
#lines(as(grid,"Spatial"), col = "blue")

# 50 km block to which each SaskSquare belongs
SaskSquares_centroids <- st_centroid(SaskSquares)
tmp <- st_intersection(SaskSquares_centroids,grid)

SaskSquares_xval <- full_join(SaskSquares,as.data.frame(tmp)[,c("SQUARE_ID","block_id","fold","grts_order")])
ggplot(SaskSquares_xval, aes(fill = fold))+
  geom_sf(col = "transparent")+
  scale_fill_gradientn(colors=viridis(5))

st_write(SaskSquares_xval,"../AOS_precision/output/SaskSquares_xval_50km.shp",append=FALSE)

