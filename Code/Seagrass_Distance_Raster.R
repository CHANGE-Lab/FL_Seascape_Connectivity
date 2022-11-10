#### DESCRIPTION ####
# This script is part of Courtney Stuart's second MSc Chapter in the Lab of 
# Dr. Stephanie Green at the University of Alberta (2019-2021). This script uses
# habitat data from Courtney's first data chapter (derived from the Unified 
# Florida Reef Map) to produce a raster data layer where each cell stores the 
# Euclidean distance from that location to the nearest continuous seagrass patch.
# This raster, along with one storing Euclidean distances to nearest mangrove 
# habitat cell created during Ch1, will be used to evaluate how proximity to 
# nursery habitats influences potential functional connectivity for sub-adult 
# gray snapper (Lutjanus griseus) and bluestriped grunt (Haemulon sciurus) occupying
# the spatially heterogeneous study seascape in the Florida Keys, USA. 


#### ASSOCIATED PUBLICATIONS ####
# For further information on the HSMs, spatial predictors, and/or reef fish
# occurrence records used in this research, read:

# Stuart, C. E., Wedding, L. M., Pittman, S. J., & Green, S. J. (2021). Habitat 
# Suitability Modeling to Inform Seascape Connectivity Conservation and Management. 
# Diversity, 13(10), 465.

# Stuart, C. E., Wedding, L. M., Pittman, S. J., Serafy, J. E., Moura, A., & 
# Green, S. J. (2022). Seascape connectivity modeling reveals potential hotspots of 
# fish-derived nutrient provisioning to restored coral reefs. Marine Ecology 
# Progress Series (in review).

#### CONTACT ####
# Courtney Stuart (courtney.e.stuart@gmail.com OR courtney.stuart@mansfield.ox.ac.uk)

#### SET-UP ####
# working directory
setwd("Z:/Courtney/Stuart_MSc_Ch2/") # main project folder for Ch2

# data directories 
spatial_wd = "Z:/Courtney/Stuart_MSc_Ch1/Spatial_Predictors/" # Ch1 spatial predictors
csv_wd = "Z:/Courtney/Stuart_MSc_Ch1/GitHub/FL_Habitat_Suitability/Data/" # Ch1 tabular data
gis_wd = "Z:/Courtney/Stuart_MSc_Ch2/GIS_Files/" # Ch2 spatial data exported from R
temp_wd = "Z:/Courtney/Stuart_MSc_Ch2/Temporary/" # temporary files

# load libraries
library(raster)
library(easypackages)
libraries("raster")
libraries("conflicted", "rgdal", "gdalUtils", "raster", "sp", "sf", "tmap", "dplyr", 
          "lwgeom", "rgeos", "cleangeo", "tidyverse", "stars", "fasterize", 
          "PNWColors", "spex", "igraph", "spatialEco")

# change where large temp rasters are saved
rasterOptions(tmpdir = "Z:/Courtney/Stuart_MSc_Ch2/Temporary/")

# save PROJ.4 string for standard projection (ESPG:26958 NAD 83/Florida East) 
my_crs = CRS("+init=epsg:26958")

# read in the habitat raster produced in chapter one
habitat = raster(paste0(spatial_wd, "Habitat.asc"))
raster::crs(habitat) = my_crs # set coordinate system info

# read in the Unified Reef Map csv file from chapter one that stores the 
# numerical and categorical benthic habitat descriptions
URM = read.csv(paste0(csv_wd, "URM_ClassLv1_IDs.csv"))
head(URM, 15)
# all cells in the habitat raster storing the value 3 are continuous seagrass

# create a new raster, where the value in each cell reflects the Euclidean 
# distance (meters) from that cell to the nearest continuous seagrass cell (ID:3). 
# Save this raster to the GIS folder for chapter two.
seagrass_dist = writeRaster(raster::mask(
  raster::crop(gridDistance(habitat, origin = 3), habitat), habitat),
  file = file.path(gis_wd, "Seagrass_Dist.asc"),
  format = "ascii", overwrite = T)
