#### DESCRIPTION ####
# This script is part of Courtney Stuart's second MSc data chapter in the lab of 
# Dr. Stephanie Green at the University of Alberta (2019-2021). This script 
# leverages previously constructed (data chapter one) MaxEnt habitat suitability
# models (HSMs) for ecological connectivity modeling. These HSMs are specific to 
# sub-adult gray snapper (Lutjanus griseus) and bluestriped grunt (Haemulon 
# sciurus) in the Florida Keys, USA. HSMs were constructed using species occurrence
# records from reef fish surveys and spatially explicit data layers of benthic 
# habitat composition and configuration, bathymetry and seafloor surface morphology,
# and seasonal water conditions.These continuous HSMs will be discretized to 
# produce binary HSMs (suitable habitat vs. matrix), and converted to resistance 
# surfaces using a negative exponential transformation function, for ecological 
# connectivity modeling using a graph-theoretic approach (via Graphab).

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
setwd("Z:/Courtney/Stuart_MSc_Ch2/") # main project folder

# data directories 
temp_wd = "Z:/Courtney/Stuart_MSc_Ch2/Temporary/" # temporary files
HSMs = "Z:/Courtney/Stuart_MSc_Ch2/HSMs/MaxEnt/" # HSMs produced in Ch1
gis_wd = "Z:/Courtney/Stuart_MSc_Ch2/GIS_Files/" # spatial data exported from R

# libraries
library(easypackages)
libraries("conflicted", "raster", "sp", "sf", "rgdal", "gdalUtils", 
          "lwgeom", "rgeos", "cleangeo", "fasterize", "parallel",
          "doParallel", "tidyverse", "dplyr")
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")

# change where large temp rasters are saved during processing
rasterOptions(tmpdir = temp_wd)

# save PROJ.4 string for standard projection (ESPG:26958 NAD 83/Florida East)
my_crs = crs("+init=epsg:26958")

#### Gray Snapper ####
# read in the raster data layer of predicted habitat suitability for sub-adult
# gray snapper (average from 10-fold cross validation in MaxEnt)
require(raster)
lg_hsm = raster(paste0(HSMs, "Subadult_Gray_Snapper/LUT_GRIS_avg.asc"))

# define CRS
crs(lg_hsm) = my_crs

# extract the Max SSS (maximum sum of training sensitivity and specificity) 
# threshold from the MaxEnt results .csv file
lg_results = read.csv(paste0(HSMs, "Subadult_Gray_Snapper/maxentResults.csv"))
lg_maxSSS = lg_results %>%
  filter(Species == "LUT_GRIS (average)") %>%
  select(Maximum.training.sensitivity.plus.specificity.Cloglog.threshold)

# create a 3 column matrix for reclassifying raster cells where the first column
# is "from", the second column is "to", and the third is "becomes". reclassify 
# so that all cells < maxSSS become 0 (matrix) and all cells => maxSSS
# become 1 (suitable habitat).
lg_maxSSS$Maximum.training.sensitivity.plus.specificity.Cloglog.threshold
lg_matrix = matrix(c(0, (as.numeric(lg_maxSSS)-0.0001), 0,
                     as.numeric(lg_maxSSS), 1, 1),
                   ncol = 3, byrow = TRUE)

# reclassify raster cells to produce binary HSM
cores = detectCores() - 2 # run in parallel
beginCluster(n = cores)
lg_binary = reclassify(lg_hsm, lg_matrix)
endCluster()

# save binary HSM as tiff (integer type)
writeRaster(lg_binary, 
            filename = paste0(gis_wd, "Subadult_Gray_Snapper_Binary_HSM.tif"),
            format = "GTiff", datatype = "INT2S", overwrite = TRUE)

# now create a cost/resistance surface representing seascape impedance to fish 
# movement. following Keeley, Beier, and Gagnon Landscape Ecol (2016) DOI 
# 10.1007/s10980-016-0387-5 and Duflot et al. J. Nat. Conserv. (2018)  DOI 
# 10.1016/j.jnc.2018.08.005, we'll use a negative exponential function to 
# describe the relationship between seascape resistance and habitat suitability:

# if suitability >= MaxSSS threshold --> suitable habitat --> resistance = 1
# if suitability < MaxSSS threshold --> matrix --> resistance = 
# e^((ln(0.001)/MaxSSS)*suitability)*10^3
beginCluster(n = cores) # run in parallel
lg_cost = calc(lg_hsm, 
               fun = function(r) {
                 r[r < as.numeric(lg_maxSSS)] <-- exp((log(0.001)/as.numeric(lg_maxSSS))*r)*(10^3)
                 r[r >= as.numeric(lg_maxSSS)] <-- 1
                 return(r*-1)
               })

endCluster()

# save the cost surface
writeRaster(lg_cost, filename = paste0(gis_wd, "Subadult_Gray_Snapper_Cost.tif"),
            format = "GTiff", overwrite = TRUE)

#### Bluestriped Grunt ####
# read in the raster data layer of predicted habitat suitability for sub-adult
# bluestriped grunt (average from 10-fold cross validation in MaxEnt)
require(raster)
hs_hsm = raster(paste0(HSMs, "Subadult_Bluestriped_Grunt/HAE_SCIU_avg.asc"))

# define CRS
crs(hs_hsm) = my_crs

# extract the Max SSS (maximum sum of training sensitivity and specificity) 
# threshold from the MaxEnt results .csv file
hs_results = read.csv(paste0(HSMs, "Subadult_Bluestriped_Grunt/maxentResults.csv"))
hs_maxSSS = hs_results %>%
  filter(Species == "HAE_SCIU (average)") %>%
  select(Maximum.training.sensitivity.plus.specificity.Cloglog.threshold)

# create a 3 column matrix for reclassifying raster cells where the first column
# is "from", the second column is "to", and the third is "becomes". reclassify 
# so that all cells < maxSSS become 0 (matrix) and all cells => maxSSS
# become 1 (suitable habitat).
hs_maxSSS$Maximum.training.sensitivity.plus.specificity.Cloglog.threshold
hs_matrix = matrix(c(0, (as.numeric(hs_maxSSS)-0.0001), 0,
                     as.numeric(hs_maxSSS), 1, 1),
                   ncol = 3, byrow = TRUE)

# reclassify raster cells to produce binary HSM
beginCluster(n = cores)
hs_binary = reclassify(hs_hsm, hs_matrix)
endCluster()

# save binary HSM as tiff (integer type)
writeRaster(hs_binary, 
            filename = paste0(gis_wd, "Subadult_Bluestriped_Grunt_Binary_HSM.tif"),
            format = "GTiff", datatype = "INT2S", overwrite = TRUE)

# now create a cost/resistance surface
# if suitability >= MaxSSS threshold --> suitable habitat --> resistance = 1
# if suitability < MaxSSS threshold --> matrix --> resistance = 
# e^((ln(0.001)/MaxSSS)*suitability)*10^3
beginCluster(n = cores)
hs_cost = calc(hs_hsm, 
               fun = function(r) {
                 r[r < as.numeric(hs_maxSSS)] <-- exp((log(0.001)/as.numeric(hs_maxSSS))*r)*(10^3)
                 r[r >= as.numeric(hs_maxSSS)] <-- 1
                 return(r*-1)
               })

endCluster()

# save the cost surface
writeRaster(hs_cost, filename = paste0(gis_wd, "Subadult_Bluestriped_Grunt_Cost.tif"),
            format = "GTiff", overwrite = TRUE)

# the 5x5 m cell resolution may preclude connectivity analyses over the whole 
# study region (exceeding memory and computation limits). Re-run these data
# prep steps using aggregated 10x10 m resolution rasters.

#### Gray Snapper 10 m ####
# aggregate continuous HSM
lg_hsm_10 = aggregate(lg_hsm, fact = 2, fun = mean, na.rm = TRUE, 
                      filename = paste0(gis_wd, "Subadult_Gray_Snapper_HSM_10m.tif"),
                      format = "GTiff", overwrite = TRUE)

# aggregate binary HSM
beginCluster(n = cores)
lg_binary_10 = reclassify(lg_hsm_10, lg_matrix)
endCluster()
writeRaster(lg_binary_10, 
            filename = paste0(gis_wd, "Subadult_Gray_Snapper_Binary_HSM_10m.tif"),
            format = "GTiff", datatype = "INT2S", overwrite = TRUE)

# aggregate cost surface
beginCluster(n = cores) # run in parallel
lg_cost_10 = calc(lg_hsm_10, 
                  fun = function(r) {
                    r[r < as.numeric(lg_maxSSS)] <-- exp((log(0.001)/as.numeric(lg_maxSSS))*r)*(10^3)
                    r[r >= as.numeric(lg_maxSSS)] <-- 1
                    return(r*-1)
                  })

endCluster()

# save the cost surface
writeRaster(lg_cost_10, filename = paste0(gis_wd, "Subadult_Gray_Snapper_Cost_10m.tif"),
            format = "GTiff", overwrite = TRUE)

#### Bluestriped Grunt 10 m ####
# aggregate continuous HSM
hs_hsm_10 = aggregate(hs_hsm, fact = 2, fun = mean, na.rm = TRUE, 
                      filename = paste0(gis_wd, "Subadult_Bluestriped_Grunt_HSM_10m.tif"),
                      format = "GTiff", overwrite = TRUE)

# aggregate binary HSM
beginCluster(n = cores)
hs_binary_10 = reclassify(hs_hsm_10, hs_matrix)
endCluster()
writeRaster(hs_binary_10, 
            filename = paste0(gis_wd, "Subadult_Bluestriped_Grunt_Binary_HSM_10m.tif"),
            format = "GTiff", datatype = "INT2S", overwrite = TRUE)

# aggregate cost surface
beginCluster(n = cores) # run in parallel
hs_cost_10 = calc(hs_hsm_10, 
                  fun = function(r) {
                    r[r < as.numeric(hs_maxSSS)] <-- exp((log(0.001)/as.numeric(hs_maxSSS))*r)*(10^3)
                    r[r >= as.numeric(hs_maxSSS)] <-- 1
                    return(r*-1)
                  })

endCluster()

# save the cost surface
writeRaster(hs_cost_10, filename = paste0(gis_wd, "Subadult_Bluestriped_Grunt_Cost_10m.tif"),
            format = "GTiff", overwrite = TRUE)

