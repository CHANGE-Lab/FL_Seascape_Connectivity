#### DESCRIPTION ####
# This script is part of Courtney Stuart's second MSc Chapter in the Lab of 
# Dr. Stephanie Green at the University of Alberta (2019-2021). The Graphab
# software was used to construct species-specific minimum planar spatial graphs
# for sub-adult gray snapper (Lutjanus griseus) and bluestriped grunt (Haemulon
# sciurus) occupying a spatially heterogeneous seascape in the Florida Keys, USA.
# These spatial graphs were analyzed quantitatively to estimate seascape- and 
# node-scale potential connectivity for each species, using the probability of 
# connectivity and equivalent connectivity (global), and interaction flux (local)
# connectivity metrics, respectively. We use interaction flux (IF) estimates of 
# potential connectivity to evaluate the spatial design of a newly established 
# (2019) broad-scale coral reef restoration program in the Florida Keys, USA - 
# Mission: Iconic Reefs. Here, examine whether distances to the nearest nursery
# habitats (i.e., mangroves and continuous seagrass beds) are correlated with IF
# values at potential reef restoration sites by calculating Pearson pairwise 
# correlation coefficients (r).

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

# data directory
data_wd = "Z:/Courtney/Stuart_MSc_Ch2/GitHub/FL_Seascape_Connectivity/Data/"

# species-specific node IF values
# gray snapper
lg_node_IF = read.csv(paste0(
  data_wd, "Subadult_Gray_Snapper_Node_Data.csv"))

# correlation between node IF values and mangrove and seagrass distances
lg_node_mg = round(cor(lg_node_IF$IF_d7852, lg_node_IF$Mean_Mangrove_Dist, 
                       use = "everything", method = "pearson"), 4)

lg_node_sg = round(cor(lg_node_IF$IF_d7852, lg_node_IF$Mean_Seagrass_Dist, 
                       use = "everything", method = "pearson"), 4)

# scientific notation for reporting
format(lg_node_mg, scientific = TRUE)
format(lg_node_sg, scientific = TRUE)

# species-specific node IF values
# bluestriped grunt
hs_node_IF = read.csv(paste0(
  data_wd, "Subadult_Bluestriped_Grunt_Node_Data.csv"))

# correlation between node IF values and mangrove and seagrass distances
hs_node_mg = round(cor(hs_node_IF$IF_d5089, hs_node_IF$Mean_Mangrove_Dist, 
                       use = "everything", method = "pearson"), 4)

hs_node_sg = round(cor(hs_node_IF$IF_d5089, hs_node_IF$Mean_Seagrass_Dist, 
                       use = "everything", method = "pearson"), 4)

# scientific notation for reporting
format(hs_node_mg, scientific = TRUE)
format(hs_node_sg, scientific = TRUE)
