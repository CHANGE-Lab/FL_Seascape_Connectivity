#### DESCRIPTION ####
# This script is part of Courtney Stuart's second MSc data chapter in the lab of 
# Dr. Stephanie Green at the University of Alberta (2019-2021). This script uses
# habitat data from Courtney's first data chapter (derived from the Unified 
# Florida Reef Map) to produce a raster data layer where each cell stores the 
# Euclidean distance from that location to the nearest continuous seagrass patch.
# This raster, along with one storing Euclidean distances to nearest mangrove 
# habitat cell created during Ch1, will be used to evaluate how proximity to 
# nursery habitats influences potential functional connectivity for sub-adult 
# gray snapper (Lutjanus griseus) and bluestriped grunt (Haemulon sciurus) 
# occupying a spatially heterogeneous study seascape in the Florida Keys, USA. 

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
lg_wd = "Z:/Courtney/Stuart_MSc_Ch2/Graphab/Subadult_Gray_Snapper_10m/" # graphab results
hs_wd = "Z:/Courtney/Stuart_MSc_Ch2/Graphab/Subadult_Bluestriped_Grunt_10m/" # graphab results
git_wd = "Z:/Courtney/Stuart_MSc_Ch2/GitHub/FL_Seascape_Connectivity/" # github data
figs_wd = "Z:/Courtney/Stuart_MSc_Ch2/Figures/" # figures

# libraries
library(easypackages)
libraries("conflicted", "sp", "sf", "rgdal", "gdalUtils", "tidyverse", 
          "dplyr", "ggplot2", "ggh4x", "PNWColors", "readxl", "scales")
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")

# save species-specific Euclidean distance to cost distance equivalencies
#  sub-adult gray snapper Lutjanus griseus (lg)
lg_5km = as.numeric(25761.9) # 5 km =~ 25761.9 cost units
lg_10km = as.numeric(78526.5) # 10 km =~ 78526.5 cost units
lg_15km = as.numeric(150716.0) # 15 km =~ 150716 cost units

# sub-adult bluestriped grunt Haemulon sciurus (hs)
hs_5km = as.numeric(17164.8) # 5 km =~ 17164.8 cost units
hs_10km = as.numeric(50898.8) # 10 km =~ 50898.8 cost units
hs_15km = as.numeric(96128.0) # 15 km =~ 96128.0 cost units

####### BATCH GLOBAL CONNECTIVITY METRICS #######
# where metrics were calculated at intervals of 5,000 cost units from 0-250,000 units
 
#### gray snapper ####
# probability of connectivity (PC)
lg_pc = read.delim(paste0(lg_wd, "Connectivity_Metrics/Batch_PC_Global_Metric.txt"))

# equivalent connectivity (EC)               
lg_ec = read.delim(paste0(lg_wd, "Connectivity_Metrics/Batch_EC_Global_Metric.txt"))

#### bluestriped grunt ####
# probability of connectivity (PC)
hs_pc = read.delim(paste0(hs_wd, "Connectivity_Metrics/Batch_PC_Global_Metric.txt"))

# equivalent connectivity (EC)                
hs_ec = read.delim(paste0(hs_wd, "Connectivity_Metrics/Batch_EC_Global_Metric.txt"))

# add species identification columns for compiling
lg_pc$Species = rep("Lutjanus griseus", nrow(lg_pc))
lg_ec$Species = rep("Lutjanus griseus", nrow(lg_ec))
hs_pc$Species = rep("Haemulon sciurus", nrow(hs_pc))
hs_ec$Species = rep("Haemulon sciurus", nrow(hs_ec))

# combine global metrics from both species
pc = rbind(lg_pc, hs_pc)
ec = rbind(lg_ec, hs_ec)

# create color palette for plotting
my_pal = pnw_palette("Bay", 8)

#### PLOT GLOBAL CONNECTIVITY METRICS ####

# plot global PC
# add three vertical lines per species to indicate potential dispersal thresholds
# a lower estimate (5 km Euclidean dispersal limit)
# our selected estimate (10 km Euclidean dispersal limit), and 
# a higher estimate (15 km Euclidean dispersal limit)
pc_plot = 
  ggplot(pc, aes(y = PC, x = d, group = Species, 
                 color = Species, fill = Species)) +
  geom_line(aes(color = Species), size = 0.8) +
  scale_color_manual(values = c(my_pal[1], my_pal[5])) +
  scale_fill_manual(values = c(my_pal[1], my_pal[5])) +
  xlab("Dispersal distance (cost units)") +
  ylab("Probability of connectivity (PC)") + 
  scale_x_continuous(minor_breaks = seq(0, 250000, by = 5000),
                     breaks = seq(0, 250000, 50000),
                     guide = "axis_minor") +
  geom_vline(xintercept = lg_5km, color = my_pal[5], size = 0.5, linetype = "dashed") +
  geom_vline(xintercept = lg_10km, color = my_pal[5], size = 0.5, linetype = "dotted") +
  geom_vline(xintercept = lg_15km, color = my_pal[5], size = 0.5, linetype = "dotdash") +
  geom_vline(xintercept = hs_5km, color = my_pal[1], size = 0.5, linetype = "dashed") +
  geom_vline(xintercept = hs_10km, color = my_pal[1], size = 0.5, linetype = "dotted") +
  geom_vline(xintercept = hs_15km, color = my_pal[1], size = 0.5, linetype = "dotdash") +
  theme_bw() +  
  theme(legend.position = "none") + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))
pc_plot

# save PC plot
png(file = paste0(temp_wd, "PC_Global_No_Legend.png"),
    res = 450,  width = 6, height = 4.5, units = "in")
pc_plot
dev.off()

png(file = paste0(temp_wd, "PC_Global_No_Legend_5x3.png"),
    res = 450,  width = 5, height = 3, units = "in")
pc_plot
dev.off()

# plot global EC
# add three vertical lines per species to indicate potential dispersal thresholds
# a lower estimate (5 km Euclidean dispersal limit)
# our selected estimate (10 km Euclidean dispersal limit), and 
# a higher estimate (15 km Euclidean dispersal limit)
ec_plot = 
  ggplot(ec, aes(y = EC, x = d, group = Species, 
                 color = Species, fill = Species)) +
  geom_line(aes(color = Species), size = 0.8) +
  scale_color_manual(values = c(my_pal[1], my_pal[5])) +
  scale_fill_manual(values = c(my_pal[1], my_pal[5])) +
  xlab("Dispersal distance (cost units)") +
  ylab("Equivalent connectivity (EC)") + 
  scale_x_continuous(minor_breaks = seq(0, 250000, by = 5000),
                     breaks = seq(0, 250000, 50000),
                     guide = "axis_minor") +
  geom_vline(xintercept = lg_5km, color = my_pal[5], size = 0.5, linetype = "dashed") +
  geom_vline(xintercept = lg_10km, color = my_pal[5], size = 0.5, linetype = "dotted") +
  geom_vline(xintercept = lg_15km, color = my_pal[5], size = 0.5, linetype = "dotdash") +
  geom_vline(xintercept = hs_5km, color = my_pal[1], size = 0.5, linetype = "dashed") +
  geom_vline(xintercept = hs_10km, color = my_pal[1], size = 0.5, linetype = "dotted") +
  geom_vline(xintercept = hs_15km, color = my_pal[1], size = 0.5, linetype = "dotdash") +
  theme_bw() +  
  theme(legend.position = "none") + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))
ec_plot

# save EC plot
png(file = paste0(temp_wd, "EC_Global_No_Legend.png"),
    res = 450,  width = 6, height = 4.5, units = "in")
ec_plot
dev.off()

png(file = paste0(temp_wd, "EC_Global_No_Legend_5x3.png"),
    res = 450,  width = 5, height = 3, units = "in")
ec_plot
dev.off()

####### BATCH LOCAL CONNECTIVITY METRICS #######
# where metrics were calculated at intervals of 5,000 cost units from 0-250,000 units
# as well as at three dispersal distance thresholds (i.e., cost unit equivalencies 
# of 5 km, 10 km, & 15 km, as defined in the code above)
# these datasets include the interaction flux (IF) local connectivity estimates
# at 5 selected and 10 alternate Mission Iconic Reefs (MIR) coral restoration sites.
lg_if = read.csv(paste0(git_wd, "Data/Subadult_Gray_Snapper_IF_MIRsites.csv"))
hs_if = read.csv(paste0(git_wd, "Data/Subadult_Bluestriped_Grunt_IF_MIRsites.csv"))

# add species columns for compiling
lg_if$Species = rep("Lutjanus griseus", nrow(lg_if))
hs_if$Species = rep("Haemulon sciurus", nrow(hs_if))

# create new data frame for ranking sites by IF at 5 km, 10 km, and 15 km
lg_if_rank = lg_if %>%
  select(Species, Site, Selected_I, LON_M, LAT_M, Id, 
         IF_d2576, IF_d7852, IF_d1507) %>%
  rename(Selected_Iconic_Reef = Selected_I,
         IF_5km = IF_d2576,
         IF_10km = IF_d7852,
         IF_15km = IF_d1507) %>%
  mutate(Rank_5km = dense_rank(desc(IF_5km)),
         Rank_10km = dense_rank(desc(IF_10km)),
         Rank_15km = dense_rank(desc(IF_15km))) %>%
  select(Species, Site, Selected_Iconic_Reef, LON_M, LAT_M,
         Id, IF_5km, Rank_5km, IF_10km, Rank_10km,
         IF_15km, Rank_15km)

hs_if_rank = hs_if %>%
  select(Species, Site, Selected_I, LON_M, LAT_M, Id, 
         IF_d1716, IF_d5089, IF_d9612) %>%
  rename(Selected_Iconic_Reef = Selected_I,
         IF_5km = IF_d1716,
         IF_10km = IF_d5089,
         IF_15km = IF_d9612) %>%
  mutate(Rank_5km = dense_rank(desc(IF_5km)),
         Rank_10km = dense_rank(desc(IF_10km)),
         Rank_15km = dense_rank(desc(IF_15km))) %>%
  select(Species, Site, Selected_Iconic_Reef, LON_M, LAT_M,
         Id, IF_5km, Rank_5km, IF_10km, Rank_10km,
         IF_15km, Rank_15km)

# combine IF data from both species
if_rank = rbind(lg_if_rank, hs_if_rank)

# save IF connectivity rankings at MIR sites to github repo
write.csv(if_rank, paste0(git_wd, "Data/IF_Rankings_MIRsites.csv"), 
          row.names = FALSE)
