#### DESCRIPTION ####
# This script is part of Courtney Stuart's second MSc Chapter in the Lab of 
# Dr. Stephanie Green at the University of Alberta (2019-2021). This script
# compares the rankings of fifteen candidate restoration sites considered under
# Florida's Mission: Iconic Reefs coral restoration initiative based on in situ
# data and spatial graph models of connectivity. Sites were first ranked based on 
# empirical data, i.e., the proportion of reef fish samples on which target 
# sub-adult gray snapper (Lutjanus griseus) and bluestriped grunt (Haemulon sciurus)
# were present. Sites were then ranked according to the local Interaction Fluc (IF)
# connectivity metric calculated from species-specific spatial graphs of potential
# functional connectivity. 


#### ASSOCIATED PUBLICATIONS ####
# For further information on the HSMs, spatial predictors, and/or reef fish
# occurrence records used in this research, read:

# Stuart, C. E., Wedding, L. M., Pittman, S. J., & Green, S. J. (2021). Habitat 
# Suitability Modeling to Inform Seascape Connectivity Conservation and Management. 
# Diversity, 13(10), 465.

# Stuart, C. E., Wedding, L. M., Pittman, S. J., Serafy, J. E., Moura, A., 
# Bruckner, A.W., & Green, S. J. (2023). Seascape connectivity modeling predicts
# potential hotspots of fish-derived nutrient provisioning to restored coral reefs.
# Marine Ecology Progress Series.

#### CONTACT ####
# Courtney Stuart (courtney.e.stuart@gmail.com OR courtney.stuart@mansfield.ox.ac.uk)

# libraries
library(easypackages)
libraries("conflicted","sp", "sf", "rgdal", "dplyr", "ggplot2", "rvc")
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")


# working directory
setwd("Z:/Courtney/Stuart_MSc_Ch2/") # main project folder

# project github repo
git_wd = "Z:/Courtney/Stuart_MSc_Ch2/GitHub/FL_Seascape_Connectivity/" 

# project geodatabase
fgdb = "Z:/Courtney/Stuart_MSc_Ch2/Geodatabases/Stuart_MSc_Ch2.gdb"
subset(ogrDrivers(), grepl("GDB", name)) 

# coordinate system information for geodatabase
my_crs = CRS("+init=epsg:26958")

# coordinate system information for fish survey data
gcs = CRS("+init=epsg:4326")


#### bluestriped grunts (Haemulon sciurus) ####
# find nodes coinciding spatially with Mission Iconic Reefs sites
# all nodes
hs_nodes = st_read(dsn = fgdb, layer = "Subadult_Bluestriped_Grunt_Batch_IF")

# 15 potential MIR sites
mir = st_read(dsn = fgdb, layer = "Available_Iconic_Reefs")

# keep only nodes with mir sites
hs_mir_nodes = hs_nodes %>%
  filter(Id %in% (st_intersection(mir, hs_nodes))$Id)

# add back in column with mir site names
# directory for which reef IDs match
hs_m = c(13218, "Davis Reef", 
         13535, "Cheeca Rocks", 
         11066, "French and Molasses Reefs",
         9064, "Key Largo Dry Rocks",
         7338, "Turtle Reef",
         19368, "Tennessee Reef", 
         24237, "South of Key Colony Beach (Marker 48)",
         9295, "Horseshoe Reef",
         24429, "Coffins Patch", 
         27229, "Newfound Harbor",
         13262, "Hen and Chickens", 
         25904, "Sombrero Reef", 
         27140, "Looe Key Reef", 
         8814, "Elbow Reef")
# save directory as data frame
hs_m_df = as.data.frame(matrix(hs_m, ncol = 2, byrow = T)) %>%
  rename(Id = V1,MIR_Site = V2)
# add and place in first column potition MIR site names
hs_mir_nodes$MIR_Site = hs_m_df$MIR_Site[match(hs_mir_nodes$Id, hs_m_df$Id)]
hs_mir_nodes = hs_mir_nodes %>% relocate(MIR_Site, .before = Id )

# remove temp data
rm(list = c("hs_m_df", "hs_m", "hs_nodes"))

# use the bluestriped grunt MIR nodes to extract observation data from the
# Reef Visual Census (RVC) program [survey years 2014, 2016, 2018]
rvc = getSampleData(years = c(2014, 2016, 2018), regions = "FLA KEYS")

# how many unique sites were sampled in 2014, 2016, and 2018?
rvc_sites = rvc %>% distinct(REGION, STRAT, PRIMARY_SAMPLE_UNIT, 
                             STATION_NR, LAT_DEGREES, LON_DEGREES, YEAR, 
                             .keep_all = T)

# first convert fork length to total length using the bluestriped grunt
# length-length conversion from FishBase  TL = 0 + 1.034 x FL
rvc = rvc %>% mutate(HS_TOT_LEN = (LEN*1.034)) 

# extract the subadult bluestriped grunt stage (11.90 cm <= TOT LEN <= 25.33 cm 
# (between size at 1 YR and size at maturation))
hs_subadult = rvc %>%
  group_by(REGION, STRAT, YEAR, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES,
           LON_DEGREES, SPECIES_CD) %>% 
  filter(SPECIES_CD == "HAE SCIU" & (HS_TOT_LEN >= 11.90 & HS_TOT_LEN <= 25.33)) %>% 
  summarize(N = sum(NUM)) %>% # was a subadult bluestriped grunt ever seen at this SSU?
  ungroup() %>%
  distinct()
hs_subadult$LIFE_STAGE = rep("SUBADULT", nrow(hs_subadult)) 
hs_subadult$PRES = ifelse(hs_subadult$N > 0, 1, 0)
hs_subadult$SOURCE = rep("REEF VISUAL CENSUS", nrow(hs_subadult))

# now the inferred absences sites (sites where either no bluestriped grunt were 
# seen or only those of another age class were seen)
hs_subadult_abs = rvc %>%
  group_by(REGION, STRAT, YEAR, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES,
           LON_DEGREES, SPECIES_CD) %>%
  filter(SPECIES_CD == "HAE SCIU" & !(HS_TOT_LEN >= 11.90 & HS_TOT_LEN <= 25.33)) %>% 
  mutate(N = 0) %>% # manually assign a value of 0 b/c these are subadult absence sites
  ungroup() %>%
  distinct(REGION, STRAT, YEAR, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES,
           LON_DEGREES, SPECIES_CD, N)
hs_subadult_abs$LIFE_STAGE = rep("SUBADULT", nrow(hs_subadult_abs))
hs_subadult_abs$PRES = ifelse(hs_subadult_abs$N > 0, 1, 0)
hs_subadult_abs$SOURCE = rep("REEF VISUAL CENSUS", nrow(hs_subadult_abs))

# combining subadult absence and presence data and using distinct because some 
# of the absence sites might be shared across both the subadult stage-specific 
# data and the inferred absence data
hs_subadult_rvc = rbind(hs_subadult, hs_subadult_abs) %>%
  distinct(REGION, STRAT, YEAR, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES,
           LON_DEGREES, .keep_all = T)

# now keep only surveys that spatially coincide with bluestriped grunt MIR nodes
# convert rvc data to spatial data and reproject
hs_sf = st_as_sf(hs_subadult_rvc, coords = c(7, 6), crs = gcs) 
hs_sf = st_transform(hs_sf, crs = my_crs)

# intersect
hs_mir_surveys = st_intersection(hs_sf, hs_mir_nodes)

# calculate the proportion of positive observations (presences) per MIR site
hs_mir_prop = hs_mir_surveys %>%
  select(MIR_Site, SPECIES_CD, LIFE_STAGE, SOURCE, PRES, MIR_Site, 
         Id, IF_d5089) %>%
  group_by(MIR_Site) %>%
  mutate(N_PRES = sum(PRES),
         N_SURVEYS = length(PRES),
         PROP = N_PRES/N_SURVEYS) %>%
  distinct(MIR_Site, SPECIES_CD, LIFE_STAGE, SOURCE, MIR_Site, 
           Id, IF_d5089, N_PRES, N_SURVEYS, PROP)
cor(hs_mir_prop$PROP, hs_mir_prop$IF_d5089)

# write data to github repo
write.csv(hs_mir_prop, 
          paste0(git_wd, "Data/Subadult_Bluestriped_Grunt_Connectivity_vs_Proportion_Positive_Samples.csv"), 
          row.names = FALSE)

#### gray snapper (Lutjanus griseus) ####
# find nodes coinciding spatially with Mission Iconic Reefs sites
# all nodes
lg_nodes = st_read(dsn = fgdb, layer = "Subadult_Gray_Snapper_Batch_IF")

# 15 potential MIR sites
mir = st_read(dsn = fgdb, layer = "Available_Iconic_Reefs")

# keep only nodes with mir sites
lg_mir_nodes = lg_nodes %>%
  filter(Id %in% (st_intersection(st_buffer(mir, dist = 5), lg_nodes))$Id)

# add back in column with mir site names
# directory for which reef IDs match
lg_m = c(12871, "Davis Reef", 
         13196, "Cheeca Rocks", 
         10702, "French Reef",
         11161, "Molasses Reef",
         9500, "Key Largo Dry Rocks",
         8137, "Turtle Reef",
         15979, "Tennessee Reef", 
         18555, "South of Key Colony Beach (Marker 48)",
         9475, "Horseshoe Reef",
         18722, "Coffins Patch", 
         20483, "Newfound Harbor",
         12804, "Hen and Chickens", 
         20618, "Sombrero Reef", 
         21144, "Looe Key Reef", 
         9431, "Elbow Reef")
# save directory as data frame
lg_m_df = as.data.frame(matrix(lg_m, ncol = 2, byrow = T)) %>%
  rename(Id = V1,MIR_Site = V2)
# add and place in first column potition MIR site names
lg_mir_nodes$MIR_Site = lg_m_df$MIR_Site[match(lg_mir_nodes$Id, lg_m_df$Id)]
lg_mir_nodes = lg_mir_nodes %>% relocate(MIR_Site, .before = Id )

# remove temp data
rm(list = c("lg_m_df", "lg_m", "lg_nodes"))

# use the gray snapper MIR nodes to extract observation data from the
# Reef Visual Census (RVC) program [survey years 2014, 2016, 2018]
rvc = getSampleData(years = c(2014, 2016, 2018), regions = "FLA KEYS")

# how many unique sites were sampled in 2014, 2016, and 2018?
rvc_sites = rvc %>% distinct(REGION, STRAT, PRIMARY_SAMPLE_UNIT, 
                             STATION_NR, LAT_DEGREES, LON_DEGREES, YEAR, 
                             .keep_all = T)

# first convert fork length to total length using the gray snapper
# length-length conversion from FishBase (TL = 0 + 1.049 x FL)
rvc = rvc %>% mutate(LG_TOT_LEN = (LEN*1.049))

# extract the subadult gray snapper stage (9.51 cm <= TOT LEN <= 24.71 cm 
# (between size at 1 YR and size at maturation))
lg_subadult = rvc %>%
  group_by(REGION, STRAT, YEAR, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES,
           LON_DEGREES, SPECIES_CD) %>% 
  filter(SPECIES_CD == "LUT GRIS" & (LG_TOT_LEN >= 9.51 & LG_TOT_LEN <= 24.71)) %>% 
  summarize(N = sum(NUM)) %>% # was a subadult gray snapper ever seen at this SSU?
  ungroup() %>%
  distinct()
lg_subadult$LIFE_STAGE = rep("SUBADULT", nrow(lg_subadult)) 
lg_subadult$PRES = ifelse(lg_subadult$N > 0, 1, 0)
lg_subadult$SOURCE = rep("REEF VISUAL CENSUS", nrow(lg_subadult))

# now the inferred absences sites (sites where either no gray snapper were 
# seen or only those of another age class were seen)
lg_subadult_abs = rvc %>%
  group_by(REGION, STRAT, YEAR, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES,
           LON_DEGREES, SPECIES_CD) %>%
  filter(SPECIES_CD == "LUT GRIS" & !(LG_TOT_LEN >= 9.51 & LG_TOT_LEN <= 24.71)) %>% 
  mutate(N = 0) %>% # manually assign a value of 0 b/c these are subadult absence sites
  ungroup() %>%
  distinct(REGION, STRAT, YEAR, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES,
           LON_DEGREES, SPECIES_CD, N)
lg_subadult_abs$LIFE_STAGE = rep("SUBADULT", nrow(lg_subadult_abs))
lg_subadult_abs$PRES = ifelse(lg_subadult_abs$N > 0, 1, 0)
lg_subadult_abs$SOURCE = rep("REEF VISUAL CENSUS", nrow(lg_subadult_abs))

# combining subadult absence and presence data and using distinct because some 
# of the absence sites might be shared across both the subadult stage-specific 
# data and the inferred absence data
lg_subadult_rvc = rbind(lg_subadult, lg_subadult_abs) %>%
  distinct(REGION, STRAT, YEAR, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES,
           LON_DEGREES, .keep_all = T)

# now keep only surveys that spatially coincide with gray snapper MIR nodes
# convert rvc data to spatial data and reproject
lg_sf = st_as_sf(lg_subadult_rvc, coords = c(7, 6), crs = gcs) 
lg_sf = st_transform(lg_sf, crs = my_crs)

# intersect
lg_mir_surveys = st_intersection(lg_sf, lg_mir_nodes)

# calculate the proportion of positive observations (presences) per MIR site
lg_mir_prop = lg_mir_surveys %>%
  select(MIR_Site, SPECIES_CD, LIFE_STAGE, SOURCE, PRES, MIR_Site, 
         Id, IF_d7852) %>%
  group_by(MIR_Site) %>%
  mutate(N_PRES = sum(PRES),
         N_SURVEYS = length(PRES),
         PROP = N_PRES/N_SURVEYS) %>%
  distinct(MIR_Site, SPECIES_CD, LIFE_STAGE, SOURCE, MIR_Site, 
           Id, IF_d7852, N_PRES, N_SURVEYS, PROP)
cor(lg_mir_prop$PROP, lg_mir_prop$IF_d7852)


# write data to github repo
write.csv(lg_mir_prop, 
          paste0(git_wd, "Data/Subadult_Gray_Snapper_Connectivity_vs_Proportion_Positive_Samples.csv"), 
          row.names = FALSE)
