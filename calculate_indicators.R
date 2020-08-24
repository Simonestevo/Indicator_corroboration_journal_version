# cd C:\\Users\\ssteven\\Dropbox\\Deakin\\Chapter_2_Extinction_test\\Extinction_test_code

# Using R version 3.5.2

## This script outputs:
## 1: Biodiversity indicator values calculated by ecoregions and saved as an .rds
## and .csv files

## 2: Results of statistical analysis

## 3: Visualisation - maps and plots

#' TODO: IMPORTANT - resolve issues with multiple polygons = same ecoregion.  Either
#' give separate IDs (might be tricky for those species not assigned via spatial
#' data), OR remove duplicates before analysis (more likely to work) (can keep for
#' mapping though)
#' TODO: Figure out best way to test correlations between different indicators
#' and years because all years, all indicators is too many.

# Load packages ----

library(raster)
library(tidyverse)
library(sf)
library(grid)
library(viridis)
library(png)
library(gridExtra)
library(reshape2)
library(rlist)
library(GGally)
library(mapview)


# Set input and output locations ----

create_new_database_version <- FALSE # Only set to true if you want to create an entirely new version from scratch
date <- Sys.Date()
country <- "Australia" # If not subsetting, set as NA, e.g. country <- NA
inputs <- "N:/Quantitative-Ecology/Simone/extinction_test/inputs"
save_outputs <- "no"
parent_outputs <- "N:/Quantitative-Ecology/Simone/extinction_test/outputs"
eco_version <- "ecoregions_2017"
#eco_version <- "official_teow_wwf"

if (!is.na(country)) {
  
  location <- tolower(country)
  
} else {
  
  location <- "global"
  
}

# Set output directory

db_version <- tail(sort(list.files(parent_outputs)), 1)
db_interim <- list.files(file.path(parent_outputs,db_version))[
  grepl("interim",list.files(file.path(parent_outputs,db_version)))]
db_outputs <- list.files(file.path(parent_outputs,db_version))[
  grepl("database",list.files(file.path(parent_outputs,db_version)))]
ind_outputs <- list.files(file.path(parent_outputs,db_version))[
  grepl("indicator",list.files(file.path(parent_outputs,db_version)))]

interim_outputs <- file.path(parent_outputs, db_version, db_interim)
outputs <- file.path(parent_outputs, db_version, db_outputs)

if( (length(ind_outputs)) == 0 ) {
  
indicator_outputs <- file.path(parent_outputs, db_version, paste(date,
                               "_indicator_output_files",sep="") )

dir.create( indicator_outputs, recursive = TRUE ) # create a new directory for today's outputs


} else {

indicator_outputs <- file.path(parent_outputs, db_version, ind_outputs)

}

# Load functions ----

# Function to scale the values of a vector from 0 to 1

scale_to_1 <- function(vector){
  
  (vector-min(vector, na.rm = TRUE))/
    (max(vector, na.rm = TRUE)-min(vector, na.rm = TRUE))
}

# Funciton to calculate the Red List Index

calculate_red_list_index <- function(data, timeframe){
  
  require(tidyverse)
  
  # Remove data without RL status
  
  data$redlist_assessment_year <- as.numeric(as.character(data$redlist_assessment_year))
  
  data <- data %>%
    filter(!is.na(redlist_status)) %>%
    group_by(tsn) 
  
  ecoregion <- as.factor(data$ecoregion_code[1])
  
  # Assign category weights
  
  weighted_data <- data %>%
    dplyr::mutate(RL_weight = ifelse(redlist_status == "LC", 0,
                                     ifelse(redlist_status == "NT", 1,
                                            ifelse(redlist_status == "VU", 2,
                                                   ifelse(redlist_status == "EN", 3,
                                                          ifelse(redlist_status == "CR", 4,
                                                                 ifelse(redlist_status == "EX", 5, NA))))))) 
  
  #weighted_data$RL_weight <- as.numeric(as.character(weighted_data$RL_weight))
  
  # Filter out rows with NE and DD
  weighted_data <- filter(weighted_data, RL_weight != "NA" )
  
  # Calculate numerical weights for each species based on risk category
  # weight.data <- calcWeights(filter.data, RL_weight)
  # weight.data <- drop_na(weight.data, .data[[RL_weight]])
  
  # Group data so the index is calculated for each taxa for each year
  grouped.data <- weighted_data %>% group_by(class, redlist_assessment_year)
  
  # Sum category weights for each group, calculate number of species per group
  summed.weights <- summarise(grouped.data, 
                              total.weight = sum(RL_weight, na.rm = TRUE), # calc sum of all weights
                              total.count = n()) # calc number of species
  
  # Calculate RLI scores for each group, rounded to 3 decimal places
  
  index.scores <- summed.weights %>%
    mutate(RLI = 1 - (total.weight/(total.count * 5)), # actual RLI formula
           Criteria = "risk",
           Ecoregion_code = ecoregion)
  
  
  #index.scores <- index.scores[seq(1, nrow(index.scores), t), ]
  
  return(index.scores)
  
}

# Function that maps indicator values to ecoregions
#' TODO: Get rid of grid background
#' 

map_indicators <- function(data, indicator_variable, title, legend) {
  
indicator_map <-  ggplot(data) +
                  geom_sf(aes(fill = indicator_variable), colour = "black", 
                          size = 0.05, show.legend = 'fill') +
                  scale_fill_viridis_c(trans = "reverse", alpha = .8, 
                                       na.value = "grey70") +         
                  theme(axis.line = element_line(),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank()) +
                  labs(fill = title) +
                  theme(legend.position = legend) +
                  facet_wrap(~ year)

return(indicator_map)

}


# Get ecoregions and their attributes ----

ecoregion_map_all <- st_read(paste(inputs,eco_version, sep = "/"))

ecoregion_map <- ecoregion_map_all %>% dplyr::select(eco_code, ECO_NAME, 
                                                     OBJECTID, ECO_ID, 
                                                     geometry) 

# Pull out only required variables

ecoregion_map <- ecoregion_map_all %>% 
                 dplyr::select(ECO_ID, ECO_NAME, OBJECTID, geometry)

# names(ecoregion_map) <- c("ecoregion_code", "ecoregion_name","eco_objectid",
#                            "geometry")

# library(rgdal)
# data.shape<- st_read(dsn="N:/Quantitative-Ecology/Simone/extinction_test/inputs/official_teow_wwf",layer="terr_biomes")

# Get ecoregion countries
if (!("ecoregion_country_data.rds" %in% list.files(outputs))) { 
  
  country_map <- st_read(paste(inputs,"countries_WGS84", sep = "/")) %>%
    rename(country_objectid = OBJECTID)
  
  ecoregion_country_sf <-  st_join(ecoregion_map, country_map, 
                                   join = st_intersects)
  
  ecoregion_country_df <- as.data.frame(ecoregion_country_sf) %>%
    dplyr::select(-geometry) %>%
    arrange(country_objectid) %>%
    rename(eco_objectid = OBJECTID)
  
  saveRDS(ecoregion_country_df, file = file.path(outputs, 
                                                 "ecoregion_country_data.rds"))
  
} else {
  
  ecoregion_country_df <- readRDS(paste(outputs, "ecoregion_country_data.rds", 
                                        sep = "/"))
  
}

# # Add countries to the ecoregion map
# 
# ecoregion_map_og <- ecoregion_map %>%
#                  merge(ecoregion_country_df[c("ECO_ID", "CNTRY_NAME")],
#                        by = "ECO_ID") 

# Simplify the geometry of the ecoregion map slightly to speed up joins


# ecoregion_map_simple <- st_simplify(ecoregion_map, preserveTopology = TRUE,
#                                     dTolerance = 0.1)


# If subsetting by country, get the ecoregion IDs you want to work with 

if (!is.na(country)) {
  
  ecoregion_subset <- ecoregion_country_df %>%
    filter(CNTRY_NAME == country) %>%
    unique(.)
  
}

if (!is.na(country)) {
  
  ecoregion_map <- ecoregion_map[ecoregion_map$ECO_ID %in% 
                                   ecoregion_subset$ECO_ID,] 
  
}
  
# Read in species data

species_data <- readRDS(file.path(inputs, "deakin_species_data/species_data.rds"))

# sp_all <- species_data

# Remove cases without redlist status and/or assessment year
#' TODO: this might be better done in the build_database script? So there's no
#' wrangling of the species data in this script?

species_data <- species_data %>%
                drop_na(redlist_assessment_year, redlist_status, class) %>%
                distinct(.) %>%
                mutate(class = replace(class, class == "AMPHIBIA", "Amphibia"),
                       class = replace(class, class == "MAMMALIA", "Mammalia"),
                       class = replace(class, class == "FLORIDEOPHYCEAE", 
                                                       "Florideophycae"),
                       class = replace(class, class == "REPTILIA", "Reptilia"),
                       class = replace(class, class == "MAGNOLIOPSIDA", 
                                                       "Magnoliopsida"),
                       class = replace(class, class == "AVES", "Aves"))


# Subset by test country

if (!is.na(country)) {
  
    # species_data <- species_data %>%
    # merge(ecoregion_subset[c("ecoregion_code", "country_name")], 
    #       by = "ecoregion_code", all = FALSE) %>%
    # unique(.)
    
    species_data <- species_data[species_data$ecoregion_code %in% 
                                   ecoregion_map$ecoregion_code, ]
  
}


# Calculate indicators ----


# Extinction/Risk status ----

species_by_ecoregion <- species_data %>%
                        select(-eco_objectid) %>%
                        distinct(.) %>%
                        group_by(ecoregion_code, redlist_assessment_year) %>%
                        mutate(number_of_species = n_distinct(tsn),
                               number_extinct = n_distinct(tsn[redlist_status == "EX"|
                                                               redlist_status == "EW"]),
                               number_atrisk = n_distinct(tsn[redlist_status == "EN"|
                                                      redlist_status == "CR"|
                                                      redlist_status == "CR(PE)"|
                                                      redlist_status == "VU"]),
                               number_lowrisk = n_distinct(tsn[redlist_status == "LC"|
                                                       redlist_status == "NT"]),
                               number_datadeficient = n_distinct(tsn[redlist_status == "DD"]),
                               check = number_of_species - (number_extinct +
                                                            number_atrisk +
                                                            number_lowrisk+
                                                            number_datadeficient),
                               proportion_extinct = number_extinct/number_of_species,
                               proportion_atrisk = number_atrisk/number_of_species,
                               proportion_lowrisk = number_lowrisk/number_of_species)
                               

# Extinctions

extinction_values <- species_by_ecoregion %>%
                     mutate(indicator = "proportion extinct") %>%
                     select(indicator, redlist_assessment_year, 
                            ecoregion_code, proportion_extinct) %>%
                     rename(year = redlist_assessment_year, 
                            raw_indicator_value = proportion_extinct)

extinction_values <- as.data.frame(extinction_values)

indicator_columns <- names(extinction_values)

# Threatened status

at_risk_values <- species_by_ecoregion %>%
                 mutate(indicator = "proportion at risk") %>%
                 select(indicator, redlist_assessment_year, 
                         ecoregion_code, proportion_atrisk) %>%
                 rename(year = redlist_assessment_year, 
                         raw_indicator_value = proportion_atrisk)

at_risk_values <- as.data.frame(at_risk_values)

# # Red List Index ----

# Split into different classes

species_data_for_rli <- species_data %>%
                        filter(class == "Amphibia"| class == "Aves"|
                               class == "Mammalia")

class_list <- split(species_data_for_rli, species_data_for_rli$class)

# Split each class into different time points (output should be dataframes of 
# species red list status in a nested list with levels: Class, Time point, Ecoregion)

class_time_list_full <- list()
class_time_ecoregion_list <- list()

for (i in seq_along(class_list)) {
  
  # Split by different assessment timepoints
  
  class_timestep <- split(class_list[[i]], 
                          class_list[[i]]$redlist_assessment_year)
  
    for (j in seq_along(class_timestep)) {
      
    # Split by each ecoregion (so have a list of species for each class, timepoint,
    # ecoregion)
    
    class_time_ecoregion_list[[j]] <- split(class_timestep[[j]], 
                                            class_timestep[[j]]$ecoregion_code)
    }
    
    class_time_list_full[[i]] <- class_time_ecoregion_list
}

# Remove the empty lists with no data in them

class_time_list <- list()

for (i in seq_along(class_time_list_full)) {

class_time_list[[i]] <- class_time_list_full[[i]][lengths(class_time_list_full[[i]]) != 0]

}

# Calculate the RLI per ecoregion, per timepoint, per class (output should be
# dataframes of RLI values by ecoregion, by timepoint, split by class)

class_time_ecoregions <- list()
classes_rli <- list()


for (i in seq_along(class_time_list)) {
  
  class <- class_time_list[[i]] # Get list of timesteps for one class
  
  class_all_timepoints <- list()
  
    for (j in seq_along(class)) {
    
      time <- class[[j]] # Get list of ecoregions for one timestep, for one class
      
        for (k in seq_along(time)) {
        
          class_time_ecoregions[[k]] <- calculate_red_list_index(time[[k]]) # Calculate RLI for each ecoregion for one timestep, one class
          
          class_time_ecoregion_df <- do.call(rbind, class_time_ecoregions) # One dataframe for one time point, one class
          
        } 
    
    class_all_timepoints[[j]] <- class_time_ecoregion_df # Put the time point into a list of timepoints
    
    }
  
  class_all_timepoints_df <- do.call(rbind, class_all_timepoints)
  
  classes_rli[[i]] <- class_all_timepoints_df # Put the list of time points into a class
  
}

# Convert back into one dataframe

rli_values <- do.call(rbind, classes_rli)

# Add scaled/inverted values

rli_values <- rli_values %>%
              filter(RLI != 0) %>%
              mutate(indicator = paste("red list index", class, sep = " ")) %>%
              rename(year = redlist_assessment_year) %>%
              rename(raw_indicator_value = RLI) %>%
              rename(ecoregion_code = Ecoregion_code) %>%
              ungroup()%>%
              dplyr::select(indicator, year, ecoregion_code, raw_indicator_value) 
                  # mutate(RLI_inverted = 1 - RLI) %>%
                  # mutate(RLI_scaled_inverted =
                  #          scale_to_1(RLI_inverted)) %>%
                  # mutate(RLI_adjusted_old = ifelse(RLI == 0, NA,
                  #                           ifelse(RLI > 0 & RLI < 0.9538,
                  #                                         0.9538, RLI))) %>%
                  # mutate(RLI_adjusted = ifelse(RLI == 0, NA,
                  #                              pmin(pmax(RLI,
                  #                       quantile(RLI, .05, na.rm = TRUE))))) %>%
                  # mutate(RLI_adjusted_inverted = 1 - RLI_adjusted)


# Human Footprint Index 2017 ----

# Read in the Human Footprint Index data 

hfp_by_ecoregion_2017 <- read.csv(paste(inputs, "human_footprint_index", 
                                        "human_footprint_index_by_ecoregion.csv",
                                        sep = "/"))

hfp_by_ecoregion_2017 <- hfp_by_ecoregion_2017 %>%
                         rename(ecoregion_name = ECO_NAME) %>%
                         merge(ecoregion_map[c("ecoregion_name", "ecoregion_code")], 
                               by = "ecoregion_name") %>% # This will subset automatically if you subset by country
                         mutate(year = "2017") %>%
                         mutate(indicator = "human footprint index") %>%
                         rename(raw_indicator_value = HFP) %>%
                         dplyr::select(indicator, year, ecoregion_code, 
                                        raw_indicator_value) 
                         # mutate(HFP_original = HFP) %>%
                         # mutate(HFP_adjusted_old = ifelse(HFP_original > 33.29037,
                         #                                   33.29037, HFP_original)) %>%
                         # mutate(HFP_adjusted =
                         #           pmin(pmax(HFP_original,quantile(HFP_original,
                         #                                           .005, na.rm = TRUE)),
                         #                quantile(HFP_original, .995,
                         #                         na.rm = TRUE))) %>%
                         # mutate(HFP_scaled_adjusted = scale_to_1(HFP_adjusted)) %>%
                         # mutate(HFP_scaled_adjusted_inverted =
                         #           1 - HFP_scaled_adjusted)


# Subset the HFP data by country 

# if (!is.na(country)) {
#   
#   full_hfp_by_ecoregion_2017 <- hfp_by_ecoregion_2017
#   
#   country_ecoregions <- ecoregion_map$ecoregion_code
#   
#   hfp_by_ecoregion_2017 <- full_hfp_by_ecoregion_2017[full_hfp_by_ecoregion_2017$ecoregion_code %in% country_ecoregions, ]
#   
# }

# Biodiversity Habitat Index Plants 2015 ----

bhi_plants_2015_all <- read.csv(file.path(inputs, 
                                       "biodiversity_habitat_index\\BHI2015_PLANTS_BY_ECOREGION.csv"))

bhi_plants_2015_values <- bhi_plants_2015_all[-1,1:2]

bhi_plants_2015_values <- bhi_plants_2015_values %>%
                          rename(eco_id = BIODIVERSITY.HABITAT.INDEX.2015,
                                 raw_indicator_value = X) %>%
                          mutate(raw_indicator_value = 
                          as.numeric(levels(raw_indicator_value))[raw_indicator_value]) %>%
                          merge(ecoregion_map[c("eco_id","ecoregion_code")], 
                                by = "eco_id") %>%
                          mutate(indicator = "biodiversity habitat index plants",
                                 year = 2015) %>%
                          dplyr::select(-geometry) %>%
                          dplyr::select(indicator_columns) %>%
                          distinct(.)

# Living Planet Index ----

lpi_data <- read.csv(file.path(inputs, 
                               "living_planet_index\\LPR2018data_public.csv"))

                          

# Biodiversity Intactness Index 2005 ----

# TEMPORARY CODE - the below takes ages to run on a laptop, so currently just
# reading in prepped bii data for development, but put the below chunk back on 
# when doing the real deal analysis

if (!is.na(country)) {

bii_2005_rich_data <- raster(file.path(inputs,
                                         "biodiversity_intactness_index\\bii_rich_2005_aus.tif"))

} else {

bii_2005_rich_data <- raster(file.path(inputs,
                        "biodiversity_intactness_index\\final-rich-bii-isl-main.tif"))

}
#' 
#' #' TODO: Work out why this is producing so many NaNs
#' 
bii_rich_ecoregion_map <- ecoregion_map %>%
                          mutate(raw_indicator_value = 
                                 raster::extract(bii_2005_rich_data,
                                                  ecoregion_map,
                                                  fun = mean, 
                                                 na.rm = TRUE))

#' saveRDS(bii_ecoregion_map, file.path(outputs, "bii_2005_aus_ecoregion_map.rds"))

#' bii_ecoregion_map <- readRDS(file.path(outputs, "bii_2005_aus_ecoregion_map.rds"))
#' 
#' #' TODO: Figure out how to deal with multiple polygons of the same ecoregion - once
#' #' you remove the geometry you end up with multiple values p/ecoregion
#' 
bii_values <- bii_ecoregion_map %>%
              st_set_geometry(NULL) %>%
              mutate(indicator = "biodiversity intactness index") %>%
              mutate(year = "2005") %>%
              #rename(ecoregion_code = eco_code) %>%
              select(names(rli_values)) %>%
              drop_na()

# Wilderness Intactness Index ----

# wii_map_data <- st_read("N:\\Quantitative-Ecology\\Simone\\extinction_test\\inputs\\intactness\\intactness\\Ecoregions2017_intactness.shp")
# 
# wii_data <- as.data.frame(wii_map_data) %>%
#   dplyr::select(OBJECTID, ECO_NAME, ECO_ID, Q2009_LL, Q2009, Q2009_UL, PLOTCAT)

# Analyse indicators ----

ecoregions <- as.data.frame(ecoregion_map) %>% dplyr::select(-geometry)

# Combine indicator values into a single dataframe ----

# TODO: Melt the adjusted values into long form
# TODO: Add the BII back in when we've dealt with multiple ecosystem issue

indicator_values <- rbind(rli_values, hfp_by_ecoregion_2017, extinction_values, 
                          at_risk_values, bhi_plants_2015_values)

indicator_values <- indicator_values %>%
                    mutate(HFP_adjusted_old = ifelse(raw_indicator_value > 33.29037,
                                                     33.29037, raw_indicator_value)) %>%
                    mutate(RLI_adjusted_old = ifelse(raw_indicator_value == 0, NA,
                                   ifelse(raw_indicator_value > 0 & 
                                            raw_indicator_value < 0.9538,
                                          0.9538, raw_indicator_value))) %>%
                    mutate(year = as.numeric(year)) %>%
                    mutate(HFP_adjusted_old = ifelse(indicator != 
                                              "human footprint index", NA,
                                              HFP_adjusted_old)) %>%
                    mutate(RLI_adjusted_old = ifelse(grepl("red list index Aves",
                                                    indicator), 
                                                    RLI_adjusted_old, NA)) %>%
                    merge(ecoregion_country_df[c("ecoregion_code", 
                                                 "country_name")], 
                          by = "ecoregion_code", all.x = TRUE, all.y = FALSE) %>%
                    distinct(.) %>%
                    mutate(indicator_year = paste(indicator, year, sep = " "))
                    # %>%
                    # mutate(HFP_adjusted =
                    #          pmin(pmax(raw_indicator_value,quantile(raw_indicator_value,
                    #                                          .005, na.rm = TRUE)),
                    #               quantile(raw_indicator_value, .995,
                    #                        na.rm = TRUE))) %>%
                    # mutate(HFP_scaled_adjusted = scale_to_1(HFP_adjusted)) %>%
                    # mutate(HFP_scaled_adjusted_inverted =
                    #          1 - HFP_scaled_adjusted) %>%
                    # mutate(scaled = scale_to_1(raw_indicator_value)) %>%
                    # mutate(inverted = 1 - raw_indicator_value) %>%
                    # mutate(scaled_inverted = 1 - scaled) %>%
                    # mutate(RLI_adjusted = ifelse(raw_indicator_value == 0, NA,
                    #                              pmin(pmax(raw_indicator_value,
                    #                   quantile(raw_indicator_value, .05, 
                    #                            na.rm = TRUE))))) %>%
                    # mutate(RLI_adjusted_inverted = 1 - RLI_adjusted)

# Test for correlations between indicators ----

# All years
#' TODO: Remove the crap years with hardly any data

indicator_values_2 <- indicator_values %>%
                      dplyr::select(ecoregion_code, indicator_year, 
                                    raw_indicator_value) %>%
                      filter(indicator_year != 
                               "biodiversity intactness index 2005") %>%
                      distinct(.)


indicator_values_wide <- indicator_values_2 %>%
                         spread(key = indicator_year, 
                                value = raw_indicator_value) 


names(indicator_values_wide) <- make.names(names(indicator_values_wide), 
                                           unique = TRUE)
# TEMPORARY - select columns that have enough data

indicator_values_wide <- indicator_values_wide %>%
                         select(biodiversity.habitat.index.plants.2015,
                                 human.footprint.index.2017,
                                 proportion.at.risk.1988,
                                 proportion.at.risk.1994,
                                 proportion.at.risk.1996,
                                 proportion.at.risk.2000,
                                 proportion.at.risk.2004,
                                 proportion.at.risk.2008,
                                 proportion.at.risk.2012,
                                 proportion.at.risk.2016,
                                 red.list.index.Amphibia.2004,
                                 red.list.index.Amphibia.2008,
                                 red.list.index.Aves.1988,
                                 red.list.index.Aves.1994,
                                 red.list.index.Aves.2000,
                                 red.list.index.Aves.2004,
                                 red.list.index.Aves.2008,
                                 red.list.index.Aves.2012,
                                 red.list.index.Aves.2016,
                                 red.list.index.Mammalia.1996,
                                 red.list.index.Mammalia.2008) %>%
                          na.omit(.)


#indicator_values_matrix <- as.matrix(indicator_values_wide[, -1])

# indicator_values_wide_complete  <- as.matrix(indicator_values_matrix[
#                                    complete.cases(indicator_values_matrix),])

correlations_all_years <- as.data.frame(cor(indicator_values_wide, 
                                            method = "pearson"))


indicators_for_scatterplots <- indicator_values_wide %>%
                               #dplyr::select(- ecoregion_code) %>%
                               mutate(human.footprint.index.2017.scaled = 
                                     scale_to_1(human.footprint.index.2017)) %>%
                               dplyr::select(-human.footprint.index.2017)

summary(indicators_for_scatterplots)

indicator_scatterplots <- ggpairs(indicators_for_scatterplots)

if (save_outputs == "yes") {
  
  ggsave(file.path(outputs, "scatterplot_all_indicators_years.png"),
         indicator_scatterplots,  device = "png")
  
}

# Split by year


indicator_values_time_list <- split(indicator_values, indicator_values$year)


if (aggregate_timepoints == "yes") {
  
#' TODO: find a programmatical way to select which years get grouped together

# Aggregate years that are close to each other

indicators_1988 <- indicator_values_time_list[[1]]

indicators_1994_1996 <- rbind(indicator_values_time_list[[2]], 
                              indicator_values_time_list[[3]])

indicators_2000 <- indicator_values_time_list[[4]]

indicators_2004_2005 <- rbind(indicator_values_time_list[[5]])

indicators_2008 <- indicator_values_time_list[[6]]

indicators_2012 <- indicator_values_time_list[[7]]

indicators_2016_2017 <- rbind(indicator_values_time_list[[8]],
                           indicator_values_time_list[[9]])

indicator_values_time_list <- list(indicators_1988, indicators_1994_1996,
                                   indicators_2000, indicators_2004_2005,
                                   indicators_2008, indicators_2012, 
                                   indicators_2016_2017)

rm(indicators_1988, indicators_1994_1996,
        indicators_2000, indicators_2004_2005,
        indicators_2008, indicators_2012, 
        indicators_2016_2017)

}




# Pairwise correlations and scatterplots ----

#' TODO: Can we use ggpair with other ggplot syntax? Including attributes but
#' also transformations? https://www.r-graph-gallery.com/199-correlation-matrix-with-ggally.html

correlations <- list()
scatterplots <- list()
years <- list()

for (i in seq_along(indicator_values_time_list)) {
  
  # Get the year name
  
  year <- range(indicator_values_time_list[[i]]$year)
  
  # Remove year column from dataframe
  
  step1 <- indicator_values_time_list[[i]] %>% 
           select(ecoregion_code, 
                  raw_indicator_value, 
                  indicator_year, 
                  -year) %>%
           distinct(.)
  
  # Convert into wide format
  
  # step2 <- dcast(step1, ecoregion_code ~ indicator, 
  #                value.var = "raw_indicator_value", sum)
  
  step2 <- step1 %>%
           spread(key = indicator_year, 
           value = raw_indicator_value) 
  
  # Convert into matrix without id variables
  
  step3 <- as.matrix(step2[2:ncol(step2)])
  
  # Check if there's enough indicators to analyse correlation (ie more than 2)
  
  if(ncol(step3) == 1) {
    
    print(paste("Year", year, "does not have values for more than one indicator to test correlations between"))
    
    rm(step1, step2, step3)
    
  } else {
    
  # Analyse correlations between indicators 
  
  step4 <- step3[complete.cases(step3),]
  
  correlations[[i]] <- as.data.frame(cor(step4, method = "pearson"))
  
  scatterplots[[i]] <- ggpairs(as.data.frame(step4))
  
  years[[i]] <- year
  
  rm(step1, step2, step3, step4)
  
  }
  
}

# Remove the empty lists with no data in them

names(correlations) <- years

correlations <- list.clean(correlations, fun = is.null, recursive = FALSE)

names(scatterplots) <- years

if(save_outputs == "yes") {

for (i in seq_along(scatterplots)) {
  
    ggsave(file.path(outputs, paste(as.character(years[[i]]),"_scatterplot",
                                    ".png", sep = "")), 
                     scatterplots[[i]],  device = "png")
    
  }
}

# Map the indicators ----

indicator_map_data <- left_join(ecoregion_map, indicator_values, 
                                 by = "ecoregion_code")


# Map the BHI
bhi_map <- indicator_map_data %>%
           filter(indicator == "biodiversity habitat index plants") %>%
           map_indicators(.$raw_indicator_value,
                          "BHI plants", 
                          "right")
bhi_map

# Create a separate sf object to save as shapefile

bhi_2015_plants_sf <- indicator_map_data %>%
                      filter(indicator == "biodiversity habitat index plants") %>%
                      select(eco_id, ecoregion_name, eco_objectid,
                             raw_indicator_value, geometry) %>%
                      rename(BHI = raw_indicator_value)

if (save_outputs == "yes") {
  
  ggsave(file.path(outputs, "2020-07-27_bhi_plants_2015_map.png"), 
         bhi_map,  device = "png")
  
  st_write(bhi_2015_plants_sf, file.path(outputs, 
                                         paste(date, 
                                               "_bhi_plants_2015_map.shp")))
  saveRDS(bhi_map, file.path(outputs, "bhi_2015_ecoregion_map.rds"))
  
}

test <- st_read("N:\\Quantitative-Ecology\\Simone\\extinction_test\\outputs\\2020-07-27_indicator_output_files")



# Map species risk status

at_risk_map <- indicator_map_data %>%
               filter(indicator == "proportion at risk") %>%
               map_indicators(.$raw_indicator_value,
                             "Proportion\nof species\nat risk", 
                             "right")

at_risk_map

if(save_outputs == "yes") {
  
ggsave(file.path(outputs, "at_risk_map.png"), 
at_risk_map,  device = "png")
    

}

extinct_map <- indicator_map_data %>%
               filter(indicator == "proportion extinct") %>%
               map_indicators(.$raw_indicator_value,
                               "Proportion\nof species\nextinct", 
                               "right")

extinct_map

if(save_outputs == "yes") {
  
  ggsave(file.path(outputs, "extinct_map.png"), 
         extinct_map,  device = "png")
  
  
}
  
  
# Map the Red List Index by class

#' TODO: Put this into a loop - currently doesn't work because they have
#' different timesteps

birds_rli_map <- indicator_map_data %>%
                 filter(indicator == "red list index Aves") %>%
                 map_indicators(.$RLI_adjusted_old,
                                "Red List\nIndex (Birds)", 
                                "right")

birds_rli_map

mammals_rli_map <- indicator_map_data %>%
                   filter(indicator == "red list index Mammalia") %>%
                   map_indicators(.$raw_indicator_value,
                                   "Red List\nIndex (Mammals)", 
                                   "right")

mammals_rli_map

amphibians_rli_map <- indicator_map_data %>%
                      filter(indicator == "red list index Amphibia") %>%
                      map_indicators(.$raw_indicator_value,
                                     "Red List\nIndex (Amphibians)", 
                                     "right")

amphibians_rli_map

# Map the Human Footprint Index

hfp_map <- indicator_map_data %>%
           filter(indicator == "human footprint index") %>%
           map_indicators(.$HFP_adjusted_old,
                          "Human\nFootprint\nIndex",
                          "right")
hfp_map

#' TODO: Turn this back on when we've figures out bii
# Map the Biodiversity Intactness Index

# bii_map <- indicator_map_data %>%
#            filter(indicator == "biodiversity intactness index") %>%
#            map_indicators(.$raw_indicator_value,
#                            "Biodiversity\nIntactness\nIndex",
#                            "right")


# Look at the data distribution

hfp_hist <- hist(indicator_values$HFP, breaks = 100)
hfp_hist # Values close to 0 = good, close to 30 = bad

rli_hist <- hist(indicator_values$RLI, breaks = 100)
rli_hist # Values close to 1 = good, close to 0 = bad

# wii_hist <- hist(indicator_values$WII_inverted, breaks = 100)
# wii_hist
# 
# indicator_values$PLOTCAT <- as.numeric(as.character(indicator_values$PLOTCAT))
# 
# wii_category_hist <- hist(indicator_values$PLOTCAT, breaks = 10)
# wii_category_hist


# RLI vs HFP scatterplot ----

# Make the scatterplot with colours mapped to RLI, size mapped to HFP

rli_hfp <- ggplot(indicator_values, aes(x = HFP_adjusted,
                                        y = RLI_scaled,
                                        size = HFP_adjusted,
                                        col = RLI_scaled)) + 
  #guides(size = guide_legend(reverse = TRUE)) +
  # annotation_custom(rasterGrob(gradient_background, 
  #                              width = unit(1,"npc"), 
  #                              height = unit(1,"npc")), 
  #         -Inf, Inf, -Inf, Inf) +
  geom_point(alpha = 0.8) +
  scale_colour_viridis_c(direction = -1 ) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  # scale_y_continuous(breaks = c(0, 1),
  #                    labels = c("Low extinction risk", "High extinction risk")) +
  # scale_x_continuous(breaks = c(0, 1),
  #                    labels = c("Small footprint", "Large footprint")) +
  labs(x = "Human Footprint Index (scaled)", 
       y = "Red List Index (Birds)",
       size = "Size of human\nfootprint\n(small to large)",
       col = "Extinction risk\n(low to high)") +
  #geom_tile() +
  #scale_colour_viridis() +
  theme(plot.background = element_rect(fill = NA),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14)) + 
  scale_y_continuous(breaks = c(0.00,0.25,0.50,0.75,1.00),
                     limits = c(0.00,1.00)) +
  scale_x_continuous(breaks = c(0,5,10,15,20,25,30))

rli_hfp

if(save_outputs == "yes") {
  
  ggsave(file.path(outputs, "hfp_rli_scatterplot.png"), rli_hfp,  device = "png")
  
}

# Make the same scatterplot with no aes mapped, but gradient background instead

gradient_background <- readPNG(file.path(inputs, "gradient_inverted_cropped.png"))

rli_hfp_2 <- ggplot(indicator_values, aes(x = HFP_adjusted,
                                          y = RLI_scaled)) + 
  guides(size = guide_legend(reverse = TRUE)) +
  annotation_custom(rasterGrob(gradient_background,
                               width = unit(1,"npc"),
                               height = unit(1,"npc")),
                    -Inf, Inf, -Inf, Inf) +
  geom_point(size = 2) +
  theme(axis.line = element_line(colour = "black")) +
  # scale_y_continuous(breaks = c(0, 1),
  #                    labels = c("Low extinction risk", "High extinction risk")) +
  # scale_x_continuous(breaks = c(0, 1),
  #                    labels = c("Small footprint", "Large footprint")) +
  labs(x = "Human Footprint Index (scaled)", 
       y = "Red List Index (Birds)",
       size = "Size of human\nfootprint\n(small to large)",
       col = "Extinction risk\n(low to high)") +
  theme(plot.background = element_rect(fill = NA),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14)) + 
  scale_y_continuous(breaks = c(0.00,0.25,0.50,0.75,1.00),
                     limits = c(0.00,1.00)) +
  scale_x_continuous(breaks = c(0,5,10,15,20,25,30))

rli_hfp_2

if(save_outputs == "yes") {
  
  ggsave(file.path(outputs, "hfp_rli_scatterplot_gradient.png"), 
         rli_hfp_2, height = 4, width = 5, device = "png")
  
}

horizontal_gradient <-  readPNG(file.path(inputs, "horizontal_gradient.png"))

rli_hfp_3 <- ggplot(indicator_values, aes(x = HFP_adjusted,
                                          y = RLI_scaled,
                                          col = RLI_scaled)) + 
  annotation_custom(rasterGrob(horizontal_gradient,
                               width = unit(1,"npc"),
                               height = unit(1,"npc")),
                    -Inf, Inf, -Inf, Inf) +
  geom_point(shape = 21, colour = "white", size = 4,
             aes(fill = RLI_scaled)) +
  scale_fill_viridis_c(direction = -1 ) +
  theme(axis.line = element_line(colour = "black")) +
  labs(x = "Human Footprint Index (scaled)", 
       y = "Red List Index (Birds)",
       size = "Size of human\nfootprint\n(small to large)",
       col = "Extinction risk\n(low to high)") +
  theme(plot.background = element_rect(fill = NA),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14)) + 
  scale_y_continuous(breaks = c(0.00,0.25,0.50,0.75,1.00),
                     limits = c(0.00,1.00)) +
  scale_x_continuous(breaks = c(0,5,10,15,20,25,30)) +
  theme(legend.position = "none")

rli_hfp_3

if(save_outputs == "yes") {
  
  ggsave(file.path(outputs, "hfp_rli_scatterplot_multi_colour_scales.png"), 
         rli_hfp_3, height = 4, width = 5, device = "png")
  
}


## Previous way of calculating RLI used for NESP/WWF docs ----

# species_data <- species_data %>% 
#   filter(class == "Aves") %>%
#   filter(redlist_assessment_year == 2016)
# 
# species_data_by_ecoregion <- split(species_data, 
#                                    species_data$ecoregion_code)
# 
# bird_summary <- species_data %>%
#   group_by(ecoregion_code, redlist_status) %>%
#   summarize(n_distinct(tsn))
# 
# # Check the outliers
# 
# eco_NA0606 <- species_data %>%
#   filter(ecoregion_code == "NA0606")
# 
# 
# # Calculate the Red List Index for each group, for each timeframe, for each ecoregion
# 
# rli_by_ecoregion <- list()
# 
# for (i in seq_along(species_data_by_ecoregion)) {
#   
#   rli_by_ecoregion[[i]] <- calculate_red_list_index(species_data_by_ecoregion[[i]])
#   
# }
# 
# # Convert back into a dataframe
# 
# birds_rli_by_ecoregion_2016 <- do.call(rbind, rli_by_ecoregion)
# 
# 
# birds_rli_by_ecoregion_2016 <- birds_rli_by_ecoregion_2016 %>%
#   rename(eco_code = Ecoregion_code) %>%
#   filter(RLI != 0) %>%
#   mutate(RLI_scaled = scale_to_1(RLI)) %>%
#   mutate(RLI_inverted = 1 - RLI) %>%
#   mutate(RLI_scaled_inverted = 
#            scale_to_1(RLI_inverted)) %>%
#   mutate(RLI_adjusted_old = ifelse(RLI == 0, NA, 
#                                    ifelse(RLI > 0 & RLI < 0.9538, 
#                                           0.9538, RLI))) %>%
#   mutate(RLI_adjusted = ifelse(RLI == 0, NA, 
#                                pmin(pmax(RLI,
#                                          quantile(RLI, .05, na.rm = TRUE))))) %>% 
#   # quantile(RLI, .995, na.rm = TRUE)))) %>%
#   # mutate(RLI_scaled = scale_to_1(RLI)) %>%
#   mutate(RLI_adjusted_inverted = 1 - RLI_adjusted)

