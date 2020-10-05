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

# library(devtools)
# Install from main ZSL repository online
# install_github("Zoological-Society-of-London/rlpi", dependencies=TRUE)
library(rlpi)
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
library(psych)
library(e1071)
library(tm)



# Set input and output locations ----

create_new_database_version <- FALSE # Only set to true if you want to create an entirely new version from scratch
date <- Sys.Date()
country <- NA #"Australia" # If not subsetting, set as NA, e.g. country <- NA
inputs <- "N:/Quantitative-Ecology/Simone/extinction_test/inputs"
save_outputs <- "yes" #only applies to maps, other things will always save
parent_outputs <- "N:/Quantitative-Ecology/Simone/extinction_test/outputs"
eco_version <- "ecoregions_2017"
#eco_version <- "official_teow_wwf"
indicator_columns <- c("indicator", "year", "ecoregion_id", "raw_indicator_value")

# Set up some ecoregions that we know how they should behave

east_australia <- 168 # Decline over time
amazon <- 473 # Decline over time
cardamom <- 223 # In good shape
mascarene <- 20 #Had a lot of extinctions


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

calculate_red_list_index <- function(data){
  
  require(tidyverse)
  
  # Remove data without RL status
  
  data$redlist_assessment_year <- as.numeric(as.character(data$redlist_assessment_year))
  
  data <- data %>%
    filter(!is.na(redlist_status)) %>%
    group_by(tsn) 
  
  ecoregion <- as.factor(data$ecoregion_id[1])
  
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
                              total.count = n(),
                              .groups = "drop_last") # calc number of species
  
  # Calculate RLI scores for each group, rounded to 3 decimal places
  
  index.scores <- summed.weights %>%
    mutate(RLI = 1 - (total.weight/(total.count * 5)), # actual RLI formula
           Criteria = "risk",
           Ecoregion_id = ecoregion)
  
  
  #index.scores <- index.scores[seq(1, nrow(index.scores), t), ]
  
  return(index.scores)
  
}

produce_scatterplots <- function(indicator_values, name, save) {
  
  scatterplot_directory <- file.path(indicator_outputs, 
                                     "scatterplots")
  
  if( !dir.exists( scatterplot_directory ) ) {
    
    dir.create( scatterplot_directory, recursive = TRUE )
    
  }
  
    indicator_values_2 <- indicator_values %>%
    dplyr::select(ecoregion_id, indicator_year,
                  raw_indicator_value) %>%
    distinct(.)
  
  
  indicator_values_wide <- indicator_values_2 %>%
    spread(key = indicator_year, 
           value = raw_indicator_value) 
  
  
  names(indicator_values_wide) <- make.names(names(indicator_values_wide), 
                                             unique = TRUE)
  # TEMPORARY - select columns that have enough data
  
  indicator_values_wide <- indicator_values_wide %>%
    select(-ecoregion_id) %>%
    na.omit(.)
  
  # indicator_values_wide[,c(2:ncol(indicator_values_wide))] <- 
  #   scale(indicator_values_wide[,c(2:ncol(indicator_values_wide))])
  
  correlations_all_years <- as.data.frame(cor(indicator_values_wide, 
                                              method = "pearson"))
  
  
  indicators_for_scatterplots <- indicator_values_wide 
  
  summary(indicators_for_scatterplots)
  
  indicator_scatterplots <- ggpairs(indicators_for_scatterplots)
  
  if (save == TRUE) {
  
  ggsave(file.path(scatterplot_directory, paste(location, eco_version,
                                                name, "scatterplots.png", sep = "_")),
         indicator_scatterplots,  device = "png")
  }
  
  return(indicator_scatterplots)
  
}


# Function that maps indicator values to ecoregions
#' TODO: Get rid of grid background
#' 

map_indicators <- function(data, indicator_variable, title, legend) {

  
indicator_map <-  ggplot(data) +
                  geom_sf(aes(fill = indicator_variable), colour = "black", 
                          size = 0.05, show.legend = 'fill') +
                  scale_fill_viridis_c(trans = "reverse",
                                       alpha = .8,
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

if (("Ecoregions2017Valid.rds" %in% list.files(file.path(inputs, 
                                                          "ecoregions_2017")))){
  
ecoregion_map_all <- readRDS(paste(file.path(inputs, "ecoregions_2017"),
                                   "Ecoregions2017valid.rds"))

} else {
  
# Read in shape file

ecoregion_map_all <- st_read(paste(inputs,eco_version, sep = "/"))

# Fix non-valid geometries - takes a little while

ecoregion_map_all <- st_make_valid(ecoregion_map_all)

saveRDS(ecoregion_map_all, file.path(paste(inputs, "ecoregions_2017", sep = "/"),
                                 "Ecoregions2017Valid.rds"))

}

# Pull out only required variables

ecoregion_map <- ecoregion_map_all %>% 
                 select(ECO_ID, ECO_NAME, OBJECTID, REALM, geometry)

# Check geometry and fix if needed



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

# Check if indicators have already been calculated ----

# if ((paste(location, eco_version, "indicator_values_master.rds", sep = "_") %in% 
#      list.files(indicator_outputs))) {
#   
# 
#   
# indicator_values_master <- readRDS(file.path(indicator_outputs,
#                                        paste(location, eco_version,
#                                              "indicator_values_master.rds",
#                                            sep = "_")))
# } else {
  

# Load species data ----


species_data_all <- readRDS(file.path(interim_outputs, 
                                  "global_species_data_3.rds"))

ecoregion_subset <- ecoregion_country_df %>%
  filter(CNTRY_NAME == country) %>%
  unique(.)

#}

if (!is.na(country)) {
  
ecoregion_map <- ecoregion_map[ecoregion_map$ECO_ID %in% 
                                   ecoregion_subset$ECO_ID,]

species_data_all <- species_data[species_data$ecoregion_id %in% 
                                       ecoregion_subset$ECO_ID,]

}
                        
# species_data <- readRDS(file.path(interim_outputs, 
#                                   "version_3_species_data_v1.rds"))

# Rename classes so they match and remove country variable

species_data <- species_data_all %>%
                distinct(.) %>%
                mutate(redlist_assessment_year = as.numeric(as.character(
                       redlist_assessment_year)),
                       decade = ifelse(redlist_assessment_year < 1990, 1980, # Group assessment date into five year bins (prob a better way to do this)
                                       ifelse(redlist_assessment_year > 1989 &
                                              redlist_assessment_year < 2000, 1990,
                                       # ifelse(redlist_assessment_year > 1994 &
                                       #        redlist_assessment_year < 2000, 1995,
                                       ifelse(redlist_assessment_year > 1999 &
                                              redlist_assessment_year < 2005, 2000,
                                       ifelse(redlist_assessment_year > 2004 &
                                              redlist_assessment_year < 2010, 2005,
                                       ifelse(redlist_assessment_year > 2009 &
                                              redlist_assessment_year < 2015, 2010,
                                       ifelse(redlist_assessment_year > 2014 &
                                              redlist_assessment_year < 2021, 2015,
                                       NA))))))) %>%
               filter(redlist_assessment_year != 2020,
                      redlist_assessment_year != 2019) %>%
               drop_na(decade)

# Check the distribution of data over years

hist(species_data$redlist_assessment_year)

ggplot(species_data) +
  geom_bar(aes(x = redlist_assessment_year)) 

spp_per_timepoint <- species_data %>%
                     group_by(redlist_assessment_year) %>%
                     summarise(n_distinct(binomial))

# Remove duplicate errors (species with more than one redlist status in the same year)

duplicates <- species_data %>%
                    group_by(binomial, redlist_assessment_year) %>%
                    summarise(check = n_distinct(redlist_status)) %>%
              filter(check > 1) %>%
              select(binomial) %>%
              distinct(.)

species_data  <- species_data[! species_data$binomial %in% 
                     duplicates$binomial,] 

# Check for species that have complete data for most time points

data_years <- c(1988, 1994, 2000, 2004, 2008, 2012, 2016) # These are the years with most of the data

# Get only the columns we need to convert into wide and check for completeness across years
data_check <- species_data %>% select(binomial, class, redlist_assessment_year, 
                                redlist_status) %>% distinct(.)

# Convert into wide format

data_check_wide <- spread(data_check, key = redlist_assessment_year, 
                    value = redlist_status)

# Subset to the years with most of the data

data_check_wide_subset_years <- data_check_wide[ , c(1, 2, 5, 7, 9, 10, 13, 17, 21)]

# Subset to only species that have a complete set of redlist values for all main years

data_check_complete <- data_check_wide_subset_years[complete.cases(data_check_wide_subset_years),]

# Get the names of species with complete data

species_with_complete_data <- unique(data_check_complete$binomial)

# Subset the original species data to only those with complete redlist across all main years
## Hint - only birds have complete data

species_data_complete_all_classes <- species_data[species_data$binomial %in% 
                                        species_with_complete_data,]

# Subset by test country

if (!is.na(country)) {
  
  species_data <- species_data[species_data$ecoregion_id %in% 
                                 ecoregion_map$ECO_ID, ]
  
}

## Get complete cases, but this time use 'decade' to determine if a species
## has complete data (decade bins redlist year into 5 year breaks)
# Get only complete cases (RL status for each decade)

complete_cases <- species_data %>% 
                  group_by(binomial) %>%
                  summarise(number_timepoints = n_distinct(decade)) %>%
                  filter(number_timepoints == 6) %>%
                  distinct(.)

species_data_complete <- species_data[species_data$binomial %in% 
                                     complete_cases$binomial,]

species_by_ecoregion_complete <- species_data_complete %>%
  distinct(.) %>%
  group_by(ecoregion_id, decade) %>%
  mutate(number_of_species_year = n_distinct(tsn),
         number_extinct = n_distinct(tsn[redlist_status == "EX"|
                                           redlist_status == "EW"]),
         number_atrisk = n_distinct(tsn[redlist_status == "EN"|
                                          redlist_status == "CR"|
                                          redlist_status == "CR(PE)"|
                                          redlist_status == "VU"]),
         number_lowrisk = n_distinct(tsn[redlist_status == "LC"|
                                           redlist_status == "NT"]),
         number_datadeficient = n_distinct(tsn[redlist_status == "DD"]),
         check = number_of_species_year - (number_extinct +
                                             number_atrisk +
                                             number_lowrisk+
                                             number_datadeficient)) %>%
  group_by(ecoregion_id) %>%
  mutate(number_of_species = max(number_of_species_year),
         proportion_extinct = number_extinct/number_of_species,
         proportion_atrisk = number_atrisk/number_of_species,
         proportion_lowrisk = number_lowrisk/number_of_species)


# All species (not just those with complete timepoints)

species_by_ecoregion <- species_data %>%
                        distinct(.) %>%
                        group_by(ecoregion_id, decade) %>%
                        mutate(number_of_species_year = n_distinct(tsn),
                               number_extinct = n_distinct(tsn[redlist_status == "EX"|
                                                               redlist_status == "EW"]),
                               number_atrisk = n_distinct(tsn[redlist_status == "EN"|
                                                      redlist_status == "CR"|
                                                      redlist_status == "CR(PE)"|
                                                      redlist_status == "VU"]),
                               number_lowrisk = n_distinct(tsn[redlist_status == "LC"|
                                                       redlist_status == "NT"]),
                               number_datadeficient = n_distinct(tsn[redlist_status == "DD"]),
                               check = number_of_species_year - (number_extinct +
                                                            number_atrisk +
                                                            number_lowrisk+
                                                            number_datadeficient)) %>%
                         group_by(ecoregion_id) %>%
                         mutate(number_of_species = max(number_of_species_year),
                                proportion_extinct = number_extinct/number_of_species,
                                proportion_atrisk = number_atrisk/number_of_species,
                                proportion_lowrisk = number_lowrisk/number_of_species)

test_ecoregion <- species_by_ecoregion %>% filter(ecoregion_id == 200)

# Proportion of species extinct ----

if ((paste(location, eco_version, "proportion_extinct.rds", sep = "_") %in% 
      list.files(indicator_outputs))) {
  
  #' TODO: Work out why this is producing so many NaNs
  
extinction_values <- readRDS(file.path(indicator_outputs, 
                                              paste(location, eco_version, 
                                                    "proportion_extinct.rds",
                                                    sep = "_"))) 
} else {

extinction_values <- species_by_ecoregion %>%
                     mutate(indicator = "proportion extinct") %>%
                     dplyr::select(indicator, decade, 
                            ecoregion_id, proportion_extinct) %>%
                     rename(year = decade) %>%
                     distinct(.) %>%
                     drop_na(year) %>%
                     group_by(ecoregion_id) %>%
                     arrange(year) %>%
                     mutate(raw_indicator_value = cummax(proportion_extinct))


extinction_values <- as.data.frame(extinction_values)

# Check extinction values behave as anticipated

ECO <- ecuador
test_spp <- species_by_ecoregion %>% filter(ecoregion_id == ECO)
test <- extinction_values %>% filter(ecoregion_id == ECO)

ggplot(test) +
  geom_line(aes(x = year, y = raw_indicator_value)) +
  geom_line(aes(x = year, y = proportion_extinct), col = "red", linetype = 2)

# Remove proportion extinct after check it is working as intended

extinction_values <- extinction_values %>% dplyr::select(-proportion_extinct)

saveRDS(extinction_values, file.path(indicator_outputs, 
                                     paste(location, eco_version, 
                                           "proportion_extinct.rds",
                                           sep = "_")))

}

indicator_df_colnames <- names(extinction_values)

# Proportion of species threatened ----

if ((paste(location, eco_version, "proportion_at_risk.rds", sep = "_") %in% 
      list.files(indicator_outputs))) {
  
  #' TODO: Work out why this is producing so many NaNs
  
at_risk_values <- readRDS(file.path(indicator_outputs, 
                                         paste(location, eco_version, 
                                               "proportion_at_risk.rds",
                                               sep = "_"))) 
} else {


at_risk_values <- species_by_ecoregion_complete %>%
                 mutate(indicator = "proportion at risk") %>%
                 select(indicator, decade, 
                         ecoregion_id, proportion_atrisk) %>%
                 rename(year = decade, 
                        raw_indicator_value = proportion_atrisk) %>%
                 distinct(.) %>%
                 drop_na(year) 


at_risk_values <- as.data.frame(at_risk_values)

saveRDS(at_risk_values, file.path(indicator_outputs, 
                                  paste(location, eco_version, 
                                        "proportion_at_risk.rds",
                                        sep = "_")))

}

# Check extinction values behave as anticipated

ECO <- 7
atrisk_test <- at_risk_values %>% filter(ecoregion_id == ECO)

ggplot(atrisk_test) +
  geom_line(aes(x = year, y = raw_indicator_value)) 


# # Red List Index by class ----
# **WARNING SLOW CODE ** ## ----
# (Takes about a day and a half?)

if ((paste(location, eco_version, "red_list_index.rds", sep = "_") %in% 
      list.files(indicator_outputs))) {
  
  #' TODO: Work out why this is producing so many NaNs
  
rli_values <- readRDS(file.path(indicator_outputs, 
                                         paste(location, eco_version, 
                                               "red_list_index.rds",
                                               sep = "_"))) 
} else {
# Split into different classes

species_data_for_rli <- species_data %>%
                        filter(class == "Amphibia"| class == "Aves"|
                               class == "Mammalia"| class == "Reptilia") %>%
                        filter(tsn != 1444976) # Remove only spp with timepoint in 2015 for amphibia

class_list <- split(species_data_for_rli, species_data_for_rli$class)

#class_list <- lapply(class_list, sample_n, size = 200)


# Split each class into different time points (output should be dataframes of 
# species red list status in a nested list with levels: Class, Time point, Ecoregion)

class_time_list_full <- list()
class_time_ecoregion_list <- list()

for (i in seq_along(class_list)) {
  
  # Split by different assessment timepoints
  
  class_timestep <- split(class_list[[i]], 
                          class_list[[i]]$decade)
  
    for (j in seq_along(class_timestep)) {
      
    # Split by each ecoregion (so have a list of species for each class, timepoint,
    # ecoregion)
    
    class_time_ecoregion_list[[j]] <- split(class_timestep[[j]], 
                                            class_timestep[[j]]$ecoregion_id)
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
  
  taxa <- class[[1]][[1]][[8]][1]
  
  print(paste("Processing class", taxa, sep = " "))
  
  class_all_timepoints <- list()
  
    for (j in seq_along(class)) {
    
      time <- class[[j]] # Get list of ecoregions for one timestep, for one class
      
      timepoint <- time[[1]][[9]][1]
      
      print(paste("Processing class", taxa,
                  "at timepoint", timepoint, sep = " "))
      
        for (k in seq_along(time)) {
          
          eco <- time[[k]][[1]][1]
          
          class_time_ecoregions[[k]] <- calculate_red_list_index(time[[k]]) # Calculate RLI for each ecoregion for one timestep, one class
          
          print(paste("Processed redlist values for class", taxa,
                      "at timepoint", timepoint, "in ecoregion", eco, sep = " "))
          
          class_time_ecoregion_df <- do.call(rbind, class_time_ecoregions) # One dataframe for one time point, one class
          
        } 
    
    class_all_timepoints[[j]] <- class_time_ecoregion_df # Put the time point into a list of timepoints
    
    }
  
  class_all_timepoints_df <- do.call(rbind, class_all_timepoints)
  
  timepoints <- unique(class_all_timepoints_df$redlist_assessment_year)
  taxon <- class_all_timepoints_df$class[1]
  
  classes_rli[[i]] <- class_all_timepoints_df # Put the list of time points into a class
  
  print(paste("Finished processing timepoints", timepoints, "for class", taxon))
  
}

# Convert back into one dataframe

rli_values <- do.call(rbind, classes_rli)

rli_values <- saveRDS(file.path(indicator_outputs,
                                paste(location, eco_version,
                                      "rli_not_formatted.rds", sep = "_")))

# Check if any amphibian values have the wrong year and fix (2008 should be 1980)

rli_values <- rli_values %>%
        mutate(redlist_assessment_year = ifelse(class == "Amphibia" & 
                                                redlist_assessment_year == 2008, 1980,
                                                redlist_assessment_year))
# Format

rli_values <- rli_values %>%
              filter(RLI != 0) %>%
              mutate(indicator = paste("RLI", class, sep = " ")) %>%
              rename(raw_indicator_value = RLI) %>%
              rename(ecoregion_id = Ecoregion_id) %>%
              ungroup()%>%
              filter(ecoregion_id != 0) %>%
              distinct(.) %>%
              filter(redlist_assessment_year < 2016) %>%
              filter(!(class == "Aves" & redlist_assessment_year == 2004)) %>%
              mutate(year = ifelse(redlist_assessment_year < 1990, 1980,
                                   ifelse(redlist_assessment_year > 1989 &
                                   redlist_assessment_year < 2000, 1990,
                                   ifelse(redlist_assessment_year > 1999 &
                                          redlist_assessment_year < 2005, 2000,
                                   ifelse(redlist_assessment_year > 2004 &
                                          redlist_assessment_year < 2010, 2005,
                                   ifelse(redlist_assessment_year > 2009 &
                                          redlist_assessment_year < 2015, 2010,
                                   ifelse(redlist_assessment_year > 2014 &
                                          redlist_assessment_year < 2021, 2015,
                                           NA))))))) %>%
              dplyr::select(indicator, year, 
                            ecoregion_id, raw_indicator_value) 
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

saveRDS(rli_values, file.path(indicator_outputs, 
          paste(location, eco_version, 
                "red_list_index.rds",
                sep = "_"))) 
}

# Remove reptiles for the moment b/c data hasn't been checked and it's behaving weirdly

rli_values <- rli_values %>%
              filter(indicator != "red list index Reptilia")

# Check rli values behave as anticipated

ECO <- cardamom
rli_test <- rli_values %>% filter(ecoregion_id == ECO) %>% 
  filter(indicator == 'RLI Amphibia')

ggplot(rli_test) +
  geom_line(aes(x = year, y = raw_indicator_value))


# # Red List Index all classes ----


# Split the species data into different time points (output should be dataframes of 
# species red list status in a nested list with levels: Time point, Ecoregion)

# only 8400 spp left by now
species_data_for_rli_all_classes <- species_data_complete_all_classes %>%
                                    filter(class != "reptile") %>%
                                    filter(tsn != 1444976)
# Split the data by timepoint

rli_timestep <- split(species_data_complete_all_classes, 
                      species_data_complete_all_classes$redlist_assessment_year) 


time_ecoregion_list <- list()

for (i in seq_along(rli_timestep)) {
    
    # Split by each ecoregion as well (so have a list of species for each timepoint,
    # ecoregion)
    
    time_ecoregions <- split(rli_timestep[[i]], 
                                      rli_timestep[[i]]$ecoregion_id)
    
    # Remove the empty lists with no data in them
    
    time_ecoregions <- list.clean(time_ecoregions, 
                                  fun = is.null, recursive = FALSE)
    
    time_ecoregion_list[[i]] <- time_ecoregions
  
}
  

# Calculate the RLI per ecoregion, per timepoint (output should be
# dataframes of RLI values by ecoregion, by timepoint)

rli_ecoregions_single_year <- list()
rli_ecoregions_all_years <- list()

for (i in seq_along(time_ecoregion_list)) {
  
  timestep <- time_ecoregion_list[[i]] # Get list of species for each ecoregion in one timestep
  
  year <- timestep[[1]][[5]][1]
  
  print(paste("Processing year", year, sep = " "))
  
  #class_all_timepoints <- list()
  
  for (j in seq_along(timestep)) {
    
      eco <- timestep[[j]][[1]][1] # Pull out species list for one ecoregion
      
      rli_ecoregions_single_year[[j]] <- calculate_red_list_index(timestep[[j]]) # Calculate RLI for each ecoregion for one timestep, one class
      
      print(paste("Processed redlist values in year", year,
                  "for ecoregion", eco, sep = " "))
      
      } 
    
    rli_ecoregions_single_year_df <- do.call(rbind, rli_ecoregions_single_year) # One dataframe for one time point
    
    rli_ecoregions_all_years[[i]] <- rli_ecoregions_single_year_df # Put the time point into a list of timepoints
    
  }

# Convert back into one dataframe

rli_ecoregions_all_years_df <- do.call(rbind, rli_ecoregions_all_years)

saveRDS(rli_ecoregions_all_years_df, 
        file.path(indicator_outputs, paste(location, eco_version,
                  "rli_all_classes_not_formatted.rds", sep = "_")))

# Check if any amphibian values have the wrong year and fix (2008 should be 1980)

rli_values <- rli_values %>%
  mutate(redlist_assessment_year = ifelse(class == "Amphibia" & 
                                            redlist_assessment_year == 2008, 1980,
                                          redlist_assessment_year))
# Format

rli_values <- rli_values %>%
  filter(RLI != 0) %>%
  mutate(indicator = paste("RLI", class, sep = " ")) %>%
  rename(raw_indicator_value = RLI) %>%
  rename(ecoregion_id = Ecoregion_id) %>%
  ungroup()%>%
  filter(ecoregion_id != 0) %>%
  distinct(.) %>%
  filter(redlist_assessment_year < 2016) %>%
  filter(!(class == "Aves" & redlist_assessment_year == 2004)) %>%
  mutate(year = ifelse(redlist_assessment_year < 1990, 1980,
                       ifelse(redlist_assessment_year > 1989 &
                                redlist_assessment_year < 2000, 1990,
                              ifelse(redlist_assessment_year > 1999 &
                                       redlist_assessment_year < 2005, 2000,
                                     ifelse(redlist_assessment_year > 2004 &
                                              redlist_assessment_year < 2010, 2005,
                                            ifelse(redlist_assessment_year > 2009 &
                                                     redlist_assessment_year < 2015, 2010,
                                                   ifelse(redlist_assessment_year > 2014 &
                                                            redlist_assessment_year < 2021, 2015,
                                                          NA))))))) %>%
  dplyr::select(indicator, year, 
                ecoregion_id, raw_indicator_value) 
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

saveRDS(rli_values, file.path(indicator_outputs, 
                              paste(location, eco_version, 
                                    "red_list_index.rds",
                                    sep = "_"))) 
}

# Remove reptiles for the moment b/c data hasn't been checked and it's behaving weirdly

rli_values <- rli_values %>%
  filter(indicator != "red list index Reptilia")

# Check rli values behave as anticipated

ECO <- cardamom
rli_test <- rli_values %>% filter(ecoregion_id == ECO) %>% 
  filter(indicator == 'RLI Amphibia')

ggplot(rli_test) +
  geom_line(aes(x = year, y = raw_indicator_value))


# Human Footprint Index  ----

# Read in each raster (each raster is a time point)

hfp_maps_file_names <-   list.files(file.path(inputs,
                         "/williams_human_footprint_index/williams_human_footprint_index_WGS1984/"),
                         recursive = FALSE)

hfp_maps_file_names <- hfp_maps_file_names[grep(".tif", hfp_maps_file_names)]

hfp_maps_files <- paste(inputs,
                        "/williams_human_footprint_index/williams_human_footprint_index_WGS1984/",
                        hfp_maps_file_names,
                        sep = "")

hfp_maps <- lapply(hfp_maps_files, raster)

names(hfp_maps) <- hfp_maps_file_names

# Extract values from each raster

if ((paste(location, eco_version,
            "human_footprint_index.rds", sep = "_") %in% 
      list.files(indicator_outputs))) {

hfp_values <- readRDS(file.path(indicator_outputs,
                                paste(location, eco_version,
                                      "human_footprint_index.rds", sep = "_")))
} else {

# * HFP 2000 ----

  if (paste(location, "hfp_2000_ecoregion_values.rds", sep = "_") %in% 
      list.files(indicator_outputs)) {
    
    hfp_2000_ecoregion_values <- readRDS(file.path(indicator_outputs, 
                                                   paste(location, 
                                                         "hfp_2000_ecoregion_values.rds",
                                                         sep = "_"))) 
  } else {
    
    hfp_2000_map <- hfp_maps[[grep("2000", names(hfp_maps))]]
    
    # SLOW CODE - 123536.0s or around 34 hours
    
    system.time(hfp_2000_ecoregion_values <- ecoregion_map %>%
                mutate(hfpmean = raster::extract(hfp_2000_map, ecoregion_map, fun = mean, na.rm = TRUE),
                       hfpsd = raster::extract(hfp_2000_map, ecoregion_map, fun = sd, na.rm = TRUE),
                       hfpmax = raster::extract(hfp_2000_map, ecoregion_map, fun = max, na.rm = TRUE),
                       hfpmin = raster::extract(hfp_2000_map, ecoregion_map, fun = min, na.rm = TRUE)))
    
    saveRDS(hfp_2000_ecoregion_values, file.path(indicator_outputs, paste(location, 
                                                                          "hfp_2000_ecoregion_values.rds",
                                                                          sep = "_")))
  }
  
  hfp_2000 <- hfp_2000_ecoregion_values %>% 
              select(ECO_ID, hfpmean) %>%
              mutate(year = "2000",
                     indicator = "mean human footprint index") %>%
              rename(ecoregion_id = ECO_ID,
                     raw_indicator_values = hfpmean) %>%
              st_set_geometry(NULL) 
  
# * HFP 2005 ----
  
  #' TODO: Get min, max and sd as well
  
  if (paste(location, "hfp_2005_ecoregion_values.rds", sep = "_") %in% 
      list.files(indicator_outputs)) {
    
    hfp_2005_ecoregion_values <- readRDS(file.path(indicator_outputs, 
                                                   paste(location, 
                                                         "hfp_2005_ecoregion_values.rds",
                                                         sep = "_"))) 
  } else {
    
    hfp_2005_map <- hfp_maps[[grep("2005", names(hfp_maps))]]
    
    # SLOW CODE - 
    
    # system.time(hfp_2005_ecoregion_values <- ecoregion_map %>%
    #               mutate(hfpmean = raster::extract(hfp_2005_map, ecoregion_map, fun = mean, na.rm = TRUE),
    #                      hfpsd = raster::extract(hfp_2005_map, ecoregion_map, fun = sd, na.rm = TRUE),
    #                      hfpmax = raster::extract(hfp_2005_map, ecoregion_map, fun = max, na.rm = TRUE),
    #                      hfpmin = raster::extract(hfp_2005_map, ecoregion_map, fun = min, na.rm = TRUE)))
    
    system.time(hfp_2005_ecoregion_values <- ecoregion_map %>%
                mutate(hfpmean = raster::extract(hfp_2005_map, 
                                                   ecoregion_map, fun = mean, 
                                                   na.rm = TRUE)))
    
    
    saveRDS(hfp_2005_ecoregion_values, file.path(indicator_outputs, 
                                                 paste(location, 
                                                 "hfp_2005_ecoregion_values.rds",
                                                 sep = "_")))
  }
  
  hfp_2005 <- hfp_2005_ecoregion_values %>% 
    select(ECO_ID, hfpmean) %>%
    mutate(year = "2005",
           indicator = "mean human footprint index") %>%
    rename(ecoregion_id = ECO_ID,
           raw_indicator_values = hfpmean) %>%
    st_set_geometry(NULL) 

  # * HFP 2010 ----
  
  if (paste(location, "hfp_2010_ecoregion_values.rds", sep = "_") %in% 
      list.files(indicator_outputs)) {
    
    hfp_2010_ecoregion_values <- readRDS(file.path(indicator_outputs, 
                                                   paste(location, 
                                                         "hfp_2010_ecoregion_values.rds",
                                                         sep = "_"))) 
  } else {
    
    hfp_2010_map <- hfp_maps[[grep("2010", names(hfp_maps))]]
    
    # SLOW CODE - 
    
    # system.time(hfp_2005_ecoregion_values <- ecoregion_map %>%
    #               mutate(hfpmean = raster::extract(hfp_2005_map, ecoregion_map, fun = mean, na.rm = TRUE),
    #                      hfpsd = raster::extract(hfp_2005_map, ecoregion_map, fun = sd, na.rm = TRUE),
    #                      hfpmax = raster::extract(hfp_2005_map, ecoregion_map, fun = max, na.rm = TRUE),
    #                      hfpmin = raster::extract(hfp_2005_map, ecoregion_map, fun = min, na.rm = TRUE)))
    
    system.time(hfp_2010_ecoregion_values <- ecoregion_map %>%
                  mutate(hfpmean = raster::extract(hfp_2010_map, 
                                                   ecoregion_map, fun = mean, 
                                                   na.rm = TRUE)))
    
    
    saveRDS(hfp_2010_ecoregion_values, file.path(indicator_outputs, 
                                                 paste(location, 
                                                 "hfp_2010_ecoregion_values.rds",
                                                 sep = "_")))
  }
  
  hfp_2010 <- hfp_2010_ecoregion_values %>% 
              select(ECO_ID, hfpmean) %>%
              mutate(year = "2010",
                     indicator = "mean human footprint index") %>%
              rename(ecoregion_id = ECO_ID,
                     raw_indicator_values = hfpmean) %>%
              st_set_geometry(NULL) 
  
  # * HFP 2013 ----
  
  if (paste(location, "hfp_2013_ecoregion_values.rds", sep = "_") %in% 
      list.files(indicator_outputs)) {
    
    hfp_2013_ecoregion_values <- readRDS(file.path(indicator_outputs, 
                                                   paste(location, 
                                                         "hfp_2013_ecoregion_values.rds",
                                                         sep = "_"))) 
  } else {
    
    hfp_2013_map <- hfp_maps[[grep("2013", names(hfp_maps))]]
    
    # SLOW CODE - 
    
    # system.time(hfp_2005_ecoregion_values <- ecoregion_map %>%
    #               mutate(hfpmean = raster::extract(hfp_2005_map, ecoregion_map, fun = mean, na.rm = TRUE),
    #                      hfpsd = raster::extract(hfp_2005_map, ecoregion_map, fun = sd, na.rm = TRUE),
    #                      hfpmax = raster::extract(hfp_2005_map, ecoregion_map, fun = max, na.rm = TRUE),
    #                      hfpmin = raster::extract(hfp_2005_map, ecoregion_map, fun = min, na.rm = TRUE)))
    
    system.time(hfp_2013_ecoregion_values <- ecoregion_map %>%
                  mutate(hfpmean = raster::extract(hfp_2013_map, 
                                                   ecoregion_map, fun = mean, 
                                                   na.rm = TRUE)))
    
    
    saveRDS(hfp_2013_ecoregion_values, file.path(indicator_outputs, 
                                                 paste(location, 
                                                 "hfp_2013_ecoregion_values.rds",
                                                 sep = "_")))
  }
  
  
# * HFP 1993 ----

# **WARNING SLOW CODE ** ## ---- 
  # (time elapsed 30088.33 ~ 8 hours)

if (paste(location, "hfp_1993_ecoregion_values.rds", sep = "_") %in% 
      list.files(indicator_outputs)) {
  
  hfp_1993_ecoregion_values <- readRDS(file.path(indicator_outputs, 
                                              paste(location, 
                                                    "hfp_1993_ecoregion_values.rds",
                                                    sep = "_"))) 
} else {
  

  hfp_1993_data <- raster(file.path(inputs,
                                    "human_footprint_index\\HFP1993RP.tif"))
  
  system.time(hfp_1993_ecoregion_values <- ecoregion_map %>%
                mutate(raw_indicator_value =
                         raster::extract(hfp_1993_data,
                                         ecoregion_map,
                                         fun = mean,
                                         na.rm = TRUE)))
  
  # system.time(hfp_1993_ecoregion_values <- 
  #   ecoregion_map %>% mutate(
  #     hfpMean = raster_extract(hfp_1993_data, ecoregion_map, fun = mean, na.rm = TRUE)
  #     #,
  #     # hfpMax = raster_extract(hfp_1993_data, ecoregion_map, fun = max, na.rm = TRUE),
  #     # hfpMin = raster_extract(hfp_1993_data, ecoregion_map, fun = min, na.rm = TRUE)
  #   ))
  
 saveRDS(hfp_1993_ecoregion_values, file.path(indicator_outputs, paste(location, 
                                            "hfp_1993_ecoregion_values.rds",
                                             sep = "_")))
}

# Create indicator values dataframe

hfp_1993_values <- hfp_1993_ecoregion_values %>%
                   st_set_geometry(NULL) %>%
                   mutate(indicator = "mean human footprint index") %>%
                   mutate(year = 1990) %>%
                   rename(ecoregion_id = ECO_ID) %>%
                   select(all_of(indicator_df_colnames)) %>%
                   mutate(raw_indicator_value = ifelse(is.nan(raw_indicator_value),
                                                       NA, raw_indicator_value))

saveRDS(hfp_1993_values, file.path(indicator_outputs, 
                                   paste(location, 
                                         "hfp_1993_ecoregion_values.rds",
                                         sep = "_"))) 
# * HFP 2009 ----

# **WARNING SLOW CODE ** ## ---- 
# (time elapsed 30088.33 ~ 8 hours)

#' if (!(paste(location, "hfp_2009_ecoregion_map.rds", sep = "_") %in% 
#'       list.files(indicator_outputs))) {
#'   
#'   #' TODO: Work out why this is producing so many NaNs
#'   
#'   hfp_2009_ecoregion_map <- readRDS(file.path(indicator_outputs, 
#'                                               paste(location, 
#'                                                     "hfp_2009_ecoregion_map.rds",
#'                                                     sep = "_"))) 
#' } else {
#'   
#'   hfp_2009_data <- raster(file.path(inputs,
#'                                     "human_footprint_index\\HFP2009RP.tif"))
#'   
#'   system.time(hfp_2009_ecoregion_values <- ecoregion_map %>%
#'                 mutate(raw_indicator_value =
#'                        raster::extract(hfp_2009_data,
#'                                          ecoregion_map,
#'                                          fun = mean,
#'                                          na.rm = TRUE),
#'                        indicator_sd = raster::extract(hfp_2009_data,
#'                                                       ecoregion_map,
#'                                                       fun = sd,
#'                                                       na.rm = TRUE)))
#'   
#'   # system.time(hfp_1993_ecoregion_values <- 
#'   #   ecoregion_map %>% mutate(
#'   #     hfpMean = raster_extract(hfp_1993_data, ecoregion_map, fun = mean, na.rm = TRUE)
#'   #     #,
#'   #     # hfpMax = raster_extract(hfp_1993_data, ecoregion_map, fun = max, na.rm = TRUE),
#'   #     # hfpMin = raster_extract(hfp_1993_data, ecoregion_map, fun = min, na.rm = TRUE)
#'   #   ))
#'   
#'   saveRDS(hfp_2009_ecoregion_values, file.path(indicator_outputs, paste(location, 
#'                                          "hfp_2009_ecoregion_values.rds",
#'                                          sep = "_")))
#'   
#' 
#' }
#' 
#' # Create indicator values dataframe
#' 
#' hfp_2009_values <- hfp_2009_ecoregion_map %>%
#'                    st_set_geometry(NULL) %>%
#'                    mutate(indicator = "mean human footprint index") %>%
#'                    mutate(year = "2009") %>%
#'                    rename(ecoregion_id = ECO_ID) %>%
#'                    select(names(rli_values)) %>%
#'                    drop_na()


# * HFP 2013 ----

# Read in the Human Footprint Index data 

hfp_by_ecoregion_2013 <- read.csv(paste(inputs, "human_footprint_index", 
                                        "human_footprint_index_by_ecoregion.csv",
                                        sep = "/"))

hfp_2013_values <- hfp_by_ecoregion_2013 %>%
                         rename(ecoregion_name = ECO_NAME) %>%
                         # merge(ecoregion_map[c("ecoregion_name", "ECO_ID")], 
                         #       by = "ECO_ID") %>% # This will subset automatically if you subset by country
                         mutate(year = 2010) %>%
                         mutate(indicator = "mean human footprint index") %>%
                         rename(raw_indicator_value = HFP,
                                ecoregion_id = ECO_ID) %>%
                         dplyr::select(indicator, year, ecoregion_id, 
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


hfp_values <- rbind(hfp_1993_values, hfp_2013_values)

saveRDS(hfp_values, file.path(indicator_outputs, paste(location, eco_version,
                          "human_footprint_index.rds", sep = "_")))

}

# Biodiversity Habitat Index Plants ----

if ((paste(location, eco_version, "bhi_plants.rds", sep = "_") %in% 
     list.files(indicator_outputs))) {
  

bhi_plants_values <- readRDS(file.path(indicator_outputs, 
                                         paste(location, eco_version, 
                                               "bhi_plants.rds",
                                               sep = "_"))) 
} else {

bhi_plants_all <- read.csv(file.path(inputs,
                                       "biodiversity_habitat_index\\BHI_PLANTS_2017_ECOREGIONS.csv"))

bhi_plants_all <- bhi_plants_all %>%
                  mutate(ecoregion_id = str_sub(REGION, start = -3)) %>%
                  rename("2005" = X2005,
                         "2010" = X2010,
                         "2015" = X2015) %>%
                  select(- REGION) %>%
                  filter(ecoregion_id != "ALL") %>%
                  filter(ecoregion_id != "000") %>%
                  mutate(ecoregion_id = str_remove(ecoregion_id, "^0+"))

bhi_plants_all_long <- melt(bhi_plants_all, id.vars = "ecoregion_id")

bhi_plants_values <- bhi_plants_all_long %>%
                     rename(year = variable, 
                            raw_indicator_value = value) %>%
                     mutate(indicator = "BHI plants",
                            year = as.numeric(as.character(year)),
                            ecoregion_id = as.numeric(ecoregion_id)) %>%
                     dplyr::select(indicator_columns) %>%
                     distinct(.)

saveRDS(bhi_plants_values, file.path(indicator_outputs, 
                                     paste(location, eco_version, 
                                           "bhi_plants.rds",
                                           sep = "_")))

}

# Living Planet Index ----

if (paste(location, eco_version, "lpi.rds", sep = "_") %in% 
    list.files(indicator_outputs)) {
  
lpi_values <- readRDS(file.path(indicator_outputs, 
                                paste(location, eco_version, 
                                      "lpi.rds",
                                      sep = "_")))
} else {
  
lpi_ecoregion_directory <- file.path(indicator_outputs, 
                                   "lpi_ecoregion_directory")

if( !dir.exists( lpi_ecoregion_directory ) ) {
  
  dir.create( lpi_ecoregion_directory, recursive = TRUE )
  
}

setwd(lpi_ecoregion_directory)

if (paste(location, eco_version, "LPI_2020_input_data_with_ecoregions.rds", sep = "_") %in% 
    list.files(indicator_outputs)) {
  
lpi_inputs <- readRDS(file.path(indicator_outputs, 
                                           paste(location, eco_version,
                                                 "LPI_2020_input_data_with_ecoregions.rds", 
                                                 sep = "_")))
} else {

lpi_data <- read.csv(file.path(inputs,
                     "living_planet_index_2020_version\\LPR2020data_public.csv"))

## Assign species to ecoregions

# Get coordinates

lpi_coordinates <- lpi_data %>%
                   select(Binomial, Latitude, Longitude) %>%
                   distinct(.) %>%
                   filter(complete.cases(Latitude, Longitude))

# Convert into simple features

lpi_sf <- st_as_sf(lpi_coordinates, coords = c('Longitude',
                                                'Latitude'),
                       crs = st_crs(ecoregion_map))

# Match coordinates of population to ecoregions

lpi_species_ecoregions <- st_intersection(lpi_sf, ecoregion_map)
lpi_inputs_ecoregions <- st_drop_geometry(lpi_species_ecoregions)

# Format the input data

lpi_inputs <- lpi_data %>%
              merge(lpi_inputs_ecoregions[c("Binomial", "ECO_ID")],
                    by = "Binomial", all = FALSE) %>% # Drops around 1000 marine spp
              select(ECO_ID, everything())

# Have a look

lpi_input_summary <- lpi_inputs %>%
                     group_by(ECO_ID) %>%
                     summarise(number_of_records = n_distinct(Binomial),
                               .groups = "drop_last")

# Save input data

saveRDS(lpi_inputs, file.path(indicator_outputs, 
                              paste(location, eco_version, 
                              "LPI_2020_input_data_with_ecoregions.rds",
                               sep = "_")))

# Split LPI data by ecoregion

lpi_inputs_ecoregions <- split(lpi_inputs, lpi_inputs$ECO_ID)

eco_names <- vector()

for (i in seq_along(lpi_inputs_ecoregions)) {

  eco_names[i] <- as.character(lpi_inputs_ecoregions[[i]]$ECO_ID[1])
  
}

names(lpi_inputs_ecoregions) <- eco_names

lpi_infiles_ecoregions <- list()
lpi_values_by_ecoregion <- list()
file_name <- vector()

for (i in seq_along(lpi_inputs_ecoregions)) {
  
single_ecoregion <- lpi_inputs_ecoregions[[i]] %>%
                    select(-ECO_ID) %>%
                    distinct(.)

index_vector_ecoregion <- c(1:length(unique(single_ecoregion$ID)))

file_name[i] <- paste("ecoregion", eco_names[i], sep = "_")

lpi_infiles_ecoregions[[i]] <- create_infile(single_ecoregion, 
                                             index_vector = index_vector_ecoregion, 
                                             name = file_name[i])

ecoregion_lpi <- LPIMain(paste(file_name[i], "infile.txt", sep = "_"),
                                use_weightings = 1, VERBOSE = FALSE)

ecoregion_lpi <- ecoregion_lpi[complete.cases(ecoregion_lpi), ]

lpi_values_by_ecoregion[[i]] <- ecoregion_lpi

rm(ecoregion_lpi)
closeAllConnections()
gc()

}

out <- list.files(lpi_ecoregion_directory)[grepl("Results",list.files(lpi_ecoregion_directory))]

get_lpi <- function(data, name) {
  
  y <- str_remove(name, "ecoregion_") 
  y <- str_remove(y, "_infile_Results.txt")
  x <- data %>%
       rownames_to_column(var = "year") %>%
       mutate(ecoregion_id = as.numeric(y),
              indicator = "LPI") %>%
       rename(raw_indicator_value = LPI_final) %>%
       select(indicator, year, ecoregion_id, raw_indicator_value) %>%
       mutate(year = as.numeric(year)) %>%
       filter(raw_indicator_value != - 99.0)
  
  print(paste("retrieved ecoregion", y, sep = " "))
  
  return(x)
  
}

lpi_values_formatted <- list()

for (i in seq_along(lpi_values_by_ecoregion)) {
  
  lpi_values_formatted[[i]] <- get_lpi(lpi_values_by_ecoregion[[i]],
                                       out[[i]])
  
}

lpi_values <- do.call(rbind, lpi_values_formatted)

lpi_test <- lpi_values %>% filter(ecoregion_id == mascarene)

ggplot(lpi_test) +
  geom_line(aes(x = year, y = raw_indicator_value)) 

saveRDS(lpi_values, file.path(indicator_outputs, 
                                     paste(location, eco_version, 
                                           "lpi.rds",
                                           sep = "_")))
}

# Richness Biodiversity Intactness Index 2005 ----

if (paste(location, eco_version, "richness_BII.rds", sep = "_") %in% 
    list.files(indicator_outputs)) {
  
bii_richness_values <- readRDS(file.path(indicator_outputs, 
                                                 paste(location, eco_version,
                                                       "richness_BII.rds", 
                                                       sep = "_")))
} else {

## WARNING - SLOW CODE - bii richness time elapsed 50753.13 (~ 14 hrs)

if ((paste(location, "richness_bii_2005_ecoregion_map.rds", sep = "_") %in% 
      list.files(indicator_outputs))) {

#' TODO: Work out why this is producing so many NaNs

bii_rich_ecoregion_map <- readRDS(file.path(indicator_outputs, 
                                            paste(location, 
                                            "richness_bii_2005_ecoregion_map.rds",
                                             sep = "_"))) 
} else {
  
  if (!is.na(country)) {
    
    bii_2005_rich_data <- raster(file.path(inputs,
                          "biodiversity_intactness_index\\bii_rich_2005_aus.tif"))
    
  } else {
    
    bii_2005_rich_data <- raster(file.path(inputs,
                          "biodiversity_intactness_index\\final-rich-bii-isl-main.tif"))
    
  }
  
system.time(bii_rich_ecoregion_map <- ecoregion_map %>%
                          mutate(raw_indicator_value = 
                                 raster::extract(bii_2005_rich_data,
                                                  ecoregion_map,
                                                  fun = mean, 
                                                 na.rm = TRUE)))

saveRDS(bii_rich_ecoregion_map, file.path(indicator_outputs, paste(location, 
                                          "richness_bii_2005_ecoregion_map.rds",
                                          sep = "_")))
}

 
#' #' TODO: Figure out how to deal with multiple polygons of the same ecoregion - once
#' #' you remove the geometry you end up with multiple values p/ecoregion

bii_richness_values <- bii_rich_ecoregion_map %>%
                      st_set_geometry(NULL) %>%
                      mutate(indicator = "richness biodiversity intactness index") %>%
                      mutate(year = 2005) %>%
                      rename(ecoregion_id = ECO_ID) %>%
                      select(all_of(indicator_columns)) %>%
                      mutate(raw_indicator_value = ifelse(
                        is.nan(raw_indicator_value),
                                      NA, raw_indicator_value))

saveRDS(bii_richness_values, file.path(indicator_outputs, 
                                       paste(location, eco_version,
                                             "richness_BII.rds", sep = "_")))

}

# Abundance Biodiversity Intactness Index 2005 ----

if (paste(location, eco_version, "abundance_BII.rds", sep = "_") %in% 
    list.files(indicator_outputs)) {
  
bii_abundance_values <- readRDS(file.path(indicator_outputs, 
                                           paste(location, eco_version,
                                                 "abundance_BII.rds", 
                                                 sep = "_")))
} else {

if (paste(location, "abundance_bii_2005_ecoregion_map.rds", sep = "_") %in% 
      list.files(indicator_outputs)) {
  
  bii_abundance_ecoregion_map <- readRDS(file.path(indicator_outputs, 
                                              paste(location, 
                                                    "abundance_bii_2005_ecoregion_map.rds",
                                                    sep = "_"))) 
} else {
  
  if (!is.na(country)) {
    
  bii_2005_abundance_data <- raster(file.path(inputs,
                                           "biodiversity_intactness_index\\bii_abund_2005_aus.tif"))
    
  } else {
    
  bii_2005_abundance_data <- raster(file.path(inputs,
                                    "biodiversity_intactness_index\\final-abund-bii-isl-main.tif"))
    
  }
  
  bii_abundance_ecoregion_map <- ecoregion_map %>%
    mutate(raw_indicator_value = 
           raster::extract(bii_2005_abundance_data,
                           ecoregion_map,
                           fun = mean, 
                           na.rm = TRUE))
  
  saveRDS(bii_abundance_ecoregion_map, file.path(indicator_outputs, paste(location, 
                                       "abundance_bii_2005_ecoregion_map.rds",
                                        sep = "_")))
}


#' #' TODO: Figure out how to deal with multiple polygons of the same ecoregion - once
#' #' you remove the geometry you end up with multiple values p/ecoregion

bii_abundance_values <- bii_abundance_ecoregion_map %>%
                        st_set_geometry(NULL) %>%
                        mutate(indicator = "abundance biodiversity intactness index") %>%
                        mutate(year = 2005) %>%
                        rename(ecoregion_id = ECO_ID,) %>%
                        select(all_of(indicator_columns)) %>%
                        mutate(raw_indicator_value = ifelse(
                          is.nan(raw_indicator_value),
                          NA, raw_indicator_value))

saveRDS(bii_abundance_values, file.path(indicator_outputs, 
                                       paste(location, eco_version,
                                             "abundance_BII.rds", sep = "_")))
}

# Wilderness Intactness Index ----

# wii_map_data <- st_read("N:\\Quantitative-Ecology\\Simone\\extinction_test\\inputs\\intactness\\intactness\\Ecoregions2017_intactness.shp")
# 
# wii_data <- as.data.frame(wii_map_data) %>%
#   dplyr::select(OBJECTID, ECO_NAME, ECO_ID, Q2009_LL, Q2009, Q2009_UL, PLOTCAT)

# Analyse indicators ----

ecoregions <- as.data.frame(ecoregion_map) %>% dplyr::select(-geometry)

# Combine indicator values into a single dataframe ----

indicator_values <- rbind(extinction_values, 
                          at_risk_values, 
                          rli_values,
                          hfp_values,
                          bhi_plants_values,
                          bii_richness_values,
                          bii_abundance_values)

indicator_values_master <- indicator_values %>%
                    merge(ecoregion_country_df[c("ECO_ID", 
                                                 "CNTRY_NAME")], 
                          by.x = "ecoregion_id", by.y = "ECO_ID",
                          all.x = TRUE, all.y = FALSE) %>%
                    distinct(.) %>%
                    mutate(indicator_year = paste(indicator, year, sep = " ")) %>%
                    rename(country = CNTRY_NAME) %>%
                    merge(ecoregions[c("ECO_ID", "REALM")], 
                          by.x = "ecoregion_id",
                          by.y = "ECO_ID") %>%
                    rename(realm = REALM) %>%
                    filter(ecoregion_id != 0)

saveRDS(indicator_values_master, file.path(indicator_outputs,
                                           paste(location, eco_version,
                                                 "indicator_values_master.rds",
                                                 sep = "_")))

}

# Data transformations ----

indicator_values_2 <- indicator_values_master %>%
                      dplyr::select(ecoregion_id, indicator_year,
                                    raw_indicator_value) %>%
                      distinct(.)

# * Convert to wide ----

indicators_wide <- indicator_values_2 %>%
                                spread(key = indicator_year, 
                                       value = raw_indicator_value) 


names(indicators_wide) <- make.names(names(indicators_wide), 
                                           unique = TRUE)

summary(indicators_wide)

# * Invert ----
# For negatively valanced variables (where high values = negative outcome)

cols_to_invert <- c("human.footprint.index.1990", "human.footprint.index.2010",
                    "proportion.at.risk.1980", "proportion.at.risk.1990",
                    "proportion.at.risk.2000", "proportion.at.risk.2005",
                    "proportion.at.risk.2010","proportion.at.risk.2015", 
                    "proportion.extinct.1980", "proportion.extinct.1990",
                    "proportion.extinct.2000", "proportion.extinct.2005",
                    "proportion.extinct.2010", "proportion.extinct.2015")

col_index <- names(indicators_wide) %in% cols_to_invert
cols_min <- as.numeric(sapply(indicators_wide, min, na.rm = TRUE))
cols_max <- as.numeric(sapply(indicators_wide, max, na.rm = TRUE))


keys <- as.data.frame(col_index) %>%
        mutate(keys = ifelse(col_index == FALSE, 1, -1)) %>%
        select(keys) %>%
        pull(.)

indicators_wide$ecoregion_id <- as.numeric(indicators_wide$ecoregion_id)

indicators_wide_inv <- as.data.frame(reverse.code(keys,indicators_wide,
                          mini = cols_min, maxi = cols_max))

summary(indicators_wide_inv)

# test

ECO <- mascarene
HFP <- indicators_wide %>% filter(ecoregion_id == ECO) %>% select(human.footprint.index.1990, 
                                                                  human.footprint.index.2010)
HFP

HFP_inv <- indicators_wide_inv %>% filter(ecoregion_id == ECO) %>% 
            select(`human.footprint.index.1990-`, `human.footprint.index.2010-`)
HFP_inv

# * Centre ----

indicators_wide_inv_centred <- indicators_wide_inv %>%
  mutate_at(c(2:ncol(indicators_wide_inv)), funs(c(scale(.)))) 

summary(indicators_wide_inv_centred)

# ** Centred boxplots ----

indicators_scaled <- melt(indicators_wide_inv_centred, 
                          id.vars = 'ecoregion_id')

boxplots <- ggplot(indicators_scaled) +
            geom_boxplot(aes(x = variable, y = value)) +
            theme(axis.text.x=element_text(angle= 45,hjust=1))

boxplots

# * Transform ----
# https://www.datanovia.com/en/lessons/transform-data-to-normal-distribution-in-r/

# Check distributions

histograms <- lapply(indicators_wide_inv_centred, hist)
names(histograms) <- colnames(indicators_wide_inv_centred)
hist(indicators_wide_inv_centred$red.list.index.Amphibia.2000)

# Measure skew (nearly all negatively skewed)

indicators_skew <- sapply(indicators_wide_inv_centred, skewness, na.rm = TRUE)
indicators_skew_index <- as.data.frame(indicators_skew) %>%
                         mutate(logtransform = ifelse(indicators_skew < -0.05, 
                                                      TRUE, FALSE)) %>%
                         select(logtransform) %>%
                         pull(.)

hist(indicators_wide_inv_centred$`human.footprint.index.2010-`)
new <- log10(max(indicators_wide_inv_centred$`human.footprint.index.2010-`+1, na.rm = TRUE) - indicators_wide_inv_centred$`human.footprint.index.2010-`) 
hist(new)
skewness(new, na.rm = TRUE)

indicators_transformed <- indicators_wide_inv_centred 

# HFP 2010
indicators_transformed$`human.footprint.index.2010-` <- log10(max(indicators_transformed$`human.footprint.index.2010-`+1, 
                                                                  na.rm = TRUE) - 
                                                                indicators_transformed$`human.footprint.index.2010-`) 
skewness(indicators_wide_inv_centred$`human.footprint.index.2010-`, na.rm = TRUE)
hist(indicators_wide_inv_centred$`human.footprint.index.2010-`)
skewness(indicators_transformed$`human.footprint.index.2010-`, na.rm = TRUE)
hist(indicators_transformed$`human.footprint.index.2010-`)

#HFP 1990

skewness(indicators_wide_inv_centred$`human.footprint.index.1990-`, na.rm = TRUE)
hist(indicators_wide_inv_centred$`human.footprint.index.1990-`)

indicators_transformed$`human.footprint.index.1990-` <- log10(max(indicators_transformed$`human.footprint.index.1990-`+1, 
                                                                  na.rm = TRUE) - 
                                                                indicators_transformed$`human.footprint.index.1990-`) 
skewness(indicators_transformed$`human.footprint.index.1990-`, na.rm = TRUE)
hist(indicators_transformed$`human.footprint.index.1990-`)

#At risk 1980

skewness(indicators_wide_inv_centred$`proportion.at.risk.1980-`, na.rm = TRUE)
hist(indicators_wide_inv_centred$`proportion.at.risk.1980-`)

indicators_transformed$`proportion.at.risk.1980-` <- 1/(max(indicators_wide_inv_centred$`proportion.at.risk.1980-`+1, 
                                                                  na.rm = TRUE) - 
                                                             indicators_wide_inv_centred$`proportion.at.risk.1980-`) 
skewness(indicators_transformed$`proportion.at.risk.1980-`, na.rm = TRUE)
hist(indicators_transformed$`proportion.at.risk.1980-`)

#At risk 1990

to_transform <- indicators_wide_inv_centred$`proportion.at.risk.1990-`
transformed <- indicators_transformed$`proportion.at.risk.1990-`

skewness(to_transform, na.rm = TRUE)
hist(to_transform)

indicators_transformed$`proportion.at.risk.1990-` <- 1/(max(to_transform + 1, na.rm = TRUE) - to_transform) 
skewness(transformed, na.rm = TRUE)
hist(transformed)

#At risk 2000

to_transform <- indicators_wide_inv_centred$`proportion.at.risk.2000-`
transformed <- indicators_transformed$`proportion.at.risk.2000-`

skewness(to_transform, na.rm = TRUE)
hist(to_transform)

indicators_transformed$`proportion.at.risk.2000-` <- 1/(max(to_transform + 1, 
                                                      na.rm = TRUE) - to_transform) 
skewness(indicators_transformed$`proportion.at.risk.2000-`, na.rm = TRUE)
hist(indicators_transformed$`proportion.at.risk.2000-`)

#At risk 2005

to_transform <- indicators_wide_inv_centred$`proportion.at.risk.2005-`
transformed <- indicators_transformed$`proportion.at.risk.2005-`

skewness(to_transform, na.rm = TRUE)
hist(to_transform)

indicators_transformed$`proportion.at.risk.2005-` <- 1/(max(to_transform + 1, 
                                                            na.rm = TRUE) - to_transform) 
skewness(indicators_transformed$`proportion.at.risk.2005-`, na.rm = TRUE)
hist(indicators_transformed$`proportion.at.risk.2005-`)

#At risk 2010

to_transform <- indicators_wide_inv_centred$`proportion.at.risk.2010-`
transformed <- indicators_transformed$`proportion.at.risk.2010-`

skewness(to_transform, na.rm = TRUE)
hist(to_transform)

indicators_transformed$`proportion.at.risk.2010-` <- 1/(max(to_transform + 1, 
                                                            na.rm = TRUE) - to_transform) 
skewness(indicators_transformed$`proportion.at.risk.2010-`, na.rm = TRUE)
hist(indicators_transformed$`proportion.at.risk.2010-`)

#At risk 2015

to_transform <- indicators_wide_inv_centred$`proportion.at.risk.2015-`
transformed <- indicators_transformed$`proportion.at.risk.2015-`

skewness(to_transform, na.rm = TRUE)
hist(to_transform)

indicators_transformed$`proportion.at.risk.2015-` <- 1/(max(to_transform + 1, 
                                                            na.rm = TRUE) - to_transform) 
skewness(indicators_transformed$`proportion.at.risk.2015-`, na.rm = TRUE)
hist(indicators_transformed$`proportion.at.risk.2015-`)

# Red List Index Aves 2005

to_transform <- indicators_wide_inv_centred$red.list.index.Aves.2005
transformed <- indicators_transformed$red.list.index.Aves.2005

skewness(to_transform, na.rm = TRUE)
hist(to_transform)

indicators_transformed$red.list.index.Aves.2005 <- 1/(max(to_transform + 1, 
                                                            na.rm = TRUE) - to_transform) 
skewness(indicators_transformed$red.list.index.Aves.2005, na.rm = TRUE)
hist(indicators_transformed$red.list.index.Aves.2005)


# ** Transformed boxplots ----

indicators_transformed_long <- melt(indicators_transformed, 
                          id.vars = 'ecoregion_id')

transformed_boxplots <- ggplot(indicators_transformed_long) +
  geom_boxplot(aes(x = variable, y = value)) +
  theme(axis.text.x=element_text(angle= 45,hjust=1))

transformed_boxplots

indicator_values_transformed <- indicators_transformed_long %>%
                                merge(indicator_values_master[c("ecoregion_id",
                                                                "country",
                                                                "realm")],
                                      by = "ecoregion_id") %>%
                                rename(indicator_year = variable,
                                       raw_indicator_value = value) %>% # TODO: Change to transformed
                                mutate(year = str_sub(indicator_year, start= -4)) %>%
                                mutate(indicator_year = as.character(indicator_year)) %>%
                                mutate(year = ifelse(year == "980-", "1980",
                                                     ifelse(year == "015-",
                                                            "2015",
                                                     ifelse(year == "000-",
                                                                   "2000",
                                                     ifelse(year == "990-",
                                                                     "1990",
                                                     ifelse(year == "005-", 
                                                            "2005",
                                                     ifelse(year == "010-", 
                                                            "2010", 
                                                     ifelse(year == "2005", 
                                                            "2005", NA)))))))) %>%
                                mutate(indicator = removeNumbers(indicator_year)) %>%
                                mutate(indicator = str_replace_all(indicator,
                                                                   '[[:punct:]]',' ')) 

saveRDS(indicator_values_transformed, file.path(indicator_outputs,
                                                paste(location, eco_version,
                                                      "indicator_vals_transformed.rds",
                                                      sep = "_")))

# Correlations ----

# Non-transformed

# indicator_names <- c("Ecoregion_id", "BII_A_2005","BHI_2005", "BHI_2010", "BHI_2015",
#                      "HFP_1990", "HFP_2010","Endangered_1980", "Endangered_1990",
#                      "At_risk_2000", "At_risk_2005", "At_risk_2010", 
#                      "At_risk_2015",  "Extinction_1980",   "Extinction_1990",
#                      "Extinction_2000", "Extinction_2005",   "Extinction_2010",
#                      "Extinction_2015", "RLI_Amph_1980", "RLI_Amph_2000",
#                      "RLI_Birds_1980", "RLI_Birds_1990", "RLI_Birds_2000",
#                      "RLI_Birds_2005", "RLI_Birds_2010", "RLI_Mamm_1990",
#                      "RLI_Mamm_2005", "BII_R_2005")


indicator_names <- c("Ecoregion_id", "BII","BHI", "BHI_2010", "BHI_2015",
                     "HFP_1990", "HFP","Endangered_1980", "Endangered_1990",
                     "At_risk_2000", "Endangered", "At_risk_2010", 
                     "At_risk_2015",  "Extinction_1980",   "Extinction_1990",
                     "Extinction_2000", "Extinction_2005",   "Extinction_2010",
                     "Extinction_2015", "RLI_Amph_1980", "RLI_Amph",
                     "RLI_Birds_1980", "RLI_Birds_1990", "RLI_Birds_2000",
                     "RLI_Birds", "RLI_Birds_2010", "RLI_Mamm_1990",
                     "RLI_Mammals", "BII_R_2005")

colnames(indicators_wide_inv_centred) <- indicator_names
                                     
                   
indicator_matrix <- indicators_wide_inv_centred %>%
                    select(-Ecoregion_id) %>%
                    na.omit(.)

# Split by location

indicators_wic <- indicators_wide_inv_centred %>%
                  merge(ecoregion_map[c("ecoregion_id", "REALM")],
                        by.x = "Ecoregion_id",
                        by.y = "ecoregion_id") 

indicators_wic_realm <- split(indicators_wic, indicators_wic$REALM)

indicators_wic_oceania <- indicators_wic_realm[["Oceania"]]

indicators_wic_australasia <- indicators_wic_realm[["Australasia"]]

indicators_wic_realm[["Global"]] <- indicators_wic

realm_matrices <- list()

for (i in seq_along(indicators_wic_realm)) {
  
  realm_matrices[[i]] <- indicators_wic_realm[[i]] %>%
                         select(-Ecoregion_id, -REALM, - geometry) %>%
                         na.omit(.)

}

names(realm_matrices) <- names(indicators_wic_realm)

# Get pictoral correlation matrices

realm_correlations <- list()

for (i in seq_along(realm_matrices)) {
  
  name <- paste(names(realm_matrices)[i], "indicator correlation matrix", sep = " ")
  
  subset_matrix <- realm_matrices[[i]][, c(1,2,6,10,20,24,27)]
  
  if(nrow(subset_matrix) == 0) {
    
    print("no data for this realm")
    
  } else {
    
    correlation_matrix <- ggcorr(subset_matrix, label = FALSE, nbreaks = 9) + 
      ggtitle(name)
    
    ggsave(file.path(indicator_outputs, 
                     paste(name, "indicator_correlation_matrix.png", sep = "_")),
           correlation_matrix, device = "png")
    
    realm_correlations[[i]] <- correlation_matrix
  }
}

realm_correlations[[9]]




pearson_correlations_all_years <- as.data.frame(cor(indicator_matrix, 
                                            method = "pearson"))

spearman_correlations_all_years <- as.data.frame(cor(indicator_matrix, 
                                                    method = "spearman"))

# All indicators with 2005 data

indicator_matrix_2005 <- indicator_matrix[, c(1,2,10,16,23,27,28)]

indicator_correlations_2005 <- ggcorr(indicator_matrix_2005)

# All indicators with 2010 data

indicator_matrix_2010 <- indicator_matrix[, c(3,6,11,17,25)]

indicator_correlations_2010 <- ggcorr(indicator_matrix_2010)

# RLI all timepoints and 1990 HFP

indicator_matrix_rli <- indicator_matrix[, c(5,21,22,23,24,25)]

indicator_correlations_rli <- ggcorr(indicator_matrix_rli)

# BHI all timepoints and 1990 HFP

indicator_matrix_bhi <- indicator_matrix[, c(5,2, 3, 4)]

indicator_correlations_bhi <- ggcorr(indicator_matrix_bhi)

# Extinctions all timepoints and 1990 HFP

indicator_matrix_ex <- indicator_matrix[, c(5,13,14,15,16,17,18)]

indicator_correlations_ex <- ggcorr(indicator_matrix_ex)

# Indicators from different data sources 2005 - 2010

indicator_matrix_ind <- indicator_matrix[, c(1,2,6,10,20,24,27)]

indicator_correlations_ind <- ggcorr(indicator_matrix_ind, palette = "BuPu",
                                     label = TRUE)


# RLI Birds 1980 and following extinctions

indicator_matrix_er <- indicator_matrix[, c(21, 13, 14,15,16,17,18)]

indicator_correlations_er <- ggcorr(indicator_matrix_er)


# Transformed

indicator_matrix_t <- indicators_transformed %>%
  select(-ecoregion_id) %>%
  na.omit(.)

pearson_correlations_all_years_t <- as.data.frame(cor(indicator_matrix_t, 
                                                    method = "pearson"))

spearman_correlations_all_years_t <- as.data.frame(cor(indicator_matrix_t, 
                                                     method = "spearman"))

# mutate(HFP_adjusted_old = ifelse(raw_indicator_value > 33.29037,
#                                  33.29037, raw_indicator_value)) %>%
# mutate(RLI_adjusted_old = ifelse(raw_indicator_value == 0, NA,
#                ifelse(raw_indicator_value > 0 & 
#                         raw_indicator_value < 0.9538,
#                       0.9538, raw_indicator_value))) %>%
# mutate(year = as.numeric(year)) %>%
# mutate(HFP_adjusted_old = ifelse(indicator != 
#                           "mean human footprint index", NA,
#                           HFP_adjusted_old)) %>%
# mutate(RLI_adjusted_old = ifelse(grepl("red list index Aves",
#                                 indicator), 
#                                 RLI_adjusted_old, NA)) %>%
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

# Scatterplots ----

# * Subset by time ----

# Non-transformed

indicators_by_year <- split(indicator_values_master, indicator_values_master$year)
indicators_by_year_names <- paste(location, names(indicators_by_year), sep = "_")

scatterplots_by_year <- list()

for (i in seq_along(indicators_by_year)){
  
  scatterplots_by_year[[i]] <- produce_scatterplots(indicators_by_year[[i]],
                                                    indicators_by_year_names[[i]],
                                                    save = FALSE)
  
}

scatterplots_by_year[[4]]

# Centred and inverted 

indicators_wic <- indicators_wic %>%
  rename(ecoregion_id = Ecoregion_id) %>%
  select(-geometry, -REALM)

indicators_wic_long <- melt(indicators_wic, 
                            id.vars = 'ecoregion_id')


indicator_values_wic <- indicators_wic_long  %>%
  rename(indicator_year = variable,
         raw_indicator_value = value) %>% # TODO: Change to transformed
  mutate(year = str_sub(indicator_year, start= -4)) %>%
  mutate(indicator_year = as.character(indicator_year)) %>%
  mutate(indicator = removeNumbers(indicator_year)) %>%
  mutate(indicator = str_replace_all(indicator,
                                     '[[:punct:]]',' '))

indicators_by_year <- split(indicator_values_wic, indicator_values_wic$year)
indicators_by_year_names <- paste(location, names(indicators_by_year), sep = "_")

scatterplots_by_year <- list()

for (i in seq_along(indicators_by_year)){
  
  scatterplots_by_year[[i]] <- produce_scatterplots(indicators_by_year[[i]],
                                                    indicators_by_year_names[[i]],
                                                    save = FALSE)
  
}

scatterplots_by_year[[1]]

# Oceania ----

indicators_wic <- indicators_wic_oceania %>%
  rename(ecoregion_id = Ecoregion_id) %>%
  select(-geometry, -REALM)

indicators_wic_long <- melt(indicators_wic, 
                            id.vars = 'ecoregion_id')


indicator_values_wic <- indicators_wic_long  %>%
  rename(indicator_year = variable,
         raw_indicator_value = value) %>% # TODO: Change to transformed
  mutate(year = NA) %>%
  mutate(indicator = as.character(indicator_year)) 

oceania_comparisons <- indicator_values_wic %>%
                       filter(indicator == "BII"| 
                              indicator == "BHI"|
                              indicator == "HFP"|
                              indicator == "Endangered"|
                              indicator == "RLI_Amph"|
                              indicator == "RLI_Birds"|
                              indicator == "RLI_Mammals")

length(unique(oceania_comparisons$ecoregion_id))

oceania_bhi_amph <- oceania_comparisons %>%
  filter(indicator == "RLI_Amph"| 
           indicator == "BHI")

produce_scatterplots(oceania_bhi_amph, "oceania_bhi_amph", save = TRUE)

oceania_scatterplots <- produce_scatterplots(oceania_comparisons, "oceania", save =FALSE)

plot(indicators_wic_oceania$HFP,indicators_wic_oceania$RLI_Amph)

oceania_master <- indicator_values_master %>%
                  filter(realm == "Oceania")

oceania_master_bhi <- oceania_master %>%
                      filter(indicator == "BHI plants") %>%
                      select(-country) %>%
                      distinct(.) %>%
                      group_by(ecoregion_id, year) %>%
                      arrange(.)

ggplot(oceania_master_bhi) +
  geom_line(aes(x = year, 
                y = raw_indicator_value, 
                color = ecoregion_id))


# Australasia ----

indicators_wic <- indicators_wic_australasia %>%
  rename(ecoregion_id = Ecoregion_id) %>%
  select(-geometry, -REALM)

indicators_wic_long <- melt(indicators_wic, 
                            id.vars = 'ecoregion_id')


indicator_values_wic <- indicators_wic_long  %>%
  rename(indicator_year = variable,
         raw_indicator_value = value) %>% # TODO: Change to transformed
  mutate(year = NA) %>%
  mutate(indicator = as.character(indicator_year)) 

australasia_comparisons <- indicator_values_wic %>%
  filter(indicator == "BII"| 
           indicator == "BHI"|
           indicator == "HFP"|
           indicator == "Endangered"|
           indicator == "RLI_Amph"|
           indicator == "RLI_Birds"|
           indicator == "RLI_Mammals")

length(unique(oceania_comparisons$ecoregion_id))

aus_bhi_amph <- australasia_comparisons %>%
  filter(indicator == "RLI_Amph"| 
           indicator == "BHI")

aus_bhi_amph$RLI_Amph_t <- 1/(max(aus_bhi_amph$RLI_Amph_t + 1, na.rm = TRUE) - 
  aus_bhi_amph$RLI_Amph_t)

1/(max(indicators_wide_inv_centred$`proportion.at.risk.1980-`+1, 
       na.rm = TRUE) - 
     indicators_wide_inv_centred$`proportion.at.risk.1980-`) 

produce_scatterplots(aus_bhi_amph, "aus_bhi_amph", save = TRUE)

oceania_scatterplots <- produce_scatterplots(oceania_comparisons, "oceania", save =FALSE)

plot(indicators_wic_oceania$HFP,indicators_wic_oceania$RLI_Amph)

oceania_master <- indicator_values_master %>%
  filter(realm == "Oceania")

oceania_master_bhi <- oceania_master %>%
  filter(indicator == "BHI plants") %>%
  select(-country) %>%
  distinct(.) %>%
  group_by(ecoregion_id, year) %>%
  arrange(.)

ggplot(oceania_master_bhi) +
  geom_line(aes(x = year, 
                y = raw_indicator_value, 
                color = ecoregion_id))

# Transformed

indicators_by_year <- split(indicators_values_transformed, indicator_values_transformed$year)
indicators_by_year_names <- paste(location, "transformed", names(indicators_by_year), sep = "_")

scatterplots_by_year <- list()

for (i in seq_along(indicators_by_year)) {
  
  scatterplots_by_year[[i]] <- produce_scatterplots(indicators_by_year[[i]],
                                                    indicators_by_year_names[[i]],
                                                    save = FALSE)
  
}

scatterplots_by_year[[2]]

# * Subset by realm ----

# * Subset by time by realm ----

# Split data into nested list (level 1 = year, level 2 = realm)

indicators_by_year_by_realm <- list()

for (i in seq_along(indicators_by_year)) {
  
  indicators_timepoint <- indicators_by_year[[i]]
  
  indicators_by_year_by_realm[[i]] <- split(indicators_timepoint, 
                                            indicators_timepoint$realm)
}

indicators_by_year_by_realm <- flatten(indicators_by_year_by_realm)
realm_names <- names(indicators_by_year_by_realm)
timepoints <- names(indicators_by_year)
time_names <- rep(timepoints, each = length(unique(realm_names)))
indicators_by_year_by_realm_names <- paste(realm_names, time_names , sep = "_")

scatterplots_by_year_by_realm <- list()

for (i in seq_along(indicators_by_year_by_realm)){
  
  scatterplots_by_year_by_realm[[i]] <- produce_scatterplots(indicators_by_year_by_realm[[i]],
                                                             indicators_by_year_by_realm_names[[i]],
                                                    save = TRUE)
  
}

names(scatterplots_by_year_by_realm) <- indicators_by_year_by_realm_names

# Look at Australiasia 2005

scatterplots_by_year_by_realm[[27]]

# Look at Oceania 2005

scatterplots_by_year_by_realm[[31]]

# * Subset effect of land use on extinction risk ----

land_use_atrisk_values <- indicator_values_master %>%
                    filter(indicator_year == "human footprint index 1990"|
                           indicator_year == "proportion at risk 1990"|
                             indicator_year == "proportion at risk 2000"|
                             indicator_year == "proportion at risk 2005"|
                             indicator_year == "proportion at risk 2010"|
                             indicator_year == "proportion at risk 2015")

hfp_vs_atrisk_scatterplots <- produce_scatterplots(land_use_atrisk_values, 
                              paste(location, eco_version,
                                    "HFP1990_vs_AtRisk_time_series",
                                    sep = "_"), save = TRUE)

# * Subset effect of land use on extinction  ----

land_use_extinction_values <- indicator_values_master %>%
    filter(indicator_year == "human footprint index 1990"|
           indicator_year == "proportion extinct 1990"|
           indicator_year == "proportion extinct 2000"|
           indicator_year == "proportion extinct 2005"|
           indicator_year == "proportion extinct 2010"|
           indicator_year == "proportion extinct 2015")

hfp_vs_extinction_scatterplots <- produce_scatterplots(land_use_extinction_values, 
                                                   paste(location, eco_version,
                                                   "HFP1990_vs_Extinction_time_series",
                                                         sep = "_"), save = TRUE)


hfp_vs_extinction_scatterplots

# * Subset only unrelated indicators  ----

subset_values <- indicator_values_master %>%
    filter(indicator_year == "human footprint index 1990"|
           indicator_year == "human footprint index 2010"|
           indicator_year == "abundance biodiversity intactness index 2005"|
           indicator_year == "richness biodiversity intactness index 2005"|
           indicator_year == "proportion extinct 2005")

subset_scatterplots <- produce_scatterplots(subset_values, 
                                            paste(location, eco_version,
                                            "subset_indicators",
                                             sep = "_"), save = FALSE)


subset_scatterplots


ggplot(data = indicator_values_master_wide_scaled, 
       aes(x = abundance.biodiversity.intactness.index.2005, 
           y = human.footprint.index.1990)) +
  geom_point() + 
  coord_trans(x = "log10", y = "log10")

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
           select(ecoregion_id, 
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

ecoregion_map <- ecoregion_map %>% rename(ecoregion_id = "ECO_ID")

# Prepare the data (inverted and centred)

indicators_wic <- indicators_wic %>%
                  rename(ecoregion_id = Ecoregion_id) %>%
                  select(-geometry, -REALM)

indicators_wic_long <- melt(indicators_wic, 
                            id.vars = 'ecoregion_id')


indicator_values_wic <- indicators_wic_long  %>%
  rename(indicator_year = variable,
         raw_indicator_value = value) %>% # TODO: Change to transformed
  mutate(year = str_sub(indicator_year, start= -4)) %>%
  mutate(indicator_year = as.character(indicator_year)) %>%
  mutate(indicator = removeNumbers(indicator_year)) %>%
  mutate(indicator = str_replace_all(indicator,
                                     '[[:punct:]]',' ')) 
# %>%
#   merge(indicator_values_master[c("ecoregion_id",
#                                   "country",
#                                   "realm")],
#         by = "ecoregion_id")
# 

indicator_map_data <- left_join(ecoregion_map, indicator_values_wic,
                                by = "ecoregion_id")
indicator_map_data_countries <- indicator_map_data %>%
        merge(indicator_values_master[c("ecoregion_id", "country")], by = "ecoregion_id")

# * BII Richness ----

richness_bii_map <- indicator_map_data %>%
                    filter(indicator == 
                             "BII R ") %>%
                    map_indicators(.$raw_indicator_value,
                                   "Richness BII 2005", 
                                   "right")
richness_bii_map

if (save_outputs == "yes") {

ggsave(file.path(indicator_outputs, paste(location,
                                          "richness_bii_ecoregion_map.png", 
                                          sep = "_")), 
       richness_bii_map,  device = "png")

}

# * BII Abundance ----

abundance_bii_map <- indicator_map_data %>%
  filter(indicator == 
           "BII A ") %>%
  map_indicators(.$raw_indicator_value,
                 "Abundance BII 2005", 
                 "right")
abundance_bii_map

if (save_outputs == "yes") {
  
  ggsave(file.path(indicator_outputs, paste(location,
                                            "abundance_bii_ecoregion_map.png", 
                                            sep = "_")),  
         abundance_bii_map,  device = "png")

}
# * BHI Plants ----

bhi_map <- indicator_map_data %>%
           filter(indicator == "BHI ") %>%
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
  
  ggsave(file.path(indicator_outputs, paste(location,
                            "bhi_ecoregion_map.png", 
                                            sep = "_")), 
         bhi_map,  device = "png")
  
  # st_write(bhi_2015_plants_sf, file.path(indicator_outputs, 
  #                                        paste(location, 
  #                                              "_bhi_plants_2015_map.shp")))
  # saveRDS(bhi_map, file.path(indicator_outputs,paste(location,
  #                                                    "bhi_plants_map.rds",
  #                                                    sep = "_")))
  
}

# * Proportion at risk ----

at_risk_map <- indicator_map_data %>%
               filter(indicator == "At risk ") %>%
               map_indicators(.$raw_indicator_value,
                             "Proportion\nof species\nat risk", 
                             "right")

at_risk_map

if(save_outputs == "yes") {
  
  ggsave(file.path(indicator_outputs, paste(location,
                                            "proportion_at_risk_ecoregion_map.png", 
                                            sep = "_")),  
at_risk_map,  device = "png")
    

}

# * Proportion extinct ----

extinct_map <- indicator_map_data %>%
               filter(indicator == "Extinction ") %>%
               map_indicators(.$raw_indicator_value,
                               "Proportion\nof species\nextinct", 
                               "right")

extinct_map

if(save_outputs == "yes") {
  
  ggsave(file.path(indicator_outputs, paste(location,
                                            "proportion_extinct_ecoregion_map.png", 
                                            sep = "_")),  
         extinct_map,  device = "png")
  
  
}
  
  
# Map the Red List Index by class

# * RLI Birds ----

indicator_map_data <- indicator_map_data_all %>% filter(year == 2005)


birds_rli_map <- indicator_map_data %>%
                 filter(indicator == "RLI Birds ") %>%
                 map_indicators(.$raw_indicator_value,
                                "Red List\nIndex (Birds)", 
                                "right")

birds_rli_map

if (save_outputs == "yes") {
  
  ggsave(file.path(indicator_outputs, paste(location,
                                            "2005_birds_rli_ecoregion_map.png", 
                                            sep = "_")),  
         birds_rli_map,  device = "png")
  
}


# * RLI Mammals ----

indicator_map_data <- indicator_map_data_all %>% filter(year == 2005)

mammals_rli_map <- indicator_map_data %>%
                   filter(indicator == "RLI Mamm ") %>%
                   map_indicators(.$raw_indicator_value,
                                   "Red List\nIndex (Mammals)", 
                                   "right")

mammals_rli_map

if (save_outputs == "yes") {
  
  ggsave(file.path(indicator_outputs, paste(location,
                                            "mammals_rli_ecoregion_map.png", 
                                            sep = "_")),  
         mammals_rli_map,  device = "png")
  
}

# * RLI Amphibians ----

indicator_map_data <- indicator_map_data_countries %>%
                      mutate(country = as.character(country)) %>%
                      filter(country == "Fiji")

amphibians_rli_map <- indicator_map_data %>%
                      filter(indicator == "RLI Amph ") %>%
                      map_indicators(.$raw_indicator_value,
                                     "Red List\nIndex (Amphibians)", 
                                     "right")

amphibians_rli_map

if (save_outputs == "yes") {
  
  ggsave(file.path(indicator_outputs, paste(location,
                                            "amphibians_rli_ecoregion_map.png", 
                                            sep = "_")),  
         amphibians_rli_map,  device = "png")
  
}

# * RLI Reptiles ----

# reptiles_rli_map <- indicator_map_data %>%
#   filter(indicator == "red list index Amphibia") %>%
#   map_indicators(.$raw_indicator_value,
#                  "Red List\nIndex (Amphibians)", 
#                  "right")
# 
# reptiles_rli_map
# 
# if (save_outputs == "yes") {
#   
#   ggsave(file.path(indicator_outputs, paste(location,
#                                             "reptiles_rli_ecoregion_map.png", 
#                                             sep = "_")),  
#          reptiles_rli_map,  device = "png")
#   
# }
# * HFP ----

indicator_map_data_all <- indicator_map_data
  
ecoregion_subset <- ecoregion_country_df %>%
    filter(CNTRY_NAME == "Australia") %>%
    unique(.)
  
indicator_map_data <- indicator_map_data_all[indicator_map_data_all$ecoregion_id %in% 
                                   ecoregion_subset$ECO_ID,] 
  
save_outputs == "yes"

indicator_map_data <- indicator_map_data_all %>% filter(year == 2010)

hfp_map <- indicator_map_data %>%
           filter(indicator == "HFP ") %>%
           map_indicators(.$raw_indicator_value,
                          "Human\nFootprint\nIndex 2010",
                          "right")
hfp_map


if (save_outputs == "yes") {
  
  ggsave(file.path(indicator_outputs, paste(location, "hfp_ecoregion_map.png", 
                                            sep = "_")), 
         hfp_map,  device = "png")
  
}




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
  
  ggsave(file.path(indicator_outputs, "hfp_rli_scatterplot.png"), rli_hfp,  device = "png")
  
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
  
  ggsave(file.path(indicator_outputs, "hfp_rli_scatterplot_gradient.png"), 
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
  
  ggsave(file.path(indicator_outputs, "hfp_rli_scatterplot_multi_colour_scales.png"), 
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

