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
library(PerformanceAnalytics)



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

produce_scatterplots <- function(indicator_values, name, save, variable) {
  
  scatterplot_directory <- file.path(indicator_outputs, 
                                     "scatterplots")
  
  if( !dir.exists( scatterplot_directory ) ) {
    
    dir.create( scatterplot_directory, recursive = TRUE )
    
  }
  
  if (variable == "centred") {
    
    indicator_values_2 <- indicator_values %>%
    dplyr::select(ecoregion_id, indicator_year,
                  centred_indicator_value) %>%
    distinct(.)
    
    indicator_values_wide <- indicator_values_2 %>%
      spread(key = indicator_year, 
             value = centred_indicator_value) 
  
  } else if (variable == "raw") {
  
    indicator_values_2 <- indicator_values %>%
      dplyr::select(ecoregion_id, indicator_year,
                    raw_indicator_value) %>%
      distinct(.)
    
    indicator_values_wide <- indicator_values_2 %>%
      spread(key = indicator_year, 
             value = raw_indicator_value) 
  
  }
  

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

map_indicators <- function(data, indicator_variable, title, legend, reverse, year) {
  
if(reverse == TRUE) {

if(missing(year)) {
    
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

  } else {

data <- data %>%
        filter(year == year)

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
                  theme(legend.position = legend) 
                
return(indicator_map)
    
  }
  
} else {

  if(missing(year)) {
    
    indicator_map <-  ggplot(data) +
      geom_sf(aes(fill = indicator_variable), colour = "black", 
              size = 0.05, show.legend = 'fill') +
      scale_fill_viridis_c(alpha = .8,
                           na.value = "grey70") +
      theme(axis.line = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank()) +
      labs(fill = title) +
      theme(legend.position = legend) +
      facet_wrap(~ year)
    
    return(indicator_map)
    
  } else {
    
    data <- data %>%
      filter(year == year)
    
    indicator_map <-  ggplot(data) +
      geom_sf(aes(fill = indicator_variable), colour = "black", 
              size = 0.05, show.legend = 'fill') +
      scale_fill_viridis_c(alpha = .8,
                           na.value = "grey70") +
      theme(axis.line = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank()) +
      labs(fill = title) +
      theme(legend.position = legend) 
    
    return(indicator_map)
    
    }
  }
}


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

# As data frame

ecoregions <- as.data.frame(ecoregion_map) %>% dplyr::select(-geometry)


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

if ((paste(location, eco_version, "indicator_values_master.rds", sep = "_") %in%
     list.files(indicator_outputs))) {



indicator_values_master <- readRDS(file.path(indicator_outputs,
                                       paste(location, eco_version,
                                             "indicator_values_master.rds",
                                           sep = "_")))
} else {
  

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

# * Tidy species data ----

#' TODO: Sort this out in the build_database script

# Check for duplicates (same tsn but different binomials)

length(unique(species_data$binomial))
length(unique(species_data$tsn))

# Species that aren't duplicated

species_data1 <- species_data %>% 
  group_by_at(vars(-binomial, -source)) %>% 
  filter(n() < 2) 

# Find species that are duplicated and pick the BI record

species_data2 <- species_data %>% 
  group_by_at(vars(-binomial, -source)) %>% 
  filter(n() > 1) %>%
  filter(source == "birdlife_international")

# Combine

species_data <- rbind(species_data1, species_data2)

rm(species_data1, species_data2)

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

rm(duplicates)

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

ECO <- east_australia
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

# Check at risk values behave as anticipated

ECO <- east_australia
atrisk_test <- at_risk_values %>% filter(ecoregion_id == ECO)

ggplot(atrisk_test) +
  geom_line(aes(x = year, y = raw_indicator_value)) 

# # Red List Index all classes ----

# * 2005 only ----

if ((paste(location, eco_version, "proportion_at_risk.rds", sep = "_") %in% 
     list.files(indicator_outputs))) {
  
  #' TODO: Work out why this is producing so many NaNs
  
rli_all_2005_values <- readRDS(file.path(indicator_outputs, 
                                      paste(location, eco_version, 
                                            "RLI_all_classes.rds",
                                            sep = "_"))) 
} else {

species_data_for_rli_2005 <- species_data %>%
                             select(-source) %>%
                             distinct(.) %>%
                             filter(class == "Amphibia"| class == "Aves"|
                                    class == "Mammalia") %>%
                             filter(tsn != 1444976) %>%
                             filter(decade == 2005) %>%
                             mutate(class = "all")


species_data_ecoregions_2005 <- split(species_data_for_rli_2005, 
                                      species_data_for_rli_2005$ecoregion_id)

# Remove the empty lists with no data in them

species_data_ecoregions_2005 <- list.clean(species_data_ecoregions_2005, 
                                fun = is.null, recursive = FALSE)

# Calculate the RLI per ecoregion

rli_ecoregions_2005 <- list()

for (i in seq_along(species_data_ecoregions_2005)) {
  
  eco <- species_data_ecoregions_2005[[i]][[1]][1] # Pull out species list for one ecoregion
  
  rli_ecoregions_2005[[i]] <- calculate_red_list_index(species_data_ecoregions_2005[[i]]) # Calculate RLI for each ecoregion for one timestep, one class
  
  print(paste("Processed redlist values in year 2005 for ecoregion",
              eco, sep = " "))
  
} 

rli_all_classes_2005 <- do.call(rbind, rli_ecoregions_2005)

# Format

rli_all_2005_values <- rli_all_classes_2005 %>%
                       ungroup(.) %>%
                       select(redlist_assessment_year, Ecoregion_id, RLI) %>%
                       mutate(indicator = "RLI",
                              year = 2005,
                              ecoregion_id = Ecoregion_id) %>%
                       rename(raw_indicator_value = RLI) %>%
                       select(indicator, year, ecoregion_id, raw_indicator_value)

saveRDS(rli_all_2005_values, file.path(indicator_outputs, 
                                  paste(location, eco_version, 
                                        "RLI_all_classes.rds",
                                        sep = "_")))

}


# * Multiple time points ----
# Split the species data into different time points (output should be dataframes of 
# species red list status in a nested list with levels: Time point, Ecoregion)

# only 8400 spp left by now
# species_data_for_rli_all_classes <- species_data_complete_all_classes %>%
#                                     filter(class != "reptile") %>%
#                                     filter(tsn != 1444976)
# # Split the data by timepoint
# 
# rli_timestep <- split(species_data_complete_all_classes, 
#                       species_data_complete_all_classes$redlist_assessment_year) 
# 
# 
# time_ecoregion_list <- list()
# 
# for (i in seq_along(rli_timestep)) {
#     
#     # Split by each ecoregion as well (so have a list of species for each timepoint,
#     # ecoregion)
#     
#     time_ecoregions <- split(rli_timestep[[i]], 
#                                       rli_timestep[[i]]$ecoregion_id)
#     
#     # Remove the empty lists with no data in them
#     
#     time_ecoregions <- list.clean(time_ecoregions, 
#                                   fun = is.null, recursive = FALSE)
#     
#     time_ecoregion_list[[i]] <- time_ecoregions
#   
# }
#   
# 
# # Calculate the RLI per ecoregion, per timepoint (output should be
# # dataframes of RLI values by ecoregion, by timepoint)
# 
# rli_ecoregions_single_year <- list()
# rli_ecoregions_all_years <- list()
# 
# for (i in seq_along(time_ecoregion_list)) {
#   
#   timestep <- time_ecoregion_list[[i]] # Get list of species for each ecoregion in one timestep
#   
#   year <- timestep[[1]][[5]][1]
#   
#   print(paste("Processing year", year, sep = " "))
#   
#   #class_all_timepoints <- list()
#   
#   for (j in seq_along(timestep)) {
#     
#       eco <- timestep[[j]][[1]][1] # Pull out species list for one ecoregion
#       
#       rli_ecoregions_single_year[[j]] <- calculate_red_list_index(timestep[[j]]) # Calculate RLI for each ecoregion for one timestep, one class
#       
#       print(paste("Processed redlist values in year", year,
#                   "for ecoregion", eco, sep = " "))
#       
#       } 
#     
#     rli_ecoregions_single_year_df <- do.call(rbind, rli_ecoregions_single_year) # One dataframe for one time point
#     
#     rli_ecoregions_all_years[[i]] <- rli_ecoregions_single_year_df # Put the time point into a list of timepoints
#     
#   }
# 
# # Convert back into one dataframe
# 
# rli_ecoregions_all_years_df <- do.call(rbind, rli_ecoregions_all_years)
# 
# saveRDS(rli_ecoregions_all_years_df, 
#         file.path(indicator_outputs, paste(location, eco_version,
#                   "rli_all_classes_not_formatted.rds", sep = "_")))
# 
# # Check if any amphibian values have the wrong year and fix (2008 should be 1980)
# 
# rli_values <- rli_values %>%
#   mutate(redlist_assessment_year = ifelse(class == "Amphibia" & 
#                                             redlist_assessment_year == 2008, 1980,
#                                           redlist_assessment_year))
# # Format
# 
# rli_values <- rli_values %>%
#   filter(RLI != 0) %>%
#   mutate(indicator = paste("RLI", class, sep = " ")) %>%
#   rename(raw_indicator_value = RLI) %>%
#   rename(ecoregion_id = Ecoregion_id) %>%
#   ungroup()%>%
#   filter(ecoregion_id != 0) %>%
#   distinct(.) %>%
#   filter(redlist_assessment_year < 2016) %>%
#   filter(!(class == "Aves" & redlist_assessment_year == 2004)) %>%
#   mutate(year = ifelse(redlist_assessment_year < 1990, 1980,
#                        ifelse(redlist_assessment_year > 1989 &
#                                 redlist_assessment_year < 2000, 1990,
#                               ifelse(redlist_assessment_year > 1999 &
#                                        redlist_assessment_year < 2005, 2000,
#                                      ifelse(redlist_assessment_year > 2004 &
#                                               redlist_assessment_year < 2010, 2005,
#                                             ifelse(redlist_assessment_year > 2009 &
#                                                      redlist_assessment_year < 2015, 2010,
#                                                    ifelse(redlist_assessment_year > 2014 &
#                                                             redlist_assessment_year < 2021, 2015,
#                                                           NA))))))) %>%
#   dplyr::select(indicator, year, 
#                 ecoregion_id, raw_indicator_value) 
# # mutate(RLI_inverted = 1 - RLI) %>%
# # mutate(RLI_scaled_inverted =
# #          scale_to_1(RLI_inverted)) %>%
# # mutate(RLI_adjusted_old = ifelse(RLI == 0, NA,
# #                           ifelse(RLI > 0 & RLI < 0.9538,
# #                                         0.9538, RLI))) %>%
# # mutate(RLI_adjusted = ifelse(RLI == 0, NA,
# #                              pmin(pmax(RLI,
# #                       quantile(RLI, .05, na.rm = TRUE))))) %>%
# # mutate(RLI_adjusted_inverted = 1 - RLI_adjusted)
# 
# saveRDS(rli_values, file.path(indicator_outputs, 
#                               paste(location, eco_version, 
#                                     "red_list_index.rds",
#                                     sep = "_"))) 
# }
# 
# # Remove reptiles for the moment b/c data hasn't been checked and it's behaving weirdly
# 
# rli_values <- rli_values %>%
#   filter(indicator != "red list index Reptilia")
# 
# # Check rli values behave as anticipated
# 
# ECO <- cardamom
# rli_test <- rli_values %>% filter(ecoregion_id == ECO) %>% 
#   filter(indicator == 'RLI Amphibia')
# 
# ggplot(rli_test) +
#   geom_line(aes(x = year, y = raw_indicator_value))


# Human Footprint Index  ----


# Extract values from each raster

if ((paste(location, eco_version,
            "human_footprint_index.rds", sep = "_") %in% 
      list.files(indicator_outputs))) {

hfp_values <- readRDS(file.path(indicator_outputs,
                                paste(location, eco_version,
                                      "human_footprint_index.rds", sep = "_")))
} else {
  
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
                     raw_indicator_value = hfpmean) %>%
              st_set_geometry(NULL) %>%
    dplyr::select(indicator, year, ecoregion_id, 
                  raw_indicator_value)
  
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
           raw_indicator_value = hfpmean) %>%
    st_set_geometry(NULL) %>%
    dplyr::select(indicator, year, ecoregion_id, 
                  raw_indicator_value)
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
                     raw_indicator_value = hfpmean) %>%
              st_set_geometry(NULL) %>%
              dplyr::select(indicator, year, ecoregion_id, 
                            raw_indicator_value)
  
# * HFP 2013 ----

# Read in the Human Footprint Index data 

hfp_by_ecoregion_2013 <- read.csv(paste(inputs, "human_footprint_index", 
                                        "human_footprint_index_by_ecoregion.csv",
                                        sep = "/"))

hfp_2013 <- hfp_by_ecoregion_2013 %>%
                         rename(ecoregion_name = ECO_NAME) %>%
                         # merge(ecoregion_map[c("ecoregion_name", "ECO_ID")], 
                         #       by = "ECO_ID") %>% # This will subset automatically if you subset by country
                         mutate(year = 2013) %>%
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


hfp_values <- rbind(hfp_2000, hfp_2005, hfp_2010, hfp_2013)

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


lpi_inputs <- lpi_inputs %>%
              filter(System == "Terrestrial")


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

# Calculate the indicator

lpi_values_formatted <- list()

for (i in seq_along(lpi_values_by_ecoregion)) {
  
  lpi_values_formatted[[i]] <- get_lpi(lpi_values_by_ecoregion[[i]],
                                       out[[i]])
  
}

lpi_values <- do.call(rbind, lpi_values_formatted)

lpi_values_subset <- lpi_values %>%
                     filter(year == 1980|
                              year == 1990|
                              year == 2000|
                              year == 2005|
                              year == 2010|
                              year == 2013|
                              year == 2015) 

length(unique(lpi_values_subset$ecoregion_id))

lpi_test <- lpi_values %>% filter(ecoregion_id == 291)

ggplot(lpi_test) +
  geom_line(aes(x = year, y = raw_indicator_value)) 

saveRDS(lpi_values, file.path(indicator_outputs, 
                                     paste(location, eco_version, 
                                           "lpi.rds",
                                           sep = "_")))
  }
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

}

# ADDITIONAL VARIABLES ---

# LPI number of records ----

if (paste(location, eco_version, "LPI_record_values.rds", sep = "_") %in% 
    list.files(indicator_outputs)) {
  
lpi_record_values <- readRDS(file.path(indicator_outputs, 
                                            paste(location, eco_version,
                                                  "LPI_record_values.rds", 
                                                  sep = "_")))
} else {

lpi_inputs <- readRDS(file.path(indicator_outputs, 
                                paste(location, eco_version,
                                      "LPI_2020_input_data_with_ecoregions.rds", 
                                      sep = "_")))

lpi_input_summary <- lpi_inputs %>%
                     group_by(ECO_ID) %>%
                     summarise(number_of_records = n_distinct(Binomial),
                                .groups = "drop_last")

lpi_record_values <- lpi_input_summary %>%
                      mutate(indicator = "LPI_records",
                             year = NA) %>%
                      rename(ecoregion_id = ECO_ID,
                             raw_indicator_value = number_of_records) %>%
                      select(indicator_columns)

saveRDS(lpi_record_values, file.path(indicator_outputs, 
                                     paste(location, eco_version, 
                                           "LPI_record_values.rds",
                                        sep = "_")))
}

# Red List number of records ----

if (paste(location, eco_version, "RLI_record_values.rds", sep = "_") %in% 
    list.files(indicator_outputs)) {
  
rli_record_values <- readRDS(file.path(indicator_outputs, 
                                         paste(location, eco_version,
                                               "RLI_record_values.rds", 
                                               sep = "_")))
} else {

rli_record_values <- species_by_ecoregion %>%
                     select(ecoregion_id, decade, number_of_species) %>%
                     distinct(.) %>%
                     mutate(indicator = "RLI_records") %>%
                     rename(raw_indicator_value = number_of_species) %>%
                     select(-decade) %>%
                     mutate(year = NA) %>%
                     dplyr::select(indicator_columns) %>%
                     distinct(.)

saveRDS(rli_record_values, file.path(indicator_outputs, 
                                     paste(location, eco_version, 
                                           "RLI_record_values.rds",
                                           sep = "_")))
}

rm(species_data, species_by_ecoregion)

# Ecoregion area ----

ecoregion_area_km2 <- read.csv(paste(file.path(inputs, "ecoregions_2017",
                               "ecoregions2017_area_km2.csv", sep ="")))

ecoregion_area_values <- ecoregion_area_km2 %>%
                         rename(ecoregion_id = Eco..ID,
                                raw_indicator_value = Ecoregion.area..km2.) %>%
                         mutate(indicator = "ecoregion area km sq",
                                year = NA) %>%
                         select(indicator_columns)

# Ecoregion biome ----

ecoregion_realm_values <- ecoregion_map_all %>%
                          select(ECO_ID, BIOME_NAME) %>%
                          mutate(indicator = "Biome",
                                 year = NA) %>%
                          rename(ecoregion_id = ECO_ID,
                                 raw_indicator_value = BIOME_NAME) %>%
                          st_set_geometry(NULL) %>%
                          select(indicator_columns)

# Islands ----

# Simplify the ecoregion map so it doesn't take forever. Note - this produces
# a warning that st_simplify doesn't work on long/lat data, but it's just 
# a warning because all our maps *are* in decimal degrees

ecoregion_map_simple <- st_simplify(ecoregion_map, preserveTopology = TRUE,
                                    dTolerance = 0.1)

# Get continents

continents_wm <- st_read(file.path(inputs,"global_islands_explorer","continents.shp"))

# Remove unnecessary columns 

continents_wm <- continents_wm %>%
                 select(OBJECTID, geometry) %>%
                 mutate(classification = "continent")

# Reproject to WGS 1984

continents <- st_transform(continents_wm, crs = crs(ecoregion_map))

# Remove mw projection

rm(continents_wm)

continents <- st_simplify(continents, preserveTopology = TRUE,
              dTolerance = 0.1)

# Join continents and ecoregions

ecoregions_continents <- st_join(continents, ecoregion_map_simple, join = st_intersects)

# NB - compared this to non-simplified versions and produced the same result, although
# unsurprising given how big continents are

ecoregions_continents <- ecoregions_continents %>%
                         st_set_geometry(NULL) %>%
                         select(ECO_ID, ECO_NAME, classification) %>%
                         distinct(.)

rm(continents)

# Big islands

big_islands_wm <- st_read(file.path(inputs,"global_islands_explorer","big_islands.shp"))
head(big_islands_wm)

big_islands_wm <- big_islands_wm %>%
                  select(OBJECTID, geometry) %>%
                  mutate(classification = "big_island")

big_islands <- st_transform(big_islands_wm, crs = crs(ecoregion_map))

rm(big_islands_wm)

# Simplify the big islands file

#big_islands_all <- big_islands

big_islands <- big_islands_all

plot(big_islands$geometry)

big_islands <- st_simplify(big_islands, preserveTopology = TRUE,
                           dTolerance = 0.3)

plot(big_islands$geometry)

# Remove the continental ecoregions

`%notin%` <- Negate(`%in%`)

ecoregion_map_islands <- ecoregion_map_simple[ecoregion_map_simple$ECO_ID %notin% 
                                              ecoregions_continents$ECO_ID,] 
# Join it to the ecoregion map

system.time(ecoregions_big_islands <- st_join(big_islands, ecoregion_map_islands,
                                             join = st_intersects))

rm(big_islands, big_islands_all)

ecoregions_big_islands <- ecoregions_big_islands %>%
                          st_set_geometry(NULL) %>%
                          select(ECO_ID, ECO_NAME, classification) %>%
                          distinct(.)

saveRDS(ecoregions_big_islands, file.path(indicator_outputs, 
                                      paste(location, eco_version, 
                                            "ecoregions_big_islands.rds",
                                            sep = "_")))

# Small islands

small_islands_wm <- st_read(file.path(inputs,"global_islands_explorer","small_islands.shp"))

small_islands_wm <- small_islands_wm %>%
                    select(OBJECTID, geometry) %>%
                    mutate(classification = "small_island")

small_islands <- st_transform(small_islands_wm, crs = crs(ecoregion_map))

rm(small_islands_wm)

small_islands_all <- small_islands


# Simplify the small islands file

# small_islands <- st_simplify(small_islands, preserveTopology = TRUE,
#                            dTolerance = 0.1)

# Join it to the ecoregion map

# ecoregion_map_islands <- ecoregion_map_simple[ecoregion_map_simple$ECO_ID %notin% 
#                                                 ecoregions_continents$ECO_ID,]
# 
# ecoregion_map_islands <- ecoregion_map_islands[ecoregion_map_islands$ECO_ID %notin% 
#                                                 ecoregions_big_islands$ECO_ID,]

system.time(ecoregions_small_islands <- st_join(small_islands, ecoregion_map_islands, 
                                              join = st_intersects))

rm(small_islands)

ecoregions_small_islands <- ecoregions_small_islands %>%
                            st_set_geometry(NULL) %>%
                            select(ECO_ID, ECO_NAME, classification) %>%
                            distinct(.)

saveRDS(ecoregions_small_islands, file.path(indicator_outputs, 
                                          paste(location, eco_version, 
                                                "ecoregions_small_islands.rds",
                                                sep = "_")))

# Very small islands

very_small_islands_wm <- st_read(file.path(inputs,"global_islands_explorer","very_small_islands.shp"))

very_small_islands_wm <- very_small_islands_wm %>%
  select(OBJECTID, geometry) %>%
  mutate(classification = "very_small_island")

very_small_islands <- st_transform(very_small_islands_wm, crs = crs(ecoregion_map))

rm(very_small_islands_wm)

# Simplify the small islands file

very_small_islands <- st_simplify(very_small_islands, preserveTopology = TRUE,
                             dTolerance = 0.1)

# Join it to the ecoregion map

system.time(ecoregions_very_small_islands <- st_join(very_small_islands, ecoregion_map_islands, 
                                                join = st_intersects))

rm(very_small_islands)

ecoregions_very_small_islands <- ecoregions_very_small_islands %>%
  st_set_geometry(NULL) %>%
  select(ECO_ID, ECO_NAME, classification) %>%
  distinct(.)

saveRDS(ecoregions_very_small_islands, file.path(indicator_outputs, 
                                          paste(location, eco_version, 
                                                "ecoregions_very_small_islands.rds",
                                                sep = "_")))

x<- rbind(ecoregions_continents, ecoregions_big_islands, ecoregions_small_islands,
          ecoregions_very_small_islands)


# Threats ----

threat_files <-  list.files(interim_outputs, recursive = FALSE)[grepl("threat_data.rds",list.files(interim_outputs))]
threat_list <- lapply(file.path(interim_outputs, threat_files), readRDS)
all_threat_data <- do.call(rbind, threat_list)

threat_scheme <- read.csv(file.path(inputs, "iucn_threats\\iucn_threat_classification_scheme.csv"))
threat_scheme$code <- as.character(threat_scheme$code)

formatted_threat_data <- all_threat_data %>%
                         merge(threat_scheme[c("code", "headline", "headline_name")], by = "code") %>%
                         filter(scope == "Majority (50-90%)" | scope == "Whole (>90%)") %>%
                         filter(timing != "Future") %>%
                         filter(score == "Medium Impact: 6" |
                                 score == "Medium Impact: 7"|
                                 score == "High Impact: 8"|
                                 score == "High Impact: 9") %>%
                         distinct(.) %>%
                         merge(species_data[c("binomial", "tsn", "ecoregion_id")], by = "binomial") %>%
                         dplyr::select(-code, - invasive) %>%
                         dplyr::select(ecoregion_id, binomial, tsn, headline, headline_name, title) %>%
                         distinct(.) %>%
                         group_by(ecoregion_id, title) %>%
                         mutate(number_of_species_affected = n_distinct(tsn)) %>%
                         group_by(ecoregion_id) %>%
                         mutate(number_of_species = n_distinct(tsn)) %>%
                         select(ecoregion_id, headline_name, title, number_of_species_affected, number_of_species) %>%
                         distinct(.) %>%
                         mutate(proportion_affected = number_of_species_affected/number_of_species) %>%
                         group_by(ecoregion_id) %>%
                         filter(proportion_affected == max(proportion_affected)) %>%
                         group_by(ecoregion_id) %>%
                         mutate(threat_count = n()) %>%
                         summarise_all(~ toString(unique(.))) %>%
                         rename(predominant_threat = title)

head(formatted_threat_data)

ecoregion_threats <- formatted_threat_data %>%
                     select(ecoregion_id, predominant_threat, threat_count)

# Combine indicator values into a single dataframe ----

indicator_values <- rbind(extinction_values, 
                          at_risk_values, 
                          rli_all_2005_values,
                          hfp_values,
                          bhi_plants_values,
                          lpi_values,
                          lpi_record_values,
                          bii_richness_values,
                          bii_abundance_values)

indicator_values_master <- indicator_values %>%
                            merge(ecoregion_country_df[c("ECO_ID", 
                                                         "CNTRY_NAME")], 
                                  by.x = "ecoregion_id", by.y = "ECO_ID",
                                  all.x = TRUE, all.y = FALSE) %>%
                            distinct(.) %>%
                            rename(country = CNTRY_NAME) %>%
                            merge(ecoregions[c("ECO_ID", "REALM")], 
                                  by.x = "ecoregion_id",
                                  by.y = "ECO_ID") %>%
                            rename(realm = REALM) %>%
                            filter(ecoregion_id != 0) %>%
                            filter(indicator != "RLI Reptilia") %>%
                            mutate(indicator_abbreviated = 
                                   ifelse(indicator == "proportion at risk",
                                            "threatened",
                                   ifelse(indicator == "proportion extinct",
                                                   "extinct",
                                   ifelse(indicator == "abundance biodiversity intactness index",
                                                   "BIIab",
                                   ifelse(indicator == "richness biodiversity intactness index",
                                                    "BIIri",
                                   ifelse(indicator == "mean human footprint index",
                                                       "HFP",
                                   ifelse(indicator == "RLI Aves",
                                                 "BirdRLI",
                                   ifelse(indicator == "RLI Amphibia",
                                                 "AmphRLI",
                                   ifelse(indicator == "RLI Mammalia",
                                                 "MammRLI",
                                            indicator))))))))) %>%
                            mutate(indicator_year = paste(indicator_abbreviated, 
                                                          year, sep = " ")) %>%
                            mutate(ecoregion_id = as.numeric(ecoregion_id))

summary(indicator_values_master)

saveRDS(indicator_values_master, file.path(indicator_outputs,
                                           paste(location, eco_version,
                                                 "indicator_values_master.rds",
                                                 sep = "_")))

write.csv(indicator_values_master, file.path(indicator_outputs,
                                             paste(location, eco_version,
                                                   "indicator_values_master.csv",
                                                   sep = "_")))

# Tidy the workspace

rm(extinction_values, 
   at_risk_values, 
   rli_all_2005_values,
   hfp_values,
   bhi_plants_values,
   lpi_values,
   lpi_record_values,
   lpi_input_summary,
   lpi_inputs,
   bii_richness_values,
   bii_abundance_values,
   species_data_all,
   species_by_ecoregion,
   indicator_values)

# Map the indicators ----

ecoregion_map_data <- ecoregion_map %>%
                     rename(ecoregion_id = ECO_ID) %>%
                     select(ecoregion_id, geometry)

indicator_values_map <- indicator_values_master %>%
                        select(ecoregion_id, indicator, year, raw_indicator_value) %>%
                        distinct(.)

indicator_map_data <- left_join(ecoregion_map_data, indicator_values_map,
                                by = "ecoregion_id")

indicator_map_data <- indicator_map_data %>%
                      select(-country) %>%
                      distinct(.)


# * BII Richness ----

richness_bii_map <- indicator_map_data %>%
                    filter(indicator == 
                           "richness biodiversity intactness index") %>%
                    map_indicators(.$raw_indicator_value,
                                   "BII\n(richness)", 
                                   "right",
                                   TRUE,
                                   2005)
richness_bii_map

if (save_outputs == "yes") {

ggsave(file.path(indicator_outputs, paste(location,
                                          "richness_bii_ecoregion_map.png", 
                                          sep = "_")), 
       richness_bii_map,  device = "png")

}

rm(richness_bii_map)

# * BII Abundance ----

abundance_bii_map <- indicator_map_data %>%
                     filter(indicator == 
                               "abundance biodiversity intactness index") %>%
                     map_indicators(.$raw_indicator_value,
                                    "BII\n(abundance)", 
                                     "right",
                                    TRUE,
                                    2005)
abundance_bii_map

if (save_outputs == "yes") {
  
  ggsave(file.path(indicator_outputs, paste(location,
                                            "abundance_bii_ecoregion_map.png", 
                                            sep = "_")),  
         abundance_bii_map,  device = "png")

}

rm(abundance_bii_map)

# * BHI Plants ----

bhi_map <- indicator_map_data %>%
           filter(indicator == "BHI plants") %>%
           map_indicators(.$raw_indicator_value,
                          "BHI plants", 
                          "right", 2005)
bhi_map

# Create a separate sf object to save as shapefile

# bhi_2015_plants_sf <- indicator_map_data %>%
#                       filter(indicator == "biodiversity habitat index plants") %>%
#                       select(eco_id, ecoregion_name, eco_objectid,
#                              raw_indicator_value, geometry) %>%
#                       rename(BHI = raw_indicator_value)

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

rm(bhi_map)

# * Proportion at risk ----

at_risk_map <- indicator_map_data %>%
               filter(indicator == "proportion at risk") %>%
               map_indicators(.$raw_indicator_value,
                             "Proportion\nof species\nat risk", 
                             "right",
                             FALSE,
                             2005)

at_risk_map

if(save_outputs == "yes") {
  
  ggsave(file.path(indicator_outputs, paste(location,
                                            "proportion_at_risk_ecoregion_map.png", 
                                            sep = "_")),  
at_risk_map,  device = "png")
    

}

rm(at_risk_map)

# * Proportion extinct ----

extinct_map <- indicator_map_data %>%
               filter(indicator == "proportion extinct") %>%
               map_indicators(.$raw_indicator_value,
                               "Proportion\nof species\nextinct", 
                               "right",
                              FALSE,
                              2005)

extinct_map

if(save_outputs == "yes") {
  
  ggsave(file.path(indicator_outputs, paste(location,
                                            "proportion_extinct_ecoregion_map.png", 
                                            sep = "_")),  
         extinct_map,  device = "png")
  
  
}
  
# * RLI ----

rli_map <- indicator_map_data %>%
           filter(indicator == "RLI") %>%
           map_indicators(.$raw_indicator_value,
                                "Red List\nIndex", 
                                "right",
                          TRUE,
                          2005)

rli_map

if (save_outputs == "yes") {
  
  ggsave(file.path(indicator_outputs, paste(location,
                                            "2005_rli_ecoregion_map.png", 
                                            sep = "_")),  
         rli_map,  device = "png")
  
}



# * HFP ----

hfp_map <- indicator_map_data %>%
           filter(indicator == "HFP ") %>%
           map_indicators(.$raw_indicator_value,
                          "Human\nFootprint\nIndex",
                          "right",
                          FALSE,
                          2005)
hfp_map


if (save_outputs == "yes") {
  
  ggsave(file.path(indicator_outputs, paste(location, "hfp_ecoregion_map.png", 
                                            sep = "_")), 
         hfp_map,  device = "png")
  
}

# * LPI ----

lpi_map <- indicator_map_data %>%
  filter(indicator == 
           "LPI") %>%
  map_indicators(.$raw_indicator_value,
                 "LPI", 
                 "right",
                 TRUE,
                 2005)
lpi_map

if (save_outputs == "yes") {
  
  ggsave(file.path(indicator_outputs, paste(location,
                                            "lpi_ecoregion_map.png", 
                                            sep = "_")), 
         lpi_map,  device = "png")
  
}


# RLI vs HFP scatterplot ----

# Make the scatterplot with colours mapped to RLI, size mapped to HFP
names(indicators_wic) <- c("ecoregion_id", "BHI.plants.2005",
                           "BIIab.2005", "BIIri.2005", "extinct.2005",
                           "HFP.2005", "RLI.2005", "threatened.2005",
                           "REALM", "geometry")

rli_hfp_data <- indicators_wic %>% 
                select(ecoregion_id, HFP.2005, RLI.2005, REALM) %>%
                merge(lpi_record_values[c("ecoregion_id", "raw_indicator_value")],
                      by = "ecoregion_id") %>%
                rename(LPI.Records.2005 = raw_indicator_value)
  
rli_hfp <- ggplot(rli_hfp_data, aes(x = HFP.2005, y = RLI.2005, col = REALM)) + 
           geom_point(aes(size = LPI.Records.2005), alpha = 0.6) +
           scale_colour_viridis(discrete = TRUE ) +
           theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                panel.background = element_blank(), 
                axis.line = element_line(colour = "black")) +
          labs(x = "Human Footprint Index Ecoregion Values", 
               y = "Red List Index Ecoregion Values") +
          theme(plot.background = element_rect(fill = NA),
                axis.text = element_text(size = 14),
                axis.title = element_text(size = 14)) 
# + 
#           scale_y_continuous(breaks = c(0.00,0.25,0.50,0.75,1.00),
#                              limits = c(0.00,1.00)) +
#           scale_x_continuous(breaks = c(0,5,10,15,20,25,30))

rli_hfp


if (save_outputs == "yes") {
  
  ggsave(file.path(indicator_outputs, paste(location, "rli_hfp_scatterplot.png", 
                                            sep = "_")), 
         rli_hfp,  device = "png")
  
}

rli_hfp_bw <- ggplot(rli_hfp_data, aes(x = HFP.2005, y = RLI.2005)) + 
  geom_point(alpha = 0.6, size = 2) +
  scale_colour_viridis(discrete = TRUE ) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  labs(x = "Human Footprint Index Ecoregion Values", 
       y = "Red List Index Ecoregion Values") +
  theme(plot.background = element_rect(fill = NA),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14)) 
# + 
#           scale_y_continuous(breaks = c(0.00,0.25,0.50,0.75,1.00),
#                              limits = c(0.00,1.00)) +
#           scale_x_continuous(breaks = c(0,5,10,15,20,25,30))

rli_hfp_bw


if (save_outputs == "yes") {
  
  ggsave(file.path(indicator_outputs, paste(location, "rli_hfp_scatterplot_bw.png", 
                                            sep = "_")), 
         rli_hfp_bw,  device = "png")
  
}

# By realm

hfp_rli_data_realms <- split(rli_hfp_data, rli_hfp_data$REALM)

plots <- list()

for (i in seq_along(hfp_rli_data_realms)) {
  
  realm <- hfp_rli_data_realms[[i]]$REALM[1]
  
  plots[[i]] <- ggplot(hfp_rli_data_realms[[i]], aes(x = HFP.2005, y = RLI.2005)) + 
    geom_point(aes(size = LPI.Records.2005),alpha = 0.6) +
    scale_colour_viridis(discrete = TRUE ) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black")) +
    labs(x = "Human Footprint Index Ecoregion Values", 
         y = "Red List Index Ecoregion Values") +
    theme(plot.background = element_rect(fill = NA),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 14)) +
    labs(title = paste(realm, 
                  "Human Footprint Index x Red List Index ecoregion values", 
                  sep = " "))
  # + 
  #           scale_y_continuous(breaks = c(0.00,0.25,0.50,0.75,1.00),
  #                              limits = c(0.00,1.00)) +
  #           scale_x_continuous(breaks = c(0,5,10,15,20,25,30))

}

plots[[4]]




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

