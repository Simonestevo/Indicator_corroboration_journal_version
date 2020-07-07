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

# Set input and output locations ----

inputs <- "N:/Quantitative-Ecology/Simone/extinction_test/inputs"
outputs_parent <- "N:/Quantitative-Ecology/Simone/extinction_test/outputs"
save_outputs <- "no"
date <- Sys.Date()
country <- "Australia" # If not subsetting, set as NA, e.g. country <- NA
outputs <- "N:\\Quantitative-Ecology\\Simone\\extinction_test\\outputs\\2020-06-11_output_files\\"
aggregate_timepoints <- "yes"

# Set output directory

# if (save_outputs == "yes") {
#   
#   outputs_dir <- paste(date,"_output_files",sep="")
#   outputs <- file.path(outputs_parent, paste(date,"_output_files",sep=""))
#   
#   if( !dir.exists( outputs ) ) {
#     
#     dir.create( outputs, recursive = TRUE ) # create a new directory for today's outputs 
#     
#   }
#   
# } else {
#   
#   previous_outputs <- list.dirs(outputs_parent, recursive = FALSE)
#   
#   outputs <- previous_outputs[length(previous_outputs)] # get the most recent outputs folder to use
#   
# }

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
  grouped.data <- weighted_data %>% group_by(class.x, redlist_assessment_year)
  
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

# Calculate indicators ----

# Get ecoregions

ecoregion_map_all <- st_read(paste(inputs,"official_teow_wwf", sep = "/"))
ecoregion_map <- ecoregion_map_all %>% dplyr::select(eco_code, ECO_NAME, geometry)

names(ecoregion_map) <- c("ecoregion_code", "ECO_NAME", "geometry")

# Get ecoregion countries

ecoregion_country_df <- readRDS(paste(outputs_parent, 
                                      "ecoregion_country_data.rds", sep = "/"))

ecoregion_country_df <- ecoregion_country_df %>%
                        rename(ecoregion_code = eco_code) %>%
                        mutate(ecoregion_code = as.character(ecoregion_code)) %>%
                        na.omit(.) %>%
                        distinct(.)

if (!is.na(country)) {
  
    ecoregion_subset <- ecoregion_country_df %>%
    filter(CNTRY_NAME == country) %>%
    unique(.) 
}

if (!is.na(country)) {
  
  ecoregion_map <- inner_join(ecoregion_map, ecoregion_subset[c("ecoregion_code", 
                                                                "OBJECTID",
                                                                "CNTRY_NAME")],
                              by = "ecoregion_code")
  
  }
  
# Read in species data

species_data <- readRDS(file.path(outputs, "species_data.rds"))

# Remove cases without redlist status and/or assessment year
#' TODO: this might be better done in the build_database script? So there's no
#' wrangling of the species data in this script?

species_data <- species_data %>%
                drop_na(redlist_assessment_year, redlist_status) %>%
                distinct(.)

# Subset by test country

if (!is.na(country)) {
  
    species_data <- species_data %>%
    merge(ecoregion_subset[c("ecoregion_code", "CNTRY_NAME")], 
          by = "ecoregion_code") %>%
    unique(.)
  
}

# Extinction/Risk status ----

# TODO: Add in redlist proportion to indicator values for analysis

species_by_ecoregion <- species_data %>%
                        group_by(ecoregion_code, redlist_assessment_year) %>%
                        summarize(number_of_species = n_distinct(tsn),
                                  proportion_extinct = sum(redlist_status == "EX"| 
                                                         redlist_status =="EW")/
                                    number_of_species,
                                  proportion_endangered = sum(redlist_status == "EN")/
                                    number_of_species,
                                  proportion_critical = sum(redlist_status == "CR"|
                                                          redlist_status == "CR(PE)")/
                                    number_of_species,
                                  proportion_vulnerable = sum(redlist_status == "VU")/
                                    number_of_species,
                                  proportion_atrisk = sum(redlist_status == "EX"| 
                                                          redlist_status =="EW"|
                                                          redlist_status == "EN"|
                                                          redlist_status == "CR"|
                                                          redlist_status == "CR(PE)"|
                                                          redlist_status == "VU")/
                                                            number_of_species,
                                  proportion_lowrisk = sum(redlist_status == "LC"|
                                                            redlist_status == "NT")/
                                    number_of_species) %>%
                        mutate(check = proportion_extinct + 
                                       proportion_endangered +
                                       proportion_critical +
                                 proportion_vulnerable +
                                 proportion_lowrisk)

# Extinctions

extinction_values <- species_by_ecoregion %>%
                     mutate(indicator = "proportion extinct") %>%
                     select(indicator, redlist_assessment_year, 
                            ecoregion_code, proportion_extinct) %>%
                     rename(year = redlist_assessment_year, 
                            raw_indicator_value = proportion_extinct)

# Threatened status

at_risk_values <- species_by_ecoregion %>%
                 mutate(indicator = "proportion at risk") %>%
                 select(indicator, redlist_assessment_year, 
                         ecoregion_code, proportion_atrisk) %>%
                 rename(year = redlist_assessment_year, 
                         raw_indicator_value = proportion_atrisk)

# # Red List Index ----

# Split into different classes

class_list <- split(species_data, species_data$class.x)

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
              mutate(indicator = paste("red list index", class.x, sep = " ")) %>%
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
                         merge(ecoregion_map[c("ECO_NAME", "ecoregion_code")], 
                               by = "ECO_NAME") %>%
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

if (!is.na(country)) {
  
  full_hfp_by_ecoregion_2017 <- hfp_by_ecoregion_2017
  
  country_ecoregions <- ecoregion_subset$ecoregion_code
  
  hfp_by_ecoregion_2017 <- full_hfp_by_ecoregion_2017[full_hfp_by_ecoregion_2017$ecoregion_code %in% country_ecoregions, ]
  
}

# Biodiversity Intactness Index 2005 ----

if (!is.na(country)) {

bii_2005_raster_data <- raster(file.path(inputs, 
                                         "biodiversity_intactness_index\\asc1.tif"))

} else {

bii_2005_raster_data <- raster(file.path(inputs, 
                                         "biodiversity_intactness_index\\lbii.asc"))

}

#' TODO: Work out why this is producing so many NaNs

bii_ecoregion_map <- ecoregion_map %>% 
                     mutate(raw_indicator_value = raster::extract(bii_2005_raster_data, 
                                                  ecoregion_map, 
                                                  fun = mean, na.rm = TRUE))

#' TODO: Figure out how to deal with multiple polygons of the same ecoregion - once
#' you remove the geometry you end up with multiple values p/ecoregion

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

# Combine indicator values into a single dataframe ----

# TODO: Melt the adjusted values into long form

names(species_by_ecoregion) <- c("eco_code", "number_of_species")

ecoregions <- as.data.frame(ecoregion_map) %>% dplyr::select(-geometry)

indicator_values <- rbind(rli_values, hfp_by_ecoregion_2017, bii_values,
                          extinction_values, at_risk_values)

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
                                                 "CNTRY_NAME")], 
                          by = "ecoregion_code", all.x = TRUE, all.y = FALSE) %>%
                    distinct(.) %>%
                    rename(country = CNTRY_NAME)
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

# Split by year

indicator_values_time_list <- split(indicator_values, indicator_values$year)


if (aggregate_timepoints == "yes") {

# Aggregate years that are close to each other

indicators_1988 <- indicator_values_time_list[[1]]

indicators_1994_1996 <- rbind(indicator_values_time_list[[2]], 
                              indicator_values_time_list[[3]])

indicators_2000 <- indicator_values_time_list[[4]]

indicators_2004_2005 <- rbind(indicator_values_time_list[[5]],
                              indicator_values_time_list[[6]])

indicators_2008 <- indicator_values_time_list[[7]]

indicators_2012 <- indicator_values_time_list[[8]]

indicators_2016_2017 <- rbind(indicator_values_time_list[[9]],
                           indicator_values_time_list[[10]])

indicator_values_time_list <- list(indicators_1988, indicators_1994_1996,
                                   indicators_2000, indicators_2004_2005,
                                   indicators_2008, indicators_2012, 
                                   indicators_2016_2017)

rm(indicators_1988, indicators_1994_1996,
        indicators_2000, indicators_2004_2005,
        indicators_2008, indicators_2012, 
        indicators_2016_2017)

}

# Test for correlations between variables, for each year

correlations <- list()
years <- list()

for (i in seq_along(indicator_values_time_list)) {
  
  year <- range(indicator_values_time_list[[i]]$year)
  
  step1 <- indicator_values_time_list[[i]] %>% 
           select(ecoregion_code, 
                  raw_indicator_value, 
                  indicator, 
                  -year) 
  
  step2 <- dcast(step1, ecoregion_code ~ indicator, 
                 value.var = "raw_indicator_value", sum)
  
  step3 <- as.matrix(step2[2:ncol(step2)])
  
  if(ncol(step3) == 1) {
    
    print(paste("Year", year, "does not have values for more than one indicator to test correlations between"))
    
    rm(step1, step2, step3)
    
  } else {
  
  step4 <- step3[complete.cases(step3),]
  
  correlations[[i]] <- as.data.frame(cor(step4, method = "pearson"))
  
  years[[i]] <- year
  
  rm(step1, step2, step3, step4)
  
  }
  
}

# Remove the empty lists with no data in them

names(correlations) <- years

correlations <- list.clean(correlations, fun = is.null, recursive = FALSE)

# Map the indicators ----

indicator_map_data <- inner_join(ecoregion_map, indicator_values, 
                                 by = "ecoregion_code")

# Map the Red List Index by class

#' TODO: Put this into a loop - currently doesn't work because they have
#' different timesteps

birds_rli_map <- indicator_map_data %>%
                 filter(indicator == "red list index Aves") %>%
                 map_indicators(.$RLI_adjusted_old,
                                "Red List\nIndex (Birds)", 
                                "right")

mammals_rli_map <- indicator_map_data %>%
                   filter(indicator == "red list index Mammalia") %>%
                   map_indicators(.$raw_indicator_value,
                                   "Red List\nIndex (Mammals)", 
                                   "right")

amphibians_rli_map <- indicator_map_data %>%
                      filter(indicator == "red list index Amphibia") %>%
                      map_indicators(.$raw_indicator_value,
                                     "Red List\nIndex (Amphibians)", 
                                     "right")
# Map the Human Footprint Index

hfp_map <- indicator_map_data %>%
           filter(indicator == "human footprint index") %>%
           map_indicators(.$HFP_adjusted_old,
                          "Human\nFootprint\nIndex",
                          "right")
# Map the Biodiversity Intactness Index

bii_map <- indicator_map_data %>%
           filter(indicator == "biodiversity intactness index") %>%
           map_indicators(.$raw_indicator_value,
                           "Biodiversity\nIntactness\nIndex",
                           "right")


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

cor(indicator_values$HFP_scaled_inverted, indicator_values$RLI_scaled, 
    method = "pearson", use = "complete.obs")

cor(indicator_values$HFP_scaled_inverted, indicator_values$PLOTCAT, 
    method = "pearson", use = "complete.obs")

cor(indicator_values$RLI_scaled, indicator_values$WII_inverted, 
    method = "pearson", use = "complete.obs")


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

# Map indicators ----

indicator_map_data <- inner_join(ecoregion_map, indicator_values[
  c("eco_code", "RLI_adjusted", "RLI_adjusted_old",
    "HFP", "HFP_adjusted", "HFP_adjusted_old", "HFP_scaled_adjusted")], 
  by = "eco_code")


indicator_map_rli <- ggplot(indicator_map_data) +
  geom_sf(aes(fill = RLI_adjusted_old), colour = "black", 
          size = 0.05, show.legend = 'fill') +
  scale_fill_viridis_c(trans = "reverse", alpha = .8, 
                       na.value = "grey70") +         
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(fill = "Red List\nIndex (Birds)") +
  theme(legend.position = "none")

# trans = "log10", 

indicator_map_rli

if(save_outputs == "yes") {
  
  ggsave(file.path(outputs, "indicator_map_rli_dimensions_nl.png"), 
         indicator_map_rli, height = 4, width = 5, device = "png")
  
}


indicator_map_hfp <- ggplot(indicator_map_data) +
  geom_sf(aes(fill = HFP_adjusted_old), colour = "black", 
          size = 0.05, show.legend = 'fill') +
  scale_fill_viridis_c(trans = "log10", alpha = .8, 
                       na.value = "grey70") +         
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(fill = "Human\nFootprint\nIndex") +
  theme(legend.position = "none")

# trans = "log10", 

indicator_map_hfp

if(save_outputs == "yes") {
  
  ggsave(file.path(outputs, "indicator_map_hfp_dimensions_nl.png"), 
         height = 4, width = 5,indicator_map_hfp,  device = "png")
  
}

indicator_maps <- grid.arrange(indicator_map_rli, indicator_map_hfp)

if(save_outputs == "yes") {
  
  ggsave(file.path(outputs, "indicator_maps.png"), indicator_maps,  device = "png")
  
}

# Red List Index for mammals ----

# ecoregion_map_data_no_geometry <- as.data.frame(ecoregion_map)
# all_map_data_temp <- ecoregion_map_data_no_geometry

# ecoregion_map_data_no_geometry <- ecoregion_map_data_no_geometry %>%
#   select(eco_code, ECO_NAME) %>%
#   distinct(.)

# mammal_rli_by_ecoregion_df_2008 <- mammal_rli_by_ecoregion_df_2008 %>%
#                                    merge(ecoregion_map_data_no_geometry[c(
#                                        "eco_code", "ECO_NAME")], 
#                                             all = TRUE) %>%
#                                      distinct(.)
# 
# 
# mammal_rli_map_data <- inner_join(ecoregion_map, mammal_rli_by_ecoregion_df_2008[
#                                 c("eco_code", "RLI")], 
#                                 by = "eco_code")
# 
# test_mammal_subset <- mammal_rli_map_data[1:100,]
# 
# length(unique(mammal_rli_map_data$eco_code))
# 
# mammal_rli_map_data <- mammal_rli_map_data %>%
#                        mutate(RLI_inverted = 1 - RLI)
# 
# mammal_rli_map <- ggplot(mammal_rli_map_data) +
#                   geom_sf(aes(fill = RLI_inverted), colour = "black", 
#                           size = 0.05) +
#                   scale_fill_viridis_c(alpha = .8, na.value = "grey70") +         
#                   theme(axis.line = element_line(),
#                         panel.grid.major = element_blank(),
#                         panel.grid.minor = element_blank(),
#                         #panel.border = element_blank(),
#                         panel.background = element_blank()) +
#                   labs(fill = "Red List Index (Mammals)")
# 
# mammal_rli_map

# ggsave(paste(outputs, "rli_mammals_map_dark.png", sep = "/"), mammal_rli_map,  
#        device = "png")


# Red List Index for BIRDS ----

ecoregion_map_data_no_geometry <- as.data.frame(ecoregion_map)
all_map_data_temp <- ecoregion_map_data_no_geometry

ecoregion_map_data_no_geometry <- ecoregion_map_data_no_geometry %>%
  dplyr::select(eco_code, ECO_NAME) %>%
  distinct(.)

# Add ecoregion names
birds_rli_by_ecoregion_2016 <- birds_rli_by_ecoregion_2016 %>%
  merge(ecoregion_map_data_no_geometry[c("eco_code", 
                                         "ECO_NAME")], 
        all = TRUE) %>%
  distinct(.)

# Add geometry back in so it can be mapped

rli_map_data <- inner_join(ecoregion_map, birds_rli_by_ecoregion_2016[
  c("eco_code", "RLI")], 
  by = "eco_code")

# Check you've got all the ecoregions

length(unique(rli_map_data$eco_code))

# Add an inverted RLI variable so the colours can go either light to dark or 
# dark to light (fix this properly in the plots later?)

rli_map_data <- rli_map_data %>%
  mutate(RLI_inverted_original = 1 - RLI) %>%
  mutate(RLI_inverted = ifelse(RLI_inverted_original == 1, NA,
                               RLI_inverted_original))

# Make a map of birds RLI for 2016

rli_map <- ggplot(rli_map_data) +
  geom_sf(aes(fill = RLI_inverted), colour = "black", 
          size = 0.05) +
  scale_fill_viridis_c(trans = "log10", alpha = .8, na.value = "grey70") +         
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(fill = "Red List Index (Birds)")

# trans = "log10", 

rli_map

ggsave(paste(outputs, "rli_birds_map_dark.png", sep = "/"), rli_map,  device = "png")

# TEMPORARY CODE - smooth outliers by putting an upper limit on RLI value (see RLI_adjusted)

rli_map_data <- rli_map_data %>%
  mutate(RLI_inverted_original = 1 - RLI) %>%
  mutate(RLI_inverted = ifelse(RLI_inverted == 1, NA, RLI_inverted_original)) %>%
  mutate(RLI_adjusted = ifelse(RLI == 0, NA,
                               ifelse(RLI > 0 & RLI < 0.9538, 
                                      0.9538, RLI)))

# Make the new map - should have a more even distribution of colours instead of
# a couple of really high values and everything else kind of the same

rli_map_2 <-  ggplot(rli_map_data) +
  geom_sf(aes(fill = RLI_adjusted), colour = "black", 
          size = 0.05) +
  scale_fill_viridis_c(trans = "sqrt",direction = -1, alpha = .8, 
                       na.value = "grey70") +         
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(fill = "Red List Index (Birds)")


rli_map_2

ggsave(paste(outputs, "rli_birds_map_adjusted_not_inverted.png", sep = "/"), 
       rli_map_2,  device = "png")


# Add ecoregion codes




hfp_map <- ggplot(hfp_map_data) +
  geom_sf(aes(fill = HFP), colour = "black", 
          size = 0.05) +
  scale_fill_viridis_c(alpha = .8, na.value = "grey70") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(fill = "Human Footprint Index")


hfp_map

ggsave(paste(outputs, "hfp_map_dark_adjusted.png", sep = "/"), hfp_map,  device = "png")

dev.off()



## Previous way of calculating RLI used for NESP/WWF docs ----

# species_data <- species_data %>% 
#   filter(class.x == "Aves") %>%
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

