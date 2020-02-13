# Using R version 3.5.2

# https://stackoverflow.com/questions/13070706/how-to-connect-r-with-access-database-in-64-bit-window

library(tidyverse)
library(readxl)
library(stringr)

# Set directory

if (Sys.info()['nodename'] == "SIMONE-PC") {
  
  basedir <- 'C:/Users/Simone/Dropbox/Deakin/Chapter_2_Extinction_test'
  setwd(basedir)
  
}  

if (Sys.info()['nodename'] == "ANALYTIX2") {
  
  basedir <- 'C:/Users/ssteven/Desktop/Chapter_2_Extinction_test'
  setwd(basedir)
}

if (Sys.info()['nodename'] == "20FMPC0C6GH9") {
  
  basedir <- 'C:/Users/ssteven/Dropbox/Deakin/Chapter_2_Extinction_test'
  setwd(basedir)
}

###############################################################################
## Functions ##
###############################################################################

## Function to subset list by name

filter_by_pattern <- function(pattern, your_list){
  
    require(purrr)
    names(your_list) %>% 
    str_detect(pattern) %>%
    keep(your_list, .)
  
}

# Function to rbind mismatched dataframes and populate columns with NA

rbind_all_columns <- function(x, y) {
  
  x_diff <- setdiff(colnames(x), colnames(y))
  y_diff <- setdiff(colnames(y), colnames(x))
  
  x[, c(as.character(y_diff))] <- NA
  
  y[, c(as.character(x_diff))] <- NA
  
  return(rbind(x, y))
}

###############################################################################
############################## 1. BUILD DATABASE ##############################
###############################################################################

# Load tables from WildFinder database (converted into .xlsx files from .mdb 
# database)

#' TODO: Figure out if there's a way to load the tables directly from the .mdb
#' TODO: Consolidate species name mismatches - refer to Stewart's code


file_path <- "N:/Quantitative-Ecology/Simone/BHI_data/wildfinder_xlsx"
file_names <- list.files(file_path)
files <- file.path(file_path, file_names)
tables <- lapply(files, read_excel) # Disregard warnings
names(tables) <- str_remove(file_names, ".xlsx")

# Change all the column names to lower case

for (i in seq_along(tables)) {

  names(tables[[i]]) <- tolower(names(tables[[i]]))
  
}

# Add ecoregion to species red list

redlist_species <- as.data.frame(do.call(cbind,tables$redlist_species))
ecoregions_species <- tables$ecoregion_species
ecoregions <- tables$ecoregions
common_names <- tables$common_names
# species <- tables$species
# wwf_species <- tables$Wwf_species

# Remove un-needed tables

# rm(tables)

# Standardise species ID colnames (BUT - check)
#' TODO: Check these species IDs are the WWF ones because redlist also have diff
#' ID numbers

colnames(ecoregions_species)[3] <- "wwf_species_id"
colnames(common_names)[3] <- "wwf_species_id"
colnames(wwf_species)[1] <- "wwf_species_id"

# Merge required tables to create a dataframe of species occurring in each ecoregion

#' TODO: Check the ecoregion endemic variable means what I think it does

species_df <- redlist_species %>% 
              merge(common_names[c("common_name_id", "common_name", "wwf_species_id")], 
                          by = "wwf_species_id") %>%
              merge(ecoregions_species[c("ecoregion_code","wwf_species_id",
                                                     "eco_endemic")], by = "wwf_species_id" ) %>%
              merge(ecoregions[c("ecoregion_code","ecoregion_name",
                                             "ecoregion_area")]) %>%
              dplyr::mutate(binomial = paste(genus, species, sep = " ")) %>%
              dplyr::mutate(source = "WildFinder")

# Clear memory by removing large objects

# rm(redlist_species, common_names, ecoregions, ecoregions_species)



# Load Cooke database - more species per ecoregion

cooke_database <- read.csv("N:/Quantitative-Ecology/Simone/BHI_data/cooke_database/species_site.csv")

## standardise columns to match species_df

cooke_database <- cooke_database %>%
                  dplyr::mutate(binomial_split = binomial) %>%
                  tidyr::separate(binomial_split, c("genus","species")) %>%
                  dplyr::mutate(source = "Cooke et al 2019") %>%
                  merge(species_df[c("binomial","wwf_species_id")], 
                  by = "binomial") %>%
                  distinct(.)


names(cooke_database) <- c("binomial", "ecoregion_code", "genus","species", 
                           "source", "wwf_species_id")

head(cooke_database)
tail(cooke_database)

## populate missing columns ready to bind with species_df

missing_cols <- setdiff(names(species_df), names(cooke_database))

cooke_wildfinder <- cooke_database %>% 
                    merge(species_df[c(missing_cols, "wwf_species_id")], 
                          by = "wwf_species_id")

# TEST - use a subset of cooke

# cooke_database <- cooke_database[1:200,]
                  
# Add cooke species to wildfinder species list

x <- rbind_all_columns(species_df, cooke_database)

head(x)








###############################################################################
########################## GET SUMMARY STATS ###########################
###############################################################################

# Not complete

extinct_species <- species_df %>%
                   filter(rlscrit == "EX")

redlist_by_ecoregion <- species_df %>%
  group_by(ecoregion_name,ecoregion_code, rlscrit) %>%
  count(ecoregion_name,ecoregion_code, rlscrit)


extinctions_by_ecoregion <- redlist_by_ecoregion %>%
  filter(rlscrit == "EX")