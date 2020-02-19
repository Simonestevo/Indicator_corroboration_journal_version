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

file_path <- "N:/Quantitative-Ecology/Simone/BHI_data/wildfinder_csv"
file_names <- list.files(file_path)
files <- file.path(file_path, file_names)
tables <- lapply(files, read_csv) # Disregard warnings
names(tables) <- str_remove(file_names, ".csv")

# Change all the column names to lower case

for (i in seq_along(tables)) {

  names(tables[[i]]) <- tolower(names(tables[[i]]))
  
}

# Get list of species per ecoregion (doesn't include binomial or common name)

ecoregions_species <- tables$ecoregion_species

# Add common names and binomial scientific names 

common_names <- tables$common_names
scientific_names <- cbind(tables$redlist_species$wwf_species_id, 
                          tables$redlist_species$genus,
                          tables$redlist_species$species)
colnames(scientific_names) <- c("wwf_species_id", "genus", "species")


species_by_ecoregion <- ecoregions_species %>%
                        merge(common_names[c("species_id", "common_name")], all = TRUE) %>%
                        dplyr::rename(wwf_species_id = species_id) %>%
                        dplyr::select(c("wwf_species_id", "common_name", 
                                        "ecoregion_code", "eco_endemic")) %>%
                        merge(scientific_names, by = "wwf_species_id", all = TRUE) %>%
                        dplyr::mutate(binomial = paste(genus, species, sep = " ")) %>%
                        dplyr::mutate(source = "WildFinder")

# Remove un-needed tables

# rm(tables)

# Clear memory by removing large objects

# rm(redlist_species, common_names, ecoregions, ecoregions_species)



# Load Cooke database - more species per ecoregion

cooke_database <- read.csv("N:/Quantitative-Ecology/Simone/BHI_data/cooke_database/species_site.csv")

## standardise columns to match species_df

cooke_database <- cooke_database %>%
                  dplyr::mutate(source2 = "Cooke et al 2019") %>%
                  dplyr::rename(ecoregion_code = eco)


# TEST - use a subset of cooke

# cooke_database <- cooke_database[1:200,]
                  
# Add cooke species to wildfinder species list then melt in to long form

test <- merge(species_by_ecoregion, cooke_database,  by = "binomial", all = TRUE)









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