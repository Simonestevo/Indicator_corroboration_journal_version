# Using R version 3.5.2

# https://stackoverflow.com/questions/13070706/how-to-connect-r-with-access-database-in-64-bit-window

library(dplyr)
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


###############################################################################
## Data wrangling ##
###############################################################################

# Load tables from WildFinder database (converted into .xlsx files from .mdb 
# database)

#' TODO: Figure out if there's a way to load the tables directly from the .mdb

file_path <- "N:/Quantitative-Ecology/Simone/BHI_data/wildfinder_xlsx"
file_names <- list.files(file_path)
files <- file.path(file_path, file_names)
tables <- lapply(files, read_excel) # Disregard warnings
names(tables) <- str_remove(file_names, ".xlsx")

# Add ecoregion to species red list

#' TODO: Just pull files we need directly from folder instead of reading all into
#' a list?

redlist_species <- tables[[14]]
ecoregions_species <- tables[[6]]
ecoregions <- tables[[7]]
common_names <- tables[[4]]

# Rename ecoregions column names to match redlist

colnames(ecoregions_species)[3] <- "WWF_SPECIES_ID"
colnames(common_names)[3] <- "WWF_SPECIES_ID"

# Merge required tables

#' TODO: Check the ecoregion endemic variable means what I think it does

species_df <- redlist_species %>% 
  merge(common_names[c("COMMON_NAME_ID", "COMMON_NAME", "WWF_SPECIES_ID")], 
        by = "WWF_SPECIES_ID") %>%
              merge(ecoregions_species[c("ECOREGION_CODE","WWF_SPECIES_ID",
                                         "ECO_ENDEMIC")], by = "WWF_SPECIES_ID" ) %>%
              merge(ecoregions[c("ECOREGION_CODE","ECOREGION_NAME",
                                 "ECOREGION_AREA")]) %>%
              dplyr::filter(ECO_ENDEMIC == "Y")

                   
                          