# Using R version 3.5.2

# https://stackoverflow.com/questions/13070706/how-to-connect-r-with-access-database-in-64-bit-window

library(tidyverse)
library(readxl)
library(stringr)
library(reshape2)
library(sf)

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
  inputs <- "N:/Quantitative-Ecology/Simone/Extinction_test/inputs"
  outputs <- "N:/Quantitative-Ecology/Simone/Extinction_test/outputs"
  
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

#' TODO: Nearly all extinct species have no ecoregion. Figure out how to fix?
#' TODO: Make sure wwf id reads in as numeric not character 
#' TODO: Add ecoregions to species_data
#' TODO: Find out what year wildfinder RL status is from
#' TODO: Figure out if there's a way to load the tables directly from the .mdb
#' TODO: Consolidate species name mismatches - refer to Stewart's code
#' TODO: List and work out how to populate missing data (incl. names)


file_path <- file.path(inputs, "wildfinder_csv")
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


wildfinder_database <-  ecoregions_species %>%
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

cooke_database <- read.csv(paste(inputs,"/cooke_database/species_site.csv", sep = ""))

## standardise columns to match species_df

cooke_database <- cooke_database %>%
                  dplyr::mutate(source = "Cooke et al 2019") %>%
                  dplyr::rename(ecoregion_code = eco)


# TEST - use a subset of cooke and wildfinder - delete this later
#' TODO: Remove the two lines below and run on entire dataset when hpc available

cooke_database <- cooke_database[sample(nrow(cooke_database), 10000), ]
wildfinder_database <- wildfinder_database[sample(nrow(wildfinder_database), 10000), ]
                  
# Fill in details so Cooke has the same columns as wildfinder

cooke_database <- cooke_database %>%
                  merge(wildfinder_database[c("wwf_species_id", "common_name", 
                                              "eco_endemic", "genus", "species",
                                              "binomial")], by = "binomial") %>%
                  dplyr::select("wwf_species_id", "common_name", "ecoregion_code",
                               "eco_endemic", "genus", "species",
                               "binomial", "source") %>%
                  dplyr::distinct(.) # line not tested yet

# Commbine the cooke and wildfinder databases
#' TODO: Work out why we're getting replicates from cooke

merged_databases <- rbind(wildfinder_database, cooke_database)

merged_databases <- distinct(merged_databases)

# To remove later when fix this up top

merged_databases$wwf_species_id <- as.numeric(merged_databases$wwf_species_id)


###############################################################################
########################## ADD RED LIST STATUS ###########################
###############################################################################

# Using Sergio's data

file_path <- file.path(inputs, "redlist_sergio")
file_names <- list.files(file_path)
files <- file.path(file_path, file_names)
tables <- lapply(files, read_csv)
names(tables) <- str_remove(file_names, ".csv")

# Add a column so we can identify the source of the data

group_tables <- list()

for (i in seq_along(tables)){
  
  group_tables[[i]] <- tables[[i]] %>% 
                       dplyr::mutate(group = names(tables[i]))
  
  
}

# Standardise column names etc. which are all inconsistent between groups

amphibians <- group_tables[[1]]

amphibians <- amphibians %>% 
              dplyr::mutate(binomial = paste(Genus, Species, sep = " ")) %>%
              dplyr::select(-c(Genus, Species)) %>%
              set_names(c("2004", "2008", "group", "binomial")) %>%
              dplyr::select("binomial","group","2004", "2008")

birds <- group_tables[[2]]

birds <- birds %>%
         set_names(c("binomial", "1988", "1994", "2000", "2004", "2008" , 
                     "2012" , "2016", "group" )) %>%
         dplyr::select("binomial", "group", everything())

mammals <- group_tables[[5]]  

mammals <- mammals %>%
           set_names(c("binomial", "1996", "2008", "group" )) %>%
           dplyr::select("binomial", "group", everything())

# Combine back into a list then dataframe

group_tables <- list(birds, mammals, amphibians)

species_redlist <- do.call(bind_rows, group_tables)

# Order columns by year

species_redlist <- species_redlist %>%
                   dplyr::select("binomial", "group", "1988", "1994", "1996", 
                                 "2000", "2004", "2008", "2012", "2016")


# Add RL status to the merged databases

species_data <- merged_databases %>%
                merge(species_redlist, by = "binomial", all = TRUE)

# Melt into long format

long_species_data <- melt(species_data, id.vars = c("binomial","wwf_species_id","common_name",
                                       "ecoregion_code", "eco_endemic", "genus",
                                       "species", "source", "group"),
                     value.name = "redlist_status")


###############################################################################
########################## GET SUMMARY STATS ###########################
###############################################################################

#' TODO: Add ecoregion name above so it is included in the species_data

# Get the number of species in each ecoregion


species_by_ecoregion <- long_species_data %>%
                        group_by(ecoregion_code) %>%
                        summarize(n_distinct(binomial))


# Get number of extinct species (can't group by ecoregion yet bc don't have them)

extinct_species <- long_species_data %>%
                   filter(redlist_status == "EX") #%>%
                   # group_by(ecoregion_code) %>%
                   # summarise(n())

# Get number of each RL status in each ecoregion

redlist_by_ecoregion <- long_species_data %>%
                        group_by(ecoregion_code, redlist_status) %>%
                        count(ecoregion_code, redlist_status)

###############################################################################
########################## MAP SUMMARY STATS ###########################
###############################################################################

library(sf)

ecoregion_map <- st_read(paste(inputs,"official_teow_wwf", sep = "/"))
