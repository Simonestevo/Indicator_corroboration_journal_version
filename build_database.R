# Using R version 3.5.2

# https://stackoverflow.com/questions/13070706/how-to-connect-r-with-access-database-in-64-bit-window

library(raster)
library(tidyverse)
library(readxl)
library(stringr)
library(reshape2)
library(sf)
library(spData)

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

inputs <- "N:/Quantitative-Ecology/Simone/Extinction_test/inputs"
outputs <- "N:/Quantitative-Ecology/Simone/Extinction_test/outputs"

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
#' TODO: Work out best way to consistently deal with variables where I've had
#' to merge multiple sources


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
#' TODO: Rename 'variable' column to "year" or similar

long_species_data <- melt(species_data, id.vars = c("binomial","wwf_species_id","common_name",
                                       "ecoregion_code", "eco_endemic", "genus",
                                       "species", "source", "group"),
                     value.name = "redlist_status")

# TEMPORARY CODE

## Save the long_species_data while we are working on it

# write_csv(long_species_data, paste(outputs, "draft_long_species_data.csv", sep = "/"))
# saveRDS(long_species_data, paste(outputs, "draft_long_species_data.rds", sep = "/"))

long_species_data <- readRDS(file.path(outputs, 'draft_long_species_data.rds'))

###############################################################################
########################## GET ECOREGIONS FOR EXTINCT SPP #####################
###############################################################################


# Get the ecoregion map & subset to required variables

ecoregion_map <- st_read(paste(inputs,"official_teow_wwf", sep = "/"))
em_small <- ecoregion_map %>% dplyr::select(eco_code, ECO_NAME, geometry)

# Get the range maps of extinct species

extinct_ranges <- st_read(paste(inputs,"redlist_extinct_species_range_maps", sep = "/"))

# Get the ecoregions for each extinct species

extinct_ranges_ecoregions <- st_join(extinct_ranges, em_small)

# Add this information to the species data.  We want to fill in the gaps
# of the missing ecoregion data for the extinct species, but also add any
# additional extinct species to our species data

extinct_species_w_ecoregions <- as.data.frame(extinct_ranges_ecoregions %>%
                                dplyr::select(BINOMIAL, eco_code)) %>%
                                dplyr::select(-geometry) %>%
                                dplyr::rename(binomial = BINOMIAL) %>%
                                dplyr::rename(ecoregion_code = eco_code) %>%
                                dplyr::mutate(source = "iucn_redlist") %>%
                                dplyr::mutate(redlist_status = "EX") %>%
                                dplyr::mutate(ecoregion_code = as.character(ecoregion_code))


#' TODO: Work out best way to fill in the missing ecoregion codes

## Current approach: Merge the two sources and create a new column which is a
## patchwork of different sources.  If this turns out to be a good approach,
## will also need to do this with the other patchwork variables as we work through
## them

long_species_data <- merge(long_species_data, extinct_species_w_ecoregions, 
                           by = "binomial", all = TRUE)

long_species_data <- long_species_data  %>%
                     dplyr::mutate(merged_ecoregion_code = 
                                   ifelse(is.na(ecoregion_code.x), 
                                   ecoregion_code.y, ecoregion_code.x))

###############################################################################
########################## GET SUMMARY STATS ###########################
###############################################################################

#' TODO: Add ecoregion name above so it is included in the species_data

# Get the number of species in each ecoregion


species_by_ecoregion <- long_species_data %>%
                        group_by(merged_ecoregion_code) %>%
                        summarize(n_distinct(binomial))

names(species_by_ecoregion) <- c("merged_ecoregion_code", "number_of_species")


# Get number of extinct species (grouping by ecoregion doesn't work well yet bc
# most extinct species haven't been assigned an ecoregion - TBD)

extinct_species <- long_species_data %>%
                   filter(redlist_status.x == "EX")  %>%
                   group_by(merged_ecoregion_code) %>%
                   summarise(n_distinct(binomial)) 

names(extinct_species) <- c("merged_ecoregion_code", "number_of_species_extinct")


proportion_extinct <- species_by_ecoregion %>%
                      merge(extinct_species, 
                            by = "merged_ecoregion_code", all = TRUE) %>%
                      dplyr::mutate(proportion_extinct = number_of_species_extinct/number_of_species) 

###############################################################################
########################## MAP SUMMARY STATS ###########################
###############################################################################

extinction_map_data <- inner_join(em_small, proportion_extinct[
  c("merged_ecoregion_code", "number_of_species_extinct")], 
                             by = c("eco_code" = "merged_ecoregion_code"))

extinction_map <- ggplot(extinction_map_data) +
                  geom_sf(aes(fill = number_of_species_extinct)) +
                  scale_fill_viridis_c(trans = "sqrt", alpha = .4)

extinction_map

