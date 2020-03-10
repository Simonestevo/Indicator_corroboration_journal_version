# https://stackoverflow.com/questions/13070706/how-to-connect-r-with-access-database-in-64-bit-window

# cd C:\\Users\\ssteven\\Dropbox\\Deakin\\Chapter_2_Extinction_test\\Extinction_test_code

# Using R version 3.5.2

## This script outputs:
## 1: A database of the WWF terrestrial ecoregions, species found
## in each ecoregion (this will not be a comprehensive list of every single species
## but aims to be as complete as possible), and their red list status, and
## 2. Summary statistics of number of species and number of extinct species

#' TODO: Remove objects when they are no longer needed

# Packages ----

library(raster)
library(tidyverse)
library(stringr)
library(reshape2)
library(sf)
library(spData)
library(taxize)
library(functionaltraits)

# Input and output locations ----

inputs <- "N:/Quantitative-Ecology/Simone/Extinction_test/inputs"
outputs <- "N:/Quantitative-Ecology/Simone/Extinction_test/outputs"

# Functions ----

## get_ecoregions: Function to load species range maps and join with the wwf ecoregion map
#' @return a dataframe with four columns (species binomial scientific name,
#' ecoregion, data source and red list status)
#' @param species directory name - a string that denotes the name of the directory
#' where the species range maps are saved
#' @param map sf object - the global wwf_teow map

get_ecoregions <- function(species_directory_name, map) {
  
  range_map <- st_read(paste(inputs,species_directory_name, sep = "/"))
  
  # Get the ecoregions for each species
  
  ranges_ecoregions <- st_join(range_map, map)
  
  # Standardise the data so we can add it to the species data easily
  
  species_w_ecoregions <- as.data.frame(ranges_ecoregions %>%
                                          dplyr::select(BINOMIAL, eco_code)) %>%
    dplyr::select(-geometry) %>%
    dplyr::rename(binomial = BINOMIAL) %>%
    dplyr::rename(ecoregion_code = eco_code) %>%
    dplyr::mutate(source = "iucn_redlist") %>%
    dplyr::mutate(redlist_status = "TBC") %>%
    dplyr::mutate(ecoregion_code = as.character(ecoregion_code))
  
  return(species_w_ecoregions)
  
}

#' Find traits and taxonomic information for a list of species
#' 
#' This function takes a list of species, collects taxonomic information for them, 
#' and then searches the databases for the valid scientific names, as well as any 
#' synonyms that were found. 
#' 
#'
#' @param species A vector of species names in the format "Genus species". Note: if the same name 
#'                is given more than once, the extra occurances are discarded. Consequently, you
#'                cannot be guaranteed that the number of rows in the resulting dataframe has the same
#'                length as the number of species that were given.
#'                
#'                Note that it also removes duplicate species names that are given. It only takes
#'                species names, NOT genus or family names. If genus or other taxonomic ranks are given,
#'                they are ignored.
#'
#' @return A dataframe with 12 columns, including 'species' (string) that includes all
#' synonyms of the binomial scientific name for each species supplied, and 'tsn'
#' (integer) which is the NCBI taxonomic idientifier (so two synonyms for the
#' same species will have the same single tsn), along with columns of other 
#' taxonomic info (kingdom, phylum etc)

#' @examples
#' find_synonyms( "Betta splendens" )
#' find_synonyms(c("Alectoris chukar", "Alectoris rufa", "Alle alle" , "Allactodipus bobrinskii" ))

find_synonyms <- function( species ) {
  
  require(taxizedb)
  
  # Get taxonomic information 
  
  sql_integer_list <- function(x){
    
    if (any(is.na(x))) {
      stop("Cannot pass NA into SQL query")
    }
    x <- as.character(x)
    if (!all(grepl('^[0-9]+$', x, perl = TRUE))) {
      stop("Found non-integer where integer required in SQL input")
    }
    paste(x, collapse = ", ")
  }
  
  taxizedb::db_download_ncbi()
  
  src_ncbi <- taxizedb::src_ncbi()
  
  species <- unique( species )
  
  # create the empty final results dataframe
  
  results <- data.frame(
    species = species,
    found = rep( FALSE, length( species ) ),
    tsn = rep( NA, length( species ) ),
    accepted_name = rep( NA, length( species) ),
    synonyms = rep( NA, length( species) ),
    common_name = rep( NA, length( species ) ),
    kingdom = rep( NA, length( species) ),
    phylum = rep( NA, length( species) ),
    class = rep( NA, length( species) ),
    order = rep( NA, length( species) ),
    family = rep( NA, length( species) ),
    genus = rep( NA, length( species) ),
    stringsAsFactors = FALSE
  )
  
  # Get the taxonomic IDs 
  
  # get NCBI IDs where available (NA if not available)
  
  species <- as.data.frame(species)
  species$species <- as.character((species$species))
  species_taxid <- sapply( species, taxizedb::name2taxid, out_type = "summary" )
  species_taxid <- as.data.frame(do.call(cbind, species_taxid))
  colnames(species_taxid) <- c("species", "tsn")
  
  # Add tsn ID to results dataframe
  
  results <- results %>%
    dplyr::select(-tsn) %>%
    merge(species_taxid, by = "species", all = TRUE) %>%
    dplyr::select(species, found, tsn, everything())
  
  # load accepted names into the results dataframe
  
  for( tsn in na.omit( unique( results$tsn ) ) ) {
    
    classification <- taxizedb::classification( tsn )[[tsn]]
    relevant_rows <- which( results$tsn == tsn )
    
    # set taxonomic information
    if( 
      is.data.frame( classification )
      && length( relevant_rows ) != 0
      && length(which(classification$rank == "kingdom")) != 0
      && length(which(classification$rank == "phylum")) != 0
      && length(which(classification$rank == "class")) != 0
      && length(which(classification$rank == "order")) != 0
      && length(which(classification$rank == "family")) != 0
      && length(which(classification$rank == "genus")) != 0
    ) {
      results[relevant_rows, "kingdom" ] <- classification[[which(classification$rank == "kingdom"), "name"]]
      results[relevant_rows, "phylum" ] <- classification[[which(classification$rank == "phylum"), "name"]]
      results[relevant_rows, "class" ] <- classification[[which(classification$rank == "class"), "name"]]
      results[relevant_rows, "order" ] <- classification[[which(classification$rank == "order"), "name"]]
      results[relevant_rows, "family" ] <- classification[[which(classification$rank == "family"), "name"]]
      results[relevant_rows, "genus" ] <- classification[[which(classification$rank == "genus"), "name"]]
    }
    
    # set accepted name
    if( 
      is.data.frame( classification )
      && length( relevant_rows ) != 0
      && length(which(classification$rank == "species")) != 0
    ) {
      results[relevant_rows, "accepted_name" ] <- classification[[which(classification$rank == "species"), "name"]]
      results[relevant_rows, "found" ] <- TRUE
    }
    
    # set common name
    common_names <- taxizedb::sql_collect(src_ncbi, paste0("SELECT * FROM names WHERE tax_id=", tsn, " AND name_class='common name'" ) )
    results[which(results$tsn == tsn), "common_name" ] <- paste0( common_names$name_txt, collapse = ", " )
    
    # set synonyms
    common_names <- taxizedb::sql_collect(src_ncbi, paste0("SELECT * FROM names WHERE tax_id=", tsn, " AND name_class='common name'" ) )
    results[which(results$tsn == tsn), "common_name" ] <- paste0( common_names$name_txt, collapse = ", " )
  }
  
  # Find synonyms
  
  # get the list of relevant ids
  relevant_tsns <- na.omit( unique( results$tsn ) )
  
  # make a list of all the names available and tsns to match
  query <- "SELECT * FROM names WHERE tax_id IN(%s) AND( name_class='scientific name' OR name_class='synonym')"
  query <- sprintf(query, sql_integer_list( relevant_tsns ))
  search_names <- taxizedb::sql_collect(src_ncbi, query)
  
  synonyms <- data.frame(
    tsn = search_names$tax_id,
    binomial = search_names$name_txt
  )
  
  synonyms$binomial <- as.character(synonyms$binomial)
  
  synonyms <- synonyms %>%
    merge(results, by = "tsn") %>%
    dplyr::select(-species) %>%
    dplyr::select(binomial, tsn, everything())
  
  # return( synonyms )
  return( synonyms )
  
}

# Standardise species databases ----

# Load tables from WildFinder database (converted into .xlsx files from .mdb 
# database)

#' TODO: Add ecoregion names to species_data (NB - this makes object size too large)
#' TODO: Find out what year wildfinder RL status is from
#' TODO: Figure out if there's a way to load the tables directly from the .mdb


# Load the wildfinder data

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

# Add common names, binomial scientific names, synonyms and tsn IDs

common_names <- tables$common_names
scientific_names <- as.data.frame(cbind(tables$redlist_species$wwf_species_id, 
                          tables$redlist_species$genus,
                          tables$redlist_species$species))
colnames(scientific_names) <- c("wwf_species_id", "genus", "species")

wildfinder_species <- scientific_names %>%
                      mutate(binomial = paste(genus, species, sep = " ")) %>%
                      dplyr::select(binomial)

wildfinder_species <- wildfinder_species[,1]

# Get synonyms

wildfinder_species <- find_synonyms( wildfinder_species)


wildfinder_database <-  ecoregions_species %>%
                        merge(common_names[c("species_id", "common_name")],
                              all = TRUE) %>%
                        rename(wwf_species_id = species_id) %>%
                        dplyr::select(c("wwf_species_id", "common_name",
                                        "ecoregion_code")) %>%
                        merge(scientific_names, by = "wwf_species_id",
                              all = TRUE) %>%
                        mutate(binomial = paste(genus, species, sep = " ")) %>%
                        mutate(source = "WildFinder") %>%
                        mutate(wwf_species_id = as.numeric(wwf_species_id)) %>%
                        merge(wildfinder_species[c("binomial", "tsn", "accepted_name")])


# Remove un-needed tables

# rm(tables)

# Clear memory by removing large objects

# rm(redlist_species, common_names, ecoregions, ecoregions_species)

# Load Cooke database - more species per ecoregion

cooke_database <- read.csv(paste(inputs,"/cooke_database/species_site.csv", sep = ""))

cooke_species <- cooke_database %>%
                 dplyr::select(binomial) 

cooke_species$binomial <- as.character(cooke_species$binomial)
cooke_species <- cooke_species[,1]

cooke_species <- find_synonyms(cooke_species)

## standardise columns to match species_df

cooke_database <- cooke_database %>%
                  dplyr::mutate(source = "Cooke et al 2019") %>%
                  dplyr::rename(ecoregion_code = eco) %>%
                  merge(cooke_species[c("binomial", "tsn", "accepted_name")])


# TEMPORARY CODE - subset data to make it manageable ----

#' TODO: Remove the two lines below and run on entire dataset when hpc available
#' 

## Save full databases

cooke_database_all <- cooke_database
wildfinder_database_all <- wildfinder_database

## Take a random sample subset

wildfinder_database <- wildfinder_database_all[sample(nrow(wildfinder_database_all), 10000), ]
cooke_database <- cooke_database_all[sample(nrow(cooke_database_all), 10000), ]



# Fill in details so Cooke has the same columns and data as wildfinder

cooke_database <- cooke_database %>%
                  merge(wildfinder_database[c("wwf_species_id","common_name","genus", 
                                              "species", "tsn")], by = "tsn", all.x = TRUE) %>%
                  dplyr::select("binomial", "wwf_species_id", "common_name", "ecoregion_code",
                                "genus", "species", "source", "tsn", "accepted_name" ) %>%
                  distinct(.)

# Load the red list data from Henriques et al 2020

file_path <- file.path(inputs, "redlist_sergio")
file_names <- list.files(file_path)
files <- file.path(file_path, file_names)
tables <- lapply(files, read_csv)
names(tables) <- str_remove(file_names, ".csv")

# Add a column so we can identify what taxa each table contains + the source
# reference

group_tables <- list()

for (i in seq_along(tables)) {
  
  group_tables[[i]] <- tables[[i]] %>% 
                       mutate(group = names(tables[i])) %>%
                       mutate(redlist_source = "Henriques etal 2020")
  
  
}

# Standardise column names etc. because they are inconsistent between the 
# different taxa tables

amphibians <- group_tables[[1]]

amphibians <- amphibians %>% 
              dplyr::mutate(binomial = paste(Genus, Species, sep = " ")) %>%
              dplyr::select(-c(Genus, Species)) %>%
              set_names(c("2004", "2008", "group", "redlist_source", "binomial")) %>%
              dplyr::select("binomial","group","2004", "2008", "redlist_source")

birds <- group_tables[[2]]

birds <- birds %>%
         set_names(c("binomial", "1988", "1994", "2000", "2004", "2008" , 
                     "2012" , "2016", "group", "redlist_source" )) %>%
         dplyr::select("binomial", "group", everything())

mammals <- group_tables[[5]]  

mammals <- mammals %>%
           set_names(c("binomial", "1996", "2008", "group", "redlist_source" )) %>%
           dplyr::select("binomial", "group", everything())

# Combine back into a list then dataframe

group_tables <- list(birds, mammals, amphibians)

#' TODO: Change name of species_redlist to henrique database?

species_redlist <- do.call(bind_rows, group_tables)

# Get the synonyms and taxonomic ids

redlist_species <- species_redlist %>%
                   dplyr::select(binomial) 

redlist_species$binomial <- as.character(redlist_species$binomial)
redlist_species <- unname(unlist(redlist_species[,1]))
redlist_species <- redlist_species[!is.na(redlist_species)]

redlist_species <- find_synonyms(redlist_species)

# Order columns, drop group variable and melt in to long format

species_redlist <- species_redlist %>%
                   dplyr::select("binomial", "1988", "1994", "1996", "2000", 
                                "2004", "2008", "2012", "2016", "redlist_source") %>%
                   melt(.,id.vars = c("binomial", "redlist_source"),
                        value.name = "redlist_status") %>%
                   rename(redlist_assessment_year = variable) %>% 
                   merge(redlist_species[c("binomial", "tsn", "accepted_name", 
                                           "common_name")])


# Merge databases ----

# Combine Cooke and Wildfinder

merged_databases <- rbind(wildfinder_database, cooke_database)

# Add red list status and year to the merged databases.

species_data_with_sources <- merged_databases %>%
                             merge(species_redlist, by = "tsn", all = TRUE) %>%
                             mutate(source = coalesce(source, redlist_source)) %>%
                             mutate(binomial = coalesce(binomial.x, binomial.y)) %>%
                             mutate(accepted_name = coalesce(accepted_name.x, 
                                                             accepted_name.y)) %>%
                             mutate(common_name = coalesce(common_name.x, 
                                                           common_name.y)) %>%
                             dplyr::select(-c(binomial.x, binomial.y, 
                                              accepted_name.x, accepted_name.y, 
                                               common_name.x, common_name.y))

# Consolidate duplicates ----


species_data <- species_data_with_sources %>%
                dplyr::select(-c("source", "redlist_source", "binomial")) %>%
                distinct(.) %>%
                dplyr::select(tsn, accepted_name, everything()) %>%
                rename(accepted_binomial = accepted_name)



# TEMPORARY CODE - checkpoint - save processed data ----

## Save the species_data while we are working on it

# write_csv(species_data, paste(outputs, "draft_species_data.csv", sep = "/"))
# saveRDS(species_data, paste(outputs, "draft_species_data.rds", sep = "/"))

# species_data <- readRDS(file.path(outputs, 'draft_species_data.rds'))


# (SLOW CODE) Get ecoregions for species missing them ----

#' TODO: Check this is working correctly - how does it cope with one species
#' in multiple ecoregions?

# Get the ecoregion map & subset to required variables

ecoregion_map_all <- st_read(paste(inputs,"official_teow_wwf", sep = "/"))
ecoregion_map <- ecoregion_map_all %>% dplyr::select(eco_code, ECO_NAME, geometry)

## TEMPORARY CODE - just using extinct ranges for now as laptop can't manage
## all species maps yet

species_ecoregions <- get_ecoregions("redlist_extinct_species_range_maps", 
                                     ecoregion_map)

species_ecoregions <- list(species_ecoregions)

# Use this code when on server, doesn't work on laptop

# range_directories <- c("redlist_extinct_species_range_maps",
#                        "redlist_reptile_range_maps", 
#                        "redlist_amphibian_range_maps",
#                        "redlist_mammal_range_maps" )
# 
# species_ecoregions <- list()
# 
# for (i in seq_along(range_directories)) {
#   
#   species_ecoregions[[i]] <- get_ecoregions(range_directories[i], em_small)
#   
# }
#
# species_ecoregions <- do.call(rbind, species_ecoregions)


# Get TSNs for the species

species_ecoregions_names <- species_ecoregions %>%
                            dplyr::select(binomial) 

species_ecoregions_names$binomial <- as.character(species_ecoregions_names$binomial)

species_ecoregions_names <- unname(unlist(species_ecoregions_names[,1]))

species_ecoregions_names <- species_ecoregions_names[!is.na(species_ecoregions_names)]

species_ecoregions_names <- find_synonyms(species_ecoregions_names)

species_ecoregions <- species_ecoregions %>%
                      merge(species_ecoregions_names[c("binomial", "tsn", 
                                                       "found")], 
                            by = "binomial")

## Add the new ecoregions to our species data

species_data <- species_data %>%
                merge(species_ecoregions[c( "tsn", "ecoregion_code")], 
                      by = "tsn", all = TRUE) %>%
                mutate(ecoregion_code = coalesce(ecoregion_code.x, 
                                                 ecoregion_code.y)) %>%
                dplyr::select(-c(ecoregion_code.x, ecoregion_code.y))

# TEMPORARY CODE - Run diagnostics on gaps in our database

# What species do we not have ecoregions for?

species_without_ecoregions <- species_data %>%
  filter(is.na(ecoregion_code)) %>%
  dplyr::select(tsn, accepted_binomial) %>%
  distinct(.)

# What species do we not have a redlist status in any year for?

species_without_redlist_status <- species_data %>%
  dplyr::select(-redlist_assessment_year, 
                genus, species) %>%
  distinct(.) %>%
  group_by(tsn) %>%
  filter(all(is.na(redlist_status)))

# What species do we have at least one redlist status for?

species_with_redlist_status <- species_data %>%
  dplyr::select(-redlist_assessment_year, 
                genus, species) %>%
  distinct(.) %>%
  group_by(tsn) %>%
  filter(!is.na(redlist_status)) %>%
  dplyr::select(-redlist_status) %>%
  distinct(.)

# Double check we haven't incorrectly detected no redlist status (should be no
# overlap between the two groups)

overlap <- species_without_redlist_status$tsn %in% species_with_redlist_status$tsn
any(overlap == TRUE) # The correct output to console should be FALSE (no overlap)

# Calculate summary statistics ----

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


# Visualise summary stats ----

extinction_map_data <- inner_join(em_small, proportion_extinct[
  c("merged_ecoregion_code", "number_of_species_extinct")], 
                             by = c("eco_code" = "merged_ecoregion_code"))

extinction_map <- ggplot(extinction_map_data) +
                  geom_sf(aes(fill = number_of_species_extinct)) +
                  scale_fill_viridis_c(trans = "sqrt", alpha = .4)

extinction_map

