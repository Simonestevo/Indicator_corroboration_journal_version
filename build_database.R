# https://stackoverflow.com/questions/13070706/how-to-connect-r-with-access-database-in-64-bit-window

# cd C:\\Users\\ssteven\\Dropbox\\Deakin\\Chapter_2_Extinction_test\\Extinction_test_code

# Using R version 3.5.2

## This script outputs:
## 1: A database of the WWF terrestrial ecoregions, species found
## in each ecoregion (this will not be a comprehensive list of every single species
## but aims to be as complete as possible), and their red list status, and
## 2. Summary statistics of number of species and number of extinct species

#' TODO: Remove objects when they are no longer needed
#' TODO: Important - use packrat or something to save version of taxizedb
#'
# install.packages("raster") 
# install.packages("tidyverse")
# install.packages("reshape2")
# install.packages("sf")
# install.packages("spData")
# install.packages("devtools")
# devtools::install_github("ropensci/taxizedb") # version 0.1.9.9130
# install.packages("rredlist")
# install.packages("rgbif")
# install.packages("rlist")
# install.packages("pryr")

# Load packages ----

library(raster)
library(tidyverse)
library(stringr)
library(stringi)
library(reshape2)
library(sf)
library(spData)
library(taxizedb)
library(functionaltraits)
library(broom)
library(rredlist)
library(rgbif)
library(rlist)
library(pryr) # Can probably remove this when finished
library(data.table)

# Set input and output locations ----

create_new_database_version <- FALSE # Only set to true if you want to create an entirely new version from scratch
date <- Sys.Date()
country <- NA #"Australia" # If not subsetting, set as NA, e.g. country <- NA
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


if (create_new_database_version == FALSE) {

# Find the most recent version of database output files
  
db_version <- tail(sort(list.files(parent_outputs)), 1)
db_interim <- list.files(file.path(parent_outputs,db_version))[
         grepl("interim",list.files(file.path(parent_outputs,db_version)))]
db_outputs <- list.files(file.path(parent_outputs,db_version))[
  grepl("database",list.files(file.path(parent_outputs,db_version)))]

interim_outputs <- file.path(parent_outputs, db_version, db_interim)
outputs <- file.path(parent_outputs, db_version, db_outputs)

} else if (create_new_database_version == TRUE) { # Dev_mode creates brand new folders for your outputs

previous_version <- tail(sort(list.files(parent_outputs)), 1)
new_version <- as.numeric(substring(previous_version, nchar(previous_version))) + 1
db_version <- paste("version", new_version, sep = "_")
new_int_dir <- file.path(parent_outputs, db_version,
                        paste(date,"_interim_output_files",sep = ""))

if( !dir.exists( file.path(new_int_dir) ) ) {
  
  dir.create( file.path(new_int_dir), recursive = TRUE )
  
  print("creating new interim directory")
  
}

new_db_dir <- file.path(parent_outputs, db_version,
                         paste(date,"_database_output_files",sep = ""))

if (!dir.exists( file.path(new_db_dir))) {
  
  dir.create( file.path(new_db_dir), recursive = TRUE )
  
  print("creating new database directory")
  
  }
}


# Functions ----

## get_ecoregions: Function to load species range maps and join with the wwf ecoregion map
#' @return a dataframe with four columns (species binomial scientific name,
#' ecoregion, data source and red list status)
#' @param rangemap_directory_path - a string that denotes the name of the directory
#' where the group (eg mammals) range maps are saved (should be one directory per map)
#' @param map sf object - the map of ecoregions, will work with either 2001 or
#' 2017 version

# map <- ecoregion_map
# rangemap_directory_path <- range_directories[[1]]

get_ecoregions <- function(range_map, rangemap_directory_path, map, 
                           interim_outputs, location, class_name) {
  
  names(range_map) <- c(toupper(names(range_map)[1:3]), names(range_map[4])) # Make column names consistent
  
  #range_map <- rename(range_map, geometry = GEOMETRY) # Convert geometry back tho
  
  # Get the ecoregions for each species
  
  ranges_ecoregions <- st_join(range_map, map, join = st_intersects)
  
  if (!is.null(ranges_ecoregions)) {
    
    print(paste("st_join complete for", 
                basename(rangemap_directory_path), sep = " "))
    
  }
  
  # Standardise the data so we can add it to the species data easily
  
  if (class_name == "bird") {
    
    species_w_ecoregions <- as.data.frame(ranges_ecoregions %>%
                                            dplyr::select(SCINAME, 
                                                          ECO_ID)) %>%
      dplyr::select(-Shape) %>%
      rename(binomial = SCINAME) %>%
      rename(ecoregion_id = ECO_ID) %>%
      mutate(eco_objectid = NA) %>%
      mutate(source = "birdlife_international") %>%
      mutate(redlist_status = "TBC") 
    
  } else {
  
  species_w_ecoregions <- as.data.frame(ranges_ecoregions %>%
                                          dplyr::select(BINOMIAL, 
                                                        ECO_ID, 
                                                        OBJECTID)) %>%
    dplyr::select(-geometry) %>%
    rename(binomial = BINOMIAL) %>%
    rename(ecoregion_id = ECO_ID) %>%
    rename(eco_objectid = OBJECTID) %>%
    mutate(source = "iucn_redlist") %>%
    mutate(redlist_status = "TBC") 

  }
  
  if(!is.null(species_w_ecoregions)) {
    
    print(paste("finished adding ecoregions to", 
                basename(rangemap_directory_path), sep = " "))
    
  }
  
  saveRDS(species_w_ecoregions, file.path(interim_outputs, 
                                          paste(location,class_name, 
                                                "ecoregions.rds", sep = "_")))
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
#' TODO: Check if need to download other species too

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

# Function to scrape data about species redlist from the IUCN website

get_redlist_data <- function(species) {
  
  library("jsonlite")
  
  output <- jsonlite::fromJSON(rl_search_(species))
  
  print(paste("Data for", species, "found", sep = " "))
  
  return(output)
}


#' This function finds locations of species occurrences and returns the ecoregion
#' they were observed in.
#' 
#' @param species A vector of strings of species accepted binomial scientific
#' name
#' 
#' @param observation An integer limiting how many observations you want to
#' pull from gbif.  Tradeoff - lower values will improve processing requirements
#' but reduce the likelihood of getting the species full range.  Larger values
#' will take longer but improve range coverage. You can begin with a small number
#' then run again for all species not found with a larger observation number, 
#' or to improve accuracy.
#' 
#' @param polygon_map An sf polygon map of terrestrial ecoregions. Not sure
#' if this would work with other polygons but probably would.
#' 
#' @return A list of two dataframes - 'extinct_coordinates' is a dataframe of 
#' species names for which coordinates were found and all the ecoregions they 
#' are found in, based on the observations you pulled. And 'no_coordinates',
#' species for which no coordinates were found.  You can then re-run the function
#' on these species with a higher observation number.
#' 

# species <- extinct_species_names
# observations <- 50
# polygon_map <- ecoregion_map

get_gbif_data <- function(species, observations, polygon_map) {
  
  require(rgbif)
  require(rlist)
  require(sf)
  require(tidyverse)
  
  #keys <- sapply(species, function(x) name_suggest(x)$key[1], USE.NAMES = FALSE)
  
  gbif_keys <- sapply(species, name_suggest, USE.NAMES = FALSE)
  
  keys <- sapply(gbif_keys, function(x) x$key[1])
  
  keys <- list.clean(keys, fun = is.null, recursive = FALSE)
  
  extinct_gbif_data <- occ_search(taxonKey = keys, limit = observations)
  
  extinct_coordinates <- list()
  no_coordinates <- list()
  
  for (i in seq_along(extinct_gbif_data)) {
    
    single_spp <- extinct_gbif_data[[i]][[3]]
    
    if ( "decimalLatitude" %in% names(single_spp) == TRUE) {
      
      print(paste("co-ordinates found for ", single_spp$species[1] ,", ", 
                  i, " of ", length(extinct_gbif_data)," ", "species", sep = ""))
     
      extinct_coordinates[[i]] <- single_spp %>%
        select(genericName, specificEpithet, decimalLatitude, decimalLongitude ) %>%
        distinct(.) %>%
        filter(complete.cases(decimalLatitude, decimalLongitude)) %>%
        mutate(species = paste(genericName, specificEpithet, sep = " ")) %>%
        select(- genericName, - specificEpithet)
      
      } else {
      
      no_coordinates[[i]] <- single_spp$species[1]
      
      print(paste("no co-ordinates found for", single_spp$species[1], "among", 
                  observations, "observations", "try increasing observation argument to a larger number"), sep = " ")
      
    }
  }
  
  extinct_coordinates <- list.clean(extinct_coordinates, 
                                    fun = is.null, 
                                    recursive = FALSE)
  extinct_coordinates <- do.call(rbind, extinct_coordinates)
  
  extinct_sf <- st_as_sf(extinct_coordinates, coords = c('decimalLongitude', 
                                                         'decimalLatitude'), 
                         crs = st_crs(polygon_map))
  
  extinct_ecoregions <- st_intersection(extinct_sf, polygon_map)
  
  extinct_ecoregions <- st_drop_geometry(extinct_ecoregions) 
  
  no_coordinates <- do.call(rbind,no_coordinates)
  
  output <- list(extinct_ecoregions, no_coordinates)
  
  names(output) <- c("extinct_ecoregions","no_coordinates")
  
  return(output)
  
}


scale_to_1 <- function(vector){
  
  (vector - min(vector, na.rm = TRUE))/
    (max(vector, na.rm = TRUE) - min(vector, na.rm = TRUE))
}

get_binomial_list <- function(data) {
  
  data$binomial <- as.character(data$binomial)
  
  data <- unname(unlist(data[,1]))
  
  data <- data[!is.na(data)]
  
  binomial_list <- unique(data)
  
}

library(futile.logger)
library(utils)

retry <- function(expr, isError=function(x) "try-error" %in% class(x), 
                  maxErrors=5, sleep=0) {
  attempts = 0
  retval = try(eval(expr))
  while (isError(retval)) {
    attempts = attempts + 1
    if (attempts >= maxErrors) {
      msg = sprintf("retry: too many retries [[%s]]", capture.output(str(retval)))
      flog.fatal(msg)
      stop(msg)
    } else {
      msg = sprintf("retry: error in attempt %i/%i [[%s]]", attempts, maxErrors, 
                    capture.output(str(retval)))
      flog.error(msg)
      warning(msg)
    }
    if (sleep > 0) Sys.sleep(sleep)
    retval = try(eval(expr))
  }
  return(retval)
}


# Function to get redlist history, because iucn limit how many calls you can
# make, so have to break the full list of reptiles up into chunks

get_redlist_history <- function(names, class_name) {
  
  redlist_history_list <- list()
  
  for (i in seq_along(names)) {
    
    species_name <- names[i]
    
    species_redlist_history <- rl_history(species_name, parse = TRUE)[[2]] 
    
    Sys.sleep(4) # Make loop pause before next, otherwise iucn web access cuts out
    
    if (length(species_redlist_history) == 0) { # if there are no results create dummy dataframe
      
      species_redlist_history <- data.frame(year = NA, code = NA, binomial = species_name,
                                            category = NA,
                                            class = class_name)
      
      print(paste("red list history not available for", species_name, sep = " "))
      
    } else {# otherwise create results dataframe
      
      species_redlist_history <- species_redlist_history %>%
        mutate(binomial = species_name,
               class = class_name)
      
      print(paste("red list history retrieved for", species_name, sep = " "))
      
    }
    
    redlist_history_list[[i]] <- species_redlist_history
    
  }
  
  return(redlist_history_list)
  
}



summarise_species_data <- function(data, number, class_name) {
    
    out1 <- data %>%
      group_by(ecoregion_id, redlist_status) %>%
      summarise(spp_number_w_status = n_distinct())
    
    out2 <- data %>%
      group_by(ecoregion_id) %>%
      summarise(number_spp_in_ecoregion =
                  n_distinct(tsn))
    
    class_ecoregion_summary <- out1 %>%
      merge(out2, by = "ecoregion_id") %>%
      mutate(class = data$class[1])
    
    class_redlist_summary <- data %>%
      group_by(redlist_status) %>%
      summarise(redlist_count = n_distinct(tsn)) %>%
      mutate(class = data$class[1]) %>%
      mutate(total_spp_in_class = length(unique(data$tsn)))
    
    out <- list(class_ecoregion_summary, class_redlist_summary)
    
    saveRDS(class_ecoregion_summary, 
            file.path(interim_outputs, 
                      paste(location, class_name, number,
                            "redlist_category_breakdown_by_ecoregion.rds",
                            sep = "_")))
    
    saveRDS(class_redlist_summary, 
            file.path(interim_outputs, 
                      paste(location, class_name, number,
                            "redlist_category_breakdown_global.rds",
                            sep = "_")))
    
    return(out)
    
}

# Load ecoregion data ----

## This will allow us to subset the data by country for development or analysis
#' TODO: Why are so many ecoregions missing countries? figure out best join method
#' TODO: Important - Add countries to species_data 

# Get the ecoregion map 

ecoregion_map_all <- st_read(paste(inputs,eco_version, sep = "/"))

# Pull out only required variables

ecoregion_map <- ecoregion_map_all %>% 
                 dplyr::select(ECO_ID, ECO_NAME, geometry, OBJECTID)

ecoregion_map <- st_make_valid(ecoregion_map)

#rm(ecoregion_map_all)

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


ecoregion_map_simple <- st_simplify(ecoregion_map, preserveTopology = TRUE,
                                    dTolerance = 0.1)


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

# Assign species to ecoregions using range maps (SLOW CODE) ----

# Process by each class of rangemaps.  

# * Amphibians ----

## elapsed time for below 9559.55 s (~ 3 hrs)

if ((paste(location, "amphibian", "ecoregions.rds", 
            sep = "_") %in% list.files(interim_outputs))) {

amphibian_ecoregions <- readRDS(file.path(interim_outputs, 
                                          paste(location, 
                                                "amphibian", "ecoregions.rds", 
                                                sep = "_" )))
} else {

amphibian_rangemap_dir <- file.path(inputs, "redlist_amphibian_range_maps")

# Read in the rangemap

amphibian_ranges <- st_read(amphibian_rangemap_dir)

# Simplify geometry slightly to improve processing time

amphibian_ranges_simple <- st_simplify(amphibian_ranges, 
                                       preserveTopology = TRUE,
                           dTolerance = 0.1)

rm(amphibian_ranges)

# Remove unneccessary columns as well

amphibian_ranges_simple <- amphibian_ranges_simple %>%
                           select(id_no, binomial, presence, geometry)


# Match amphibian ranges to ecoregions

amphibian_ecoregions <- get_ecoregions(amphibian_ranges_simple, 
                                       amphibian_rangemap_dir, 
                                       ecoregion_map_simple,
                                       interim_outputs,
                                       location,
                                       "amphibian")
rm(amphibian_ranges_simple)

}

# * Mammals ----

## elapsed time for below 14264.58 s (~ 4 hrs)

if ((paste(location, "mammal", "ecoregions.rds", 
           sep = "_") %in% list.files(interim_outputs))) {
  
  mammal_ecoregions <- readRDS(file.path(interim_outputs, 
                                            paste(location, 
                                                  "mammal", "ecoregions.rds", 
                                                  sep = "_" )))
} else {

mammal_rangemap_dir <- file.path(inputs, "redlist_mammal_range_maps")

# Read in the rangemap

mammal_ranges <- st_read(mammal_rangemap_dir)

# Simplify geometry slightly to improve processing time

mammal_ranges_simple <- st_simplify(mammal_ranges, 
                                    preserveTopology = TRUE,
                                    dTolerance = 0.1)

rm(mammal_ranges)

# Remove unneccessary columns as well

mammal_ranges_simple <- mammal_ranges_simple %>%
                        select(id_no, 
                               binomial, 
                               presence, 
                               geometry)


# Match mammal ranges to ecoregions

mammal_ecoregions <- get_ecoregions(mammal_ranges_simple, 
                                                mammal_rangemap_dir, 
                                                ecoregion_map_simple,
                                                interim_outputs,
                                                location,
                                                "mammal")

rm(mammal_ranges_simple)

}

# * Birds ----

## elapsed time for below 58943.19 s (~16 hours)

## Note that bird maps come from birdlife international, not iucn, so
## they are stored in a geodatabase with slightly different geometry types, so
## require a couple of extra steps to process

if ((paste(location, "bird", "ecoregions.rds", 
           sep = "_") %in% list.files(interim_outputs))) {
  
bird_ecoregions <- readRDS(file.path(interim_outputs, 
                                     paste(location, 
                                     "bird", "ecoregions.rds", 
                                     sep = "_" )))
} else {

bird_rangemap_dir <- file.path(inputs, "birdlife_avian_range_maps","BOTW.gdb")

# Read in the rangemap

bird_ranges <- st_read(bird_rangemap_dir, layer = "All_Species")

# Remove unneccessary columns 

bird_ranges <- bird_ranges %>%
               select(SISID,
                      SCINAME, 
                      PRESENCE, 
                      Shape)

# bird_ranges contains both multipolygon and multisurface
# geometry types, which means other st functions won't work.

# Get only multisurface geoms
#' TODO: Figure out how to recast these to multipolygons

bird_ranges_ms <- bird_ranges %>%
                  filter(st_geometry_type(Shape) == "MULTISURFACE")

# Get only multipolygon geoms

bird_ranges <- bird_ranges %>%
               filter(st_geometry_type(Shape) == "MULTIPOLYGON")


# Simplify geometry slightly to improve processing time

bird_ranges_simple <- st_simplify(bird_ranges, 
                                  preserveTopology = TRUE,
                                  dTolerance = 0.1)

rm(bird_ranges)

# saveRDS(bird_ranges_simple, file.path(interim_outputs, "all_bird_ranges_simplified.rds"))
# bird_ranges_simple <- readRDS(file.path(interim_outputs, "all_bird_ranges_simplified.rds"))

# Match bird ranges to ecoregions

bird_ecoregions <- get_ecoregions(bird_ranges_simple, 
                                              bird_rangemap_dir, 
                                                 ecoregion_map_simple,
                                                 interim_outputs,
                                                 location,
                                                 "bird")
rm(bird_ranges_simple)

}

# * Reptiles ----

## elapsed time for below 15089.72  s (~ 4 hrs)

if ((paste(location, "reptile", "ecoregions.rds", 
           sep = "_") %in% list.files(interim_outputs))) {
  
  reptile_ecoregions <- readRDS(file.path(interim_outputs, 
                                          paste(location, 
                                                "reptile", "ecoregions.rds", 
                                                sep = "_" )))
} else {
  
  reptile_rangemap_dir <- file.path(inputs, "redlist_reptile_range_maps")
  
  # Read in the rangemap
  
  reptile_ranges <- st_read(reptile_rangemap_dir)
  
  # Subset to terrestrial only
  
  reptile_ranges <- reptile_ranges %>%
    filter(terrestial == "true") #spelling error is in the sf column
  
  # Simplify geometry slightly to improve processing time
  
  reptile_ranges_simple <- st_simplify(reptile_ranges, 
                                       preserveTopology = TRUE,
                                       dTolerance = 0.1)
  
  rm(reptile_ranges)
  
  # Remove unneccessary columns as well
  
  reptile_ranges_simple <- reptile_ranges_simple %>%
    select(id_no, 
           binomial, 
           presence, 
           geometry)
  
  
  # Match reptile ranges to ecoregions
  
  reptile_ecoregions <- get_ecoregions(reptile_ranges_simple, 
                                       reptile_rangemap_dir, 
                                       ecoregion_map_simple,
                                       interim_outputs,
                                       location,
                                       "reptile")
  rm(reptile_ranges_simple)
  
  # Add in the point data
  
  reptile_points <- read.csv(file.path(inputs, "redlist_reptile_range_maps",
                                       "REPTILES_points.csv"))
  # Select necessary columns
  
  reptile_points <- reptile_points %>%
    dplyr::select(binomial, latitude, longitude, category) 
  
  # Convert to sf object and set crs
  
  reptile_points_sf <- st_as_sf(reptile_points, coords = c('longitude', 
                                                           'latitude'), 
                                crs = st_crs(ecoregion_map_simple))
  
  # Get ecoregions the points fall within
  
  reptile_point_ecoregions <- st_intersection(reptile_points_sf, 
                                              ecoregion_map_simple)
  # Format to match polygon outputs
  
  reptile_point_ecoregions <- reptile_point_ecoregions %>%
    mutate(source = "iucn_redlist_point_data") %>%
    dplyr::select(binomial, ECO_ID, 
                  OBJECTID, source, category) %>%
    rename(ecoregion_id = ECO_ID,
           eco_objectid = OBJECTID,
           redlist_status = category) %>%
    st_drop_geometry(.)
  
  # Bind the two
  
  reptile_ecoregions <- rbind(reptile_ecoregions, reptile_point_ecoregions)
  
  # Save new version
  
  saveRDS(reptile_ecoregions, file.path(interim_outputs, 
                                        paste(location,"reptile", 
                                              "ecoregions.rds", sep = "_")))
  
}


# Get rangemap synonyms ----

# * Amphibians  ----

if ((paste(location,"amphibian", 
                 "rangemap_synonyms.rds", sep = "_") %in% 
     list.files(interim_outputs))) {

amphibian_rangemap_binomials <- get_binomial_list(amphibian_ecoregions)
  
amphibian_rangemap_synonyms <- readRDS(file.path(interim_outputs, 
                                          paste(location,"amphibian", 
                                                "rangemap_synonyms.rds", sep = "_")))
} else {

amphibian_rangemap_binomials <- get_binomial_list(amphibian_ecoregions)

system.time(amphibian_rangemap_synonyms <- find_synonyms(amphibian_rangemap_binomials))

saveRDS(amphibian_rangemap_synonyms, file.path(interim_outputs, 
                                      paste(location,"amphibian", 
                                            "rangemap_synonyms.rds", sep = "_")))
}

# Join the synonyms and tsn to the rangemap data

amphibian_ecoregions <- amphibian_ecoregions %>%
                        merge(amphibian_rangemap_synonyms[c("tsn", "binomial", 
                                                            "accepted_name")],
                        by = "binomial") %>%
                        dplyr::select(binomial, tsn, ecoregion_id, 
                                      eco_objectid, source, redlist_status)

# * Mammals  ----

if ((paste(location,"mammal", 
           "rangemap_synonyms.rds", sep = "_") %in% 
     list.files(interim_outputs))) {
  
  mammal_rangemap_binomials <- get_binomial_list(mammal_ecoregions)
  
  mammal_rangemap_synonyms <- readRDS(file.path(interim_outputs, 
                                                   paste(location,"mammal", 
                                                         "rangemap_synonyms.rds", sep = "_")))
} else {

mammal_rangemap_binomials <- get_binomial_list(mammal_ecoregions)

system.time(mammal_rangemap_synonyms <- find_synonyms(mammal_rangemap_binomials))

saveRDS(mammal_rangemap_synonyms, file.path(interim_outputs, 
                                      paste(location,"mammal", 
                                            "rangemap_synonyms.rds", sep = "_")))
}

# Join the synonyms and tsn to the rangemap data

mammal_ecoregions <- mammal_ecoregions %>%
                     merge(mammal_rangemap_synonyms[c("tsn", "binomial", 
                                                    "accepted_name")],
                            by = "binomial") %>%
                     dplyr::select(binomial, tsn, ecoregion_id, 
                                   eco_objectid, source, redlist_status)

# * Reptiles ----

if ((paste(location,"reptile", 
           "rangemap_synonyms.rds", sep = "_") %in% 
     list.files(interim_outputs))) {
  
  reptile_rangemap_binomials <- get_binomial_list(reptile_ecoregions)
  
  reptile_rangemap_synonyms <- readRDS(file.path(interim_outputs, 
                                                   paste(location,"reptile", 
                                                         "rangemap_synonyms.rds", sep = "_")))
} else {

reptile_rangemap_binomials <- get_binomial_list(reptile_ecoregions)

system.time(reptile_rangemap_synonyms <- find_synonyms(reptile_rangemap_binomials))

saveRDS(reptile_rangemap_synonyms, file.path(interim_outputs, 
                                      paste(location,"reptile", 
                                            "rangemap_synonyms.rds", sep = "_")))
}
# Join the synonyms and tsn to the rangemap data

reptile_ecoregions <- reptile_ecoregions %>%
                      merge(reptile_rangemap_synonyms[c("tsn", "binomial", 
                                                            "accepted_name")],
                              by = "binomial") %>%
                      dplyr::select(binomial, tsn, ecoregion_id, 
                                      eco_objectid, source, redlist_status)

# * Birds ----

if ((paste(location,"bird", 
           "rangemap_synonyms.rds", sep = "_") %in% 
     list.files(interim_outputs))) {
  
  bird_rangemap_binomials <- get_binomial_list(bird_ecoregions)
  
  bird_rangemap_synonyms <- readRDS(file.path(interim_outputs, 
                                                   paste(location,"bird", 
                                                         "rangemap_synonyms.rds", sep = "_")))
} else {

bird_rangemap_binomials <- get_binomial_list(bird_ecoregions)

system.time(bird_rangemap_synonyms <- find_synonyms(bird_rangemap_binomials))

saveRDS(bird_rangemap_synonyms, file.path(interim_outputs, 
                                      paste(location,"bird", 
                                            "rangemap_synonyms.rds", sep = "_")))
}

# Join the synonyms and tsn to the rangemap data

bird_ecoregions <- bird_ecoregions %>%
                   merge(bird_rangemap_synonyms[c("tsn", "binomial", 
                                                        "accepted_name")],
                         by = "binomial", all = TRUE) %>%
                   dplyr::select(binomial, tsn, ecoregion_id, 
                                 eco_objectid, source, redlist_status)

# Get red list status' ----

# * Amphibians ----

# Read in Henriques data which gives history of species' red list status

if ((paste(location, "amphibian", "ecoregion_redlist.rds", 
           sep = "_") %in% list.files(interim_outputs))) {
  
amphibian_ecoregion_redlist <- readRDS(file.path(interim_outputs, 
                                            paste(location, 
                                                  "amphibian", "ecoregion_redlist.rds", 
                                                  sep = "_" )))
} else {
  
amphibian_redlist_data <- read_csv(file.path(inputs, 
                          "henriques_redlist_history",
                           "RLTS_amphibian_data_organised.csv"))

# Format and melt

amphibian_redlist_data <- amphibian_redlist_data %>% 
                          mutate(class = "amphibia") %>%
                          mutate(redlist_source = "Henriques etal 2020") %>% 
                          dplyr::mutate(binomial = paste(Genus, Species, 
                                                         sep = " ")) %>%
                          dplyr::select(-c(Genus, Species)) %>%
                          set_names(c("2004", "1980", "class", "redlist_source", 
                                      "binomial")) %>%
                          dplyr::select("binomial","class","2004", "1980", 
                                        "redlist_source") %>%
                          reshape2::melt(.,id.vars = c("binomial", "class", 
                                             "redlist_source"), 
                               value.name = "redlist_status",
                               variable.name = "redlist_assessment_year") 
                          

# Get the synonyms for this data source

if ((paste(location, "amphibian", "redlist_synonyms.rds", 
           sep = "_") %in% list.files(interim_outputs))) {

amphibian_redlist_synonyms <- readRDS(file.path(interim_outputs,
                                        paste(location, "amphibian", "redlist_synonyms.rds", 
                                        sep = "_")))

amphibian_redlist_binomials <- get_binomial_list(amphibian_redlist_data)
  
} else {
  
amphibian_redlist_binomials <- get_binomial_list(amphibian_redlist_data)

system.time(amphibian_redlist_synonyms <- find_synonyms(amphibian_redlist_binomials))

saveRDS(amphibian_redlist_synonyms, file.path(interim_outputs, 
                                              paste(location,"amphibian", 
                                                    "redlist_synonyms.rds", 
                                                    sep = "_")))
}
# Join the synonyms and tsn to the redlist data

# Because we want to delete any false NAs (where a species actually does have a 
# redlist status, but because it has several synonyms it shows up more than once
# sometimes with redlist status, and also without). We want to remove those NA
# values because they aren't true NAs and will inflate the NA count and send you
# on a goose chase looking for the redlist data that you actually already have

amphibians_multi_tsn <- amphibian_redlist_data %>%
                        merge(amphibian_redlist_synonyms[c("tsn", "binomial", 
                                                           "accepted_name")],
                                                    by = "binomial",
                                                    all = TRUE) %>%
                        dplyr::select(binomial, tsn, redlist_assessment_year, 
                                      redlist_status, redlist_source) %>%
                        distinct(tsn, redlist_assessment_year, redlist_status,  
                                  redlist_source, .keep_all = TRUE) %>%
                                  group_by(tsn) %>%
                                  filter(n() > 1) %>%
                                  ungroup(.) %>%
                        filter(!is.na(redlist_status)) 

# But if there's only one name for that species and it is listed as NA, then
# it's truly NA and we want to keep it.  Annoying I can't work out how to 
# satisfy both conditions in one pipe

amphibians_single_tsn <- amphibian_redlist_data %>%
                         merge(amphibian_redlist_synonyms[c("tsn", "binomial", 
                                                            "accepted_name")],
                                by = "binomial",
                                all = TRUE) %>%
                         dplyr::select(binomial, tsn, redlist_assessment_year, 
                                        redlist_status, redlist_source) %>%
                         distinct(tsn, redlist_assessment_year, redlist_status,  
                                   redlist_source, .keep_all = TRUE) %>%
                         group_by(tsn) %>%
                         filter(n() == 1) %>%
                         ungroup(.)

amphibian_redlist_data <- as.data.frame(rbind(amphibians_single_tsn, 
                                              amphibians_multi_tsn))

# Join ecoregion and redlist data together

amphibian_ecoregion_redlist <- amphibian_ecoregions %>%
  dplyr::select(-redlist_status) %>%
  merge(amphibian_redlist_data[c("tsn", 
                              "binomial",
                              "redlist_assessment_year",
                              "redlist_status", 
                              "redlist_source")], by = "tsn", 
        all = TRUE) %>%
  mutate(binomial = coalesce(binomial.x, binomial.y),
         class = "Amphibia") %>%
  select(ecoregion_id, tsn, binomial, source, 
         redlist_assessment_year,
         redlist_status, redlist_source, class)

saveRDS(amphibian_ecoregion_redlist, file.path(interim_outputs, 
                                               paste(location,"amphibian", 
                                                     "ecoregion_redlist.rds", 
                                                     sep = "_")))
}


# Summarise to check for data gaps

amphibian_summaries <- summarise_species_data(amphibian_ecoregion_redlist, "1", 
                                              "amphibian")
amphibian_redlist_by_ecoregion <- amphibian_summaries[[1]]
amphibian_redlist_global <- amphibian_summaries[[2]]

# Check data

# Check the data

test <- amphibian_redlist_data %>% filter(redlist_status == "EX") # only one ex mammal?
newtest <- amphibian_ecoregion_redlist %>% filter(binomial == "Atelopus vogli")


# Check what species we still don't have red list data for

# amphibians_missing_data <- amphibian_rangemap_synonyms[!
#                            amphibian_rangemap_synonyms$tsn %in%
#                            amphibian_redlist_synonyms$tsn,]
# 
# amphibian_binomials <- unique(get_binomial_list(amphibians_missing_data))
# 
# amphibian_binomials <- amphibian_binomials[-grep("'", 
#                                amphibian_binomials)]

# Create a sub-directory for the outputs so they will be easy to collate later

# if( !dir.exists( file.path(interim_outputs, "amphibian_redlist_history") ) ) {
#   
#   dir.create( file.path(interim_outputs, "amphibian_redlist_history"), 
#               recursive = TRUE )
#   
#   amphibian_redlist_directory <- file.path(interim_outputs, 
#                                          "amphibian_redlist_history")
# }
# 
# # Loop through the names and get redlist history data, saving the output in
# # sections
# 
# amphibian_binomials <- unique(c(amphibian_rangemap_binomials, 
#                                 amphibian_redlist_binomials))
# 
# amphibian_binomial_list <- split(amphibian_binomials, 
#                                  ceiling(seq_along(amphibian_binomials)/50))
# 
# out <- list()
# 
# for (i in seq_along(amphibian_binomial_list)) {
#   
#   section_list <- amphibian_binomial_list[[i]]
#   
#   section_out <- retry(get_redlist_history(section_list, "amphibian"),
#                        maxErrors = 1000, sleep = 300) # retry function pauses the loop and tries again later when IUCN server isn't accessible
#   
#   df <- do.call(rbind, section_out)
#   
#   saveRDS(df, file.path(amphibian_redlist_directory, paste("section", i, 
#                                                          "amphibian_redlist_history.rds", 
#                                                          sep = "_")))
#   out[[i]] <- df
# 
# }
# 
# amphibian_redlist_history <- do.call(rbind, out)
# 
# saveRDS(amphibian_redlist_history, file.path(interim_outputs, paste(location,
#   "amphibian_redlist_history.rds", sep = "_")))
# 
# # Compare
# frog <- 125270 
# amphibian_redlist_data %>% filter(tsn == frog)
# amphibian_redlist_history %>% filter(tsn == frog)
# all_amphibian_synonyms %>% filter(tsn == frog)
# 
# amphibian_binomials[grep("robinsoni", amphibian_binomials)]
# 
# # Because we want to delete any false NAs (where a species actually does have a 
# # redlist status, but because it has several synonyms it shows up more than once
# # sometimes with redlist status, and also without). We want to remove those NA
# # values because they aren't true NAs and will inflate the NA count and send you
# # on a goose chase looking for the redlist data that you actually already have
# 
# all_amphibian_synonyms <- rbind(amphibian_rangemap_synonyms, amphibian_redlist_synonyms)
# 
# all_amphibian_synonyms <- distinct(all_amphibian_synonyms)
# 
# amphibians_multi_tsn_iucn <- amphibian_redlist_history %>%
#                              merge(all_amphibian_synonyms[c("tsn", "binomial", 
#                                                                  "accepted_name")],
#                                     by = "binomial",
#                                     all = TRUE) %>%
#                              rename(redlist_assessment_year = year,
#                                     redlist_status = code) %>%
#                              mutate(redlist_source = "IUCN_API") %>%
#                              dplyr::select(binomial, tsn, redlist_assessment_year, 
#                                             redlist_status, redlist_source) %>%
#                              distinct(tsn, redlist_assessment_year, redlist_status,  
#                                        redlist_source, .keep_all = TRUE) %>%
#                              group_by(tsn) %>%
#                              filter(n() > 1) %>%
#                              ungroup(.) %>%
#                              filter(!is.na(redlist_status)) 
# 
# # But if there's only one name for that species and it is listed as NA, then
# # it's truly NA and we want to keep it.  Annoying I can't work out how to 
# # satisfy both conditions in one pipe
# 
# amphibians_single_tsn_iucn <- amphibian_redlist_history %>%
#                               merge(all_amphibian_synonyms[c("tsn", "binomial", 
#                                                                  "accepted_name")],
#                                     by = "binomial",
#                                     all = TRUE) %>%
#                               rename(redlist_assessment_year = year,
#                                      redlist_status = code) %>%
#                               mutate(redlist_source = "IUCN_API") %>%
#                               dplyr::select(binomial, tsn, redlist_assessment_year, 
#                                             redlist_status, redlist_source) %>%
#                               distinct(tsn, redlist_assessment_year, redlist_status,  
#                                        redlist_source, .keep_all = TRUE) %>%
#                               group_by(tsn) %>%
#                               filter(n() == 1) %>%
#                               ungroup(.)
# 
# amphibian_redlist_data_iucn <- rbind(amphibians_single_tsn_iucn, amphibians_multi_tsn_iucn)
# 
# 
# amphibian_ecoregion_redlist_iucn <- amphibian_ecoregions %>%
#   select(-redlist_status) %>%
#   merge(amphibian_redlist_data_iucn[c("tsn", 
#                                  "redlist_assessment_year",
#                                  "redlist_status", 
#                                  "redlist_source")],
#         by = "tsn",
#         all = TRUE) %>%
#   rename(location_source = source) %>%
#   mutate(class = "Amphibia") %>%
#   merge(ecoregion_country_df[c("ECO_ID", 
#                                "CNTRY_NAME")],
#         by.x = "ecoregion_id",
#         by.y = "ECO_ID") %>%
#   filter(ecoregion_id != 0)
# 
# 
# saveRDS(amphibian_ecoregion_redlist_iucn, file.path(interim_outputs, 
#                                                paste(location,"amphibian", 
#                                                      "ecoregion_redlist_iucn.rds", 
#                                                      sep = "_")))
# 
# # Compare
# frog <- 8313 
# amphibian_redlist_data %>% filter(tsn == frog)
# amphibian_redlist_data_iucn %>% filter(tsn == frog)
# all_amphibian_synonyms %>% filter(tsn == frog)
# 
# output <- jsonlite::fromJSON(rl_search_("Amphiuma tridactylum"))

## UP TO HERE ....



# * Mammals ----

# Read in Henriques data which gives history of species' red list status

if ((paste(location, "mammal", "ecoregion_redlist.rds", 
           sep = "_") %in% list.files(interim_outputs))) {
  
mammal_ecoregion_redlist <- readRDS(file.path(interim_outputs, 
                                            paste(location, 
                                                  "mammal", "ecoregion_redlist.rds", 
                                                  sep = "_" )))
} else {

mammal_redlist_data <- read_csv(file.path(inputs, 
                                             "henriques_redlist_history",
                                             "RLTS_mammal_data_organised.csv"))

# Format and melt

mammal_redlist_data <-  mammal_redlist_data %>% 
                        mutate(class = "Mammalia") %>%
                        mutate(redlist_source = "Henriques etal 2020") %>% 
                        rename(binomial = Name) %>%
                        set_names(c("binomial","1996", "2008", "class", "redlist_source")) %>%
                        dplyr::select("binomial","class","1996", "2008", "redlist_source") %>%
                        reshape2::melt(.,id.vars = c("binomial", "class", 
                                           "redlist_source"), 
                             value.name = "redlist_status",
                             variable.name = "redlist_assessment_year") %>%
                        mutate(redlist_status = ifelse(redlist_status == "CR(PE)",
                                                       "CR", redlist_status))

# Get the synonyms for this data source

if ((paste(location, "mammal", "redlist_synonyms.rds", 
           sep = "_") %in% list.files(interim_outputs))) {
  
  mammal_redlist_synonyms <- readRDS(file.path(interim_outputs,
                                                  paste(location, "mammal", 
                                                        "redlist_synonyms.rds", 
                                                        sep = "_")))
  
  mammal_redlist_binomials <- get_binomial_list(mammal_redlist_data)
  
} else {

mammal_redlist_binomials <- get_binomial_list(mammal_redlist_data)

system.time(mammal_redlist_synonyms <- find_synonyms(mammal_redlist_binomials))

saveRDS(mammal_redlist_synonyms, file.path(interim_outputs, 
                                              paste(location,"mammal", 
                                                    "redlist_synonyms.rds", 
                                                sep = "_")))
}

# Join the synonyms and tsn to the redlist data

# Because we want to delete any false NAs (where a species actually does have a 
# redlist status, but because it has several synonyms it shows up more than once
# sometimes with redlist status, and also without). We want to remove those NA
# values because they aren't true NAs and will inflate the NA count and send you
# on a goose chase looking for the redlist data that you actually already have

mammal_multi_tsn <- mammal_redlist_data %>%
                    merge(mammal_redlist_synonyms[c("tsn", "binomial", 
                                                       "accepted_name")],
                          by = "binomial",
                          all = TRUE) %>%
                    dplyr::select(binomial, tsn, redlist_assessment_year, 
                                  redlist_status, redlist_source) %>%
                    distinct(tsn, redlist_assessment_year, redlist_status,  
                             redlist_source, .keep_all = TRUE) %>%
                    group_by(tsn) %>%
                    filter(n() > 1) %>%
                    ungroup(.) %>%
                    filter(!is.na(redlist_status)) 

# But if there's only one name for that species and it is listed as NA, then
# it's truly NA and we want to keep it.  Annoying I can't work out how to 
# satisfy both conditions in one pipe

mammal_single_tsn <- mammal_redlist_data %>%
                     merge(mammal_redlist_synonyms[c("tsn", "binomial", 
                                                       "accepted_name")],
                          by = "binomial",
                          all = TRUE) %>%
                     dplyr::select(binomial, tsn, redlist_assessment_year, 
                                  redlist_status, redlist_source) %>%
                     distinct(tsn, redlist_assessment_year, redlist_status,  
                             redlist_source, .keep_all = TRUE) %>%
                     group_by(tsn) %>%
                     filter(n() == 1) %>%
                     ungroup(.)

mammal_redlist_data <- as.data.frame(rbind(mammal_single_tsn, mammal_multi_tsn))

# Join ecoregion and redlist data together

mammal_ecoregion_redlist <- mammal_ecoregions %>%
                            dplyr::select(-redlist_status) %>%
                            merge(mammal_redlist_data[c("tsn", 
                                            "binomial",
                                            "redlist_assessment_year",
                                            "redlist_status", 
                                            "redlist_source")], by = "tsn", 
                                  all = TRUE) %>%
                            mutate(binomial = coalesce(binomial.x, binomial.y),
                                   class = "Mammalia") %>%
                            select(ecoregion_id, tsn, binomial, source, 
                                   redlist_assessment_year,
                                   redlist_status, redlist_source, class)


saveRDS(mammal_ecoregion_redlist, file.path(interim_outputs, 
                                               paste(location,"mammal", 
                                                     "ecoregion_redlist.rds", 
                                                     sep = "_")))
}

# Summarise to check for data gaps

mammal_summaries <- summarise_species_data(mammal_ecoregion_redlist, "1", "mammal")
mammal_redlist_by_ecoregion <- mammal_summaries[[1]]
mammal_redlist_global <- mammal_summaries[[2]]

# Check the data

test <- mammal_ecoregion_redlist %>% filter(redlist_status == "EX") # only one ex mammal?
rangemap <- mammal_ecoregions %>% filter(binomial == "Thylacinus cynocephalus")
redlist <- mammal_redlist_data %>% filter(binomial == "Thylacinus cynocephalus")
synonyms <- mammal_redlist_synonyms %>% filter(binomial == "Thylacinus cynocephalus")
newtest <- mammal_ecoregion_redlist %>% filter(tsn == 9275)
ecoregion <- mammal_ecoregions[c(1,80),]
x <- test_output_all %>% filter(binomial == "Thylacinus cynocephalus")

# * Birds ----

# Read in Henriques data which gives history of species' red list status

if ((paste(location, "bird", "ecoregion_redlist.rds", 
           sep = "_") %in% list.files(interim_outputs))) {
  
bird_ecoregion_redlist <- readRDS(file.path(interim_outputs, 
                                         paste(location, 
                                               "bird", "ecoregion_redlist.rds", 
                                               sep = "_" )))
} else {

bird_redlist_data <- read_csv(file.path(inputs, 
                                          "henriques_redlist_history",
                                          "RLTS_bird_data_organised.csv"))

# Format and melt

bird_redlist_data <-  bird_redlist_data %>% 
                      mutate(class = "Aves") %>%
                      mutate(redlist_source = "Henriques etal 2020") %>% 
                      set_names(c("binomial", "1988", "1994", "2000", "2004",
                                  "2008", "2012", "2016", "class", "redlist_source")) %>%
                      dplyr::select("binomial","class","1988", "1994", "2000", "2004",
                                    "2008", "2012", "2016", "redlist_source") %>%
                      melt(.,id.vars = c("binomial", "class", 
                                         "redlist_source"), 
                           value.name = "redlist_status",
                           variable.name = "redlist_assessment_year") %>%
                      mutate(redlist_status = ifelse(redlist_status == "CR(PE)"|
                                                     redlist_status == "CR (PE)"|
                                                     redlist_status == "CR(PEW)",
                                 "CR", redlist_status))

# Get the synonyms for this data source

if ((paste(location, "bird", "redlist_synonyms.rds", 
           sep = "_") %in% list.files(interim_outputs))) {
  
  bird_redlist_synonyms <- readRDS(file.path(interim_outputs,
                                               paste(location, "bird", 
                                                     "redlist_synonyms.rds", 
                                                     sep = "_")))
  
  bird_redlist_binomials <- get_binomial_list(bird_redlist_data)
  
} else {

bird_redlist_binomials <- get_binomial_list(bird_redlist_data)

system.time(bird_redlist_synonyms <- find_synonyms(bird_redlist_binomials))

saveRDS(bird_redlist_synonyms, file.path(interim_outputs, 
                                           paste(location,"bird", 
                                                 "redlist_synonyms.rds", 
                                                 sep = "_")))
}

# Join the synonyms and tsn to the redlist data

# Because we want to delete any false NAs (where a species actually does have a 
# redlist status, but because it has several synonyms it shows up more than once
# sometimes with redlist status, and also without). We want to remove those NA
# values because they aren't true NAs and will inflate the NA count and send you
# on a goose chase looking for the redlist data that you actually already have

bird_multi_tsn <- bird_redlist_data %>%
                  merge(bird_redlist_synonyms[c("tsn", "binomial", 
                                                  "accepted_name")],
                        by = "binomial",
                        all = TRUE) %>%
                  dplyr::select(binomial, tsn, redlist_assessment_year, 
                                redlist_status, redlist_source) %>%
                  distinct(tsn, redlist_assessment_year, redlist_status,  
                           redlist_source, .keep_all = TRUE) %>%
                  group_by(tsn) %>%
                  filter(n() > 1) %>%
                  ungroup(.) %>%
                  filter(!is.na(redlist_status)) 

# But if there's only one name for that species and it is listed as NA, then
# it's truly NA and we want to keep it.  Annoying I can't work out how to 
# satisfy both conditions in one pipe

bird_single_tsn <- bird_redlist_data %>%
  merge(bird_redlist_synonyms[c("tsn", "binomial", 
                                  "accepted_name")],
        by = "binomial",
        all = TRUE) %>%
  dplyr::select(binomial, tsn, redlist_assessment_year, 
                redlist_status, redlist_source) %>%
  distinct(tsn, redlist_assessment_year, redlist_status,  
           redlist_source, .keep_all = TRUE) %>%
  group_by(tsn) %>%
  filter(n() == 1) %>%
  ungroup(.)


# Merge only the species that we have a tsn for

bird_redlist_data <- as.data.frame(rbind(bird_single_tsn, 
                                              bird_multi_tsn))

# Join ecoregion and redlist data together

bird_ecoregion_redlist <- bird_ecoregions %>%
                          dplyr::select(-redlist_status) %>%
                          merge(bird_redlist_data[c("tsn", 
                                                         "binomial",
                                                         "redlist_assessment_year",
                                                         "redlist_status", 
                                                         "redlist_source")], by = "tsn", 
                                all = TRUE) %>%
                          mutate(binomial = coalesce(binomial.x, binomial.y),
                                 class = "Aves") %>%
                          select(ecoregion_id, tsn, binomial, source, 
                                 redlist_assessment_year,
                                 redlist_status, redlist_source, class)


saveRDS(bird_ecoregion_redlist, file.path(interim_outputs, 
                                            paste(location,"bird", 
                                                  "ecoregion_redlist.rds", 
                                                  sep = "_")))


}

# Summarise to check for data gaps

bird_summaries <- summarise_species_data(bird_ecoregion_redlist, "1", "bird")
bird_redlist_by_ecoregion <- bird_summaries[[1]]
bird_redlist_global <- bird_summaries[[2]]

# Bird data checks!
# bird_rangemaps_no_tsn <- bird_rangemap_synonyms %>% filter(is.na(tsn))
# bird_ecoregions_no_tsn <- bird_ecoregions %>% filter(is.na(tsn))
# bird_redlist_no_tsn <- bird_redlist_data %>% filter(is.na(tsn))
# bird_ER_RL_no_tsn <- bird_ecoregion_redlist %>% filter(is.na(tsn))
# 
# australian_birds <- bird_ecoregion_redlist[bird_ecoregion_redlist$ecoregion_id %in% 
#                                     ecoregion_subset$ECO_ID,] 
# extinct_birds <- bird_ecoregion_redlist %>% filter(redlist_status == "EX")

# * Reptiles ----

if ((paste(location, "reptile", "ecoregion_redlist.rds", 
           sep = "_") %in% list.files(interim_outputs))) {
  
reptile_ecoregion_redlist <- readRDS(file.path(interim_outputs, 
                                              paste(location, 
                                                    "reptile", "ecoregion_redlist.rds", 
                                                    sep = "_" )))
} else {

## Get list of names to search (search all synonyms)
# all <- reptile_binomials_all

reptile_binomials_all <- unique(reptile_rangemap_synonyms$binomial)

# Remove names that have extra quotation marks because they stop the function working

reptile_binomials_all <- reptile_binomials_all[-grep("'", reptile_binomials_all)]

# Because the IUCN server is a bit unpredictable, split the list into sections 
# so we can get the redlist history data in a loop and save progress iteratively

reptile_binomial_list <- split(reptile_binomials_all, 
                               ceiling(seq_along(reptile_binomials_all)/50))

# reptile_binomial_list_all <- reptile_binomial_list

# Create a sub-directory for the outputs so they will be easy to collate later

if ((paste(location, "reptile_redlist_data.rds", 
           sep = "_") %in% list.files(interim_outputs))) {
  
reptile_redlist_data <- readRDS(file.path(interim_outputs, 
                                   paste(location, "reptile_redlist_data.rds", 
                                                       sep = "_" )))
} else {


reptile_redlist_directory <- file.path(interim_outputs, 
                                       "reptile_redlist_history")

if( !dir.exists( file.path(interim_outputs, "reptile_redlist_history") ) ) {
  
  dir.create( file.path(interim_outputs, "reptile_redlist_history"), 
              recursive = TRUE )
  
}

# Loop through the names and get redlist history data, saving the output in
# sections

# WARNING - SLOW CODE, TAKES AROUND 12 HOURS

out <- list()

for (i in seq_along(reptile_binomial_list)) {
  
  section_list <- reptile_binomial_list[[i]]
  
  section_out <- retry(get_redlist_history(section_list, "reptile"),
                                           maxErrors = 1000, sleep = 300) # retry function pauses the loop and tries again later when IUCN server isn't accessible
  df <- do.call(rbind, section_out)
  
  saveRDS(df, file.path(reptile_redlist_directory, paste("section", i ,
                                                "reptile_redlist_history.rds", 
                                                sep = "_")))
  out[[i]] <- df
}

# Convert list of redlist data back into a nice dataframe 

all_reptile_out <- list.files(reptile_redlist_directory)
out <- lapply(file.path(reptile_redlist_directory, all_reptile_out), readRDS)

reptile_redlist_data <- do.call(rbind, out)

saveRDS(reptile_redlist_data, 
        file.path(interim_outputs, paste(location, "reptile_redlist_data.rds",
                                         sep = "_")))

}

status_list <- reptile_redlist_data %>%
               dplyr::select(code, category) %>%
               distinct(.)

#' TODO: Important - check how 'equivalent' the older Red List categories 
#' are with the new ones. Might have to turf the older time points
#' https://www.animalinfo.org/notes.htm

reptile_redlist_data<- reptile_redlist_data %>%
                        rename(redlist_assessment_year = year) %>%
                        mutate(redlist_source = "IUCN API") %>%
                        mutate(redlist_status = ifelse(category == "Data Deficient"|
                                                       category == "Indeterminate"|
                                                       category == "Insufficiently known"|
                                                       category == "Rare",
                                                       "DD",
                                                ifelse(category == "Least Concern"|
                                                       category == "Lower risk/least concern",
                                                       "LC",
                                                ifelse(category == "Near Threatened"|
                                                       category == "Lower Risk/near threatened",
                                                       "NT",
                                                ifelse(category == "Vulnerable"|
                                                       category == "Lower Risk/near threatened",
                                                       "VU",
                                                ifelse(category == "Endangered",
                                                       "EN",
                                                ifelse(category == "Critically endangered",
                                                       "CR",
                                                ifelse(category == "Extinct in the wild",
                                                       "EW",
                                                ifelse(category == "Extinct",
                                                              "EX", NA ))))))))) %>%
                          dplyr::select(-code,-category)
    
# We don't need to get the synonyms for this data source because it's the same
# as the rangemap data source

reptile_redlist_synonyms <- reptile_rangemap_synonyms

# Join the synonyms and tsn to the redlist data

# Because we want to delete any false NAs (where a species actually does have a 
# redlist status, but because it has several synonyms it shows up more than once
# sometimes with redlist status, and also without). We want to remove those NA
# values because they aren't true NAs and will inflate the NA count and send you
# on a goose chase looking for the redlist data that you actually already have

reptile_multi_tsn <- reptile_redlist_data %>%
  merge(reptile_redlist_synonyms[c("tsn", "binomial", 
                                "accepted_name")],
        by = "binomial",
        all = TRUE) %>%
  dplyr::select(binomial, tsn, redlist_assessment_year, 
                redlist_status, redlist_source) %>%
  distinct(tsn, redlist_assessment_year, redlist_status,  
           redlist_source, .keep_all = TRUE) %>%
  group_by(tsn) %>%
  filter(n() > 1) %>%
  ungroup(.) %>%
  filter(!is.na(redlist_status)) 

# But if there's only one name for that species and it is listed as NA, then
# it's truly NA and we want to keep it.  Annoying I can't work out how to 
# satisfy both conditions in one pipe

reptile_single_tsn <- reptile_redlist_data %>%
  merge(reptile_redlist_synonyms[c("tsn", "binomial", 
                                "accepted_name")],
        by = "binomial",
        all = TRUE) %>%
  dplyr::select(binomial, tsn, redlist_assessment_year, 
                redlist_status, redlist_source) %>%
  distinct(tsn, redlist_assessment_year, redlist_status,  
           redlist_source, .keep_all = TRUE) %>%
  group_by(tsn) %>%
  filter(n() == 1) %>%
  ungroup(.)

# Merge only the species that we have a tsn for

reptile_redlist_data <- as.data.frame(rbind(reptile_single_tsn, 
                                            reptile_multi_tsn))

# Join ecoregion and redlist data together

reptile_ecoregion_redlist <- reptile_ecoregions %>%
  dplyr::select(-redlist_status) %>%
  merge(reptile_redlist_data[c("tsn", 
                            "binomial",
                            "redlist_assessment_year",
                            "redlist_status", 
                            "redlist_source")], by = "tsn", 
        all = TRUE) %>%
  mutate(binomial = coalesce(binomial.x, binomial.y),
         class = "Reptilia") %>%
  select(ecoregion_id, tsn, binomial, source, 
         redlist_assessment_year,
         redlist_status, redlist_source, class)



saveRDS(reptile_ecoregion_redlist, file.path(interim_outputs, 
                                          paste(location,"reptile", 
                                                "ecoregion_redlist.rds", 
                                                sep = "_")))


}


# Summarise to check for data gaps

reptile_summaries <- summarise_species_data(reptile_ecoregion_redlist, "1", "reptile")
reptile_redlist_by_ecoregion <- reptile_summaries[[1]]
reptile_redlist_global <- reptile_summaries[[2]]

# Species Data 1: Combine rangemaps and redlist data ----

species_data <- rbind(amphibian_ecoregion_redlist,
                      mammal_ecoregion_redlist,
                      bird_ecoregion_redlist,
                      reptile_ecoregion_redlist)

saveRDS(species_data, file.path(interim_outputs, 
                                             paste(location,"species_data_1.rds", 
                                                   sep = "_")))

# Species Data 2: Find and add missing ecoregions to species data 1 ----

if ((paste(location, "species_data_2.rds", 
                      sep = "_") %in% list.files(interim_outputs))) {
  
species_data_2 <- readRDS(file.path(interim_outputs,paste(location,
                                                          "species_data_2.rds", 
                                                          sep = "_") ))

} else {
  
species_missing_ecoregions <- species_data %>% filter(is.na(ecoregion_id))
species_missing_ecoregions_binomials <- unique(species_missing_ecoregions$binomial)
species_found_ecoregions <- get_gbif_data(species_missing_ecoregions_binomials, 
                                          100, ecoregion_map)

species_found_ecoregions <- species_found_ecoregions[[1]] %>%
                            distinct(.) %>%
                            rename(binomial = species) 

saveRDS(species_found_ecoregions, file.path(interim_outputs, "global_species_found_ecoregions.rds"))

species_data_2 <- species_data %>%
                  merge(species_found_ecoregions[c("ECO_ID", "binomial")],
                        by ="binomial", all = TRUE) %>%
                  mutate(eco_source = "gbif",
                         ecoregion_id = coalesce(ecoregion_id, ECO_ID),
                         source = coalesce(source, eco_source)) %>%
                  select(ecoregion_id, tsn, binomial, source, 
                         redlist_assessment_year,
                         redlist_status, redlist_source, class)

saveRDS(species_data_2, file.path(interim_outputs, 
                                paste(location,"species_data_2.rds", 
                                      sep = "_")))
}

# Species Data 3: Find and add extinct species to species data 2 ----

# Get the names of all extinct species in IUCN Red List
  
  iucn_scraped_extinct_species <- retry(rl_sp_category("EX", key = NULL, 
                                                       parse = TRUE),
                                        maxErrors = 1000, sleep = 300)
  
  iucn_scraped_ew_species <- retry(rl_sp_category("EW", key = NULL, 
                                                  parse = TRUE),
                                     maxErrors = 1000, sleep = 300)
  
  # Get their Red List Assessment history
  
  ## Clean names and put into a list so can get and save the data in chunks from RLI API
 
  extinct_species_names <- unique(c(iucn_scraped_extinct_species[[3]]$scientific_name,
                                    iucn_scraped_ew_species[[3]]$scientific_name))
  
  ## Remove the subspecies
  
  extinct_species_names <- extinct_species_names[stri_count(extinct_species_names,
                                                 regex="\\S+") < 3]
    
  extinct_species_names <- unique(extinct_species_names)
  
   
  ## Create a sub-directory for the outputs so they will be easy to collate later
  
  if ((paste(location, "extinct_species_redlist_data.rds", 
             sep = "_") %in% list.files(interim_outputs))) {
    
   extinct_species_redlist_data <- readRDS(file.path(interim_outputs, 
                                              paste(location, 
                                                    "extinct_species_redlist_data.rds", 
                                                    sep = "_" )))
  } else {
    
    
  extinct_species_redlist_directory <- file.path(interim_outputs, 
                                           "extinct_species_redlist_history")
    
  if( !dir.exists( file.path(interim_outputs, "extinct_species_redlist_history") ) ) {
      
      dir.create( file.path(interim_outputs, "extinct_species_redlist_history"), 
                  recursive = TRUE )
      
  }
  
  extinct_species_binomial_list <- split(extinct_species_names, 
                                           ceiling(seq_along(extinct_species_names)/50))
    
  
  out <- list()
  
  for (i in seq_along(extinct_species_binomial_list)) {
    
    section_list <- extinct_species_binomial_list[[i]]
    
    section_out <- retry(get_redlist_history(section_list, "extinct_species"),
                         maxErrors = 1000, sleep = 300) # retry function pauses the loop and tries again later when IUCN server isn't accessible
    df <- do.call(rbind, section_out)
    
    saveRDS(df, file.path(extinct_species_redlist_directory, paste("section", i ,
                                          "extinct_species_redlist_history.rds", 
                                                           sep = "_")))
    out[[i]] <- df
  }
  
  ## Convert list of redlist history back into a nice dataframe 
  
  all_extinctions_out <- list.files(extinct_species_redlist_directory)
  out <- lapply(file.path(extinct_species_redlist_directory, 
                          all_extinctions_out), readRDS)
  
  extinct_species_redlist_data <- do.call(rbind, out)
  
 # Older Red List categories https://www.animalinfo.org/notes.htm
  
  extinct_species_redlist_data <- extinct_species_redlist_data %>%
                                  select(-class) %>%
                                  merge(extinct_species_synonyms[c('tsn', 
                                                                   'class', 
                                                                   'binomial')],
                                        by = 'binomial',
                                        all = TRUE) %>%
                                  mutate(redlist_status = 
                                  ifelse(
                                         category == "Data Deficient"|
                                         category == "Indeterminate"|
                                         category == "Insufficiently known"|
                                         category == "Insufficiently Known"|
                                         category == "Not Recognized"|
                                         category == "Status inadequately known-survey required or data sought",
                                         "DD",
                                  ifelse(category == "Least Concern"|
                                         category == "Lower risk/least concern",
                                         "LC",
                                  ifelse(category == "Near Threatened"|
                                         category == "Lower Risk/near threatened"|
                                         category == "Rare"|
                                         category == "Very rare but believed to be stable or increasing",
                                         "NT",
                                  ifelse(category == "Vulnerable"|
                                         category == "Very rare and believed to be decreasing in numbers",
                                         "VU",
                                  ifelse(category == "Endangered"|
                                         category == "Threatened",
                                         "EN",
                                  ifelse(category == "Critically endangered"|
                                         category == "Critically Endangered",
                                                                      "CR",
                                  ifelse(category == "Extinct in the wild"|
                                         category == "Extinct in the Wild",
                                                      "EW",
                                  ifelse(category == "Extinct"|
                                         category == "Extinct?"|
                                         category == "Extinct/Endangered",
                                                     "EX", NA ))))))))) %>%
    mutate(redlist_status = ifelse(str_detect(category, "inadequately"), "DD",
                            ifelse(str_detect(category, "decreasing"), "VU",
                            ifelse(str_detect(category, "increasing"), "NT",
                                   redlist_status)))) %>%
    rename(redlist_assessment_year = year) %>%
    mutate(redlist_source = "IUCN API")
  
  
  saveRDS(extinct_species_redlist_data, file.path(interim_outputs, 
                                        paste(location, 
                                        "extinct_species_redlist_data.rds", 
                                                        sep = "_" )))
  
  ## Get ecoregions for the extinct species
  
  ## Get extinct species synonyms
  
  #' TODO: IMPORTANT - why are there only around 200 species for which synonyms were found?
  
  extinct_species_synonyms <- find_synonyms(extinct_species_names)
  
  extinct_species_names_all <- unique(c(extinct_species_synonyms$binomial, 
                               extinct_species_synonyms$accepted_name,
                               extinct_species_names))
  
  extinct_species_ecoregions_gbif <- get_gbif_data(extinct_species_names_all, 
                                            1000, ecoregion_map)
  
  extinct_species_ecoregions <- extinct_species_ecoregions_gbif[[1]] %>%
                                distinct(.) %>%
                                rename(binomial = species) 
  
  extinct_species_ecoregions <- extinct_species_ecoregions %>%
                                merge(extinct_species_synonyms[c("binomial", "tsn")], 
                                      by = "binomial", all = TRUE) %>%
                                mutate(source = "gbif")
  
  saveRDS(extinct_species_ecoregions, file.path(interim_outputs, 
                                                paste(location, 
                                                "extinct_species_ecoregions.rds", 
                                                                       sep = "_" )))
  
  extinct_species_names_df <- as.data.frame(extinct_species_names)
  
  extinct_species_without_ecoregions <- extinct_species_names_df[!extinct_species_names_df$extinct_species_names %in% 
                                        extinct_species_ecoregions$binomial,] 
 
  #### UP TO HERE!!! THE BELOW WON'T WORK PROPERLY ####
  
   # Join ecoregion and redlist data
  
 extinct_species_ecoregion_redlist <- extinct_species_ecoregions %>%
                                      merge(extinct_species_redlist_data[c( 
                                   "binomial",
                                   "redlist_assessment_year",
                                   "redlist_status", 
                                   "redlist_source",
                                   "class")], by = "binomial", 
          all = TRUE) %>%
    rename(ecoregion_id = ECO_ID) %>% # This part needs fixing
    select(ecoregion_id, tsn, binomial, source, 
           redlist_assessment_year,
           redlist_status, redlist_source, class)
  
  saveRDS(extinct_species_ecoregion_redlist, file.path(interim_outputs, 
                                                 paste(location,"extinct_species", 
                                                       "ecoregion_redlist.rds", 
                                                       sep = "_")))
  
  
  
  species_data_3 <- species_data_2 %>%
    merge(extinct_species_ecoregion_redlist[c("ecoregion_id", "binomial")],
          by ="binomial", all = TRUE) %>%
    mutate(eco_source = "gbif",
           ecoregion_id = coalesce(ecoregion_id.x, ecoregion_id.y),
           source = coalesce(source, eco_source)) %>%
    select(ecoregion_id, tsn, binomial, source, 
           redlist_assessment_year,
           redlist_status, redlist_source, class)
  

  saveRDS(species_data_3, file.path(interim_outputs, 
                                    paste(location,"species_data_3.rds", 
                                          sep = "_")))
  
  
  }



# Check data for gaps ----

species_redlist_global <- rbind(amphibian_redlist_global,
                                mammal_redlist_global,
                                bird_redlist_global) %>%
                          group_by(redlist_status) %>%
                          summarise(redlist_status = sum(redlist_status))



#' intersect_ranges_w_ecoregions <- function(class_name, location, shapefile, data) {
#'   
#' 
#' if (!(paste(class_name, location, "iucn_range_map_species_with_ecoregions.rds", 
#'             sep = "_") %in% list.files(interim_outputs))) {
#'   
#'   range_directories <- list(file.path(inputs, paste("redlist", class_name,  
#'                                                    "range_maps", sep = "_")))
#'   
#'   iucn_rangemap_database <- list()
#'   
#'   for (i in seq_along(range_directories)) {
#'     
#'     iucn_rangemap_database[[i]] <- get_ecoregions(range_directories[[i]], 
#'                                                   shapefile)
#'     
#'     print(paste("processed", i, "of" , length(range_directories), 
#'                 "group range maps", sep = " "))
#'     
#'   }
#'   
#'   iucn_rangemap_database <- do.call(rbind, iucn_rangemap_database)
#'   iucn_rangemap_database <- distinct(iucn_rangemap_database)
#'   
#'   saveRDS(iucn_rangemap_database, 
#'           file = file.path(interim_outputs, paste(class_name,
#'           location,"iucn_range_map_species_with_ecoregions.rds",
#'           sep = "_")))
#' 
#' } else {
#'   
#'   iucn_rangemap_database <- readRDS(file.path(interim_outputs, 
#'                                     paste(class_name,
#'                                           location,
#'                                     "iucn_range_map_species_with_ecoregions.rds",
#'                                     sep = "_")))
#' }
#' 
#' if (!is.na(country)) {
#'   
#'   iucn_rangemap_database <- iucn_rangemap_database[
#'                             iucn_rangemap_database$ecoregion_id %in% 
#'                             ecoregion_subset$ECO_ID, ]
#' }
#' 
#' # Get TSNs for the species from the IUCN range maps
#' 
#' if (!(file.path(outputs, paste(class_name, location, "iucn_rangemap_synonyms.rds",
#'                       sep = "_")) %in% list.files(outputs))) {
#'   
#'   # pull out unique species names only
#'   
#'   iucn_rangemap_species <- iucn_rangemap_database %>%
#'     dplyr::select(binomial) 
#'   
#'   iucn_rangemap_species$binomial <- as.character(iucn_rangemap_species$binomial)
#'   
#'   iucn_rangemap_species <- unname(unlist(iucn_rangemap_species[,1]))
#'   
#'   iucn_rangemap_species <- iucn_rangemap_species[!is.na(iucn_rangemap_species)]
#'   
#'   iucn_rangemap_species <- unique(iucn_rangemap_species)
#'   
#'   # search for synonyms and add taxonomic serial number (slow)
#'   
#'   iucn_rangemap_species <- find_synonyms(iucn_rangemap_species)
#'   
#'   # save synonyms for this groups of species
#'   
#'   saveRDS(iucn_rangemap_species, file.path(outputs, 
#'                                            paste(class_name,
#'                                                  location,
#'                                                  "iucn_rangemap_synonyms.rds",
#'                                                  sep = "_")))
#'   
#'   # add the synonyms and taxonomic info to the ecoregion-species data
#'   
#'   iucn_rangemap_database <-  iucn_rangemap_database %>%
#'     merge(iucn_rangemap_species[c("binomial", "tsn", 
#'                                   "found","class")], 
#'           by = "binomial")
#'   
#'   # save the species with ecoregions data
#'   
#'   saveRDS(iucn_rangemap_database, file.path(interim_outputs, 
#'                                             paste(class_name,
#'                                                   location,
#'                                             "iucn_rangemap_database.rds", 
#'                                             sep = "_")))
#'   
#' } else {
#'   
#'   # Read in the synonyms previously stored
#'   
#'   iucn_rangemap_species <- readRDS(file.path(outputs, 
#'                                              paste(class_name,
#'                                                    location,
#'                                                    "iucn_rangemap_synonyms.rds",
#'                                                     sep = "_")))
#'   
#'   # Merge with ecoregion data
#'   
#'   iucn_rangemap_database <-  iucn_rangemap_database %>%
#'     merge(iucn_rangemap_species[c("binomial", "tsn", 
#'                                   "found", "class")], 
#'           by = "binomial")
#'   
#'   saveRDS(iucn_rangemap_database, file.path(interim_outputs, paste(class_name,
#'                                                            location,
#'                                             "iucn_rangemap_database.rds", 
#'                                             sep = "_")))
#'   
#' }
#' 
#' # Load red list status data ----
#' 
#' file_path <- file.path(inputs, "redlist_sergio")
#' file_names <- list.files(file_path)
#' file <- file_names[str_detect(file_names, class_name)]
#' files <- file.path(file_path, file)
#' table <- read_csv(files)
#' #names(tables) <- str_remove(file_names, ".csv")
#' 
#' # Add a column so we can identify what taxa each table contains + the source
#' # reference
#' 
#' if (class_name == "amphibian") {
#'   
#' matched_class <- "Amphibia"
#' 
#' table <- table %>% 
#'          mutate(class = matched_class) %>%
#'          mutate(redlist_source = "Henriques etal 2020") %>% 
#'          dplyr::mutate(binomial = paste(Genus, Species, sep = " ")) %>%
#'          dplyr::select(-c(Genus, Species)) %>%
#'          set_names(c("2004", "2008", "class", "redlist_source", "binomial")) %>%
#'          dplyr::select("binomial","class","2004", "2008", "redlist_source")
#' 
#' }
#' 
#' if (class_name == "mammal") {
#'   
#'   matched_class <- "Mammalia"
#'   
#'   table <- table %>% 
#'            mutate(class = matched_class) %>%
#'            mutate(redlist_source = "Henriques etal 2020") %>% 
#'            rename(binomial = Name) %>%
#'            set_names(c("binomial","1996", "2008", "class", "redlist_source")) %>%
#'            dplyr::select("binomial","class","1996", "2008", "redlist_source")
#'   
#' }
#' 
#' if (class_name == "reptile") {
#'   
#'   matched_class <- "Reptilia"
#'   
#'   table <- table %>% 
#'            mutate(class = matched_class) %>%
#'            mutate(redlist_source = "Henriques etal 2020") %>% 
#'            rename(binomial = Name) %>%
#'            set_names(c("binomial","1996", "2008", "class", "redlist_source")) %>%
#'            dplyr::select("binomial","class","1996", "2008", "redlist_source")
#'   
#' }
#' 
#' # Get synonyms for henriques species ----
#' 
#' henriques_species <- table %>%
#'                      dplyr::select(binomial) 
#' 
#' henriques_species$binomial <- as.character(henriques_species$binomial)
#' henriques_species <- unname(unlist(henriques_species[,1]))
#' henriques_species <- henriques_species[!is.na(henriques_species)]
#' henriques_species <- unique(henriques_species)
#' 
#' if (!(paste(class_name,
#'             "henriques_synonyms.rds",
#'             sep = "_") %in% list.files(interim_outputs))) { 
#'   
#'   henriques_species <- find_synonyms(henriques_species)
#'   
#'   saveRDS(henriques_species, file.path(outputs, 
#'                                        paste(class_name,
#'                                              "henriques_synonyms.rds",
#'                                              sep = "_")))
#'   
#' } else {
#'   
#'   
#'   henriques_species <- readRDS(file.path(outputs, 
#'                                          paste(class_name,
#'                                                "henriques_synonyms.rds",
#'                                                sep = "_")))
#'   
#' }
#' 
#' # Order columns, drop group variable and melt in to long format
#' 
#' redlist_status_data <- table %>%
#'                        select(-class) %>%
#'                        melt(.,id.vars = c("binomial",  
#'                                           "redlist_source"),
#'                              value.name = "redlist_status") %>%
#'                        rename(redlist_assessment_year = variable) %>% 
#'                        merge(henriques_species[c("binomial", "tsn", 
#'                                                  "accepted_name", 
#'                                                  "common_name", "class")]) %>%
#'                        rename(accepted_binomial = accepted_name)
#'     
#' 
#' # Merge redlist and ecoregion data ----
#' 
#' ecoregion_redlist_data <- redlist_status_data %>%
#'                           merge(iucn_rangemap_database[c("tsn",
#'                                                          "ecoregion_id",
#'                                                          "eco_objectid")], 
#'                                 by = "tsn", all = TRUE) %>%
#'                           select(tsn, accepted_binomial, 
#'                                  class, redlist_assessment_year,
#'                                  redlist_status, common_name, eco_objectid,
#'                                  ecoregion_id) %>% 
#'                           merge(data[c("ECO_ID", "CNTRY_NAME")], 
#'                                 by.x = "ecoregion_id", by.y =  "ECO_ID" ) %>%
#'                           rename(country = "CNTRY_NAME")
#' 
#' 
#' saveRDS(ecoregion_redlist_data, file.path(interim_outputs, 
#'                                      paste(class_name,
#'                                            location,
#'                                            "species_data_1.rds",
#'                                            sep = "_")))
#' 
#' return(ecoregion_redlist_data)
#' 
#' }
#' 
#' # Intersect ranges and ecoregions for each class
#' # (can put this in a loop but takes a long time to run each,
#' # probably more pragmatic to split them up or parallelise)
#' 
#' amphibians <- intersect_ranges_w_ecoregions("amphibian", location, 
#'                                             ecoregion_map, ecoregion_country_df)
#' 
#' mammals <- intersect_ranges_w_ecoregions("mammal", location, 
#'                                          ecoregion_map, ecoregion_country_df)
#' 
#' #reptiles <- intersect_ranges_w_ecoregions("reptile", ecoregion_map)
#' 
#' # Combine into one file and save
#' 
#' species_data <- rbind(amphibians, mammals)
#' 
#' saveRDS(species_data, file.path(interim_outputs, 
#'                                 paste(location, "species_data_1.rds",
#'                                       sep = "_")))
#' 

#' 
#' amphibian_summary <- summarise_species_data(amphibians, 1)
#' mammal_summary <- summarise_species_data(mammals, 1)
#' 
#' 
#' # old_sp_data <- readRDS(file.path(inputs, "deakin_species_data", "species_data.rds"))
#' # group_tables <- list()
#' # 
#' # for (i in seq_along(tables)) {
#' #   
#' #   group_tables[[i]] <- tables[[i]] %>% 
#' #     mutate(group = names(tables[i])) %>%
#' #     mutate(redlist_source = "Henriques etal 2020")
#' #   
#' #   
#' # }
#' 
#' # Standardise column names etc. because they are inconsistent between the 
#' # different taxa tables
#' 
#' # amphibians <- group_tables[[1]]
#' # 
#' # amphibians <- amphibians %>% 
#' #   dplyr::mutate(binomial = paste(Genus, Species, sep = " ")) %>%
#' #   dplyr::select(-c(Genus, Species)) %>%
#' #   set_names(c("2004", "2008", "group", "redlist_source", "binomial")) %>%
#' #   dplyr::select("binomial","group","2004", "2008", "redlist_source") 
#' 
#' birds <- group_tables[[2]]
#' 
#' birds <- birds %>%
#'   set_names(c("binomial", "1988", "1994", "2000", "2004", "2008" , 
#'               "2012" , "2016", "group", "redlist_source" )) %>%
#'   dplyr::select("binomial", "group", everything())
#' 
#' mammals <- group_tables[[5]]  
#' 
#' mammals <- mammals %>%
#'   set_names(c("binomial", "1996", "2008", "group", "redlist_source" )) %>%
#'   dplyr::select("binomial", "group", everything())
#' 
#' # Combine back into a list then dataframe
#' 
#' group_tables <- list(birds, mammals, amphibians)
#' 
#' #' TODO: Change name of species_redlist to henrique database?
#' 
#' henriques_redlist_database <- do.call(bind_rows, group_tables)
#' 
#' # Get Henriques synonyms ----
#' 
#' henriques_species <- henriques_redlist_database %>%
#'   dplyr::select(binomial) 
#' 
#' henriques_species$binomial <- as.character(henriques_species$binomial)
#' henriques_species <- unname(unlist(henriques_species[,1]))
#' henriques_species <- henriques_species[!is.na(henriques_species)]
#' 
#' if (!("henrique_species_synonyms.rds" %in% list.files(interim_outputs))) { 
#'   
#'   henriques_species <- find_synonyms(henriques_species)
#'   
#'   saveRDS(henriques_species, file.path(interim_outputs, "henrique_species_synonyms.rds"))
#'   
#' } else {
#'   
#'   
#'   henriques_species <- readRDS(file.path(interim_outputs, "henrique_species_synonyms.rds"))
#'   
#' }
#' 
#' # Order columns, drop group variable and melt in to long format
#' 
#' henriques_redlist_database <- henriques_redlist_database %>%
#'   dplyr::select("binomial", "1988", "1994", "1996", "2000", 
#'                 "2004", "2008", "2012", "2016", "redlist_source") %>%
#'   melt(.,id.vars = c("binomial", "redlist_source"),
#'        value.name = "redlist_status") %>%
#'   rename(redlist_assessment_year = variable) %>% 
#'   merge(henriques_species[c("binomial", "tsn", "accepted_name", 
#'                             "common_name", "class")])
#' 
#' 
#' # Load Wildfinder database ----
#' 
#' # Load tables from WildFinder database (converted into .xlsx files from .mdb 
#' # database)
#' 
#' #' TODO: Find out what year wildfinder RL status is from
#' #' TODO: Figure out if there's a way to load the tables directly from the .mdb
#' 
#' # Load the wildfinder data
#' 
#' file_path <- file.path(inputs, "wildfinder_csv")
#' file_names <- list.files(file_path)
#' files <- file.path(file_path, file_names)
#' tables <- lapply(files, read_csv) # Disregard warnings
#' names(tables) <- str_remove(file_names, ".csv")
#' 
#' # Change all the column names to lower case
#' 
#' for (i in seq_along(tables)) {
#' 
#'   names(tables[[i]]) <- tolower(names(tables[[i]]))
#'   
#' }
#' 
#' # Get list of species per ecoregion (doesn't include binomial or common name)
#' 
#' ecoregions_species <- tables$ecoregion_species
#' 
#' # Add common names, binomial scientific names, synonyms and tsn IDs
#' 
#' common_names <- tables$common_names
#' 
#' scientific_names <- as.data.frame(cbind(tables$redlist_species$wwf_species_id, 
#'                           tables$redlist_species$genus,
#'                           tables$redlist_species$species))
#' colnames(scientific_names) <- c("wwf_species_id", "genus", "species")
#' 
#' wildfinder_species <- scientific_names %>%
#'                       mutate(binomial = paste(genus, species, sep = " ")) %>%
#'                       dplyr::select(binomial)
#' 
#' wildfinder_species <- wildfinder_species[,1]
#' 
#' # Get Wildfinder synonyms (SLOW code) ----
#' 
#' if(!("wildfinder_species_synonyms.rds" %in% list.files(interim_outputs))) { 
#'   
#' wildfinder_species <- find_synonyms( wildfinder_species)
#' 
#' saveRDS(wildfinder_species, file.path(interim_outputs,"wildfinder_species_synonyms.rds"))
#' 
#' } else {
#' 
#' 
#' wildfinder_species <- readRDS(file.path(interim_outputs, 
#'                                         "wildfinder_species_synonyms.rds"))
#' 
#' }
#' 
#' wildfinder_database <-  ecoregions_species %>%
#'                         merge(common_names[c("species_id", "common_name")],
#'                               all = TRUE) %>%
#'                         rename(wwf_species_id = species_id) %>%
#'                         dplyr::select(c("wwf_species_id", "common_name",
#'                                         "ecoregion_code")) %>%
#'                         merge(scientific_names, by = "wwf_species_id",
#'                               all = TRUE) %>%
#'                         mutate(binomial = paste(genus, species, sep = " ")) %>%
#'                         mutate(source = "WildFinder") %>%
#'                         mutate(wwf_species_id = as.numeric(wwf_species_id)) %>%
#'                         merge(wildfinder_species[c("binomial", "tsn", 
#'                                                    "accepted_name", "class")])
#' 
#' 
#' # Remove un-needed tables
#' 
#' rm(tables)
#' 
#' # Clear memory by removing large objects
#' 
#' rm(common_names)
#' 
#' # Subset the wildfinder database by country if working with subset
#' 
#' if (!is.na(country)) {
#'   
#'   full_WF_database <- wildfinder_database
#'   
#'   country_ecoregions <- ecoregion_subset$ECO_ID
#'   
#'   wildfinder_database <- full_WF_database[full_WF_database$ecoregion_code %in% country_ecoregions, ]
#'   
#' }
#' 
#' # Load Cooke database ----
#' 
#' cooke_database <- read.csv(paste(inputs,"/cooke_database/species_site.csv", sep = ""))
#' 
#' cooke_species <- cooke_database %>%
#'                  dplyr::select(binomial) 
#' 
#' cooke_species$binomial <- as.character(cooke_species$binomial)
#' cooke_species <- cooke_species[,1]
#' 
#' # Get Cooke synonyms (SLOW code) ----
#' 
#' if(!("cooke_species_synonyms.rds" %in% list.files(interim_outputs))) { 
#'   
#'   cooke_species <- find_synonyms(cooke_species)
#'   
#'   saveRDS(cooke_species, file.path(interim_outputs, "cooke_species_synonyms.rds"))
#'   
#' } else {
#'   
#'   cooke_species <- readRDS(file.path(interim_outputs, "cooke_species_synonyms.rds"))
#' 
#' }
#' 
#' 
#' 
#' ## standardise columns to match species_df
#' 
#' cooke_database <- cooke_database %>%
#'                   dplyr::mutate(source = "Cooke et al 2019") %>%
#'                   dplyr::rename(ecoregion_code = eco) %>%
#'                   merge(cooke_species[c("binomial", "tsn", "accepted_name", "class")])
#' 
#' 
#' # Subset the cooke database by country if subsetting
#' 
#' if (!is.na(country)) {
#'   
#'  full_cooke_database <- cooke_database
#'   
#'  country_ecoregions <- ecoregion_subset$ECO_ID
#'   
#'  cooke_database <- full_cooke_database[full_cooke_database$ecoregion_code %in% country_ecoregions, ]
#'   
#' }
#' 
#' 
#' # TEMPORARY CODE - subset data to make it manageable ----
#' 
#' #' TODO: Remove the two lines below and run on entire dataset when hpc available
#' #' 
#' 
#' ## Save full databases
#' 
#' # cooke_database_all <- cooke_database
#' # wildfinder_database_all <- wildfinder_database
#' # 
#' # ## Take a random sample subset
#' # 
#' # wildfinder_database <- wildfinder_database_all[sample(nrow(wildfinder_database_all), 10000), ]
#' # cooke_database <- cooke_database_all[sample(nrow(cooke_database_all), 10000), ]
#' 
#' 
#' # Fill in details so Cooke has the same columns and data as wildfinder
#' 
#' cooke_database <- cooke_database %>%
#'                   merge(wildfinder_database[c("wwf_species_id","common_name","genus", 
#'                                               "species", "tsn")], by = "tsn", all.x = TRUE) %>%
#'                   dplyr::select("binomial", "wwf_species_id", "common_name", "ecoregion_code",
#'                                 "class", "genus", "species", "source", "tsn", "accepted_name" ) %>%
#'                   distinct(.) 
#' 
#' 
#' # Merge databases ----
#' 
#' # Combine Cooke and Wildfinder
#' 
#' merged_databases <- rbind(wildfinder_database, cooke_database)
#' 
#' # Check how many species
#' 
#' message("There are ", length(unique(merged_databases$tsn)), 
#'         " species in the merged database")
#' 
#' # Clean up large objects no longer needed ----
#' 
#' # Add red list status and year to the merged databases.
#' 
#' if (!is.na(country)) {
#' 
#' species_data_with_sources <- merged_databases %>%
#'                              merge(henriques_redlist_database, by = "tsn", all = FALSE) %>% 
#'                              mutate(source = coalesce(source, redlist_source)) %>%
#'                              mutate(binomial = coalesce(binomial.x, binomial.y)) %>%
#'                              mutate(accepted_name = coalesce(accepted_name.x, 
#'                                                              accepted_name.y)) %>%
#'                              mutate(common_name = coalesce(common_name.x, 
#'                                                            common_name.y)) %>%
#'                              dplyr::select(-c(binomial.x, binomial.y, 
#'                                               accepted_name.x, accepted_name.y, 
#'                                                common_name.x, common_name.y, class.y))
#' 
#' } else {
#'   
#' species_data_with_sources <- merged_databases %>%
#'                              merge(henriques_redlist_database, by = "tsn", all = TRUE) %>% 
#'                              mutate(source = coalesce(source, redlist_source)) %>%
#'                              mutate(binomial = coalesce(binomial.x, binomial.y)) %>%
#'                              mutate(accepted_name = coalesce(accepted_name.x, 
#'                                                                 accepted_name.y)) %>%
#'                              mutate(common_name = coalesce(common_name.x, 
#'                                                               common_name.y)) %>%
#'                              dplyr::select(-c(binomial.x, binomial.y, 
#'                                                  accepted_name.x, accepted_name.y, 
#'                                                  common_name.x, common_name.y, class.y))
#'   
#' }
#' 
#' message("There are ", length(unique(species_data_with_sources$tsn)), 
#'         " species in the species data with sources dataframe")
#' 
#' # Consolidate duplicates 
#' 
#' species_data <- species_data_with_sources %>%
#'                 dplyr::select(-c("source", "redlist_source", "binomial")) %>%
#'                 distinct(.) %>%
#'                 dplyr::select(tsn, accepted_name, everything()) %>%
#'                 rename(accepted_binomial = accepted_name) %>%
#'                 rename(class = class.x)
#' 
#' message("There are ", length(unique(species_data$tsn)), 
#'         " species in the species data dataframe")
#' 
#' if(save_outputs == "yes"){
#'   
#'   saveRDS(species_data, file.path(outputs, "species_data_incomplete.rds"))
#' 
#' }
#' 
#' # Clean up large objects no longer needed ----
#' 
#' rm(merged_databases, tables, wildfinder_database, wildfinder_species, scientific_names, 
#'    henriques_redlist_database, henriques_species, group_tables, amphibians,
#'    birds, mammals, cooke_database, cooke_species)
#' 
#' 
#' 
#' # Species data checkpoint ----
#' 
#' ## Add the new ecoregions to our species data
#' 
#' if (("species_data.rds" %in% list.files(outputs))) { 
#'   
#' species_data <- readRDS(file.path(outputs, "species_data.rds"))
#' 
#' } else {
#' 
#' species_data <- species_data %>%
#'                 merge(iucn_rangemap_database[c( "tsn", "ecoregion_code", 
#'                                                 "eco_objectid")], 
#'                       by = "tsn", all = TRUE) %>%
#'                 mutate(ecoregion_code = coalesce(ecoregion_code.x, 
#'                                                  ecoregion_code.y)) %>%
#'                 dplyr::select(-c(ecoregion_code.x, ecoregion_code.y)) %>%
#'                 distinct(.)
#' 
#' if (save_outputs == "yes") {
#' 
#'   saveRDS(species_data, file.path(outputs,"species_data_2.rds", sep = ""))
#' 
#'   }
#' 
#' }
#' 
#' # Check species_data for gaps ----
#' 
#' # What species are extinct, or extinct in the wild?
#' 
#' extinct_species <- species_data %>%
#'                    dplyr::select(tsn, accepted_binomial, redlist_status) %>%
#'                    filter(redlist_status == "EX" | redlist_status == "EW")  %>%
#'                    distinct(.) %>%
#'                    mutate(extinct = 1)
#' 
#' # What species do we not have accepted binomials for?
#' 
#' #' TODO: Is there a way - or do we need - to pull the other details for these species?
#' 
#' species_without_binomials <- species_data %>%
#'                              filter(is.na(accepted_binomial)) %>%
#'                              select(tsn) %>%
#'                              distinct(.)
#' 
#' message("There are ", length(unique(species_without_binomials$tsn)), 
#'         " species without binomials in the species data dataframe")
#' 
#' # Remove all the no binomial species for the time being
#' 
#' species_data <- species_data %>%
#'                 filter(accepted_binomial != is.na(accepted_binomial))
#' 
#' message("After removing missing binomials, there are ", 
#'         length(unique(species_data$tsn)), 
#'         " species with binomials in the species data dataframe")
#' 
#' # What species do we not have ecoregions for? Note this will produce nothing
#' # if subsetting by country
#' 
#' 
#' species_without_ecoregions <- species_data %>%
#'                               filter(is.na(ecoregion_code)) %>%
#'                               dplyr::select(tsn, accepted_binomial, ecoregion_code) %>%
#'                               distinct(.) %>%
#'                               mutate(no_ecoregion = 1)
#' 
#' message("There are ", 
#'         length(unique(species_without_ecoregions$tsn)), 
#'         " species without ecoregions in the species data dataframe")
#' 
#' missing_ecoregions_species_data <- species_data %>%
#'                                    filter(is.na(ecoregion_code))
#' 
#' # What species do we not have ecoregion object id for?
#' 
#' species_without_eco_objectid <- species_data %>%
#'                                 filter(is.na(eco_objectid)) %>%
#'                                 filter(!is.na(redlist_status)) %>%
#'                                 dplyr::select(tsn, accepted_binomial, 
#'                                               eco_objectid) %>%
#'                                 distinct(.) %>%
#'                                 mutate(no_objectid = 1)
#' 
#' message("There are ", 
#'         length(unique(species_without_eco_objectid$tsn)), 
#'         " species without ecoregions in the species data dataframe")
#' 
#' 
#' # What species do we not have a redlist status in any year for?
#' 
#' species_without_redlist_status <- species_data %>%
#'                                   dplyr::select(-c(redlist_assessment_year,
#'                                                 genus, species, ecoregion_code,
#'                                                 eco_objectid)) %>%
#'                                   distinct(.) %>%
#'                                   group_by(tsn) %>%
#'                                   filter(all(is.na(redlist_status))) %>%
#'                                   mutate(no_redlist_status = 1) 
#' 
#' message("There are ", 
#'         length(unique(species_without_redlist_status$tsn)), 
#'         " species without a redlist status in the species data dataframe")
#' 
#' 
#' missing_redlist_species_data <- species_data %>%
#'                                 group_by(tsn) %>%
#'                                 filter(all(is.na(redlist_status)))
#' 
#'                               
#' 
#' # if(save_outputs == "yes") {
#' # 
#' #   write_csv(species_without_redlist_status, paste(outputs, "/", date, "_species_without_redlist_status.csv", sep = ""))
#' #   saveRDS(species_without_redlist_status, paste(outputs,"/", date, "_species_without_redlist_status.rds", sep = ""))
#' # 
#' # }
#' # 
#' # # What species do we have at least one redlist status for?
#' # 
#' # species_with_redlist_status <- species_data %>%
#' #                                dplyr::select(-redlist_assessment_year,
#' #                                               genus, species) %>%
#' #                                distinct(.) %>%
#' #                                group_by(tsn) %>%
#' #                                filter(!is.na(redlist_status)) %>%
#' #                                dplyr::select(-redlist_status) %>%
#' #                                distinct(.)
#' 
#' # Double check we haven't incorrectly detected no redlist status (should be no
#' # overlap between the two groups)
#' 
#' # overlap <- species_without_redlist_status$tsn %in% species_with_redlist_status$tsn
#' # any(overlap == TRUE) # The correct output to console should be FALSE (no overlap)
#' 
#' 
#' # What species are included in our database, and do we have the necessary data for them?
#' # prev_species_in_database <- species_in_database
#' 
#' # species_in_database <- species_data %>%
#' #                        select(tsn, accepted_binomial) %>%
#' #                        distinct(.) %>%
#' #                        merge(extinct_species[
#' #                        c("tsn", "extinct")], by = "tsn", all = TRUE) 
#' # 
#' species_with_incomplete_data <- species_without_redlist_status %>%
#'                                 merge(species_without_ecoregions[
#'                                 c("tsn", "no_ecoregion")], by = "tsn", all = TRUE) %>%
#'                                 merge(species_without_eco_objectid[
#'                                 c("tsn", "no_objectid")], by = "tsn", all = TRUE) %>%
#'                                 mutate(incomplete_data = ifelse(no_redlist_status == 1|
#'                                                                 no_ecoregion == 1|
#'                                                                 no_objectid == 1|
#'                                                                 is.na(accepted_binomial),
#'                                                                       1, 0))
#' 
#'                               
#' 
#' # write.csv(species_with_incomplete_data, file.path(outputs, 
#' #                                                  "200520_species_with_incomplete_data.csv"))
#' 
#' # species_with_complete_data <- species_in_database %>%
#' #                               filter(is.na(incomplete_data))
#' 
#' print(paste("The database contains", length(unique(species_data$tsn)),
#'             "species", sep = " "))
#' 
#' print(paste(length(unique(species_with_incomplete_data$tsn)),
#'             "of the species have incomplete information", sep = " "))
#' 
#' print(paste(sum(species_with_incomplete_data$no_redlist_status, na.rm = TRUE), 
#'             " have no redlist status, ", 
#'             sum(species_with_incomplete_data$no_ecoregion, na.rm = TRUE),
#'             " have no ecoregion, and ",
#'             sum(species_with_incomplete_data$no_objectid , na.rm = TRUE),
#'             " have no ecoregion object id", sep = ""))
#' 
#' 
#' # Find missing Red List Status ----
#' 
#' # Find red list status for the species we are missing info for
#' 
#' #species_with_incomplete_data <- read.csv("N:\\Quantitative-Ecology\\Simone\\extinction_test\\outputs\\2020-05-19_output_files\\200520_species_with_incomplete_data.csv")
#' #species_with_incomplete_data <- species_with_incomplete_data %>%
#'                                 #select(-X)
#' # Add your iucn API token to you environment
#' 
#' # file.edit("~/.Renviron")
#' 
#' # In the file that opens, add the text below, save that file, then close and reopen R
#' 
#' # IUCN_REDLIST_KEY='6d46333255cbb1843dd8f6984e45f23a49e8011585b1f128c3a46d17c6a8b5ae'
#' 
#' rl_citation()
#' 
#' 
#' # test <- get_redlist_data('Thylacinus cynocephalus')
#' 
#' # Create list of search names
#' 
#' # species_list <- c('Thylacinus cynocephalus',"Fratercula arctica")
#' 
#' if (!("iucn_scraped_data.rds" %in% list.files(outputs))) {
#'   
#' species_list <- as.character(species_without_redlist_status$accepted_binomial)
#' 
#' species_list <- unique(species_list)
#' 
#' # Scrape the data in batches because otherwise the iucn website throws an error
#' 
#' species_list_1 <- species_list[1:970]
#' 
#' iucn_scraped_data_1 <- lapply(species_list_1, get_redlist_data)
#' 
#' species_list_2 <- species_list[971:1940]
#' 
#' iucn_scraped_data_2 <- lapply(species_list_2, get_redlist_data)
#' 
#' species_list_3 <- species_list[1941:length(species_list)]
#' 
#' iucn_scraped_data_3 <- lapply(species_list_3, get_redlist_data)
#'  
#' # species_list_4 <- species_list[2911:length(species_list)]
#' 
#' #iucn_scraped_data_4 <- lapply(species_list_4, get_redlist_data)
#' 
#' iucn_scraped_data <- c(iucn_scraped_data_1, iucn_scraped_data_2,
#'                        iucn_scraped_data_3)
#' 
#' saveRDS(iucn_scraped_data, file.path(outputs,"iucn_scraped_data.rds"))
#' 
#' rm(iucn_scraped_data_1, iucn_scraped_data_2, iucn_scraped_data_3)
#' 
#' } else {
#' 
#' iucn_scraped_data <- readRDS(file.path(outputs,"iucn_scraped_data.rds"))
#' 
#' }
#' 
#' # Turn it into a nice data frame (website downloads as a list)
#' 
#' iucn_scraped_data_list <- list()
#' no_data <- list()
#' 
#' for(i in seq_along(iucn_scraped_data)) {
#'   
#' if (length(iucn_scraped_data[[i]][[2]]) == 30) {
#'   
#' accepted_binomial <- iucn_scraped_data[[i]][[1]]
#' redlist_status <- as.character(iucn_scraped_data[[i]][[2]][13][[1]])
#' class <- as.character(iucn_scraped_data[[i]][[2]][5][[1]])
#' redlist_assessment_year <- iucn_scraped_data[[i]][[2]][11][[1]]
#' source <- "IUCN"
#' iucn_scraped_data_list[[i]] <- as.data.frame(cbind(accepted_binomial, redlist_status, 
#'                                              published_year, class, source))
#' 
#' } else if (length(iucn_scraped_data[[i]][[2]]) == 1) {
#'   
#' accepted_binomial <- iucn_scraped_data[[i]][2][[1]]
#' redlist_status <- NA
#' redlist_assessment_year <- NA
#' class <- NA
#' source <- "IUCN"
#' no_data[[i]] <- as.data.frame(cbind(accepted_binomial, redlist_status, 
#'                                              published_year, class, source))
#'   }
#' }
#' 
#' iucn_scraped_dataframe <- do.call(rbind, iucn_scraped_data_list)
#' 
#' iucn_scraped_dataframe$redlist_status <- as.character(iucn_scraped_dataframe$redlist_status)
#' iucn_scraped_dataframe$accepted_binomial <- as.character(iucn_scraped_dataframe$accepted_binomial)
#' 
#' names(iucn_scraped_dataframe) <- c("accepted_binomial", 
#'                                    "redlist_status",
#'                                    "redlist_assessment_year",
#'                                    "class",
#'                                    "source")
#' # species we had no luck with
#' 
#' no_data <- do.call(rbind, no_data)
#' 
#' # Add new red list status back into the species data
#' 
#' correct_order <- names(species_data)
#' 
#' found_redlist_species_data <- missing_redlist_species_data %>%
#'                               merge(iucn_scraped_dataframe[c("accepted_binomial", 
#'                                                              "redlist_status", 
#'                                                              "redlist_assessment_year")], 
#'                                     by = "accepted_binomial", all = FALSE) %>%
#'                               dplyr::select(-redlist_assessment_year.x, redlist_status.x) %>%
#'                               rename(redlist_status = redlist_status.y,
#'                                      redlist_assessment_year = redlist_assessment_year.y) %>%
#'                               select(correct_order) %>%
#'                               distinct(.)
#' 
#' 
#' # Save a version of the old species data just in case
#' 
#' species_data_all <- species_data
#' 
#' # Add in the new redlist data
#' 
#' species_data <- rbind(species_data, found_redlist_species_data)
#' 
#' # Remove the rest of the missing redlist species
#' 
#' species_data <- species_data %>%
#'                 filter(!is.na(redlist_status))
#' 
#' 
#' print(paste(length(unique(species_data$tsn)),
#'             "species with redlist status", sep = " "))
#' 
#' 
#' 

#' 
#' # Add in the new extinct species data
#' 
#' species_data <- rbind(species_data, extinct_species_data)
#' 
#' # Find missing ecoregions for species without them ----
#' 
#' if (!("species_ecoregions_found.rds" %in% list.files(outputs))) {
#' 
#' species_no_ecoregion_names <- species_without_ecoregions %>%
#'                               select(accepted_binomial) %>%
#'                               unique(.) %>%
#'                               pull(.)
#' 
#' species_no_ecoregion_names <- word(species_no_ecoregion_names, 1, 2)
#' 
#' species_ecoregions_found <- get_gbif_data(species_no_ecoregion_names,
#'                                           20, ecoregion_map)
#' 
#' if (save_outputs == "yes") {
#'   
#'   saveRDS(species_ecoregions_found, 
#'           file.path(outputs, "species_ecoregions_found.rds"))
#'   
#'   } 
#' 
#' } else {
#'   
#'   species_ecoregions_found <- readRDS(file.path(outputs, 
#'                                                 'species_ecoregions_found.rds'))
#' }
#' 
#' species_with_new_ecoregions <- unique(species_ecoregions_found[[1]])
#' 
#' names(species_with_new_ecoregions) <- c("accepted_binomial", "ecoregion_code",
#'                                         "ecoregion_name", "eco_objectid")
#' 
#' # Add new red list status back into the species data
#' 
#' correct_order <- names(species_data)
#' 
#' head(species_with_new_ecoregions)
#' 
#' found_ecoregion_species_data <- missing_ecoregions_species_data %>%
#'                                 merge(species_with_new_ecoregions[c("accepted_binomial", 
#'                                                                "ecoregion_code", 
#'                                                                "eco_objectid")], 
#'                                       by = "accepted_binomial", all = FALSE) %>%
#'                                 dplyr::select(- eco_objectid.x, ecoregion_code.x) %>%
#'                                 rename(eco_objectid = eco_objectid.y,
#'                                        ecoregion_code = ecoregion_code.y) %>%
#'                                 select(correct_order) %>%
#'                                 distinct(.)
#' 
#' 
#' # Save a version of the old species data just in case
#' 
#' species_data_all <- species_data
#' 
#' # Add in the new redlist data
#' 
#' species_data <- rbind(species_data, found_ecoregion_species_data)
#' 
#' # Remove the rest of the missing redlist species
#' 
#' species_data <- species_data %>%
#'                 filter(!is.na(ecoregion_code))
#' 
#' 
#' print(paste(length(unique(species_data$tsn)),
#'             "species with complete redlist and ecoregion data", sep = " "))
#' 
#' if (save_outputs == "yes") {
#'   
#'   saveRDS(species_data, 
#'           file.path(outputs, "species_data.rds"))
#'   
#' } 
#' 
#' # Calculate summary statistics ----
#' 
#' # Get the number of species in each ecoregion
#' 
#' 
#' species_by_ecoregion <- species_data %>%
#'                         group_by(ecoregion_code) %>%
#'                         summarize(n_distinct(tsn))
#' 
#' 
#' names(species_by_ecoregion) <- c("ecoregion_code", "number_of_species")
#' 
#' # Get the number of species in each redlist category
#' 
#' species_by_redlist_status <- species_data %>%
#'                              group_by(redlist_status) %>%
#'                              summarize(n_distinct(tsn))
#' 
#' 
#' 
#' # Get number of extinct species (grouping by ecoregion doesn't work well yet bc
#' # most extinct species haven't been assigned an ecoregion - TBD)
#' 
#' number_of_extinct_species <- species_data %>%
#'                    filter(redlist_status == "EX")  %>%
#'                    group_by(ecoregion_code) %>%
#'                    summarise(n_distinct(tsn))
#' 
#' names(number_of_extinct_species) <- c("ecoregion_code", "number_of_species_extinct")
#' 
#' 
#' number_of_extinct_wild_species <- species_data %>%
#'                         filter(redlist_status == "EX" | redlist_status == "EW")  %>%
#'                         group_by(ecoregion_code) %>%
#'                         summarise(n_distinct(tsn))
#' 
#' names(number_of_extinct_wild_species) <- c("ecoregion_code", "number_of_species_extinct")
#' 
#' if(save_outputs == "yes") {
#' 
#'   write_csv(number_of_extinct_species, paste(outputs, "/", date, "_extinct_species.csv", sep = ""))
#'   saveRDS(number_of_extinct_species, paste(outputs,"/", date, "_extinct_species.rds", sep = ""))
#' 
#' }
#' 
#' proportion_extinct <- species_by_ecoregion %>%
#'                       merge(number_of_extinct_species,
#'                             by = "ecoregion_code", all = TRUE) %>%
#'                       dplyr::mutate(proportion_extinct =
#'                                     number_of_species_extinct/number_of_species)
#' 
#' 
#' # Visualise summary stats ----
#' 
#' extinction_map_data <- inner_join(ecoregion_map, proportion_extinct[
#'                         c("ecoregion_code", "number_of_species_extinct")],
#'                         by = c("ECO_ID" = "ecoregion_code"))
#' 
#' extinction_map <- ggplot(extinction_map_data) +
#'                   geom_sf(aes(fill = number_of_species_extinct)) +
#'                   scale_fill_viridis_c(trans = "sqrt", alpha = .4)
#' 
#' extinction_map
#' 
#' if(save_outputs == "yes") {
#'   
#'   ggsave(file.path(outputs, "extinctions_map.png", extinction_map))
#'   
#' }
#' 
#' dev.off()
#' 
#' 
#' 
#' species_map_data <- inner_join(ecoregion_map, species_by_ecoregion[
#'                     c("ecoregion_code", "number_of_species")],
#'                     by = c("ECO_ID" = "ecoregion_code"))
#' 
#' species_map <- ggplot(species_map_data) +
#'                geom_sf(aes(fill = number_of_species)) +
#'                scale_fill_viridis_c(trans = "sqrt", alpha = .4)
#' 
#' species_map
#' 
#' if(save_outputs == "yes") {
#'   
#'   ggsave(file.path(outputs, "species_richness_map.png", species_map))
#'   
#' }
#' 
