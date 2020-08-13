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

# Load packages ----

library(raster)
library(tidyverse)
library(stringr)
library(reshape2)
library(sf)
library(spData)
#library(taxize)
# devtools::install_github("ropensci/taxizedb") # version 0.1.9.9130
library(taxizedb)
library(functionaltraits)
library(broom)
library(rredlist)
library(rgbif)
library(rlist)

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

get_ecoregions <- function(rangemap, rangemap_directory_path, map) {
  
  names(range_map) <- c(toupper(names(range_map)[1:27]), names(range_map[28])) # Make column names consistent
  
  #range_map <- rename(range_map, geometry = GEOMETRY) # Convert geometry back tho
  
  # Get the ecoregions for each species
  
  ranges_ecoregions <- st_join(range_map, map, join = st_intersects)
  
  if (!is.null(ranges_ecoregions)) {
    
    print(paste("st_join complete for", 
                basename(rangemap_directory_path), sep = " "))
    
  }
  
  # Standardise the data so we can add it to the species data easily
  
  
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
  # %>%
  #   mutate(ecoregion_id = as.character(ecoregion_id))
  
  
  if(!is.null(species_w_ecoregions)) {
    
    print(paste("finished adding ecoregions to", 
                basename(rangemap_directory_path), sep = " "))
    
  }
  
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
        select(species, decimalLatitude, decimalLongitude ) %>%
        distinct(.) %>%
        filter(complete.cases(decimalLatitude, decimalLongitude))
      
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

# Load ecoregion data ----

## This will allow us to subset the data by country for development or analysis
#' TODO: Why are so many ecoregions missing countries? figure out best join method
#' TODO: Important - Add countries to species_data 

# Get the ecoregion map 

ecoregion_map_all <- st_read(paste(inputs,eco_version, sep = "/"))

# Pull out only required variables

ecoregion_map <- ecoregion_map_all %>% 
                 dplyr::select(ECO_ID, ECO_NAME, geometry, OBJECTID)

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

# Add countries to the ecoregion map

ecoregion_map <- ecoregion_map %>%
                 merge(ecoregion_country_df[c("ECO_ID", "CNTRY_NAME")],
                       by = "ECO_ID") 
#%>%
                 #filter(ECO_ID != 0) # Remove rock and ice because it gets allocated to multiple ocuntries

# If subsetting by country, get the ecoregion IDs you want to work with 

if (!is.na(country)) {
  
ecoregion_subset <- ecoregion_country_df %>%
                    filter(CNTRY_NAME == country) %>%
                    unique(.)

}

if (!is.na(country)) {
  
ecoregion_map <- ecoregion_map %>%
                 filter(CNTRY_NAME == country)
                   
  
}

# Intersect range maps with ecoregions (SLOW CODE) ----


# To fill in missing ecoregions

# If starting from here, load species data

# species_data <- readRDS(file.path(outputs, "species_data_incomplete.rds"))

intersect_ranges_w_ecoregions <- function(class_name, location, shapefile, data) {
  

if (!(paste(class_name, location, "iucn_range_map_species_with_ecoregions.rds", 
            sep = "_") %in% list.files(interim_outputs))) {
  
  range_directories <- list(file.path(inputs, paste("redlist", class_name,  
                                                   "range_maps", sep = "_")))
  
  iucn_rangemap_database <- list()
  
  for (i in seq_along(range_directories)) {
    
    iucn_rangemap_database[[i]] <- get_ecoregions(range_directories[[i]], 
                                                  shapefile)
    
    print(paste("processed", i, "of" , length(range_directories), 
                "group range maps", sep = " "))
    
  }
  
  iucn_rangemap_database <- do.call(rbind, iucn_rangemap_database)
  iucn_rangemap_database <- distinct(iucn_rangemap_database)
  
  saveRDS(iucn_rangemap_database, 
          file = file.path(interim_outputs, paste(class_name,
          location,"iucn_range_map_species_with_ecoregions.rds",
          sep = "_")))

} else {
  
  iucn_rangemap_database <- readRDS(file.path(interim_outputs, 
                                    paste(class_name,
                                          location,
                                    "iucn_range_map_species_with_ecoregions.rds",
                                    sep = "_")))
}

if (!is.na(country)) {
  
  iucn_rangemap_database <- iucn_rangemap_database[
                            iucn_rangemap_database$ecoregion_id %in% 
                            ecoregion_subset$ECO_ID, ]
}

# Get TSNs for the species from the IUCN range maps

if (!(file.path(outputs, paste(class_name, location, "iucn_rangemap_synonyms.rds",
                      sep = "_")) %in% list.files(outputs))) {
  
  # pull out unique species names only
  
  iucn_rangemap_species <- iucn_rangemap_database %>%
    dplyr::select(binomial) 
  
  iucn_rangemap_species$binomial <- as.character(iucn_rangemap_species$binomial)
  
  iucn_rangemap_species <- unname(unlist(iucn_rangemap_species[,1]))
  
  iucn_rangemap_species <- iucn_rangemap_species[!is.na(iucn_rangemap_species)]
  
  iucn_rangemap_species <- unique(iucn_rangemap_species)
  
  # search for synonyms and add taxonomic serial number (slow)
  
  iucn_rangemap_species <- find_synonyms(iucn_rangemap_species)
  
  # save synonyms for this groups of species
  
  saveRDS(iucn_rangemap_species, file.path(outputs, 
                                           paste(class_name,
                                                 location,
                                                 "iucn_rangemap_synonyms.rds",
                                                 sep = "_")))
  
  # add the synonyms and taxonomic info to the ecoregion-species data
  
  iucn_rangemap_database <-  iucn_rangemap_database %>%
    merge(iucn_rangemap_species[c("binomial", "tsn", 
                                  "found","class")], 
          by = "binomial")
  
  # save the species with ecoregions data
  
  saveRDS(iucn_rangemap_database, file.path(interim_outputs, 
                                            paste(class_name,
                                                  location,
                                            "iucn_rangemap_database.rds", 
                                            sep = "_")))
  
} else {
  
  # Read in the synonyms previously stored
  
  iucn_rangemap_species <- readRDS(file.path(outputs, 
                                             paste(class_name,
                                                   location,
                                                   "iucn_rangemap_synonyms.rds",
                                                    sep = "_")))
  
  # Merge with ecoregion data
  
  iucn_rangemap_database <-  iucn_rangemap_database %>%
    merge(iucn_rangemap_species[c("binomial", "tsn", 
                                  "found", "class")], 
          by = "binomial")
  
  saveRDS(iucn_rangemap_database, file.path(interim_outputs, paste(class_name,
                                                           location,
                                            "iucn_rangemap_database.rds", 
                                            sep = "_")))
  
}

# Load red list status data ----

file_path <- file.path(inputs, "redlist_sergio")
file_names <- list.files(file_path)
file <- file_names[str_detect(file_names, class_name)]
files <- file.path(file_path, file)
table <- read_csv(files)
#names(tables) <- str_remove(file_names, ".csv")

# Add a column so we can identify what taxa each table contains + the source
# reference

if (class_name == "amphibian") {
  
matched_class <- "Amphibia"

table <- table %>% 
         mutate(class = matched_class) %>%
         mutate(redlist_source = "Henriques etal 2020") %>% 
         dplyr::mutate(binomial = paste(Genus, Species, sep = " ")) %>%
         dplyr::select(-c(Genus, Species)) %>%
         set_names(c("2004", "2008", "class", "redlist_source", "binomial")) %>%
         dplyr::select("binomial","class","2004", "2008", "redlist_source")

}

if (class_name == "mammal") {
  
  matched_class <- "Mammalia"
  
  table <- table %>% 
           mutate(class = matched_class) %>%
           mutate(redlist_source = "Henriques etal 2020") %>% 
           rename(binomial = Name) %>%
           set_names(c("binomial","1996", "2008", "class", "redlist_source")) %>%
           dplyr::select("binomial","class","1996", "2008", "redlist_source")
  
}

if (class_name == "reptile") {
  
  matched_class <- "Reptilia"
  
  table <- table %>% 
           mutate(class = matched_class) %>%
           mutate(redlist_source = "Henriques etal 2020") %>% 
           rename(binomial = Name) %>%
           set_names(c("binomial","1996", "2008", "class", "redlist_source")) %>%
           dplyr::select("binomial","class","1996", "2008", "redlist_source")
  
}

# Get synonyms for henriques species ----

henriques_species <- table %>%
                     dplyr::select(binomial) 

henriques_species$binomial <- as.character(henriques_species$binomial)
henriques_species <- unname(unlist(henriques_species[,1]))
henriques_species <- henriques_species[!is.na(henriques_species)]
henriques_species <- unique(henriques_species)

if (!(paste(class_name,
            "henriques_synonyms.rds",
            sep = "_") %in% list.files(interim_outputs))) { 
  
  henriques_species <- find_synonyms(henriques_species)
  
  saveRDS(henriques_species, file.path(outputs, 
                                       paste(class_name,
                                             "henriques_synonyms.rds",
                                             sep = "_")))
  
} else {
  
  
  henriques_species <- readRDS(file.path(outputs, 
                                         paste(class_name,
                                               "henriques_synonyms.rds",
                                               sep = "_")))
  
}

# Order columns, drop group variable and melt in to long format

redlist_status_data <- table %>%
                       select(-class) %>%
                       melt(.,id.vars = c("binomial",  
                                          "redlist_source"),
                             value.name = "redlist_status") %>%
                       rename(redlist_assessment_year = variable) %>% 
                       merge(henriques_species[c("binomial", "tsn", 
                                                 "accepted_name", 
                                                 "common_name", "class")]) %>%
                       rename(accepted_binomial = accepted_name)
    

# Merge redlist and ecoregion data ----

ecoregion_redlist_data <- redlist_status_data %>%
                          merge(iucn_rangemap_database[c("tsn",
                                                         "ecoregion_id",
                                                         "eco_objectid")], 
                                by = "tsn", all = TRUE) %>%
                          select(tsn, accepted_binomial, 
                                 class, redlist_assessment_year,
                                 redlist_status, common_name, eco_objectid,
                                 ecoregion_id) %>% 
                          merge(data[c("ECO_ID", "CNTRY_NAME")], 
                                by.x = "ecoregion_id", by.y =  "ECO_ID" ) %>%
                          rename(country = "CNTRY_NAME")


saveRDS(ecoregion_redlist_data, file.path(interim_outputs, 
                                     paste(class_name,
                                           location,
                                           "species_data_1.rds",
                                           sep = "_")))

return(ecoregion_redlist_data)

}

# Intersect ranges and ecoregions for each class
# (can put this in a loop but takes a long time to run each,
# probably more pragmatic to split them up or parallelise)

amphibians <- intersect_ranges_w_ecoregions("amphibian", location, 
                                            ecoregion_map, ecoregion_country_df)

mammals <- intersect_ranges_w_ecoregions("mammal", location, 
                                         ecoregion_map, ecoregion_country_df)

#reptiles <- intersect_ranges_w_ecoregions("reptile", ecoregion_map)

# Combine into one file and save

species_data <- rbind(amphibians, mammals)

saveRDS(species_data, file.path(interim_outputs, 
                                paste(location, "species_data_1.rds",
                                      sep = "_")))

summarise_species_data <- function(data, number) {
  
    out1 <- data %>%
            group_by(ecoregion_id, redlist_status) %>%
            summarise(spp_number_w_status = n())
    
    out2 <- data %>%
            group_by(ecoregion_id) %>%
            summarise(number_spp_in_ecoregion = 
                        n_distinct(tsn))
    
    class_ecoregion_summary <- out1 %>%
           merge(out2, by = "ecoregion_id") %>%
           mutate(class = data$class[1])
    
    class_redlist_summary <- data %>%
                             group_by(redlist_status) %>%
                             summarise(redlist_count = n()) %>%
                             mutate(class = data$class[1])
    
    out <- list(class_ecoregion_summary, class_redlist_summary)
    
    saveRDS(out, file.path(outputs, paste(number,
                                          data$class[1],location,
                                          "data_summary.rds",
                                          sep = "_")))
    
}

amphibian_summary <- summarise_species_data(amphibians, 1)
mammal_summary <- summarise_species_data(mammals, 1)


# old_sp_data <- readRDS(file.path(inputs, "deakin_species_data", "species_data.rds"))
# group_tables <- list()
# 
# for (i in seq_along(tables)) {
#   
#   group_tables[[i]] <- tables[[i]] %>% 
#     mutate(group = names(tables[i])) %>%
#     mutate(redlist_source = "Henriques etal 2020")
#   
#   
# }

# Standardise column names etc. because they are inconsistent between the 
# different taxa tables

# amphibians <- group_tables[[1]]
# 
# amphibians <- amphibians %>% 
#   dplyr::mutate(binomial = paste(Genus, Species, sep = " ")) %>%
#   dplyr::select(-c(Genus, Species)) %>%
#   set_names(c("2004", "2008", "group", "redlist_source", "binomial")) %>%
#   dplyr::select("binomial","group","2004", "2008", "redlist_source") 

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

henriques_redlist_database <- do.call(bind_rows, group_tables)

# Get Henriques synonyms ----

henriques_species <- henriques_redlist_database %>%
  dplyr::select(binomial) 

henriques_species$binomial <- as.character(henriques_species$binomial)
henriques_species <- unname(unlist(henriques_species[,1]))
henriques_species <- henriques_species[!is.na(henriques_species)]

if (!("henrique_species_synonyms.rds" %in% list.files(interim_outputs))) { 
  
  henriques_species <- find_synonyms(henriques_species)
  
  saveRDS(henriques_species, file.path(interim_outputs, "henrique_species_synonyms.rds"))
  
} else {
  
  
  henriques_species <- readRDS(file.path(interim_outputs, "henrique_species_synonyms.rds"))
  
}

# Order columns, drop group variable and melt in to long format

henriques_redlist_database <- henriques_redlist_database %>%
  dplyr::select("binomial", "1988", "1994", "1996", "2000", 
                "2004", "2008", "2012", "2016", "redlist_source") %>%
  melt(.,id.vars = c("binomial", "redlist_source"),
       value.name = "redlist_status") %>%
  rename(redlist_assessment_year = variable) %>% 
  merge(henriques_species[c("binomial", "tsn", "accepted_name", 
                            "common_name", "class")])


# Load Wildfinder database ----

# Load tables from WildFinder database (converted into .xlsx files from .mdb 
# database)

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

# Get Wildfinder synonyms (SLOW code) ----

if(!("wildfinder_species_synonyms.rds" %in% list.files(interim_outputs))) { 
  
wildfinder_species <- find_synonyms( wildfinder_species)

saveRDS(wildfinder_species, file.path(interim_outputs,"wildfinder_species_synonyms.rds"))

} else {


wildfinder_species <- readRDS(file.path(interim_outputs, 
                                        "wildfinder_species_synonyms.rds"))

}

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
                        merge(wildfinder_species[c("binomial", "tsn", 
                                                   "accepted_name", "class")])


# Remove un-needed tables

rm(tables)

# Clear memory by removing large objects

rm(common_names)

# Subset the wildfinder database by country if working with subset

if (!is.na(country)) {
  
  full_WF_database <- wildfinder_database
  
  country_ecoregions <- ecoregion_subset$ECO_ID
  
  wildfinder_database <- full_WF_database[full_WF_database$ecoregion_code %in% country_ecoregions, ]
  
}

# Load Cooke database ----

cooke_database <- read.csv(paste(inputs,"/cooke_database/species_site.csv", sep = ""))

cooke_species <- cooke_database %>%
                 dplyr::select(binomial) 

cooke_species$binomial <- as.character(cooke_species$binomial)
cooke_species <- cooke_species[,1]

# Get Cooke synonyms (SLOW code) ----

if(!("cooke_species_synonyms.rds" %in% list.files(interim_outputs))) { 
  
  cooke_species <- find_synonyms(cooke_species)
  
  saveRDS(cooke_species, file.path(interim_outputs, "cooke_species_synonyms.rds"))
  
} else {
  
  cooke_species <- readRDS(file.path(interim_outputs, "cooke_species_synonyms.rds"))

}



## standardise columns to match species_df

cooke_database <- cooke_database %>%
                  dplyr::mutate(source = "Cooke et al 2019") %>%
                  dplyr::rename(ecoregion_code = eco) %>%
                  merge(cooke_species[c("binomial", "tsn", "accepted_name", "class")])


# Subset the cooke database by country if subsetting

if (!is.na(country)) {
  
 full_cooke_database <- cooke_database
  
 country_ecoregions <- ecoregion_subset$ECO_ID
  
 cooke_database <- full_cooke_database[full_cooke_database$ecoregion_code %in% country_ecoregions, ]
  
}


# TEMPORARY CODE - subset data to make it manageable ----

#' TODO: Remove the two lines below and run on entire dataset when hpc available
#' 

## Save full databases

# cooke_database_all <- cooke_database
# wildfinder_database_all <- wildfinder_database
# 
# ## Take a random sample subset
# 
# wildfinder_database <- wildfinder_database_all[sample(nrow(wildfinder_database_all), 10000), ]
# cooke_database <- cooke_database_all[sample(nrow(cooke_database_all), 10000), ]


# Fill in details so Cooke has the same columns and data as wildfinder

cooke_database <- cooke_database %>%
                  merge(wildfinder_database[c("wwf_species_id","common_name","genus", 
                                              "species", "tsn")], by = "tsn", all.x = TRUE) %>%
                  dplyr::select("binomial", "wwf_species_id", "common_name", "ecoregion_code",
                                "class", "genus", "species", "source", "tsn", "accepted_name" ) %>%
                  distinct(.) 


# Merge databases ----

# Combine Cooke and Wildfinder

merged_databases <- rbind(wildfinder_database, cooke_database)

# Check how many species

message("There are ", length(unique(merged_databases$tsn)), 
        " species in the merged database")

# Clean up large objects no longer needed ----

# Add red list status and year to the merged databases.

if (!is.na(country)) {

species_data_with_sources <- merged_databases %>%
                             merge(henriques_redlist_database, by = "tsn", all = FALSE) %>% 
                             mutate(source = coalesce(source, redlist_source)) %>%
                             mutate(binomial = coalesce(binomial.x, binomial.y)) %>%
                             mutate(accepted_name = coalesce(accepted_name.x, 
                                                             accepted_name.y)) %>%
                             mutate(common_name = coalesce(common_name.x, 
                                                           common_name.y)) %>%
                             dplyr::select(-c(binomial.x, binomial.y, 
                                              accepted_name.x, accepted_name.y, 
                                               common_name.x, common_name.y, class.y))

} else {
  
species_data_with_sources <- merged_databases %>%
                             merge(henriques_redlist_database, by = "tsn", all = TRUE) %>% 
                             mutate(source = coalesce(source, redlist_source)) %>%
                             mutate(binomial = coalesce(binomial.x, binomial.y)) %>%
                             mutate(accepted_name = coalesce(accepted_name.x, 
                                                                accepted_name.y)) %>%
                             mutate(common_name = coalesce(common_name.x, 
                                                              common_name.y)) %>%
                             dplyr::select(-c(binomial.x, binomial.y, 
                                                 accepted_name.x, accepted_name.y, 
                                                 common_name.x, common_name.y, class.y))
  
}

message("There are ", length(unique(species_data_with_sources$tsn)), 
        " species in the species data with sources dataframe")

# Consolidate duplicates 

species_data <- species_data_with_sources %>%
                dplyr::select(-c("source", "redlist_source", "binomial")) %>%
                distinct(.) %>%
                dplyr::select(tsn, accepted_name, everything()) %>%
                rename(accepted_binomial = accepted_name) %>%
                rename(class = class.x)

message("There are ", length(unique(species_data$tsn)), 
        " species in the species data dataframe")

if(save_outputs == "yes"){
  
  saveRDS(species_data, file.path(outputs, "species_data_incomplete.rds"))

}

# Clean up large objects no longer needed ----

rm(merged_databases, tables, wildfinder_database, wildfinder_species, scientific_names, 
   henriques_redlist_database, henriques_species, group_tables, amphibians,
   birds, mammals, cooke_database, cooke_species)



# Species data checkpoint ----

## Add the new ecoregions to our species data

if (("species_data.rds" %in% list.files(outputs))) { 
  
species_data <- readRDS(file.path(outputs, "species_data.rds"))

} else {

species_data <- species_data %>%
                merge(iucn_rangemap_database[c( "tsn", "ecoregion_code", 
                                                "eco_objectid")], 
                      by = "tsn", all = TRUE) %>%
                mutate(ecoregion_code = coalesce(ecoregion_code.x, 
                                                 ecoregion_code.y)) %>%
                dplyr::select(-c(ecoregion_code.x, ecoregion_code.y)) %>%
                distinct(.)

if (save_outputs == "yes") {

  saveRDS(species_data, file.path(outputs,"species_data_2.rds", sep = ""))

  }

}

# Check species_data for gaps ----

# What species are extinct, or extinct in the wild?

extinct_species <- species_data %>%
                   dplyr::select(tsn, accepted_binomial, redlist_status) %>%
                   filter(redlist_status == "EX" | redlist_status == "EW")  %>%
                   distinct(.) %>%
                   mutate(extinct = 1)

# What species do we not have accepted binomials for?

#' TODO: Is there a way - or do we need - to pull the other details for these species?

species_without_binomials <- species_data %>%
                             filter(is.na(accepted_binomial)) %>%
                             select(tsn) %>%
                             distinct(.)

message("There are ", length(unique(species_without_binomials$tsn)), 
        " species without binomials in the species data dataframe")

# Remove all the no binomial species for the time being

species_data <- species_data %>%
                filter(accepted_binomial != is.na(accepted_binomial))

message("After removing missing binomials, there are ", 
        length(unique(species_data$tsn)), 
        " species with binomials in the species data dataframe")

# What species do we not have ecoregions for? Note this will produce nothing
# if subsetting by country


species_without_ecoregions <- species_data %>%
                              filter(is.na(ecoregion_code)) %>%
                              dplyr::select(tsn, accepted_binomial, ecoregion_code) %>%
                              distinct(.) %>%
                              mutate(no_ecoregion = 1)

message("There are ", 
        length(unique(species_without_ecoregions$tsn)), 
        " species without ecoregions in the species data dataframe")

missing_ecoregions_species_data <- species_data %>%
                                   filter(is.na(ecoregion_code))

# What species do we not have ecoregion object id for?

species_without_eco_objectid <- species_data %>%
                                filter(is.na(eco_objectid)) %>%
                                filter(!is.na(redlist_status)) %>%
                                dplyr::select(tsn, accepted_binomial, 
                                              eco_objectid) %>%
                                distinct(.) %>%
                                mutate(no_objectid = 1)

message("There are ", 
        length(unique(species_without_eco_objectid$tsn)), 
        " species without ecoregions in the species data dataframe")


# What species do we not have a redlist status in any year for?

species_without_redlist_status <- species_data %>%
                                  dplyr::select(-c(redlist_assessment_year,
                                                genus, species, ecoregion_code,
                                                eco_objectid)) %>%
                                  distinct(.) %>%
                                  group_by(tsn) %>%
                                  filter(all(is.na(redlist_status))) %>%
                                  mutate(no_redlist_status = 1) 

message("There are ", 
        length(unique(species_without_redlist_status$tsn)), 
        " species without a redlist status in the species data dataframe")


missing_redlist_species_data <- species_data %>%
                                group_by(tsn) %>%
                                filter(all(is.na(redlist_status)))

                              

# if(save_outputs == "yes") {
# 
#   write_csv(species_without_redlist_status, paste(outputs, "/", date, "_species_without_redlist_status.csv", sep = ""))
#   saveRDS(species_without_redlist_status, paste(outputs,"/", date, "_species_without_redlist_status.rds", sep = ""))
# 
# }
# 
# # What species do we have at least one redlist status for?
# 
# species_with_redlist_status <- species_data %>%
#                                dplyr::select(-redlist_assessment_year,
#                                               genus, species) %>%
#                                distinct(.) %>%
#                                group_by(tsn) %>%
#                                filter(!is.na(redlist_status)) %>%
#                                dplyr::select(-redlist_status) %>%
#                                distinct(.)

# Double check we haven't incorrectly detected no redlist status (should be no
# overlap between the two groups)

# overlap <- species_without_redlist_status$tsn %in% species_with_redlist_status$tsn
# any(overlap == TRUE) # The correct output to console should be FALSE (no overlap)


# What species are included in our database, and do we have the necessary data for them?
# prev_species_in_database <- species_in_database

# species_in_database <- species_data %>%
#                        select(tsn, accepted_binomial) %>%
#                        distinct(.) %>%
#                        merge(extinct_species[
#                        c("tsn", "extinct")], by = "tsn", all = TRUE) 
# 
species_with_incomplete_data <- species_without_redlist_status %>%
                                merge(species_without_ecoregions[
                                c("tsn", "no_ecoregion")], by = "tsn", all = TRUE) %>%
                                merge(species_without_eco_objectid[
                                c("tsn", "no_objectid")], by = "tsn", all = TRUE) %>%
                                mutate(incomplete_data = ifelse(no_redlist_status == 1|
                                                                no_ecoregion == 1|
                                                                no_objectid == 1|
                                                                is.na(accepted_binomial),
                                                                      1, 0))

                              

# write.csv(species_with_incomplete_data, file.path(outputs, 
#                                                  "200520_species_with_incomplete_data.csv"))

# species_with_complete_data <- species_in_database %>%
#                               filter(is.na(incomplete_data))

print(paste("The database contains", length(unique(species_data$tsn)),
            "species", sep = " "))

print(paste(length(unique(species_with_incomplete_data$tsn)),
            "of the species have incomplete information", sep = " "))

print(paste(sum(species_with_incomplete_data$no_redlist_status, na.rm = TRUE), 
            " have no redlist status, ", 
            sum(species_with_incomplete_data$no_ecoregion, na.rm = TRUE),
            " have no ecoregion, and ",
            sum(species_with_incomplete_data$no_objectid , na.rm = TRUE),
            " have no ecoregion object id", sep = ""))


# Find missing Red List Status ----

# Find red list status for the species we are missing info for

#species_with_incomplete_data <- read.csv("N:\\Quantitative-Ecology\\Simone\\extinction_test\\outputs\\2020-05-19_output_files\\200520_species_with_incomplete_data.csv")
#species_with_incomplete_data <- species_with_incomplete_data %>%
                                #select(-X)
# Add your iucn API token to you environment

# file.edit("~/.Renviron")

# In the file that opens, add the text below, save that file, then close and reopen R

# IUCN_REDLIST_KEY='6d46333255cbb1843dd8f6984e45f23a49e8011585b1f128c3a46d17c6a8b5ae'

rl_citation()


# test <- get_redlist_data('Thylacinus cynocephalus')

# Create list of search names

# species_list <- c('Thylacinus cynocephalus',"Fratercula arctica")

if (!("iucn_scraped_data.rds" %in% list.files(outputs))) {
  
species_list <- as.character(species_without_redlist_status$accepted_binomial)

species_list <- unique(species_list)

# Scrape the data in batches because otherwise the iucn website throws an error

species_list_1 <- species_list[1:970]

iucn_scraped_data_1 <- lapply(species_list_1, get_redlist_data)

species_list_2 <- species_list[971:1940]

iucn_scraped_data_2 <- lapply(species_list_2, get_redlist_data)

species_list_3 <- species_list[1941:length(species_list)]

iucn_scraped_data_3 <- lapply(species_list_3, get_redlist_data)
 
# species_list_4 <- species_list[2911:length(species_list)]

#iucn_scraped_data_4 <- lapply(species_list_4, get_redlist_data)

iucn_scraped_data <- c(iucn_scraped_data_1, iucn_scraped_data_2,
                       iucn_scraped_data_3)

saveRDS(iucn_scraped_data, file.path(outputs,"iucn_scraped_data.rds"))

rm(iucn_scraped_data_1, iucn_scraped_data_2, iucn_scraped_data_3)

} else {

iucn_scraped_data <- readRDS(file.path(outputs,"iucn_scraped_data.rds"))

}

# Turn it into a nice data frame (website downloads as a list)

iucn_scraped_data_list <- list()
no_data <- list()

for(i in seq_along(iucn_scraped_data)) {
  
if (length(iucn_scraped_data[[i]][[2]]) == 30) {
  
accepted_binomial <- iucn_scraped_data[[i]][[1]]
redlist_status <- as.character(iucn_scraped_data[[i]][[2]][13][[1]])
class <- as.character(iucn_scraped_data[[i]][[2]][5][[1]])
redlist_assessment_year <- iucn_scraped_data[[i]][[2]][11][[1]]
source <- "IUCN"
iucn_scraped_data_list[[i]] <- as.data.frame(cbind(accepted_binomial, redlist_status, 
                                             published_year, class, source))

} else if (length(iucn_scraped_data[[i]][[2]]) == 1) {
  
accepted_binomial <- iucn_scraped_data[[i]][2][[1]]
redlist_status <- NA
redlist_assessment_year <- NA
class <- NA
source <- "IUCN"
no_data[[i]] <- as.data.frame(cbind(accepted_binomial, redlist_status, 
                                             published_year, class, source))
  }
}

iucn_scraped_dataframe <- do.call(rbind, iucn_scraped_data_list)

iucn_scraped_dataframe$redlist_status <- as.character(iucn_scraped_dataframe$redlist_status)
iucn_scraped_dataframe$accepted_binomial <- as.character(iucn_scraped_dataframe$accepted_binomial)

names(iucn_scraped_dataframe) <- c("accepted_binomial", 
                                   "redlist_status",
                                   "redlist_assessment_year",
                                   "class",
                                   "source")
# species we had no luck with

no_data <- do.call(rbind, no_data)

# Add new red list status back into the species data

correct_order <- names(species_data)

found_redlist_species_data <- missing_redlist_species_data %>%
                              merge(iucn_scraped_dataframe[c("accepted_binomial", 
                                                             "redlist_status", 
                                                             "redlist_assessment_year")], 
                                    by = "accepted_binomial", all = FALSE) %>%
                              dplyr::select(-redlist_assessment_year.x, redlist_status.x) %>%
                              rename(redlist_status = redlist_status.y,
                                     redlist_assessment_year = redlist_assessment_year.y) %>%
                              select(correct_order) %>%
                              distinct(.)


# Save a version of the old species data just in case

species_data_all <- species_data

# Add in the new redlist data

species_data <- rbind(species_data, found_redlist_species_data)

# Remove the rest of the missing redlist species

species_data <- species_data %>%
                filter(!is.na(redlist_status))


print(paste(length(unique(species_data$tsn)),
            "species with redlist status", sep = " "))



# Get extinct species data from IUCN website ----

#' TODO: IMPORTANT - check if the 'tsn' is actually tsn or some other id
#' 

if (!("extinct_species_ecoregions.rds" %in% list.files(outputs))) {

iucn_scraped_extinct_species <- rl_sp_category("EX", key = NULL, parse = TRUE)
iucn_scraped_ew_species <- rl_sp_category("EW", key = NULL, parse = TRUE)

# Put names into a vector and clean it (names that don't follow the binomial
# 'genus species' format will break the next loop)

extinct_species_names <- unique(c(iucn_scraped_extinct_species[[3]]$scientific_name, 
                                  iucn_scraped_ew_species[[3]]$scientific_name))

extinct_species_names <- word(extinct_species_names, 1, 2)
extinct_species_names <- unique(extinct_species_names)

# Feed the names back to the iucn website to get the other variables we need

extinct_species_data_list <- list()
leftovers <- list()

for (i in seq_along(extinct_species_names)) {
  
  temp <- rl_search(extinct_species_names[i], parse = FALSE)
  
  if (is.list(temp[[2]]) && (length(temp[[2]]) != 0)) {
  
  tsn <- temp[[2]][[1]]$taxonid
  accepted_binomial <- temp[[2]][[1]]$scientific_name
  wwf_species_id <- NA
  genus <- temp[[2]][[1]]$genus
  species <- tail(accepted_binomial, n = 1)
  class <- temp[[2]][[1]]$class
  redlist_assessment_year <- temp[[2]][[1]]$published_year
  redlist_status <- "EX"
  common_name <- NA
  ecoregion_code <- NA
  df <- cbind(tsn, accepted_binomial, wwf_species_id, genus, species,
              class, redlist_assessment_year, redlist_status, common_name, ecoregion_code)
  
  extinct_species_data_list[[i]] <- df
  
  } else {
    
  leftovers[[i]] <- temp[[1]]
  
  }
  
}

extinct_species_data <- as.data.frame(do.call(rbind, extinct_species_data_list))

extinct_species_data$class <- as.character(extinct_species_data$class)
extinct_species_data$accepted_binomial <- as.character(extinct_species_data$accepted_binomial)

# Subset to the classes we are looking for

# extinct_species_data <- extinct_species_data %>%
#                         filter(class == c("REPTILIA", "AVES", "MAMMALIA", 
#                                          "AMPHIBIA"))

# Get ecoregions for the extinct species

extinct_species_names <- extinct_species_data %>%
                         dplyr::select(accepted_binomial) %>%
                         pull(.) %>%
                         unique(.)

# Add ecoregions to extinct species via gbif (SLOW CODE) ----

extinct_species_ecoregions <- get_gbif_data(extinct_species_names, 20, ecoregion_map)

no_ecoregions <- as.vector(extinct_species_ecoregions[[2]])

all_extinct_species_ecoregions <- extinct_species_ecoregions

extinct_species_ecoregions <- extinct_species_ecoregions[[1]] %>%
                              distinct(.) %>%
                              rename(accepted_binomial = species)

# Add ecoregions to our other data

extinct_species_data <- extinct_species_data %>%
                        merge(extinct_species_ecoregions[c("accepted_binomial",
                                  "ECO_ID", "OBJECTID")], 
                                  by = "accepted_binomial", all = TRUE) %>%
                        dplyr::select(-ecoregion_code) %>%
                        rename(ecoregion_code = ECO_ID,
                               eco_objectid = OBJECTID) %>%
                        mutate(redlist_status = "EX") %>%
                        select(correct_order)

} else {

  extinct_species_data <- readRDS(file.path(outputs,"extinct_species_ecoregions.rds"))

}


if (save_outputs == "yes") {
  
  saveRDS(extinct_species_data, 
          file.path(outputs, "extinct_species_ecoregions.rds"))
  
} 

# Add in the new extinct species data

species_data <- rbind(species_data, extinct_species_data)

# Find missing ecoregions for species without them ----

if (!("species_ecoregions_found.rds" %in% list.files(outputs))) {

species_no_ecoregion_names <- species_without_ecoregions %>%
                              select(accepted_binomial) %>%
                              unique(.) %>%
                              pull(.)

species_no_ecoregion_names <- word(species_no_ecoregion_names, 1, 2)

species_ecoregions_found <- get_gbif_data(species_no_ecoregion_names,
                                          20, ecoregion_map)

if (save_outputs == "yes") {
  
  saveRDS(species_ecoregions_found, 
          file.path(outputs, "species_ecoregions_found.rds"))
  
  } 

} else {
  
  species_ecoregions_found <- readRDS(file.path(outputs, 
                                                'species_ecoregions_found.rds'))
}

species_with_new_ecoregions <- unique(species_ecoregions_found[[1]])

names(species_with_new_ecoregions) <- c("accepted_binomial", "ecoregion_code",
                                        "ecoregion_name", "eco_objectid")

# Add new red list status back into the species data

correct_order <- names(species_data)

head(species_with_new_ecoregions)

found_ecoregion_species_data <- missing_ecoregions_species_data %>%
                                merge(species_with_new_ecoregions[c("accepted_binomial", 
                                                               "ecoregion_code", 
                                                               "eco_objectid")], 
                                      by = "accepted_binomial", all = FALSE) %>%
                                dplyr::select(- eco_objectid.x, ecoregion_code.x) %>%
                                rename(eco_objectid = eco_objectid.y,
                                       ecoregion_code = ecoregion_code.y) %>%
                                select(correct_order) %>%
                                distinct(.)


# Save a version of the old species data just in case

species_data_all <- species_data

# Add in the new redlist data

species_data <- rbind(species_data, found_ecoregion_species_data)

# Remove the rest of the missing redlist species

species_data <- species_data %>%
                filter(!is.na(ecoregion_code))


print(paste(length(unique(species_data$tsn)),
            "species with complete redlist and ecoregion data", sep = " "))

if (save_outputs == "yes") {
  
  saveRDS(species_data, 
          file.path(outputs, "species_data.rds"))
  
} 

# Calculate summary statistics ----

# Get the number of species in each ecoregion


species_by_ecoregion <- species_data %>%
                        group_by(ecoregion_code) %>%
                        summarize(n_distinct(tsn))


names(species_by_ecoregion) <- c("ecoregion_code", "number_of_species")

# Get the number of species in each redlist category

species_by_redlist_status <- species_data %>%
                             group_by(redlist_status) %>%
                             summarize(n_distinct(tsn))



# Get number of extinct species (grouping by ecoregion doesn't work well yet bc
# most extinct species haven't been assigned an ecoregion - TBD)

number_of_extinct_species <- species_data %>%
                   filter(redlist_status == "EX")  %>%
                   group_by(ecoregion_code) %>%
                   summarise(n_distinct(tsn))

names(number_of_extinct_species) <- c("ecoregion_code", "number_of_species_extinct")


number_of_extinct_wild_species <- species_data %>%
                        filter(redlist_status == "EX" | redlist_status == "EW")  %>%
                        group_by(ecoregion_code) %>%
                        summarise(n_distinct(tsn))

names(number_of_extinct_wild_species) <- c("ecoregion_code", "number_of_species_extinct")

if(save_outputs == "yes") {

  write_csv(number_of_extinct_species, paste(outputs, "/", date, "_extinct_species.csv", sep = ""))
  saveRDS(number_of_extinct_species, paste(outputs,"/", date, "_extinct_species.rds", sep = ""))

}

proportion_extinct <- species_by_ecoregion %>%
                      merge(number_of_extinct_species,
                            by = "ecoregion_code", all = TRUE) %>%
                      dplyr::mutate(proportion_extinct =
                                    number_of_species_extinct/number_of_species)


# Visualise summary stats ----

extinction_map_data <- inner_join(ecoregion_map, proportion_extinct[
                        c("ecoregion_code", "number_of_species_extinct")],
                        by = c("ECO_ID" = "ecoregion_code"))

extinction_map <- ggplot(extinction_map_data) +
                  geom_sf(aes(fill = number_of_species_extinct)) +
                  scale_fill_viridis_c(trans = "sqrt", alpha = .4)

extinction_map

if(save_outputs == "yes") {
  
  ggsave(file.path(outputs, "extinctions_map.png", extinction_map))
  
}

dev.off()



species_map_data <- inner_join(ecoregion_map, species_by_ecoregion[
                    c("ecoregion_code", "number_of_species")],
                    by = c("ECO_ID" = "ecoregion_code"))

species_map <- ggplot(species_map_data) +
               geom_sf(aes(fill = number_of_species)) +
               scale_fill_viridis_c(trans = "sqrt", alpha = .4)

species_map

if(save_outputs == "yes") {
  
  ggsave(file.path(outputs, "species_richness_map.png", species_map))
  
}

