
get_ecoregions <- function(rangemap_directory_path, map) {
  
  range_map <- st_read(rangemap_directory_path)
  
  names(range_map) <- toupper(names(range_map)) # Make column names consistent
  
  range_map <- rename(range_map, geometry = GEOMETRY) # Convert geometry back tho
  
  # Get the ecoregions for each species
  
  ranges_ecoregions <- st_join(range_map, map, join = st_intersects)
  
  if (!is.null(ranges_ecoregions)) {
    
    print(paste("st_join complete for", 
                basename(rangemap_directory_path), sep = " "))
    
  }
  
  # Standardise the data so we can add it to the species data easily
  

  species_w_ecoregions <- as.data.frame(ranges_ecoregions %>%
                                          dplyr::select(BINOMIAL, 
                                                        eco_code, 
                                                        OBJECTID)) %>%
    dplyr::select(-geometry) %>%
    rename(binomial = BINOMIAL) %>%
    rename(ecoregion_code = eco_code) %>%
    rename(eco_objectid = OBJECTID) %>%
    mutate(source = "iucn_redlist") %>%
    mutate(redlist_status = "TBC") %>%
    mutate(ecoregion_code = as.character(ecoregion_code))
  

  if(!is.null(species_w_ecoregions)) {
    
    print(paste("finished adding ecoregions to", 
                basename(rangemap_directory_path), sep = " "))
    
  }
  
  return(species_w_ecoregions)
  
}


get_redlist_data <- function(species) {
  
  library("jsonlite")
  
  output <- jsonlite::fromJSON(rl_search_(species))
  
  print(paste("Data for", species, "found", sep = " "))
  
  return(output)
}


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

# Specify file locations ----

inputs <- "N:/Quantitative-Ecology/Simone/extinction_test/inputs"

#version <- "ecoregions_2017"
version <- "official_teow_wwf"

# Load new ecoregion map ----

ecoregion_map_all <- st_read(paste(inputs,version, sep = "/"))

# Subset ecoregion map for testing ----

if (version == "ecoregions_2017")  { 
  
  ecoregion_map_all <- rename(ecoregion_map_all, eco_code = ECO_ID)
  
  realm_name <- "Australasia"

}

if (version == "official_teow_wwf") {
  
  realm_name <- "AA"
  
}
  

ecoregion_map <- ecoregion_map_all %>% 
                 filter(REALM == realm_name) %>% # Subset to only one realm
                 select(OBJECTID, ECO_NAME, eco_code, geometry) # select only essential columns


# Sort the species into ecoregions ----

range_directory <- file.path(inputs, "redlist_amphibian_range_maps")
species_by_ecoregion <- get_ecoregions(range_directory, ecoregion_map)

# Look at the results (Species allocated NA are not within the 
# ecoregions? So if you run this globally you shouldn't get any/many NAs )

table2 <- species_by_ecoregion %>%
          group_by(ecoregion_code) %>%
          summarize(n_distinct(binomial))

# Strip down to just species by ecoregion

species_by_ecoregion_aust <- species_by_ecoregion %>%
                             filter(!is.na(ecoregion_code)) %>%
                             dplyr::select(-eco_objectid, -source, 
                                           -redlist_status) %>%
                             distinct(.) %>%
                             rename(accepted_binomial = binomial)

# Test the new method against existing ----

# Load in existing species data 

species_data <- readRDS(file.path(inputs, "deakin_species_data/species_data.rds"))

# and subset it to relevant area (this won't work bc eco codes don't match)
species_data <- species_data[species_data$ecoregion_code %in% 
                               ecoregion_map$eco_code , ]

# Subset it by class and strip down to species by ecoregion

species_data <- species_data %>%
                filter(class == "Amphibia") %>%
                dplyr::select(accepted_binomial, 
                              ecoregion_code) %>%
                distinct(.)
                
# Merge the old version and new so we can compare

names(species_data) <- paste(names(species_data), "old", sep = "_")
names(species_by_ecoregion_aust) <- paste(names(species_by_ecoregion_aust), 
                                          "new", sep = "_")
ecoregion <- "AA0701"

test_ecoregion_old <- species_data %>%
                      filter(ecoregion_code_old == ecoregion)

test_ecoregion_new <- species_by_ecoregion_aust %>%
                      filter(ecoregion_code_new == ecoregion)

species_in_new_but_not_old <- test_ecoregion_new[!test_ecoregion_new$accepted_binomial_new %in% 
                        test_ecoregion_old$accepted_binomial_old , ]

species_in_old_but_not_new <- test_ecoregion_old[!test_ecoregion_new$accepted_binomial_old %in% 
                                                   test_ecoregion_new$accepted_binomial_new , ]

# Add  in the red list information

species_data_RL <- species_by_ecoregion_aust %>%
                   merge(henriques_redlist_database,
                         by.x = "accepted_binomial_new",
                         by.y = "binomial")

test_ecoregion_new_RL <- species_data_RL %>%
                         filter(ecoregion_code_new == ecoregion)



# Ecoregion testing ----

# Get the list of changes to the ecoregions from dinerstein supp info

eco_changes_all <- read.csv(file.path(inputs, "ecoregions_2017", 
                                      "ecoregion_changes.csv"))

eco_changes <- eco_changes_all %>%
  dplyr::select(Realm, Ecoregion.Name,
                Biome, 
                Ecoregion.Change.from.Olson.et.al...2001..Terrestrial.Ecoregions.of.the.World) %>%
  rename(true_change = Ecoregion.Change.from.Olson.et.al...2001..Terrestrial.Ecoregions.of.the.World) %>%
  filter(Realm == "Australasia")


ecoregion_map_all_2017 <- st_read(file.path(inputs, "ecoregions_2017"))

# Prepare both ecoregion versions to join

eco_map_2017 <- ecoregion_map_all_2017 %>% 
  filter(REALM == "Australasia") %>% # Subset to only one realm
  select(OBJECTID, ECO_NAME, ECO_ID, geometry) %>% # select only essential columns
  # st_simplify(., 
  #             preserveTopology = TRUE, 
  #             dTolerance = 100) %>%  # Simplify geometry to improve processing times
  rename(OBJECTID_2017 = OBJECTID,
         ECO_NAME_2017 = ECO_NAME,
         ECO_ID_2017 = ECO_ID) # Simplify geometry

object_size(eco_map_2017) # Check object size

# Format 2001 map reading for joining

eco_map_2001 <- ecoregion_map_all %>% 
  filter(REALM == "AA") %>% # Subset to only one realm
  select(OBJECTID, ECO_NAME, ECO_ID, geometry) %>%  # select only essential columns
  # st_simplify(., 
  #             preserveTopology = TRUE, 
  #             dTolerance = 100) %>%  # Simplify geometry to improve processing times
  rename(OBJECTID_2001 = OBJECTID,
         ECO_NAME_2001 = ECO_NAME,
         ECO_ID_2001 = ECO_ID)

# Join the two maps

eco_maps <- st_join(eco_map_2017, eco_map_2001)

# Add the correct changes dataframe

eco_maps <- eco_maps %>%
  mutate(ECO_NAME_2001 = as.character(ECO_NAME_2001),
         ECO_NAME_2017 = as.character(ECO_NAME_2017)) %>%
  mutate(found_change = ifelse(ECO_NAME_2001 == ECO_NAME_2017, 0, 1)) %>%
  merge(eco_changes[c("Ecoregion.Name", "true_change")],
        by.x = "ECO_NAME_2017",
        by.y = "Ecoregion.Name") %>%
  st_set_geometry(., NULL)

write.csv(eco_maps, file.path(interim_outputs, "ecoregion_join_changes.csv"))

# Subset the species data to only species found in the area we are looking at
# included in the ecoregion 2001 map

species_data <- species_data[species_data$ecoregion_code %in% 
                               eco_map_2001$eco_code , ]

sp_data <- species_data %>%
  select(tsn, accepted_binomial, eco_objectid, ecoregion_code) %>%
  distinct(.)

# sp_data <- species_data[1:100,]
# 
# sp_data <- sp_data %>% 
#   select(tsn, accepted_binomial, ecoregion_code) %>%
#   distinct(.)
# 
# sp_map_2001 <- right_join(ecoregion_map, sp_data, by = "ecoregion_code")

# sp_map_2017 <- st_join(sp_map_2001, eco_map_2017)





