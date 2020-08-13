
library(pryr)
library(dplyr)
library(sf)

get_ecoregions <- function(rangemap_directory_path, map) {
  
  range_map <- st_read(rangemap_directory_path)
  
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

get_simple_ecoregions <- function(rangemap, rangemap_directory_path, map) {
  
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

# Get_ecoregions performance testing ----

inputs <- "N:/Quantitative-Ecology/Simone/extinction_test/inputs"
parent_outputs <- "N:/Quantitative-Ecology/Simone/extinction_test/outputs"
eco_version <- "ecoregions_2017"
country <- "Australia"

db_version <- tail(sort(list.files(parent_outputs)), 1)
db_interim <- list.files(file.path(parent_outputs,db_version))[
  grepl("interim",list.files(file.path(parent_outputs,db_version)))]
db_outputs <- list.files(file.path(parent_outputs,db_version))[
  grepl("database",list.files(file.path(parent_outputs,db_version)))]

interim_outputs <- file.path(parent_outputs, db_version, db_interim)
outputs <- file.path(parent_outputs, db_version, db_outputs)

if (!is.na(country)) {
  
  location <- tolower(country)
  
} else {
  
  location <- "global"
  
}

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

# Review object size of subset (Australian) ecoregion_map

object_size(ecoregion_map)

# Check number of vertices

pts_ecoregion_map <- st_cast(ecoregion_map$geometry, "MULTIPOINT")
cnt_ecoregion_map <- sapply(pts_ecoregion_map, length)
sum(cnt_ecoregion_map)

# Simplify the map

ecoregion_map_simple <- st_simplify(ecoregion_map, preserveTopology = TRUE,
                                    dTolerance = 0.1)

# Check it's not too simple

plot(st_geometry(ecoregion_map))

plot(st_geometry(ecoregion_map_simple), border = "red")

# Review object size of simplifed ecoregion map

object_size(ecoregion_map)

object_size(ecoregion_map_simple)

# Check number of vertices

pts_ecoregion_map_simple <- st_cast(ecoregion_map_simple$geometry, "MULTIPOINT")
cnt_ecoregion_map_simple <- sapply(pts_ecoregion_map_simple, length)
sum(cnt_ecoregion_map_simple)

# Simplify the rangemap

range_map <- st_read("N:\\Quantitative-Ecology\\Simone\\extinction_test\\inputs\\redlist_amphibian_range_maps")

range_map <- st_simplify(range_map, preserveTopology = TRUE,
                         dTolerance = 0.1)

# Compare processing time for normal and simplified ecoregion maps

normal_time <- system.time(amphibians <- get_ecoregions("N:\\Quantitative-Ecology\\Simone\\extinction_test\\inputs\\redlist_amphibian_range_maps",
                                         ecoregion_map))

simple_time <- system.time(amphibians_simple <- get_ecoregions("N:\\Quantitative-Ecology\\Simone\\extinction_test\\inputs\\redlist_amphibian_range_maps",
                                         ecoregion_map_simple))

simple_range_time <- system.time(amphibians_simple_ranges <- 
                                 get_simple_ecoregions(rangemap,
                                                       "N:\\Quantitative-Ecology\\Simone\\extinction_test\\inputs\\redlist_amphibian_range_maps",
                                 ecoregion_map_simple))

# combine into named vector

times <- c(normal = normal_time[3], simple = simple_time[3],
           simple_ranges = simple_range_time[3])

# Get quickest

which.min(times)

# Get slowest
which.max(times)

# Compare the outputs

## Remove non aus spp

amphibians <- amphibians %>%
              filter(!is.na(ecoregion_id))

amphibians_simple <- amphibians_simple %>%
                     filter(!is.na(ecoregion_id))

amphibians_simple_ranges <- amphibians_simple_ranges %>%
                            filter(!is.na(ecoregion_id))
# Compare number of individual species

spp_no <- c(normal = length(unique(amphibians$binomial)), 
            simple = length(unique(amphibians_simple$binomial)),
            simple_ranges = length(unique(amphibians_simple_ranges$binomial)))

# Check we have same number of spp (should be 0)

spp_no[1] - spp_no[2]

spp_no[1] - spp_no[3]

# Check we have same spp in each (should be TRUE)

identical(unique(amphibians$binomial), 
          unique(amphibians_simple$binomial),
          unique(amphibians_simple_ranges$binomial))

normal_ecos <- split(amphibians, amphibians$ecoregion_id)
simple_ecos <- split(amphibians_simple, amphibians_simple$ecoregion_id)
ranges_ecos <- split(amphibians_simple_ranges, amphibians_simple_ranges$ecoregion_id)

quantities <- list()
problems <- list()

for (i in seq_along(normal_ecos)) {
  
  outcome <- identical(unique(normal_ecos[[i]]$binomial), 
                       unique(simple_ecos[[i]]$binomial))
  
  # Get logical vector of whether spp are the same or not
  
  x <- normal_ecos[[i]]$binomial %in% simple_ecos[[i]]$binomial
 
  # Number of different spp
  y <- length(x[x == FALSE]) 
  
  # Number of same spp
  z <- length(x[x = TRUE])
  
  # Number of species total
  
  total <- length(x)
  
  eco_num <- normal_ecos[[i]]$ecoregion_id[1]
  
  df <- data.frame(ecoregion_id = eco_num, num_spp = total,
                   num_diff = y, num_same = z)
  
  quantities[[i]] <- df
  
  
  if (outcome == FALSE) {
    
    different <- list(normal_ecos[[i]], simple_ecos[[i]])
    
    problems[[i]] <- different
    
  } else {
    
    eco_num <- normal_ecos[[i]]$ecoregion_id[1]
    print(paste("normal ecoregion", eco_num, 
          "contains the same species as simplified ecoregion", 
          eco_num, sep = " "))
  }
  
  problems <- list.clean(problems)
  
}

simple_quantities_df <- do.call(rbind, quantities)

sum(simple_quantities_df$num_diff)/sum(simple_quantities_df$num_spp)

ranges_quantities <- list()
ranges_problems <- list()

for (i in seq_along(normal_ecos)) {
  
  outcome <- identical(unique(normal_ecos[[i]]$binomial), 
                       unique(ranges_ecos[[i]]$binomial))
  
  # Get logical vector of whether spp are the same or not
  
  x <- normal_ecos[[i]]$binomial %in% ranges_ecos[[i]]$binomial
  
  # Number of different spp
  y <- length(x[x == FALSE]) 
  
  # Number of same spp
  z <- length(x[x = TRUE])
  
  # Number of species total
  
  total <- length(x)
  
  eco_num <- normal_ecos[[i]]$ecoregion_id[1]
  
  df <- data.frame(ecoregion_id = eco_num, num_spp = total,
                   num_diff = y, num_same = z)
  
  ranges_quantities[[i]] <- df
  
  
  if (outcome == FALSE) {
    
    different <- list(normal_ecos[[i]], ranges_ecos[[i]])
    
    ranges_problems[[i]] <- different
    
  } else {
    
    eco_num <- normal_ecos[[i]]$ecoregion_id[1]
    print(paste("normal ecoregion", eco_num, 
                "contains the same species as simplified ranges and ecoregions", 
                eco_num, sep = " "))
  }
  
  ranges_problems <- list.clean(ranges_problems)
  
}

ranges_quantities_df <- do.call(rbind, ranges_quantities)

# Calculate percentage difference between two

sum(ranges_quantities_df$num_diff)/sum(ranges_quantities_df$num_spp)










