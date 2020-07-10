# Get GBIF data

#Tutorial - https://ropensci.org/tutorials/rgbif_tutorial/

# library(rgbif)
# library(sf)
# library(rgdal)
# library(tidyverse)
# 
# extinct_species <- readRDS('N:/Quantitative-Ecology/Simone/extinction_test/outputs/2020-06-11_output_files/iucn_scraped_extinct.rds')
# extinct_species_full <- extinct_species
# extinct_species <- as.character(extinct_species_full[1:5,2])
# 
# ecoregion_map_all <- st_read("N:/Quantitative-Ecology/Simone/extinction_test/inputs/official_teow_wwf/wwf_terr_ecos.shp")
# ecoregion_map <- ecoregion_map_all %>% dplyr::select(eco_code, ECO_NAME, geometry)
# names(ecoregion_map) <- c("ecoregion_code", "ECO_NAME", "geometry")


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

get_gbif_data <- function(species, observations, polygon_map) {
  
  require(rgbif)
  require(rlist)
  require(sf)
  require(tidyverse)
  
  keys <- sapply(species, function(x) name_suggest(x)$key[1], USE.NAMES=FALSE)
  
  keys <- list.clean(keys, fun = is.null, recursive = FALSE)
  
  extinct_gbif_data <- occ_search(taxonKey = keys, limit = observations)
  
  extinct_coordinates <- list()
  no_coordinates <- list()
  
  for (i in seq_along(extinct_gbif_data)) {
    
    single_spp <- extinct_gbif_data[[i]][[3]]
    
    if ( "decimalLatitude" %in% names(single_spp) == TRUE) {
      
      extinct_coordinates[[i]] <- single_spp %>%
        select(name, decimalLatitude, decimalLongitude ) %>%
        distinct(.) %>%
        filter(complete.cases(decimalLatitude, decimalLongitude))
      
    } else {
      
      no_coordinates[[i]] <- single_spp[1,1]
      
      print(paste("no co-ordinates found for", single_spp[1,1], "among", 
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
  
  output <- list(extinct_coordinates, no_coordinates)
  
  names(output) <- c("extinct_coordinates","no_coordinates")
  
  return(output)
  
}

# # Examples of use
# 
# # Get ecoregions for list of extinct species according to 5 observations per spp
# 
# extinct_ecoregions_1 <- get_gbif_data(extinct_species,5, ecoregion_map)
# 
# # Make list of species with no coordinates found
# 
# not_found <- pull(as.vector(extinct_ecoregions_1[[2]]))
# 
# # Get ecoregions for species not found by increasing observations searched to 500
# 
# extinct_ecoregions_2 <- get_gbif_data(not_found,500, ecoregion_map)

