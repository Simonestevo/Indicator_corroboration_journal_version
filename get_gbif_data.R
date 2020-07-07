# Get GBIF data

#Tutorial - https://ropensci.org/tutorials/rgbif_tutorial/

library(rgbif)

extinct_species <- readRDS('N:/Quantitative-Ecology/Simone/extinction_test/outputs/2020-06-11_output_files/iucn_scraped_extinct.rds')

extinct_species_full <- extinct_species

extinct_species <- as.character(extinct_species_full[1:10,2])



get_gbif_data <- function(species, observations) {
  
  require(rgbif)
  require(rlist)
  
  keys <- sapply(species, function(x) name_suggest(x)$key[1], USE.NAMES=FALSE)
  
  keys <- list.clean(keys, fun = is.null, recursive = FALSE)
  
  extinct_gbif_data <- occ_search(taxonKey=keys, limit= observations)
  
  extinct_coordinates <- list()
  no_coordinates <- list()
  
  for (i in seq_along(extinct_gbif_data)) {
    
    single_spp <- extinct_gbif_data[[i]][[3]]
    
    if( "decimalLatitude" %in% names(single_spp) == TRUE) {
      
      extinct_coordinates[[i]]  <- single_spp %>%
        select(name, year, decimalLatitude, decimalLongitude ) %>%
        distinct(.) %>%
        filter(complete.cases(decimalLatitude, decimalLongitude))
      
    } else {
      
      # no_coordinates[[i]] <- single_spp[1,1]
      
      print(paste("no co-ordinates found for", single_spp[1,1], "among", 
                  observations, "observations"), sep = " ")
      
    }
  }
  
  extinct_coordinates <- list.clean(extinct_coordinates, 
                                    fun = is.null, 
                                    recursive = FALSE)
  return(extinct_coordinates)
  
}