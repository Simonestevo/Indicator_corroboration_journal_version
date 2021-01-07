
# Using R version 4.0.3

## Clear the space
rm(list = ls()) # clear memory

# Load packages ----

# Data handling and table reshaping
library(tidyverse)
library(reshape2)
library(data.table)
library(rlist)
library(ALA4R)

# Plotting
library(ggplot2)
library(viridis)

# Maps
library(sf)
library(leaflet)

# Functions ----

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



# Set input and output locations ----

date <- Sys.Date()
country <- "Australia" # If not subsetting, set as NA, e.g. country <- NA
inputs <- "N:/Quantitative-Ecology/Simone/extinction_test/outputs/version_3/2020-08-10_interim_output_files"
parent_outputs <- "N:/Quantitative-Ecology/Simone/extinction_test/outputs"
  
output_directory <- file.path(parent_outputs, 
                              paste(date, 
                                    "species_data_check_output_files", 
                                    sep = "_"))

dir.create(output_directory)

# Get input data ----

# Load species data

species_data <- readRDS(file.path(inputs, "global_species_data_3.rds"))


# Load ecoregion map

ecoregion_map_all <- readRDS(paste(
  file.path("N:/Quantitative-Ecology/Simone/extinction_test/inputs",
            "ecoregions_2017"),
            "Ecoregions2017valid.rds"))

ecoregion_map <- ecoregion_map_all %>%
                 dplyr::select(ECO_ID, ECO_NAME, OBJECTID, REALM, geometry)

rm(ecoregion_map_all)

world_map <- st_read("N:/Quantitative-Ecology/Simone/extinction_test/inputs/countries_WGS84")

australia <- world_map %>%
             filter(CNTRY_NAME == "Australia")

# Subset by country

if (!is.na(country)) {

ecoregion_country_df <- readRDS(file.path("N:/Quantitative-Ecology/Simone/extinction_test/outputs/version_3/2020-08-10_database_output_files",
                                          "ecoregion_country_data.rds"))

ecoregion_subset <- ecoregion_country_df %>%
                    filter(CNTRY_NAME == country) %>%
                    unique(.)
  
ecoregion_map <- ecoregion_map[ecoregion_map$ECO_ID %in% 
                                   ecoregion_subset$ECO_ID,]

species_data <- species_data[species_data$ecoregion_id %in% 
                               ecoregion_subset$ECO_ID,]
  
}


species_data_unique <- species_data %>%
  group_by(ecoregion_id) %>%
  mutate(spp_number = match(tsn, unique(tsn))) %>%
  select(ecoregion_id, binomial, class, spp_number) %>%
  distinct(.)

# Australia comparison ----

# source - https://www.environment.gov.au/cgi-bin/sprat/public/sprat.pl

v2_gov_species_data_all <- read.csv("N:/Quantitative-Ecology/Simone/extinction_test/inputs/australian_threatened_species_gov/07012021-024331-report.csv", head = TRUE)

# Subset to relevant columns

aus_species_data <- v2_gov_species_data_all %>%
                    select(Scientific.Name, EPBC.Threat.Status,
                           Class, IUCN.Red.List.Listed.Names,
                           IUCN.Red.List) %>%
                    filter(Class == c("Aves", "Mammalia", "Amphibia"))

# Get taxonomic serial numbers so can be matched to our data

aus_binomials <- unique(aus_species_data$Scientific.Name)
aus_synonyms <- find_synonyms(aus_binomials)

aus_species_data <- aus_species_data %>%
                    merge(aus_synonyms[c("binomial", "tsn")],
                          by.x = "Scientific.Name",
                          by.y = "binomial")

head(aus_species_data)

# Check how many species on australian list
length(unique(aus_species_data$Scientific.Name))

# EPBC status breakdown
table(aus_species_data$EPBC.Threat.Status)

# Australian red list status breakdown
table(aus_species_data$IUCN.Red.List)

# Our red list status breakdown
table(species_data$redlist_status)

# Join the australian government data to ours so we can compare

species_comparison <- species_data %>%
                      group_by(tsn) %>%
                      mutate(redlist_assessment_year = 
                             as.numeric(as.character(redlist_assessment_year))) %>%
                      filter(redlist_assessment_year == max(redlist_assessment_year)) %>%
                      select(binomial, tsn, redlist_status) %>%
                      distinct(.) %>%
                      merge(aus_species_data, by = "tsn", all = TRUE) %>%
                      mutate(epbc_status = ifelse(EPBC.Threat.Status == 
                                                 "Vulnerable", "VU",
                                                 ifelse(EPBC.Threat.Status == 
                                                 "Endangered","EN",
                                                 ifelse(EPBC.Threat.Status == 
                                                 "Critically Endangered", "CR",
                                                 ifelse(EPBC.Threat.Status ==
                                                 "Extinct", "EX",
                                                 ifelse(EPBC.Threat.Status ==
                                                 "Extinct in the wild", "EX",
                                                 ifelse(EPBC.Threat.Status ==
                                                 "Conservation Dependent", "CD",
                                                 NA))))))) %>%
                      select(tsn, binomial,redlist_status, epbc_status,
                             IUCN.Red.List,Scientific.Name, Class) %>%
                      mutate(aus_redlist_status = ifelse(IUCN.Red.List == 
                             "Vulnerable", "VU",
                             ifelse(IUCN.Red.List == 
                             "Endangered","EN",
                             ifelse(IUCN.Red.List == 
                             "Critically Endangered", "CR",
                             ifelse(IUCN.Red.List ==
                             "Extinct", "EX",
                             ifelse(IUCN.Red.List ==
                             "Extinct in the wild", "EX",
                             ifelse(IUCN.Red.List ==
                             "Least Concern", "LC",
                             ifelse(IUCN.Red.List ==
                             "Data Deficient", "DD",
                             ifelse(IUCN.Red.List ==
                             "Near Threatened", "NT",
                             NA))))))))) %>%
                      mutate(match_status = 
                             ifelse(redlist_status == aus_redlist_status,
                                      TRUE, FALSE)) 

head(species_comparison)

# Check our red list data matches theirs 

table(species_comparison$match_status)

par(mfrow=c(1,2))
# Breakdown of our red list categories

barplot(prop.table(table(species_comparison$redlist_status)), main = "Our data")

# Breakdown of their red list categories

barplot(prop.table(table(species_comparison$aus_redlist_status)),
                main = "Australian\ngovernment\ndata")

### Note - our numbers are different, but proportions mostly the same, 
### with the exceptions of extinctions - they have the same number we do, 
### even though overall they have about a quarter of the number of species







# Create map of species present ----

# Add an index


head(species_data_unique)

max_sp <- max(species_data_unique$spp_number)

sp_cols <- 1:max_sp

sp_cols <- rep(sp_cols, times = length(unique(species_data_unique$ecoregion_id)))

ecoregion_id <- rep(unique(species_data_unique$ecoregion_id), each = max_sp)

x <- species_data_unique %>%
     group_by(ecoregion_id) %>%
     add_row()

# Spread the species data

species_data_wide <- species_data_unique %>%
                     spread(key = ecoregion_id,
                            value = binomial)

species_data_wide <- species_data_wide[1:5, 1:5]

for (i in seq(nrow(species_data_wide))) {
  
for (j in seq(ncol(species_data_wide[2:length(species_data_wide)]))) {

if (!is.na(species_data_wide[i,j])) {

species_data_wide[,j] <- as.character(species_data_wide[,j])
  
species_data_wide[i,j] <- colnames(species_data_wide[,j])

print(paste("replacing cell value with", colnames(species_data_wide[,j]), sep = " "))
  
    }
  }
}

species_data_wide <- species_data_wide %>%
                     mutate(number_species = rowSums(!is.na(.))))
        
species_map_data <- left_join(ecoregion_map[c("ECO_ID", "ECO_NAME", "geometry")], 
                              species_data_wide,
                              by.x = "ECO_ID",
                              by.y = "ecoregion_id") 



# Create a colour palette for proportion of extinctions

ext_pal <- colorNumeric("PuBu", domain = indicator_map_input_data$number_species)
hfp_pal <- colorNumeric("PuBu", domain = indicator_map_input_data$HFP_2005)

#' TODO: Add real, not scaled indicator values
#' TODO: Add ecoregion labels?
#' TODO: Try adding all the indicator values into one label

indicator_map_input_data %>% 
  leaflet() %>% 
  addTiles() %>% 
  addPolygons(weight = 1, color = ~ext_pal(extinct_2005),
              fillOpacity = 0.8, group = "Proportion of species extinct",
              # add labels that display indicator value and ecoregion id
              label = ~paste(ecoregion_id, "proportion extinct =", extinct_2005, 
                             sep = " "),
              # highlight polygons on hover
              highlightOptions = highlightOptions(weight = 5, color = "white",
                                                  bringToFront = TRUE)) %>%
  addPolygons(weight = 1, color = ~hfp_pal(HFP_2005),
              fillOpacity = 0.8, group = "Human footprint index",
              # add labels that display indicator value and ecoregion id
              label = ~paste(ecoregion_id, "HFP =", HFP_2005, sep = " "),
              # highlight polygons on hover
              highlightOptions = highlightOptions(weight = 5, color = "white",
                                                  bringToFront = TRUE)) %>%
  addLayersControl(baseGroups = c("OSM", "Carto", "Esri"), 
                   overlayGroups = c("Proportion of species extinct", 
                                     "Human footprint index"))


  
  
  
  