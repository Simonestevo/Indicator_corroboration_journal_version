
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
external_inputs <- "N:/Quantitative-Ecology/Simone/extinction_test/inputs/species_data_checks"

  
output_directory <- file.path(parent_outputs, 
                              paste(date, 
                                    "species_data_check_output_files", 
                                    sep = "_"))

dir.create(output_directory)

# Get input data ----

# Load species data

species_data <- readRDS(file.path(inputs, "global_species_data_3.rds"))

species_data_unique <- species_data %>%
                       group_by(ecoregion_id) %>%
                       mutate(spp_number = match(tsn, unique(tsn))) %>%
                       select(ecoregion_id, binomial, class, spp_number) %>%
                       distinct(.)

species_data_2008 <- species_data %>%
                     filter(redlist_assessment_year == 2008,
                            ecoregion_id != 0,
                            class != "Reptilia") %>%
                     select(ecoregion_id, binomial, 
                            redlist_status, class) %>%
                     distinct(.)


# Load ecoregion map

ecoregion_map_all <- readRDS(paste(
  file.path("N:/Quantitative-Ecology/Simone/extinction_test/inputs",
            "ecoregions_2017"),
            "Ecoregions2017valid.rds"))

ecoregion_map <- ecoregion_map_all %>%
                 dplyr::select(ECO_ID, ECO_NAME, OBJECTID, REALM, BIOME_NAME, geometry)

rm(ecoregion_map_all)

# Add Biome to species data

ecoregion_biomes <- ecoregion_map %>%
                    select(ECO_ID, BIOME_NAME) %>%
                    st_set_geometry(NULL) %>%
                    rename(ecoregion_id = ECO_ID,
                           biome = BIOME_NAME)

species_data <- species_data %>%
                merge(ecoregion_biomes,
                      by = "ecoregion_id") 
                

world_map <- st_read("N:/Quantitative-Ecology/Simone/extinction_test/inputs/countries_WGS84")

australia <- world_map %>%
             filter(CNTRY_NAME == "Australia")

ecoregion_country_df <- readRDS(file.path("N:/Quantitative-Ecology/Simone/extinction_test/outputs/version_3/2020-08-10_database_output_files",
                                          "ecoregion_country_data.rds"))


# Australia comparison ----

ecoregion_subset <- ecoregion_country_df %>%
                    filter(CNTRY_NAME == "Australia") %>%
                    unique(.)

aus_ecoregion_map <- ecoregion_map[ecoregion_map$ECO_ID %in% 
                                 ecoregion_subset$ECO_ID,]

aus_species_data_ours <- species_data[species_data$ecoregion_id %in% 
                               ecoregion_subset$ECO_ID,]

aus_species_data_ours <- aus_species_data_ours %>%
                         filter(class != "Reptilia") %>%
                          group_by(tsn) %>%
                          mutate(redlist_assessment_year = 
                                   as.numeric(as.character(redlist_assessment_year))) %>%
                          filter(redlist_assessment_year == max(redlist_assessment_year)) %>%
                          select(binomial, tsn, redlist_status) %>%
                          distinct(.)

# source - https://www.environment.gov.au/cgi-bin/sprat/public/sprat.pl

v2_gov_species_data_all <- read.csv(file.path(external_inputs, 
                                              "aus_07012021-024331-report.csv"))

# Subset to relevant columns

aus_species_data_gov <- v2_gov_species_data_all %>%
                    select(Scientific.Name, EPBC.Threat.Status,
                           Class, IUCN.Red.List.Listed.Names,
                           IUCN.Red.List) %>%
                    filter(Class == c("Aves", "Mammalia", "Amphibia"))

# Get taxonomic serial numbers so can be matched to our data

aus_binomials <- unique(aus_species_data_gov$Scientific.Name)
aus_synonyms <- find_synonyms(aus_binomials)

aus_species_data_gov <- aus_species_data_gov %>%
                    merge(aus_synonyms[c("binomial", "tsn")],
                          by.x = "Scientific.Name",
                          by.y = "binomial")

head(aus_species_data_gov)

# Check how many species on australian list
length(unique(aus_species_data_gov$Scientific.Name))

length(unique(aus_species_data_ours$tsn))

# EPBC status breakdown
table(aus_species_data_gov$EPBC.Threat.Status)

# Australian red list status breakdown
table(aus_species_data_gov$IUCN.Red.List)

# Our red list status breakdown
table(species_data$redlist_status)

# Join the australian government data to ours so we can compare

species_comparison <- aus_species_data_ours %>%
                      merge(aus_species_data_gov, by = "tsn", all = TRUE) %>%
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

table(aus_species_data_gov$IUCN.Red.List)
table(aus_species_data_ours$redlist_status)

### Note - our numbers are different, but proportions mostly the same, 
### with the exceptions of extinctions - they have the same number we do, 
### even though overall they have about a quarter of the number of species

spp_on_aus_list_only <- species_comparison %>%
                        filter(is.na(redlist_status))

length(unique(spp_on_aus_list_only$tsn))

spp_on_our_list_only <- species_comparison %>%
                        filter(is.na(aus_redlist_status))

length(unique(spp_on_our_list_only$tsn))

### Note - there are quite a few extinct species not on our list, pretty sure
### because there are no range maps available for them

# Serengeti comparison ----

serengeti_species_external <- read.csv(file.path(external_inputs, "serengeti_map_of_life_2.csv"))

serengeti_species_external <- serengeti_species_external %>%
                              dplyr::filter(Taxonomic.Group == "Birds"|
                                     Taxonomic.Group == "Mammals"|
                                       Taxonomic.Group == "Amphibians")

unique(serengeti_species_external$Taxonomic.Group)

# Subset to ecoregion 57 Southern Acacia-Commiphora bushlands and thickets

serengeti_species <- species_data %>%
                     filter(ecoregion_id == 57) %>%
                     filter(class != "Reptilia") %>%
                     group_by(tsn) %>%
                     mutate(redlist_assessment_year = 
                               as.numeric(as.character(redlist_assessment_year))) %>%
                     filter(redlist_assessment_year == max(redlist_assessment_year)) %>%
                     select(binomial, tsn, redlist_status) %>%
                     distinct(.) 

head(serengeti_species)

# Get tsn for external list so we can match properly

serengeti_binomials <- unique(serengeti_species_external$Scientific.Name)
serengeti_synonyms <- find_synonyms(serengeti_binomials)

serengeti_synonyms <- serengeti_synonyms %>%
                      rename(external_species_list = accepted_name)


# Check how many species in each data set

length(unique(serengeti_species$tsn))
length(unique(serengeti_synonyms$tsn))

serengeti_comparison <- serengeti_species %>%
                        merge(serengeti_synonyms[c("external_species_list", 
                                                           "tsn")],
                              by = "tsn", all = TRUE) %>%
                        rename(internal_species_list = binomial) %>%
                        distinct(.)

head(serengeti_comparison)

serengeti_internal_only <- serengeti_comparison %>%
                           filter(is.na(external_species_list))

length(unique(serengeti_internal_only$tsn))

serengeti_external_only <- serengeti_comparison %>%
                           filter(is.na(internal_species_list))

length(unique(serengeti_external_only$tsn))

# Xeric shrublands check ----

xeric_species_data <- species_data %>%
                      filter(biome == "Deserts & Xeric Shrublands") %>%
                      mutate(redlist_assessment_year = 
                               as.numeric(as.character(redlist_assessment_year))) %>%
                      select(ecoregion_id, tsn, binomial, redlist_assessment_year,
                             redlist_status, biome) %>%
                      merge(ecoregion_country_df[c("ECO_ID", "CNTRY_NAME")],
                            by.x = "ecoregion_id",
                            by.y = "ECO_ID") %>%
                      filter(redlist_assessment_year == 2008) %>%
                      distinct(.) %>%
                      rename(country = CNTRY_NAME)

xeric_species_data_simple <- xeric_species_data %>%
                             select(tsn, binomial, redlist_status) %>%
                             distinct(.)

# IUCN breakdown
table(xeric_species_data_simple$redlist_status)

# Total spp number

xeric_tot <- length(unique(xeric_species_data_simple$tsn))
xeric_tot

# Proportion extinct

xeric_rl <- xeric_species_data_simple %>%
  group_by(redlist_status) %>%
  mutate(rl = n_distinct(tsn)) %>%
  select(redlist_status, rl) %>%
  distinct(.)

xeric_rl

14/xeric_tot

# Proportion threatened

(180+269+14+62+1)/xeric_tot


# Look at a biome positively correlated - Mangroves

mangrove_species_data <- species_data %>%
                          filter(biome == "Mangroves") %>%
                          mutate(redlist_assessment_year = 
                                   as.numeric(as.character(redlist_assessment_year))) %>%
                          select(ecoregion_id, tsn, binomial, redlist_assessment_year,
                                 redlist_status, biome) %>%
                          merge(ecoregion_country_df[c("ECO_ID", "CNTRY_NAME")],
                                by.x = "ecoregion_id",
                                by.y = "ECO_ID") %>%
                          filter(redlist_assessment_year == 2008) %>%
                          distinct(.) %>%
                          rename(country = CNTRY_NAME)

mangrove_species_data_simple <- mangrove_species_data %>%
  select(tsn, binomial, redlist_status) %>%
  distinct(.)


# IUCN breakdown
table(mangrove_species_data_simple$redlist_status)

# Total spp number

mangrove_tot <- length(unique(mangrove_species_data_simple$tsn))
mangrove_tot

# Proportion extinct

mangrove_rl <- mangrove_species_data_simple %>%
                    group_by(redlist_status) %>%
                    mutate(rl = n_distinct(tsn)) %>%
                select(redlist_status, rl) %>%
                distinct(.)

5/mangrove_tot

# Proportion threatened

(180+49+322+5)/mangrove_tot

# Look at Madagascar

mad_xeric <- xeric_species_data %>%
  filter(country == "Madagascar") %>%
  select(binomial, redlist_status) %>%
  distinct(.)

mad_xeric_rl <- mad_xeric %>%
  group_by(redlist_status) %>%
  mutate(rl = n_distinct(binomial)) %>%
  select(redlist_status, rl) %>%
  distinct(.)

mad_xeric_rl

length(unique(mad_xeric$binomial))

# Proportion extinct

0/329

# Proportion threatened

(2+14+25+24)/329


## Australia

au_xeric <- xeric_species_data %>%
            filter(country == "Australia") %>%
            select(binomial, redlist_status) %>%
            distinct(.)

au_xeric_rl <- au_xeric %>%
               group_by(redlist_status) %>%
               mutate(rl = n_distinct(binomial)) %>%
               select(redlist_status, rl) %>%
               distinct(.)

au_xeric_rl

length(unique(au_xeric$binomial))

4/556

table(au_xeric$redlist_status)

# TNC comparison ----

# Compare our species count to TNC species count (noting not all ecoregions match)

ecoregion_differences <- read.csv(file.path(external_inputs, 
                                            "ecoregion_join_changes.csv"))

tnc_data_all <- st_read(file.path(external_inputs, "tnc_wwf_spp_list/data0"))

tnc_key <- tnc_data_all %>%
           dplyr::select(ECO_NAME, ECO_ID) %>%
           st_set_geometry(NULL) %>%
           distinct(.)

names(tnc_key)

ecoregions <- ecoregion_map %>% st_set_geometry(NULL)

ecoregion_key <- ecoregions %>%
                  rename(ecoregion_id = ECO_ID) %>%
                  merge(tnc_key,
                        by = "ECO_NAME") %>%
                  rename(ecoregion_id_tnc = ECO_ID) %>%
                  select(ECO_NAME, ecoregion_id, ecoregion_id_tnc) %>%
                  distinct(.)
  
head(ecoregion_key)

tnc_data <- tnc_data_all %>%
            st_set_geometry(NULL) %>%
            rename(external_threatened_n = thrtnd_cnt) %>%
            dplyr::select(ECO_ID, 
                   external_threatened_n) %>%
            merge(ecoregion_key, 
                  by.x = "ECO_ID",
                  by.y = "ecoregion_id_tnc") %>%
            distinct(.) %>%
            filter(ECO_ID != -9999) %>%
            rename(ecoregion_id_tnc = ECO_ID) %>%
            dplyr::select(ECO_NAME, ecoregion_id, ecoregion_id_tnc, everything())

head(tnc_data)

internal_cnt <- species_data %>%
                filter(redlist_assessment_year == 2008) %>%
                group_by(ecoregion_id) %>%
                mutate(internal_threatened_n = n_distinct(binomial[redlist_status == "EN"|
                                                      redlist_status == "CR"|
                                                      redlist_status == "CR(PE)"|
                                                      redlist_status == "VU"])) %>%
                select(ecoregion_id, internal_threatened_n) %>%
                distinct(.)
  
head(internal_cnt)

count_comparison <- internal_cnt %>%
                    merge(tnc_data,
                          by = "ecoregion_id") %>%
                    mutate(difference = internal_threatened_n - external_threatened_n)

head(count_comparison)

hist(count_comparison$difference, breaks = 40)

## Look at species in biggest difference ecoregions

# We have 122 fewer species in borneo rainforest

borneo_rforest <- species_data_2008 %>%
                  filter(ecoregion_id == 219)

# We have around 70 species more in north andean paramo

n_andean_paramo <- species_data_2008 %>%
                   filter(ecoregion_id == 593)

# Look at area - are the biggest discrepancies in the small ecoregions?

ecoregion_area <- read.csv("N:/Quantitative-Ecology/Simone/extinction_test/inputs/ecoregions_2017/ecoregions2017_area_km2.csv")
ecoregion_area <- rename(ecoregion_area, ecoregion_id = Eco..ID,
                         area = Ecoregion.area..km2.)
count_comparison <- count_comparison %>%
                    merge(ecoregion_area[c("ecoregion_id","area")],
                          by = "ecoregion_id") %>%
                    mutate(rel_area = area/max(area))

diff_mean <- mean(count_comparison$difference)
diff_sd <- sd(count_comparison$difference)

confidence_interval <- function(vector, interval) {
  # Standard deviation of sample
  vec_sd <- sd(vector)
  # Sample size
  n <- length(vector)
  # Mean of sample
  vec_mean <- mean(vector)
  # Error according to t distribution
  error <- qt((interval + 1)/2, df = n - 1) * vec_sd / sqrt(n)
  # Confidence interval as a vector
  result <- c("lower" = vec_mean - error, "upper" = vec_mean + error)
  return(result)
}

confidence_interval(count_comparison$difference, 0.95)

# Find extinct species without an ecoregion ----

extinct_species_redlist_data <- readRDS("N:\\Quantitative-Ecology\\Simone\\extinction_test\\outputs\\version_3\\2020-08-10_interim_output_files\\global_extinct_species_redlist_data.rds")

extinct_w_ecoregions <- species_data %>%
                        filter(redlist_status == "EX") %>%
                        select(ecoregion_id, tsn, binomial, 
                               redlist_status, class) %>%
                        distinct(.)

extinct_stocktake <- extinct_species_redlist_data %>%
                     merge(extinct_w_ecoregions[c("binomial", "ecoregion_id")],
                           by = "binomial",
                           all = TRUE) %>%
                     select(-redlist_assessment_year) %>%
                     dplyr::select(binomial, tsn, ecoregion_id) %>%
                     distinct(.)

extinct_with_no_ecoregion <- extinct_stocktake %>%
                             filter(is.na(ecoregion_id))

write.csv(extinct_with_no_ecoregion, file.path(output_directory, 
                                     "extinct_species_needing_ecoregions.csv"))

extinct_with_ecoregions <- extinct_stocktake %>%
                           filter(!is.na(ecoregion_id))

write.csv(extinct_with_ecoregions, file.path(output_directory, 
                                               "extinct_species_location_known.csv"))

length(unique(extinct_with_no_ecoregion$binomial))


# Create map of species present ----

species_data_for_map <- species_data_unique[complete.cases(species_data_unique),]
species_data_for_map <- ungroup(species_data_for_map)

country <- "Australia"

if (!is.na(country)) {
  
  ecoregion_subset <- ecoregion_country_df %>%
    filter(CNTRY_NAME == country) %>%
    unique(.)
  
}

if (!is.na(country)) {
  
species_data_for_map <- species_data_for_map[species_data_for_map$ecoregion_id %in% 
                          ecoregion_subset$ECO_ID,] 

ecoregion_map_subset <- ecoregion_map[ecoregion_map$ECO_ID %in% 
                                 ecoregion_subset$ECO_ID,] 
  
}

# Create a species list per ecoregion

ecos <- unique(species_data_for_map$ecoregion_id)

eco_spp <- list()

for (i in seq_along(ecos)){

id <- ecos[i]  

spp <- species_data_for_map %>%
       filter(ecoregion_id == id) %>%
       select(binomial) %>%
       pull(.)

ecoregion_id <- as.numeric(id)

number_species <- length(spp)

species_list <- paste(spp, collapse = ",  ")

eco_spp[[i]] <- as.data.frame(cbind(ecoregion_id, number_species, species_list))

print(paste("Completed loop for ecoregion", ecoregion_id, sep = " "))

}

ecoregion_sp_lists <- do.call(rbind,eco_spp)
ecoregion_sp_lists$ecoregion_id <- as.numeric(ecoregion_sp_lists$ecoregion_id)
class(ecoregion_sp_lists$ecoregion_id)

ecoregion_map_subset <- rename(ecoregion_map_subset, ecoregion_id = ECO_ID)
names(ecoregion_map_subset)

View(ecoregion_sp_lists)
        
species_map_data <- left_join(ecoregion_map_subset,
                              ecoregion_sp_lists, by = "ecoregion_id") 

names(species_map_data)

# Create a colour palette for proportion of extinctions

species_map_data$number_species <- as.numeric(species_map_data$number_species)

num_pal <- colorNumeric("PuBu", domain = species_map_data$number_species)
spp_pal <- colorNumeric("PuBu", domain = species_map_data$HFP_2005)

# Make map

species_map_data %>% 
  leaflet() %>% 
  addTiles() %>% 
  addPolygons(weight = 1, color = ~num_pal(number_species),
              fillOpacity = 0.8, group = "Number of species in ecoregion",
              # add labels that display indicator value and ecoregion id
              label = ~paste("ecoregion", ecoregion_id, "number of species =", 
                             number_species, 
                             sep = " "),
              # highlight polygons on hover
              highlightOptions = highlightOptions(weight = 5, color = "white",
                                                  bringToFront = TRUE)) %>%
  addPolygons(weight = 1, #color = ~num_pal(species_list),
              #fillOpacity = 0.8, 
              group = "Species in ecoregion",
              # add labels that display indicator value and ecoregion id
              label = ~paste("ecoregion", ecoregion_id, "Species list:", 
                             species_list, sep = " "),
              # highlight polygons on hover
              highlightOptions = highlightOptions(weight = 5, color = "white",
                                                  bringToFront = TRUE)) %>%
  addPolygons(weight = 1, #color = ~num_pal(species_list),
              #fillOpacity = 0.8, 
              group = "Species in ecoregion",
              # add labels that display indicator value and ecoregion id
              label = ~paste("ecoregion", ecoregion_id, "Species list:", 
                             species_list, sep = " "),
              # highlight polygons on hover
              highlightOptions = highlightOptions(weight = 5, color = "white",
                                                  bringToFront = TRUE)) %>%
  addLayersControl(#baseGroups = c("OSM", "Carto", "Esri"), 
                   overlayGroups = c("Number of species in ecoregion", 
                                     "Species identities"))


  
  
  
  