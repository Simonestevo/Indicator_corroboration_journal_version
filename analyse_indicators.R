
# Using R version 4.0.3

# Load packages ----

library(tidyverse)
library(ggplot2)
library(grid)
library(viridis)
library(png)
library(gridExtra)
library(reshape2)
library(GGally)
library(mapview)
library(psych)
library(e1071)
library(tm)
library(PerformanceAnalytics)

# Set input and output locations ----

create_new_database_version <- FALSE # Only set to true if you want to create an entirely new version from scratch
date <- Sys.Date()
country <- NA #"Australia" # If not subsetting, set as NA, e.g. country <- NA
inputs <- "N:/Quantitative-Ecology/Simone/extinction_test/outputs/version_3/2020-08-25_indicator_output_files"
save_outputs <- "yes" #only applies to maps, other things will always save
eco_version <- "ecoregions_2017"
parent_outputs <- "N:/Quantitative-Ecology/Simone/extinction_test/outputs"
#eco_version <- "official_teow_wwf"
indicator_columns <- c("indicator", "year", "ecoregion_id", "raw_indicator_value")

# Set up some ecoregions that we know how they should behave

east_australia <- 168 # Decline over time
amazon <- 473 # Decline over time
cardamom <- 223 # In good shape
mascarene <- 20 #Had a lot of extinctions


if (!is.na(country)) {
  
  location <- tolower(country)
  
} else {
  
  location <- "global"
  
}

# Set output directory

db_version <- tail(sort(list.files(parent_outputs)), 1)
db_interim <- list.files(file.path(parent_outputs,db_version))[
  grepl("interim",list.files(file.path(parent_outputs,db_version)))]
db_outputs <- list.files(file.path(parent_outputs,db_version))[
  grepl("database",list.files(file.path(parent_outputs,db_version)))]
ind_outputs <- list.files(file.path(parent_outputs,db_version))[
  grepl("indicator",list.files(file.path(parent_outputs,db_version)))]
analysis_outputs <- list.files(file.path(parent_outputs,db_version))[
  grepl("analysis",list.files(file.path(parent_outputs,db_version)))]

interim_outputs <- file.path(parent_outputs, db_version, db_interim)
outputs <- file.path(parent_outputs, db_version, db_outputs)

if( (length(analysis_outputs)) == 0 ) {
  
  analysis_outputs <- file.path(parent_outputs, db_version, paste(date,
                                 "_analysis_output_files",sep="") )
  
  dir.create(analysis_outputs, recursive = TRUE ) # create a new directory for today's outputs
  
  
} else {
  
  analysis_outputs <- file.path(parent_outputs, db_version, analysis_outputs)
  
}


# Load functions ----

# Read in data ----




# Remove some outliers

# lpi_values_subset <- lpi_values_subset %>%
#                      filter(raw_indicator_value < 2)
# 
# hfp_values_subset <- hfp_values %>%
#                      filter(raw_indicator_value < 40)
# 
# extinction_values_subset <- extinction_values %>%
#                             filter(ecoregion_id != 20) %>%
#                             filter(ecoregion_id != 637)



# Analyse indicators ----

# Data transformations ----

indicator_values_2 <- indicator_values_master %>%
  dplyr::select(ecoregion_id, indicator_year,
                raw_indicator_value) %>%
  distinct(.)

# * Convert to wide ----

indicators_wide <- indicator_values_2 %>%
  spread(key = indicator_year, 
         value = raw_indicator_value) 


names(indicators_wide) <- make.names(names(indicators_wide), 
                                     unique = TRUE)

summary(indicators_wide)

# * Invert ----
# For negatively valanced variables (where high values = negative outcome)

cols_to_invert <- c("HFP.2000","HFP.2005","HFP.2010", "HFP.2013",
                    "threatened.1980", "threatened.1990",
                    "threatened.2000", "threatened.2005",
                    "threatened.2010","threatened.2015", 
                    "extinct.1980", "extinct.1990",
                    "extinct.2000", "extinct.2005",
                    "extinct.2010", "extinct.2015")

col_index <- names(indicators_wide) %in% cols_to_invert
cols_min <- as.numeric(sapply(indicators_wide, min, na.rm = TRUE))
cols_max <- as.numeric(sapply(indicators_wide, max, na.rm = TRUE))


keys <- as.data.frame(col_index) %>%
  mutate(keys = ifelse(col_index == FALSE, 1, -1)) %>%
  select(keys) %>%
  pull(.)

indicators_wide$ecoregion_id <- as.numeric(indicators_wide$ecoregion_id)

indicators_wide_inv <- as.data.frame(reverse.code(keys,indicators_wide,
                                                  mini = cols_min, maxi = cols_max))

summary(indicators_wide_inv)

# test

ECO <- east_australia
HFP <- indicators_wide %>% filter(ecoregion_id == ECO) %>% select(HFP.2000, 
                                                                  HFP.2010)
HFP

HFP_inv <- indicators_wide_inv %>% filter(ecoregion_id == ECO) %>% 
  select(`HFP.2000-`, `HFP.2010-`)
HFP_inv

# * Centre ----

indicators_wide_inv_centred <- indicators_wide_inv %>%
  mutate_at(c(2:ncol(indicators_wide_inv)), funs(c(scale(.)))) 

summary(indicators_wide_inv_centred)

# ** Centred boxplots ----

indicators_scaled <- melt(indicators_wide_inv_centred, 
                          id.vars = 'ecoregion_id')

boxplots <- ggplot(indicators_scaled) +
  geom_boxplot(aes(x = variable, y = value)) +
  theme(axis.text.x=element_text(angle= 45,hjust=1))

boxplots

# * Transform ----
# https://www.datanovia.com/en/lessons/transform-data-to-normal-distribution-in-r/

# Check distributions

histograms <- lapply(indicators_wide_inv_centred, hist)
names(histograms) <- colnames(indicators_wide_inv_centred)
hist(indicators_wide_inv_centred$red.list.index.Amphibia.2000)

# Measure skew (nearly all negatively skewed)

indicators_skew <- sapply(indicators_wide_inv_centred, skewness, na.rm = TRUE)
indicators_skew_index <- as.data.frame(indicators_skew) %>%
  mutate(logtransform = ifelse(indicators_skew < -0.05, 
                               TRUE, FALSE)) %>%
  select(logtransform) %>%
  pull(.)

hist(indicators_wide_inv_centred$`human.footprint.index.2010-`)
new <- log10(max(indicators_wide_inv_centred$`human.footprint.index.2010-`+1, na.rm = TRUE) - indicators_wide_inv_centred$`human.footprint.index.2010-`) 
hist(new)
skewness(new, na.rm = TRUE)

indicators_transformed <- indicators_wide_inv_centred 

# HFP 2010
indicators_transformed$`human.footprint.index.2010-` <- log10(max(indicators_transformed$`human.footprint.index.2010-`+1, 
                                                                  na.rm = TRUE) - 
                                                                indicators_transformed$`human.footprint.index.2010-`) 
skewness(indicators_wide_inv_centred$`human.footprint.index.2010-`, na.rm = TRUE)
hist(indicators_wide_inv_centred$`human.footprint.index.2010-`)
skewness(indicators_transformed$`human.footprint.index.2010-`, na.rm = TRUE)
hist(indicators_transformed$`human.footprint.index.2010-`)

#HFP 1990

skewness(indicators_wide_inv_centred$`human.footprint.index.1990-`, na.rm = TRUE)
hist(indicators_wide_inv_centred$`human.footprint.index.1990-`)

indicators_transformed$`human.footprint.index.1990-` <- log10(max(indicators_transformed$`human.footprint.index.1990-`+1, 
                                                                  na.rm = TRUE) - 
                                                                indicators_transformed$`human.footprint.index.1990-`) 
skewness(indicators_transformed$`human.footprint.index.1990-`, na.rm = TRUE)
hist(indicators_transformed$`human.footprint.index.1990-`)

#At risk 1980

skewness(indicators_wide_inv_centred$`proportion.at.risk.1980-`, na.rm = TRUE)
hist(indicators_wide_inv_centred$`proportion.at.risk.1980-`)

indicators_transformed$`proportion.at.risk.1980-` <- 1/(max(indicators_wide_inv_centred$`proportion.at.risk.1980-`+1, 
                                                            na.rm = TRUE) - 
                                                          indicators_wide_inv_centred$`proportion.at.risk.1980-`) 
skewness(indicators_transformed$`proportion.at.risk.1980-`, na.rm = TRUE)
hist(indicators_transformed$`proportion.at.risk.1980-`)

#At risk 1990

to_transform <- indicators_wide_inv_centred$`proportion.at.risk.1990-`
transformed <- indicators_transformed$`proportion.at.risk.1990-`

skewness(to_transform, na.rm = TRUE)
hist(to_transform)

indicators_transformed$`proportion.at.risk.1990-` <- 1/(max(to_transform + 1, na.rm = TRUE) - to_transform) 
skewness(transformed, na.rm = TRUE)
hist(transformed)

#At risk 2000

to_transform <- indicators_wide_inv_centred$`proportion.at.risk.2000-`
transformed <- indicators_transformed$`proportion.at.risk.2000-`

skewness(to_transform, na.rm = TRUE)
hist(to_transform)

indicators_transformed$`proportion.at.risk.2000-` <- 1/(max(to_transform + 1, 
                                                            na.rm = TRUE) - to_transform) 
skewness(indicators_transformed$`proportion.at.risk.2000-`, na.rm = TRUE)
hist(indicators_transformed$`proportion.at.risk.2000-`)

#At risk 2005

to_transform <- indicators_wide_inv_centred$`proportion.at.risk.2005-`
transformed <- indicators_transformed$`proportion.at.risk.2005-`

skewness(to_transform, na.rm = TRUE)
hist(to_transform)

indicators_transformed$`proportion.at.risk.2005-` <- 1/(max(to_transform + 1, 
                                                            na.rm = TRUE) - to_transform) 
skewness(indicators_transformed$`proportion.at.risk.2005-`, na.rm = TRUE)
hist(indicators_transformed$`proportion.at.risk.2005-`)

#At risk 2010

to_transform <- indicators_wide_inv_centred$`proportion.at.risk.2010-`
transformed <- indicators_transformed$`proportion.at.risk.2010-`

skewness(to_transform, na.rm = TRUE)
hist(to_transform)

indicators_transformed$`proportion.at.risk.2010-` <- 1/(max(to_transform + 1, 
                                                            na.rm = TRUE) - to_transform) 
skewness(indicators_transformed$`proportion.at.risk.2010-`, na.rm = TRUE)
hist(indicators_transformed$`proportion.at.risk.2010-`)

#At risk 2015

to_transform <- indicators_wide_inv_centred$`proportion.at.risk.2015-`
transformed <- indicators_transformed$`proportion.at.risk.2015-`

skewness(to_transform, na.rm = TRUE)
hist(to_transform)

indicators_transformed$`proportion.at.risk.2015-` <- 1/(max(to_transform + 1, 
                                                            na.rm = TRUE) - to_transform) 
skewness(indicators_transformed$`proportion.at.risk.2015-`, na.rm = TRUE)
hist(indicators_transformed$`proportion.at.risk.2015-`)

# Red List Index Aves 2005

to_transform <- indicators_wide_inv_centred$red.list.index.Aves.2005
transformed <- indicators_transformed$red.list.index.Aves.2005

skewness(to_transform, na.rm = TRUE)
hist(to_transform)

indicators_transformed$red.list.index.Aves.2005 <- 1/(max(to_transform + 1, 
                                                          na.rm = TRUE) - to_transform) 
skewness(indicators_transformed$red.list.index.Aves.2005, na.rm = TRUE)
hist(indicators_transformed$red.list.index.Aves.2005)


# ** Transformed boxplots ----

indicators_transformed_long <- melt(indicators_transformed, 
                                    id.vars = 'ecoregion_id')

transformed_boxplots <- ggplot(indicators_transformed_long) +
  geom_boxplot(aes(x = variable, y = value)) +
  theme(axis.text.x=element_text(angle= 45,hjust=1))

transformed_boxplots

indicator_values_transformed <- indicators_transformed_long %>%
  merge(indicator_values_master[c("ecoregion_id",
                                  "country",
                                  "realm")],
        by = "ecoregion_id") %>%
  rename(indicator_year = variable,
         raw_indicator_value = value) %>% # TODO: Change to transformed
  mutate(year = str_sub(indicator_year, start= -4)) %>%
  mutate(indicator_year = as.character(indicator_year)) %>%
  mutate(year = ifelse(year == "980-", "1980",
                       ifelse(year == "015-",
                              "2015",
                              ifelse(year == "000-",
                                     "2000",
                                     ifelse(year == "990-",
                                            "1990",
                                            ifelse(year == "005-", 
                                                   "2005",
                                                   ifelse(year == "010-", 
                                                          "2010", 
                                                          ifelse(year == "2005", 
                                                                 "2005", NA)))))))) %>%
  mutate(indicator = removeNumbers(indicator_year)) %>%
  mutate(indicator = str_replace_all(indicator,
                                     '[[:punct:]]',' ')) 

saveRDS(indicator_values_transformed, file.path(indicator_outputs,
                                                paste(location, eco_version,
                                                      "indicator_vals_transformed.rds",
                                                      sep = "_")))

# Correlations ----

# Non-transformed

# indicator_names <- c("Ecoregion_id", "BII_A_2005","BHI_2005", "BHI_2010", "BHI_2015",
#                      "HFP_1990", "HFP_2010","Endangered_1980", "Endangered_1990",
#                      "At_risk_2000", "At_risk_2005", "At_risk_2010", 
#                      "At_risk_2015",  "Extinction_1980",   "Extinction_1990",
#                      "Extinction_2000", "Extinction_2005",   "Extinction_2010",
#                      "Extinction_2015", "RLI_Amph_1980", "RLI_Amph_2000",
#                      "RLI_Birds_1980", "RLI_Birds_1990", "RLI_Birds_2000",
#                      "RLI_Birds_2005", "RLI_Birds_2010", "RLI_Mamm_1990",
#                      "RLI_Mamm_2005", "BII_R_2005")


# indicator_names <- c("Ecoregion_id", "BII","BHI", "BHI_2010", "BHI_2015",
#                      "HFP_1990", "HFP","Endangered_1980", "Endangered_1990",
#                      "At_risk_2000", "Endangered", "At_risk_2010", 
#                      "At_risk_2015",  "Extinction_1980",   "Extinction_1990",
#                      "Extinction_2000", "Extinction_2005",   "Extinction_2010",
#                      "Extinction_2015", "RLI_Amph_1980", "RLI_Amph",
#                      "RLI_Birds_1980", "RLI_Birds_1990", "RLI_Birds_2000",
#                      "RLI_Birds", "RLI_Birds_2010", "RLI_Mamm_1990",
#                      "RLI_Mammals", "BII_R_2005")
# 
# colnames(indicators_wide_inv_centred) <- indicator_names


indicator_matrix <- indicators_wide_inv_centred %>%
  select(-ecoregion_id) %>%
  na.omit(.)

# Split by location

# Get 2005 only
indicators_wide_inv_centred_2005 <- indicators_wide_inv_centred[,c(1,2,5,6,10,14,17,21)]


indicators_wic <- indicators_wide_inv_centred_2005 %>%
  merge(ecoregion_map[c("ECO_ID", "REALM")],
        by.x = "ecoregion_id",
        by.y = "ECO_ID") 

indicators_wic_realm <- split(indicators_wic, indicators_wic$REALM)

indicators_wic_oceania <- indicators_wic_realm[["Oceania"]]

indicators_wic_australasia <- indicators_wic_realm[["Australasia"]]

indicators_wic_realm[["Global"]] <- indicators_wic

realm_matrices <- list()

for (i in seq_along(indicators_wic_realm)) {
  
  realm_matrices[[i]] <- indicators_wic_realm[[i]] %>%
    select(-ecoregion_id, -REALM, - geometry) %>%
    na.omit(.)
  
}

names(realm_matrices) <- names(indicators_wic_realm)

# Get pictoral correlation matrices

realm_correlations <- list()

for (i in seq_along(realm_matrices)) {
  
  name <- paste(names(realm_matrices)[i], "indicator correlation matrix", sep = " ")
  
  subset_matrix <- realm_matrices[[i]]#[, c(1,2,6,10,20,24,27)]
  
  if(nrow(subset_matrix) == 0) {
    
    print("no data for this realm")
    
  } else {
    
    correlation_matrix <- ggcorr(subset_matrix, label = FALSE, nbreaks = 9) + 
      ggtitle(name)
    
    ggsave(file.path(indicator_outputs, 
                     paste(name, "indicator_correlation_matrix.png", sep = "_")),
           correlation_matrix, device = "png")
    
    realm_correlations[[i]] <- correlation_matrix
  }
}

realm_correlations[[9]]




pearson_correlations_all_years <- as.data.frame(cor(indicator_matrix, 
                                                    method = "pearson"))

spearman_correlations_all_years <- as.data.frame(cor(indicator_matrix, 
                                                     method = "spearman"))

# All indicators with 2005 data

indicator_matrix_2005 <- indicator_matrix[, c(1,2,10,16,23,27,28)]

indicator_correlations_2005 <- ggcorr(indicator_matrix_2005)

# All indicators with 2010 data

indicator_matrix_2010 <- indicator_matrix[, c(3,6,11,17,25)]

indicator_correlations_2010 <- ggcorr(indicator_matrix_2010)

# RLI all timepoints and 1990 HFP

indicator_matrix_rli <- indicator_matrix[, c(5,21,22,23,24,25)]

indicator_correlations_rli <- ggcorr(indicator_matrix_rli)

# BHI all timepoints and 1990 HFP

indicator_matrix_bhi <- indicator_matrix[, c(5,2, 3, 4)]

indicator_correlations_bhi <- ggcorr(indicator_matrix_bhi)

# Extinctions all timepoints and 1990 HFP

indicator_matrix_ex <- indicator_matrix[, c(5,13,14,15,16,17,18)]

indicator_correlations_ex <- ggcorr(indicator_matrix_ex)

# Indicators from different data sources 2005 - 2010

indicator_matrix_ind <- indicator_matrix[, c(1,2,6,10,20,24,27)]

indicator_correlations_ind <- ggcorr(indicator_matrix_ind, palette = "BuPu",
                                     label = TRUE)


# RLI Birds 1980 and following extinctions

indicator_matrix_er <- indicator_matrix[, c(21, 13, 14,15,16,17,18)]

indicator_correlations_er <- ggcorr(indicator_matrix_er)


# Transformed

indicator_matrix_t <- indicators_transformed %>%
  select(-ecoregion_id) %>%
  na.omit(.)

pearson_correlations_all_years_t <- as.data.frame(cor(indicator_matrix_t, 
                                                      method = "pearson"))

spearman_correlations_all_years_t <- as.data.frame(cor(indicator_matrix_t, 
                                                       method = "spearman"))

# mutate(HFP_adjusted_old = ifelse(raw_indicator_value > 33.29037,
#                                  33.29037, raw_indicator_value)) %>%
# mutate(RLI_adjusted_old = ifelse(raw_indicator_value == 0, NA,
#                ifelse(raw_indicator_value > 0 & 
#                         raw_indicator_value < 0.9538,
#                       0.9538, raw_indicator_value))) %>%
# mutate(year = as.numeric(year)) %>%
# mutate(HFP_adjusted_old = ifelse(indicator != 
#                           "mean human footprint index", NA,
#                           HFP_adjusted_old)) %>%
# mutate(RLI_adjusted_old = ifelse(grepl("red list index Aves",
#                                 indicator), 
#                                 RLI_adjusted_old, NA)) %>%
# mutate(HFP_adjusted =
#          pmin(pmax(raw_indicator_value,quantile(raw_indicator_value,
#                                          .005, na.rm = TRUE)),
#               quantile(raw_indicator_value, .995,
#                        na.rm = TRUE))) %>%
# mutate(HFP_scaled_adjusted = scale_to_1(HFP_adjusted)) %>%
# mutate(HFP_scaled_adjusted_inverted =
#          1 - HFP_scaled_adjusted) %>%
# mutate(scaled = scale_to_1(raw_indicator_value)) %>%
# mutate(inverted = 1 - raw_indicator_value) %>%
# mutate(scaled_inverted = 1 - scaled) %>%
# mutate(RLI_adjusted = ifelse(raw_indicator_value == 0, NA,
#                              pmin(pmax(raw_indicator_value,
#                   quantile(raw_indicator_value, .05, 
#                            na.rm = TRUE))))) %>%
# mutate(RLI_adjusted_inverted = 1 - RLI_adjusted)

# Scatterplots ----

indicator_values_master_1 <- indicator_values_master

indicator_values_master_2 <- indicator_values_master_1 %>%
  mutate(key2 = paste(ecoregion_id, indicator_year, sep = " ")) %>%
  mutate(key = make.names(key2)) %>%
  mutate(key = str_remove(key, "X")) %>%
  select(-key2) 

indicators_scaled <- indicators_scaled %>%
  mutate(variable = str_remove(variable, "-"))

indicator_values_scaled <- indicators_scaled %>%
  mutate(key = paste(ecoregion_id, variable, sep = ".")) %>%
  rename(centred_indicator_value = value)

unique(indicator_values_scaled$variable)
unique(indicator_values_master$indicator_year)

indicator_values_master <- indicator_values_master_2 %>%
  merge(indicator_values_scaled[c("key", 
                                  "centred_indicator_value")], by = "key") 

unique(indicator_values_master$indicator_year)

# * Subset by indicator ----

indicator_values_master_list <- split(indicator_values_master, 
                                      indicator_values_master$indicator)

indicator_values_master_wide <- list()
indicator_by_time_plots <- list()

for (i in seq_along(indicator_values_master_list)) {
  
  wide_indicator_data <- indicator_values_master_list[[i]] %>%
    select(ecoregion_id, indicator_year, 
           centred_indicator_value) %>%
    distinct(.) %>%
    spread(key = indicator_year, 
           value = centred_indicator_value) %>%
    select(-ecoregion_id)
  
  indicator_plot <- ggpairs(wide_indicator_data)
  
  indicator_name <- indicator_values_master_list[[i]][1,2]
  
  indicator_values_master_wide[[i]] <- wide_indicator_data
  indicator_by_time_plots[[i]] <- indicator_plot
  
  scatterplot_directory <- file.path(indicator_outputs, 
                                     "scatterplots")
  
  if( !dir.exists( scatterplot_directory ) ) {
    
    dir.create( scatterplot_directory, recursive = TRUE )
    
  }
  
  ggsave(file.path(scatterplot_directory, paste(location, eco_version,
                                                indicator_name, 
                                                "time_scatterplots.png", 
                                                sep = "_")),
         indicator_plot,  device = "png")
  
}


# * Subset by time ----


# Non-transformed

indicators_by_year <- split(indicator_values_master, indicator_values_master$year)
indicators_by_year_names <- paste(location, names(indicators_by_year), sep = "_")

scatterplots_by_year <- list()

for (i in seq_along(indicators_by_year)){
  
  scatterplots_by_year[[i]] <- produce_scatterplots(indicators_by_year[[i]],
                                                    indicators_by_year_names[[i]],
                                                    save = TRUE, "centred")
  
}

scatterplots_by_year[[4]]


# Centred and inverted 

indicators_wic <- indicators_wic %>%
  rename(ecoregion_id = Ecoregion_id) %>%
  select(-geometry, -REALM)

indicators_wic_long <- melt(indicators_wic, 
                            id.vars = 'ecoregion_id')


indicator_values_wic <- indicators_wic_long  %>%
  rename(indicator_year = variable,
         raw_indicator_value = value) %>% # TODO: Change to transformed
  mutate(year = str_sub(indicator_year, start= -4)) %>%
  mutate(indicator_year = as.character(indicator_year)) %>%
  mutate(indicator = removeNumbers(indicator_year)) %>%
  mutate(indicator = str_replace_all(indicator,
                                     '[[:punct:]]',' '))

indicators_by_year <- split(indicator_values_wic, indicator_values_wic$year)
indicators_by_year_names <- paste(location, names(indicators_by_year), sep = "_")

scatterplots_by_year <- list()

for (i in seq_along(indicators_by_year)){
  
  scatterplots_by_year[[i]] <- produce_scatterplots(indicators_by_year[[i]],
                                                    indicators_by_year_names[[i]],
                                                    save = FALSE)
  
}

scatterplots_by_year[[1]]


# RLI and HFP 2005 only ---

indicator_values_rli_hfp <- indicator_values_wic %>%
  filter(indicator_year == "HFP.2005-" |
           indicator_year ==  "RLI.2005") 

rli_hfp_corr <- produce_scatterplots(indicator_values_rli_hfp, "test", save = FALSE,"raw")

indicators_by_year <- split(indicator_values_wic, indicator_values_wic$year)
indicators_by_year_names <- paste(location, names(indicators_by_year), sep = "_")

scatterplots_by_year <- list()

for (i in seq_along(indicators_by_year)){
  
  scatterplots_by_year[[i]] <- produce_scatterplots(indicators_by_year[[i]],
                                                    indicators_by_year_names[[i]],
                                                    save = FALSE)
  
}

scatterplots_by_year[[1]]

# Oceania ----

indicators_wic <- indicators_wic_oceania %>%
  rename(ecoregion_id = Ecoregion_id) %>%
  select(-geometry, -REALM)

indicators_wic_long <- melt(indicators_wic, 
                            id.vars = 'ecoregion_id')


indicator_values_wic <- indicators_wic_long  %>%
  rename(indicator_year = variable,
         raw_indicator_value = value) %>% # TODO: Change to transformed
  mutate(year = NA) %>%
  mutate(indicator = as.character(indicator_year)) 

oceania_comparisons <- indicator_values_wic %>%
  filter(indicator == "BII"| 
           indicator == "BHI"|
           indicator == "HFP"|
           indicator == "Endangered"|
           indicator == "RLI_Amph"|
           indicator == "RLI_Birds"|
           indicator == "RLI_Mammals")

length(unique(oceania_comparisons$ecoregion_id))

oceania_bhi_amph <- oceania_comparisons %>%
  filter(indicator == "RLI_Amph"| 
           indicator == "BHI")

produce_scatterplots(oceania_bhi_amph, "oceania_bhi_amph", save = TRUE)

oceania_scatterplots <- produce_scatterplots(oceania_comparisons, "oceania", save =FALSE)

plot(indicators_wic_oceania$HFP,indicators_wic_oceania$RLI_Amph)

oceania_master <- indicator_values_master %>%
  filter(realm == "Oceania")

oceania_master_bhi <- oceania_master %>%
  filter(indicator == "BHI plants") %>%
  select(-country) %>%
  distinct(.) %>%
  group_by(ecoregion_id, year) %>%
  arrange(.)

ggplot(oceania_master_bhi) +
  geom_line(aes(x = year, 
                y = raw_indicator_value, 
                color = ecoregion_id))


# Australasia ----

indicators_wic <- indicators_wic_australasia %>%
  rename(ecoregion_id = Ecoregion_id) %>%
  select(-geometry, -REALM)

indicators_wic_long <- melt(indicators_wic, 
                            id.vars = 'ecoregion_id')


indicator_values_wic <- indicators_wic_long  %>%
  rename(indicator_year = variable,
         raw_indicator_value = value) %>% # TODO: Change to transformed
  mutate(year = NA) %>%
  mutate(indicator = as.character(indicator_year)) 

australasia_comparisons <- indicator_values_wic %>%
  filter(indicator == "BII"| 
           indicator == "BHI"|
           indicator == "HFP"|
           indicator == "Endangered"|
           indicator == "RLI_Amph"|
           indicator == "RLI_Birds"|
           indicator == "RLI_Mammals")

length(unique(oceania_comparisons$ecoregion_id))

aus_bhi_amph <- australasia_comparisons %>%
  filter(indicator == "RLI_Amph"| 
           indicator == "BHI")

aus_bhi_amph$RLI_Amph_t <- 1/(max(aus_bhi_amph$RLI_Amph_t + 1, na.rm = TRUE) - 
                                aus_bhi_amph$RLI_Amph_t)

1/(max(indicators_wide_inv_centred$`proportion.at.risk.1980-`+1, 
       na.rm = TRUE) - 
     indicators_wide_inv_centred$`proportion.at.risk.1980-`) 

produce_scatterplots(aus_bhi_amph, "aus_bhi_amph", save = TRUE)

oceania_scatterplots <- produce_scatterplots(oceania_comparisons, "oceania", save =FALSE)

plot(indicators_wic_oceania$HFP,indicators_wic_oceania$RLI_Amph)

oceania_master <- indicator_values_master %>%
  filter(realm == "Oceania")

oceania_master_bhi <- oceania_master %>%
  filter(indicator == "BHI plants") %>%
  select(-country) %>%
  distinct(.) %>%
  group_by(ecoregion_id, year) %>%
  arrange(.)

ggplot(oceania_master_bhi) +
  geom_line(aes(x = year, 
                y = raw_indicator_value, 
                color = ecoregion_id))

# Transformed

indicators_by_year <- split(indicators_values_transformed, indicator_values_transformed$year)
indicators_by_year_names <- paste(location, "transformed", names(indicators_by_year), sep = "_")

scatterplots_by_year <- list()

for (i in seq_along(indicators_by_year)) {
  
  scatterplots_by_year[[i]] <- produce_scatterplots(indicators_by_year[[i]],
                                                    indicators_by_year_names[[i]],
                                                    save = FALSE)
  
}

scatterplots_by_year[[2]]

# * Subset by realm ----

# * Subset by time by realm ----

# Split data into nested list (level 1 = year, level 2 = realm)

indicators_by_year_by_realm <- list()

for (i in seq_along(indicators_by_year)) {
  
  indicators_timepoint <- indicators_by_year[[i]]
  
  indicators_by_year_by_realm[[i]] <- split(indicators_timepoint, 
                                            indicators_timepoint$realm)
}

indicators_by_year_by_realm <- flatten(indicators_by_year_by_realm)
realm_names <- names(indicators_by_year_by_realm)
timepoints <- names(indicators_by_year)
time_names <- rep(timepoints, each = length(unique(realm_names)))
indicators_by_year_by_realm_names <- paste(realm_names, time_names , sep = "_")

scatterplots_by_year_by_realm <- list()

for (i in seq_along(indicators_by_year_by_realm)){
  
  scatterplots_by_year_by_realm[[i]] <- produce_scatterplots(indicators_by_year_by_realm[[i]],
                                                             indicators_by_year_by_realm_names[[i]],
                                                             save = TRUE)
  
}

names(scatterplots_by_year_by_realm) <- indicators_by_year_by_realm_names

# Look at Australiasia 2005

scatterplots_by_year_by_realm[[27]]

# Look at Oceania 2005

scatterplots_by_year_by_realm[[31]]

# * Subset effect of land use on extinction risk ----

land_use_atrisk_values <- indicator_values_master %>%
  filter(indicator_year == "human footprint index 1990"|
           indicator_year == "proportion at risk 1990"|
           indicator_year == "proportion at risk 2000"|
           indicator_year == "proportion at risk 2005"|
           indicator_year == "proportion at risk 2010"|
           indicator_year == "proportion at risk 2015")

hfp_vs_atrisk_scatterplots <- produce_scatterplots(land_use_atrisk_values, 
                                                   paste(location, eco_version,
                                                         "HFP1990_vs_AtRisk_time_series",
                                                         sep = "_"), save = TRUE)

# * Subset effect of land use on extinction  ----

land_use_extinction_values <- indicator_values_master %>%
  filter(indicator_year == "human footprint index 1990"|
           indicator_year == "proportion extinct 1990"|
           indicator_year == "proportion extinct 2000"|
           indicator_year == "proportion extinct 2005"|
           indicator_year == "proportion extinct 2010"|
           indicator_year == "proportion extinct 2015")

hfp_vs_extinction_scatterplots <- produce_scatterplots(land_use_extinction_values, 
                                                       paste(location, eco_version,
                                                             "HFP1990_vs_Extinction_time_series",
                                                             sep = "_"), save = TRUE)


hfp_vs_extinction_scatterplots

# * Subset only unrelated indicators  ----

subset_values <- indicator_values_master %>%
  filter(indicator_year == "human footprint index 1990"|
           indicator_year == "human footprint index 2010"|
           indicator_year == "abundance biodiversity intactness index 2005"|
           indicator_year == "richness biodiversity intactness index 2005"|
           indicator_year == "proportion extinct 2005")

subset_scatterplots <- produce_scatterplots(subset_values, 
                                            paste(location, eco_version,
                                                  "subset_indicators",
                                                  sep = "_"), save = FALSE)


subset_scatterplots


ggplot(data = indicator_values_master_wide_scaled, 
       aes(x = abundance.biodiversity.intactness.index.2005, 
           y = human.footprint.index.1990)) +
  geom_point() + 
  coord_trans(x = "log10", y = "log10")

# Pairwise correlations and scatterplots ----

#' TODO: Can we use ggpair with other ggplot syntax? Including attributes but
#' also transformations? https://www.r-graph-gallery.com/199-correlation-matrix-with-ggally.html

correlations <- list()
scatterplots <- list()
years <- list()

for (i in seq_along(indicator_values_time_list)) {
  
  # Get the year name
  
  year <- range(indicator_values_time_list[[i]]$year)
  
  # Remove year column from dataframe
  
  step1 <- indicator_values_time_list[[i]] %>% 
    select(ecoregion_id, 
           raw_indicator_value, 
           indicator_year, 
           -year) %>%
    distinct(.)
  
  # Convert into wide format
  
  # step2 <- dcast(step1, ecoregion_code ~ indicator, 
  #                value.var = "raw_indicator_value", sum)
  
  step2 <- step1 %>%
    spread(key = indicator_year, 
           value = raw_indicator_value) 
  
  # Convert into matrix without id variables
  
  step3 <- as.matrix(step2[2:ncol(step2)])
  
  # Check if there's enough indicators to analyse correlation (ie more than 2)
  
  if(ncol(step3) == 1) {
    
    print(paste("Year", year, "does not have values for more than one indicator to test correlations between"))
    
    rm(step1, step2, step3)
    
  } else {
    
    # Analyse correlations between indicators 
    
    step4 <- step3[complete.cases(step3),]
    
    correlations[[i]] <- as.data.frame(cor(step4, method = "pearson"))
    
    scatterplots[[i]] <- ggpairs(as.data.frame(step4))
    
    years[[i]] <- year
    
    rm(step1, step2, step3, step4)
    
  }
  
}

# Remove the empty lists with no data in them

names(correlations) <- years

correlations <- list.clean(correlations, fun = is.null, recursive = FALSE)

names(scatterplots) <- years

if(save_outputs == "yes") {
  
  for (i in seq_along(scatterplots)) {
    
    ggsave(file.path(outputs, paste(as.character(years[[i]]),"_scatterplot",
                                    ".png", sep = "")), 
           scatterplots[[i]],  device = "png")
    
  }
}
