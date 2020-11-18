
# Using R version 4.0.3

# Load packages ----

# install.packages("tidyverse", dependencies = TRUE)
# install.packages("ggplot2", dependencies = TRUE)
# #install.packages("grid")
# install.packages("viridis", dependencies = TRUE)
# install.packages("png", dependencies = TRUE)
# install.packages("gridExtra", dependencies = TRUE)
# install.packages("reshape2", dependencies = TRUE)
# install.packages("GGally", dependencies = TRUE)
# install.packages("mapview", dependencies = TRUE)
# install.packages("psych", dependencies = TRUE)
# install.packages("e1071", dependencies = TRUE)
# install.packages("tm", dependencies = TRUE)
# install.packages("PerformanceAnalytics", dependencies = TRUE)

library(tidyverse)
library(ggplot2)
#library(grid)
library(viridis)
library(png)
library(gridExtra)
library(reshape2)
library(GGally)
library(mapview)
library(psych) #not loaded
library(e1071)
library(tm) #notloaded
library(PerformanceAnalytics) #not loaded
library(data.table)
library(rlist)
library(factoextra)
library(ggcorrplot)
library(MASS)

# Set input and output locations ----

create_new_database_version <- FALSE # Only set to true if you want to create an entirely new version from scratch
date <- Sys.Date()
country <- NA #"Australia" # If not subsetting, set as NA, e.g. country <- NA
analysis_inputs <- "N:/Quantitative-Ecology/Simone/extinction_test/outputs/version_3/2020-08-25_indicator_output_files"
save_outputs <- "yes" #only applies to maps, other things will always save
eco_version <- "ecoregions_2017"
parent_outputs <- "N:/Quantitative-Ecology/Simone/extinction_test/outputs"
#eco_version <- "official_teow_wwf"
indicator_columns <- c("indicator", "year", "ecoregion_id", "raw_indicator_value")
timepoint <- 2005


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

# Create a folder for the day because you'll probably make lots of version


current_analysis_outputs <- file.path(analysis_outputs,paste(date,
                                            "_analysis_output_files",sep=""))

dir.create(current_analysis_outputs)


# Load functions ----

# Function to scale the values of a vector from 0 to 1

scale_to_1 <- function(vector){
  
  (vector-min(vector, na.rm = TRUE))/
    (max(vector, na.rm = TRUE)-min(vector, na.rm = TRUE))
}

add_grouping_variable <- function(variable, indicator_data, ecoregion_data) {
  
  out <- indicator_data %>%
    merge(ecoregion_data[c("ecoregion_id", variable)],
          by = "ecoregion_id") 
  
  return(out)
  
}

# Read in data ----

# Indicator data

raw_indicators_long <- readRDS(file.path(analysis_inputs,
                          "global_ecoregions_2017_indicator_values_master.rds"))

raw_indicators_long$indicator_year <- str_replace(raw_indicators_long$indicator_year, " ", "_")

indicators <- unique(raw_indicators_long$indicator_year)

raw_indicators_wide <- readRDS(file.path(analysis_inputs,
                          "global_ecoregions_2017_indicator_values_master_wide.rds"))

# Ecoregion data

raw_ecoregions_long <- readRDS(file.path(analysis_inputs,
                         "global_ecoregions_2017_ecoregion_values_master.rds")) 
raw_ecoregions_wide <- readRDS(file.path(analysis_inputs,
                         "global_ecoregions_2017_ecoregion_values_master_wide.rds"))

# Threat scheme

threat_scheme <- read.csv(file.path("N:/Quantitative-Ecology/Simone/extinction_test/inputs/iucn_threats\\iucn_threat_classification_scheme.csv"))

headline_threats <- threat_scheme %>%
  dplyr::select(headline_name) %>%
  distinct(.) %>%
  pull()

headline_threats <- c(headline_threats, "All threats")

# Ecoregion map

# ecoregion_map_all <- readRDS(paste(
#                      file.path("N:/Quantitative-Ecology/Simone/extinction_test/inputs",
#                                "ecoregions_2017"),
#                                "Ecoregions2017valid.rds"))
# 
# ecoregion_map <- ecoregion_map_all %>%
#   dplyr::select(ECO_ID, ECO_NAME, OBJECTID, REALM, geometry)
# 
# ecoregion_map_renamed <- ecoregion_map %>% dplyr::rename(ecoregion_id = ECO_ID)
# 
# rm(ecoregion_map_all, ecoregion_map)
# 
# # Countries
# 
# ecoregion_countries <- readRDS(file.path(parent_outputs, "version_3", 
#                         "2020-08-10_database_output_files",
#                         "ecoregion_country_data.rds"))

### TEMPORARY CODE - REMOVE LPI BC IT LOOKS WEIRD AND LOSES ALOT OF DATA ###

raw_indicators_long <- raw_indicators_long %>%
                       filter(indicator != "LPI")

if (!is.na(timepoint)) {

raw_indicators_long <- raw_indicators_long %>%
                         filter(year == timepoint) %>%
                         dplyr::select(ecoregion_id, indicator_year,
                                      raw_indicator_value) %>%
                         distinct(.)

raw_indicators_wide <- raw_indicators_long %>%
                         spread(key = indicator_year, 
                                value = raw_indicator_value) 

indicators <- unique(raw_indicators_long$indicator_year)

}

# Ecoregion value cleaning ----

# Convert numeric to factors

ecoregions_wide <- raw_ecoregions_wide %>%
                   mutate(LPI_records = as.numeric(ifelse(is.na(LPI_records), 
                                               0, LPI_records)))

ecoregions_wide$scientific.publications.factor <- cut(ecoregions_wide$mean.scientific.publications, breaks = 3, 
                            labels = c("Few_publications", 
                                       "Moderate_publications",
                                       "Many_publications"))

ecoregions_wide$area.factor <- cut(ecoregions_wide$ecoregion.area.km.sq, breaks = 5, 
                                  labels = c("Very_small_ecoregions", 
                                             "Small_ecoregions", 
                                             "Medium_ecoregions",
                                             "Large_ecoregions", 
                                             "Very_large_ecoregions"))

lpi_breaks <- c(-1, 50, 150, 400)
ecoregions_wide$lpi.records.factor <- cut(ecoregions_wide$LPI_records, 
                                          breaks = lpi_breaks, 
                                   labels = c("Few (< 50) LPI records", 
                                              "Moderate (between 50 and 100) LPI records",
                                              "Many (> 100) LPI records"))

ecoregions_wide$rli.records.factor <- cut(ecoregions_wide$RLI_records, breaks = 3, 
                                   labels = c("Few (< 590) RLI records", 
                                              "Moderate (between 590 & 1200) RLI records",
                                              "Many (> 1200) RLI records"))


ecoregions_wide <- ecoregions_wide %>% 
                   mutate(included.in.HFP = ifelse(included.in.HFP == 1, 
                                                   "Threat related to HFP",
                                                   ifelse(included.in.HFP == 0,
                                                          "Threat external to HFP",
                                                          NA)))

ecoregions_wide <- ecoregions_wide %>%
        mutate(scenario = as.factor(paste(rli.records.factor, included.in.HFP, sep = " & ")),
               scenario.numeric = as.factor(as.numeric(scenario))) 

ecoregions_wide_countries <- ecoregions_wide %>%
                             merge(ecoregion_countries[c("ECO_ID", "CNTRY_NAME")],
                                   by.x = "ecoregion_id", by.y = "ECO_ID") %>%
                             rename(country = CNTRY_NAME)

rm(raw_ecoregions_long, raw_ecoregions_wide)

# Indicator data cleaning ----

# * Invert ----
# For negatively valanced variables (where high values = negative outcome)

negatives_index <- names(raw_indicators_wide) %like% 'HFP|extinct|threatened'
negatives <- names(raw_indicators_wide)[negatives_index]

cols_min <- as.numeric(sapply(raw_indicators_wide, min, na.rm = TRUE))
cols_max <- as.numeric(sapply(raw_indicators_wide, max, na.rm = TRUE))


keys <- as.data.frame(negatives_index) %>%
        mutate(keys = ifelse(negatives_index == FALSE, 1, -1)) %>%
        dplyr::select(keys) %>%
        pull(.)

indicators_wide <- as.data.frame(reverse.code(keys,raw_indicators_wide,
                                 mini = cols_min, maxi = cols_max))

# Remove the little negative thing at the end of reversed column names, otherwise
# they become difficult to use with dplyr

names(indicators_wide) <- str_replace(names(indicators_wide), "-", "")

# test

# ECO <- east_australia # This ecoregion has a low HFP value (bad outcome)
# HFP <- raw_indicators_wide %>% filter(ecoregion_id == ECO) %>% select(`HFP 2005`)
# HFP
# 
# # Once inverted, it should have a high ecoregion value to represent the bad outcome
# HFP_inv <- indicators_wide %>% filter(ecoregion_id == ECO) %>% 
#   select(`HFP 2005-`)
# HFP_inv

# * Centre ----

#TODO: Do we need to transform any variables? bc probably need to do so before scaling
# https://www.datanovia.com/en/lessons/transform-data-to-normal-distribution-in-r/

indicators_wide_centred <- indicators_wide %>%
  mutate_at(c(2:ncol(indicators_wide)), funs(c(scale(.)))) 

summary(indicators_wide_centred)

# ** Centred boxplots ----

#' TODO: IMPORTANT - DECIDE WHETHER TO REMOVE LPI ALTOGETHER

indicator_boxplot_data <- reshape2::melt(indicators_wide_centred, 
                          id.vars = 'ecoregion_id')

boxplots <- ggplot(indicator_boxplot_data) +
            geom_boxplot(aes(x = variable, y = value)) +
            theme(axis.text.x = element_text(angle= 45,hjust=1))

boxplots

boxplot(indicators_wide_centred$`HFP 2005-`)

# * Manage outliers ----

indicators_wide_centred_trunc <- indicators_wide_centred %>%
                           # filter(`LPI 2005` < quantile(`LPI 2005`, 
                           #                              0.99, na.rm = TRUE)) %>%
                           # filter(`extinct 2005-` < quantile(`extinct 2005-`, 
                           #                              0.99, na.rm = TRUE)) %>%
                           filter(HFP_2005 < quantile(HFP_2005, 
                                                            0.99, na.rm = TRUE)) 

summary(indicators_wide_centred_trunc)

indicator_boxplot_data_2 <- reshape2::melt(indicators_wide_centred_trunc, 
                                           id.vars = 'ecoregion_id')

boxplots_2 <- ggplot(indicator_boxplot_data_2) +
              geom_boxplot(aes(x = variable, y = value)) +
              theme(axis.text.x = element_text(angle = 45,hjust = 1))

boxplots_2

rm(raw_indicators_long, raw_indicators_wide)

# * Transform ----

# ??

# Finalise analysis data ----

analysis_input_data <- indicators_wide_centred_trunc %>%
  merge(ecoregions_wide, by = "ecoregion_id") 

# Scatterplots ----

# analysis_input_data <- analysis_input_data_all
# analysis_input_data_all <- analysis_input_data

indicator_combinations <- combn(indicators, 2)

indicator_list <- as.list(as.data.frame(indicator_combinations))

# Loop through each combination of indicators and make a scatterplot.  Note
# these are really ugly with all the labels and stuff - they're just to identify
# the weird ecoregions, we can tidy up anything that would actually go into the paper

scatterplots <- list()

for (i in seq_along(indicator_list)) {
  
scatterplot_data <- analysis_input_data[, c("ecoregion_id", 
                                            "island.status",
                                            indicator_list[[i]][1],
                                            indicator_list[[i]][2])]
labx <- names(scatterplot_data)[3]
laby <- names(scatterplot_data)[4]
names(scatterplot_data) <- c("ecoregion_id", "vargroup", "varx", "vary")

scatterplot <- ggplot(scatterplot_data,aes(x = varx, 
                                           y = vary, 
                                           color = vargroup)) +
               geom_point() +
               geom_text(label = analysis_input_data$ecoregion_id) + 
               geom_smooth(method=lm) +
               labs(x= labx,
                    y = laby)

ggsave(file.path(current_analysis_outputs, paste(labx, "_x_", laby, ".png", sep = "")),
       scatterplot, device = "png")

scatterplots[[i]] <- scatterplot

}

scatterplots[[1]]
scatterplots[[2]]
scatterplots[[3]]
scatterplots[[4]]
scatterplots[[5]]
scatterplots[[6]]
scatterplots[[7]]
scatterplots[[8]]
scatterplots[[9]]
scatterplots[[10]]
scatterplots[[11]]
scatterplots[[12]]
scatterplots[[13]]
scatterplots[[14]]
scatterplots[[15]]
scatterplots[[16]]
scatterplots[[17]]
scatterplots[[18]]
scatterplots[[19]]
scatterplots[[20]]

# Get the ecoregions with the strangest results (eg really high HFP, low RLI)

weird_ecoregions <- analysis_input_data %>%
                    filter(ecoregion_id %in% c(325, 400, 637, 20, 632, 548))


# Correlation plots ----

# * Country ----

names(raw_ecoregions_wide)

country_indicators <- add_grouping_variable("country", indicators_wide_centred_trunc,
                                          ecoregions_wide_countries )

country_indicator_list <- split(country_indicators, country_indicators$country)

country_indicator_list[["All countries"]] <- indicators_wide_centred_trunc %>%
  mutate(country = "All countrie")

country_matrices <- list()

for (i in seq_along(country_indicator_list)) {
  
  country_matrices[[i]] <- country_indicator_list[[i]] %>%
    dplyr::select(-ecoregion_id, -country) %>%
    na.omit(.)
  
}

names(country_matrices) <- names(country_indicator_list)

country_correlations <- list()

for (i in seq_along(country_matrices)) {
  
  subset_matrix <- country_matrices[[i]]#[, c(1,2,6,10,20,24,27)]
  
  n <- nrow(subset_matrix)
  
  name <- paste(names(country_matrices)[i], "n =", n, sep = " ")
  
  if (nrow(subset_matrix) < 5) {
    
    print(paste("insufficient data for", names(country_matrices)[i], sep = " "))
    
  } else {
    
    p.mat <- cor_pmat(subset_matrix)
    
    correlation_plot <- ggcorrplot(cor(subset_matrix, method = "spearman"),
                                   title = name, #names(country_matrices)[i],
                                   p.mat = p.mat, 
                                   type = "lower", insig = "blank",
                                   outline.color = "white",
                                   colors = c("#453781FF", "white", "#287D8EFF"),
                                   lab = TRUE)
    
    # ggsave(file.path(current_analysis_outputs, 
    #                  paste(names(country_matrices)[i],
    #                        "indicator_correlation_matrix.png", sep = "_")),
    #        correlation_plot, device = "png")
    
    country_correlations[[i]] <- correlation_plot
    
  }
}

names(country_correlations) <- names(country_matrices)

country_correlations <- list.clean(country_correlations)

country_correlations[[6]]

# Get the coefficients

country_coefficients <- list()

for (i in seq_along(country_correlations)){
  
  country_matrix <- country_correlations[[i]][[1]]
  
  country <- names(country_correlations)[i]
  
  hfp_rli_coefficient <- country_matrix %>%
                         filter(Var1 == "RLI 2005") %>%
                         filter(Var2 == "HFP 2005-") %>%
                         select(value)
  
  ecoregion_n <- ecoregion_countries %>%
                 filter(CNTRY_NAME == country) 
  
  ecoregion_n <- length(unique(ecoregion_n$ECO_ID))
  
  country_coefficient <- cbind(country, hfp_rli_coefficient, ecoregion_n)
               
  
  country_coefficients[[i]] <- country_coefficient 
                              
  
}

country_coefficients_hfp_rli <- do.call(rbind, country_coefficients)

head(country_coefficients_hfp_rli)

# * Biome ----

biome_indicators <- add_grouping_variable("Biome", indicators_wide_centred_trunc,
                                            ecoregions_wide_countries )

biome_indicator_list <- split(biome_indicators, biome_indicators$Biome)

# biome_indicator_list[["All countries"]] <- indicators_wide_centred_trunc %>%
#   mutate(country = "All countrie")

biome_matrices <- list()

for (i in seq_along(biome_indicator_list)) {
  
  biome_matrices[[i]] <- biome_indicator_list[[i]] %>%
    dplyr::select(-ecoregion_id, -Biome) %>%
    na.omit(.)
  
}

names(biome_matrices) <- names(biome_indicator_list)

biome_correlations <- list()

for (i in seq_along(biome_matrices)) {
  
  subset_matrix <- biome_matrices[[i]]#[, c(1,2,6,10,20,24,27)]
  
  n <- nrow(subset_matrix)
  
  name <- paste(names(biome_matrices)[i], "n =", n, sep = " ")
  
  if (nrow(subset_matrix) < 5) {
    
    print(paste("insufficient data for", names(biome_matrices)[i], sep = " "))
    
  } else {
    
    p.mat <- cor_pmat(subset_matrix)
    
    correlation_plot <- ggcorrplot(cor(subset_matrix, method = "spearman"),
                                   title = name, #names(biome_matrices)[i],
                                   p.mat = p.mat, 
                                   type = "lower", insig = "blank",
                                   outline.color = "white",
                                   colors = c("#453781FF", "white", "#287D8EFF"),
                                   lab = TRUE)
    
    # ggsave(file.path(current_analysis_outputs, 
    #                  paste(names(biome_matrices)[i],
    #                        "indicator_correlation_matrix.png", sep = "_")),
    #        correlation_plot, device = "png")
    
    biome_correlations[[i]] <- correlation_plot
    
  }
}

names(biome_correlations) <- names(biome_matrices)

biome_correlations <- list.clean(biome_correlations)

biome_correlations[[6]]

# Get the coefficients

biome_coefficients <- list()

for (i in seq_along(biome_correlations)){
  
  biome_matrix <- biome_correlations[[i]][[1]]
  
  Biome <- names(biome_correlations)[i]
  
  hfp_rli_coefficient <- biome_matrix %>%
    filter(Var1 == "RLI 2005") %>%
    filter(Var2 == "HFP 2005-") %>%
    select(value)
  
  n <- length(unique(biome_indicator_list[[i]]$ecoregion_id)) 
  
  biome_coefficient <- cbind(Biome, hfp_rli_coefficient, n)
  
  
  biome_coefficients[[i]] <- biome_coefficient 
  
  
}

biome_coefficients_hfp_rli <- do.call(rbind, biome_coefficients)

head(biome_coefficients_hfp_rli)

# * Realm ----

names(raw_ecoregions_wide)

realm_indicators <- add_grouping_variable("realm", indicators_wide_centred_trunc,
                                          ecoregions_wide)

realm_indicator_list <- split(realm_indicators, realm_indicators$realm)

realm_indicator_list[["Global"]] <- indicators_wide_centred %>%
                                    mutate(realm = "Global")

realm_matrices <- list()

for (i in seq_along(realm_indicator_list)) {
  
  realm_matrices[[i]] <- realm_indicator_list[[i]] %>%
    dplyr::select(-ecoregion_id, -realm) %>%
    na.omit(.)
  
}

names(realm_matrices) <- names(realm_indicator_list)

# Get pictoral correlation matrices

realm_correlations <- list()

for (i in seq_along(realm_matrices)) {
  
  name <- paste(names(realm_matrices)[i], "indicator correlation matrix", sep = " ")
  
  subset_matrix <- realm_matrices[[i]]#[, c(1,2,6,10,20,24,27)]
  
  if (nrow(subset_matrix) == 0) {
    
    print(paste("no data for", names(realm_matrices)[i], sep = " "))
    
  } else {
    
    correlation_matrix <- ggcorr(subset_matrix, 
                                 method = c("pairwise.complete.obs","spearman"),
                                 label = FALSE, nbreaks = 9) + 
      ggtitle(name)
    
    ggsave(file.path(current_analysis_outputs, 
                     paste(name, "indicator_correlation_matrix.png", sep = "_")),
           correlation_matrix, device = "png")
    
    realm_correlations[[i]] <- correlation_matrix
  }
}

realm_correlations[[8]]

realm_correlations_v2 <- list()

for (i in seq_along(realm_matrices)) {
  
  subset_matrix <- realm_matrices[[i]]#[, c(1,2,6,10,20,24,27)]
  
  n <- nrow(subset_matrix)
  
  name <- paste(names(realm_matrices)[i], "n =", n, sep = " ")
  
  if (nrow(subset_matrix) < 10) {
    
    print(paste("insufficient data for", names(realm_matrices)[i], sep = " "))
    
  } else {
    
    p.mat <- cor_pmat(subset_matrix)
    
    correlation_plot <- ggcorrplot(cor(subset_matrix, method = "spearman"),
                                   title = name, #names(realm_matrices)[i],
                                   p.mat = p.mat, 
                                   type = "lower", insig = "blank",
                                   outline.color = "white",
                                   colors = c("#453781FF", "white", "#287D8EFF"),
                                   lab = TRUE)
    
    ggsave(file.path(current_analysis_outputs, 
                     paste(names(realm_matrices)[i],
                           "indicator_correlation_matrix.png", sep = "_")),
           correlation_plot, device = "png")
    
    realm_correlations_v2[[i]] <- correlation_plot
    
  }
}

realm_correlations_v2[[7]]

# * Headline threat ----


headline_threat_indicators <- add_grouping_variable("headline.threat.type")

headline_threat_indicator_list <- split(headline_threat_indicators, 
                                        headline_threat_indicators$headline.threat.type)

headline_threat_indicator_list[["All threats"]] <- indicators_wide_centred %>%
  mutate(headline.threat.type = "All threats")

headline_threat_indicator_list <- headline_threat_indicator_list[headline_threats]
headline_threat_indicator_list <- list.clean(headline_threat_indicator_list)



headline_threat_matrices <- list()

for (i in seq_along(headline_threat_indicator_list)) {
  
  headline_threat_matrices[[i]] <- headline_threat_indicator_list[[i]] %>%
    dplyr::select(-ecoregion_id, -headline.threat.type) %>%
    na.omit(.)
  
}

names(headline_threat_matrices) <- names(headline_threat_indicator_list)

# Get pictoral correlation matrices

headline_threat_correlations <- list()

for (i in seq_along(headline_threat_matrices)) {
  
  subset_matrix <- headline_threat_matrices[[i]]#[, c(1,2,6,10,20,24,27)]
  
  n <- nrow(subset_matrix)
  
  name <- paste(names(headline_threat_matrices)[i], "indicator correlation matrix,", "n =", n, sep = " ")
  
  if (nrow(subset_matrix) == 0) {
    
    print(paste("no data for", names(headline_threat_matrices)[i], sep = " "))
    
  } else {
    
    correlation_matrix <- ggcorr(subset_matrix, 
                                 method = c("pairwise.complete.obs","spearman"),
                                 label = FALSE, nbreaks = 9) + 
      ggtitle(name)
    
    ggsave(file.path(current_analysis_outputs, 
                     paste(name, "indicator_correlation_matrix.png", sep = "_")),
           correlation_matrix, device = "png")
    
    headline_threat_correlations[[i]] <- correlation_matrix
  }
}

headline_threat_correlations[[7]]

# * Island status ----

names(raw_ecoregions_wide)

island_status_indicators <- add_grouping_variable("island.status")

island_status_indicator_list <- split(island_status_indicators, 
                                      island_status_indicators$island.status)

island_status_indicator_list[["Global"]] <- indicators_wide_centred %>%
  mutate(island.status = "Global")

island_status_matrices <- list()

for (i in seq_along(island_status_indicator_list)) {
  
  island_status_matrices[[i]] <- island_status_indicator_list[[i]] %>%
    dplyr::select(-ecoregion_id, -island.status) %>%
    na.omit(.)
  
}

names(island_status_matrices) <- names(island_status_indicator_list)

# Get pictoral correlation matrices

island_status_correlations <- list()

for (i in seq_along(island_status_matrices)) {
  
  name <- paste(names(island_status_matrices)[i], "indicator correlation matrix", sep = " ")
  
  subset_matrix <- island_status_matrices[[i]]#[, c(1,2,6,10,20,24,27)]
  
  if (nrow(subset_matrix) == 0) {
    
    print(paste("no data for", names(island_status_matrices)[i], sep = " "))
    
  } else {
    
    correlation_matrix <- ggcorr(subset_matrix, 
                                 method = c("pairwise.complete.obs","spearman"),
                                 label = FALSE, nbreaks = 9) + 
      ggtitle(name)
    
    ggsave(file.path(current_analysis_outputs, 
                     paste(name, "indicator_correlation_matrix.png", sep = "_")),
           correlation_matrix, device = "png")
    
    island_status_correlations[[i]] <- correlation_matrix
  }
}

island_status_correlations[[1]]

# * Scientific capacity ----

#' TODO: Is this actually a proxy for data bias?

names(ecoregions_wide)

scientific_capacity_indicators <- add_grouping_variable("scientific.publications.factor")

scientific_capacity_indicator_list <- split(scientific_capacity_indicators, 
                                      scientific_capacity_indicators$scientific.publications.factor)

scientific_capacity_indicator_list[["All capacities"]] <- indicators_wide_centred %>%
  mutate(scientific.publications.factor = "All capacities")

scientific_capacity_matrices <- list()

for (i in seq_along(scientific_capacity_indicator_list)) {
  
  scientific_capacity_matrices[[i]] <- scientific_capacity_indicator_list[[i]] %>%
    dplyr::select(-ecoregion_id, -scientific.publications.factor) %>%
    na.omit(.)
  
}

names(scientific_capacity_matrices) <- names(scientific_capacity_indicator_list)

# Get pictoral correlation matrices

scientific_capacity_correlations <- list()

for (i in seq_along(scientific_capacity_matrices)) {
  
  subset_matrix <- scientific_capacity_matrices[[i]]#[, c(1,2,6,10,20,24,27)]
  
  n <- nrow(subset_matrix)
  
  name <- paste(names(scientific_capacity_matrices)[i], 
                "indicator correlation matrix,", "n =", n, sep = " ")
  
  
  if (nrow(subset_matrix) == 0) {
    
    print(paste("no data for", names(scientific_capacity_matrices)[i], sep = " "))
    
  } else {
    
    correlation_matrix <- ggcorr(subset_matrix, 
                                 method = c("pairwise.complete.obs","spearman"),
                                 label = FALSE, nbreaks = 9) + 
      ggtitle(name)
    
    ggsave(file.path(current_analysis_outputs, 
                     paste(name, "indicator_correlation_matrix.png", sep = "_")),
           correlation_matrix, device = "png")
    
    scientific_capacity_correlations[[i]] <- correlation_matrix
  }
}

scientific_capacity_correlations[[5]]

# * Number of records ----

names(ecoregions_wide)

lpi_records_indicators <- add_grouping_variable("lpi.records.factor")

lpi_records_indicator_list <- split(lpi_records_indicators, 
                                    lpi_records_indicators$lpi.records.factor)

lpi_records_indicator_list[["All records"]] <- indicators_wide_centred %>%
  mutate(lpi.records.factor = "All records")

lpi_records_matrices <- list()

for (i in seq_along(lpi_records_indicator_list)) {
  
  lpi_records_matrices[[i]] <- lpi_records_indicator_list[[i]] %>%
    dplyr::select(-ecoregion_id, -lpi.records.factor) %>%
    na.omit(.)
  
}

names(lpi_records_matrices) <- names(lpi_records_indicator_list)

# Get pictoral correlation matrices

lpi_records_correlations <- list()

for (i in seq_along(lpi_records_matrices)) {
  
  subset_matrix <- lpi_records_matrices[[i]]#[, c(1,2,6,10,20,24,27)]
  
  n <- nrow(subset_matrix)
  
  name <- paste(names(lpi_records_matrices)[i], 
                "indicator correlation matrix,", "n =", n, sep = " ")
  
  
  if (nrow(subset_matrix) == 0) {
    
    print(paste("no data for", names(lpi_records_matrices)[i], sep = " "))
    
  } else {
    
    correlation_matrix <- ggcorr(subset_matrix, 
                                 method = c("pairwise.complete.obs","spearman"),
                          label = FALSE, nbreaks = 9) + 
                            ggtitle(name)
    
    ggsave(file.path(current_analysis_outputs, 
                     paste(name, "indicator_correlation_matrix.png", sep = "_")),
           correlation_matrix, device = "png")
    
    lpi_records_correlations[[i]] <- correlation_matrix
  }
}

lpi_records_correlations[[6]]

# * Test scenarios ----

scenario_indicators <- add_grouping_variable("scenario")

scenario_indicator_list <- split(scenario_indicators, scenario_indicators$scenario)

scenario_indicator_list[["Global"]] <- indicators_wide_centred %>%
  mutate(scenario = "Global")

scenario_matrices <- list()

for (i in seq_along(scenario_indicator_list)) {
  
  scenario_matrices[[i]] <- scenario_indicator_list[[i]] %>%
    dplyr::select(-ecoregion_id, -scenario) %>%
    na.omit(.)
  
}

names(scenario_matrices) <- names(scenario_indicator_list)

# Get pictoral correlation matrices

scenario_correlations <- list()

for (i in seq_along(scenario_matrices)) {
  
  subset_matrix <- scenario_matrices[[i]]#[, c(1,2,6,10,20,24,27)]
  
  n <- nrow(subset_matrix)
  
  name <- paste(names(scenario_matrices)[i], 
                "indicator correlation matrix,", "n =", n, sep = " ")
  
  if (nrow(subset_matrix) < 10) {
    
    print(paste("no data for", names(scenario_matrices)[i], sep = " "))
    
  } else {
    
    correlation_matrix <- ggcorr(subset_matrix, 
                                 method = c("pairwise.complete.obs","spearman"),
                                 label = FALSE, nbreaks = 9) + 
      ggtitle(name)
    
    # ggsave(file.path(current_analysis_outputs, 
    #                  paste(name, "indicator_correlation_matrix.png", sep = "_")),
    #        correlation_matrix, device = "png")
    
    scenario_correlations[[i]] <- correlation_matrix
  }
}

# The ggcorr and ggpairs outputs don't appear to always match?
n <- 5
test <- list.clean(scenario_correlations)
test[[n]]

x <- scenario_matrices[[n]]
x

y <- ggpairs(x,lower = list(continuous = wrap("points", size=0.1)),
        upper = list(continuous = wrap("cor", method = "spearman",size= 2)))
y

n <- 8

plotnames <- paste("Scenario", c(1:7), sep = "_")

scenario_correlations_v2 <- list()

for (i in seq_along(scenario_matrices)) {
  
  subset_matrix <- scenario_matrices[[i]]#[, c(1,2,6,10,20,24,27)]
  
  n <- nrow(subset_matrix)
  
  name <- paste(names(scenario_matrices)[i], "n =", n, sep = " ")
  
  if (nrow(subset_matrix) < 10) {
    
    print(paste("insufficient data for", names(scenario_matrices)[i], sep = " "))
    
  } else {
  
p.mat <- cor_pmat(subset_matrix)

correlation_plot <- ggcorrplot(cor(subset_matrix, method = "spearman"),
           title = name, #names(scenario_matrices)[i],
           p.mat = p.mat, 
           type = "lower", insig = "blank",
           outline.color = "white",
           colors = c("#453781FF", "white", "#287D8EFF"),
           lab = TRUE)

ggsave(file.path(current_analysis_outputs, 
       paste(plotnames[[i]], "indicator_correlation_matrix.png", sep = "_")),
       correlation_plot, device = "png")

scenario_correlations_v2[[i]] <- correlation_plot

  }
}

scenario_correlations_v2[[5]]

# Could use this instead, allows you to mark significant correlations

# https://rpkgs.datanovia.com/ggcorrplot/reference/ggcorrplot.html


# * Predominant threat included in HFP ----

HFP_threat_indicators <- add_grouping_variable("included.in.HFP")

HFP_threat_indicator_list <- split(HFP_threat_indicators, HFP_threat_indicators$included.in.HFP)

HFP_threat_indicator_list[["Global"]] <- indicators_wide_centred %>%
  mutate(included.in.HFP = "Global")

HFP_threat_matrices <- list()

for (i in seq_along(HFP_threat_indicator_list)) {
  
  HFP_threat_matrices[[i]] <- HFP_threat_indicator_list[[i]] %>%
    dplyr::select(-ecoregion_id, -included.in.HFP) %>%
    na.omit(.)
  
}

names(HFP_threat_matrices) <- names(HFP_threat_indicator_list)

# Get pictoral correlation matrices

HFP_threat_correlations <- list()

for (i in seq_along(HFP_threat_matrices)) {
  
  subset_matrix <- HFP_threat_matrices[[i]]#[, c(1,2,6,10,20,24,27)]
  
  n <- nrow(subset_matrix)
  
  name <- paste(names(HFP_threat_matrices)[i], 
                "indicator correlation matrix,", "n =", n, sep = " ")
  
  if (nrow(subset_matrix) < 10) {
    
    print(paste("no data for", names(HFP_threat_matrices)[i], sep = " "))
    
  } else {
    
    correlation_matrix <- ggcorr(subset_matrix, 
                                 method = c("pairwise.complete.obs","spearman"),
                                 label = FALSE, nbreaks = 9) + 
      ggtitle(name)
    
    # ggsave(file.path(current_analysis_outputs, 
    #                  paste(name, "indicator_correlation_matrix.png", sep = "_")),
    #        correlation_matrix, device = "png")
    
    HFP_threat_correlations[[i]] <- correlation_matrix
  }
}

test <- list.clean(HFP_threat_correlations)
test[[1]]


# Map scenarios ----

# Global ----

ecoregion_values_scenario <- ecoregions_wide %>%
                             dplyr::select(ecoregion_id, scenario)

ecoregion_values_map_data <- left_join(ecoregion_map_renamed, 
                                       ecoregion_values_scenario,
                                       by = "ecoregion_id")

global_scenario_map <-  ggplot(ecoregion_values_map_data) +
  geom_sf(aes(fill = scenario), colour = "black", 
          size = 0.05, show.legend = 'fill') +
  scale_fill_viridis_d(alpha = .8,
                       na.value = "grey70") +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  labs(fill = "test") +
  theme(legend.position = "right")

global_scenario_map

ggsave(file.path(current_analysis_outputs,
                 paste(location, eco_version,
                       "data_threat_scenarios.png",
                                          sep = "_")), 
       scenario_map,  device = "png")

# Realms ----

ecoregion_values_scenario_oceania <- ecoregions_wide %>%
  dplyr::select(ecoregion_id, scenario, realm) 

ecoregion_values_map_data <- left_join(ecoregion_map_renamed, 
                                       ecoregion_values_scenario_oceania,
                                       by = "ecoregion_id") %>%
  filter(realm == "Oceania")

scenario_map <-  ggplot(ecoregion_values_map_data) +
  geom_sf(aes(fill = scenario), colour = "black", 
          size = 0.05, show.legend = 'fill') +
  scale_fill_viridis_d(alpha = .8,
                       na.value = "grey70") +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  labs(fill = "test") +
  theme(legend.position = "right") +
  scale_x_continuous(limits = c(110, 300)) + 
  scale_y_continuous(limits = c(-50, 70)) 

scenario_map

ggsave(file.path(current_analysis_outputs,
                 paste(location, eco_version,
                       "data_threat_scenarios.png",
                       sep = "_")), 
       scenario_map,  device = "png")

# Correlation matrices ----

spearman_correlations_all_years <- as.data.frame(cor(indicator_matrix, 
                                                     method = "spearman"))

# PCA ----

# Prepare data

pca_data_1 <- gather(indicators_wide_centred_trunc, indicator, indicator_value, 
                     "BHI_plants 2005":"threatened_2005", factor_key=TRUE)

pca_data_2 <- pca_data_1 %>%
              merge(ecoregions_wide,by = "ecoregion_id") %>%
              dplyr::select(ecoregion_id, realm, headline.threat.type,
                            Biome, island.status, scientific.publications.factor,
                            area.factor,
                            lpi.records.factor,
                            rli.records.factor,
                            scenario,
                            included.in.HFP,
                            predominant.threat.type,
                            scenario.numeric,
                            predominant.threat.count,
                            everything()) 

pca_data_3 <- spread(pca_data_2,indicator,indicator_value) 

pca_data_4 <- pca_data_3[complete.cases(pca_data_3[,19:ncol(pca_data_3)]),]

pca_data_5 <- pca_data_4 %>%
              dplyr::select(-BIIab_2005, -threatened_2005) #[pca_data_4$headline.threat.type %in% headline_threats,]

# Conduct the PCA

pca=prcomp(pca_data_5[,19:ncol(pca_data_5)],center=TRUE,scale=TRUE)
summary(pca)
print(pca)

var <- get_pca_var(pca)
var$contrib


# look at the eigen values and cumulative variance plot

screeplot(pca, type = "l", npcs = 5, main = "Screeplot of the first 10 PCs")
abline(h = 1, col="red", lty=5)
legend("topright", legend=c("Eigenvalue = 1"),
       col=c("red"), lty=5, cex=0.6)
cumpro <- cumsum(pca$sdev^2 / sum(pca$sdev^2))
plot(cumpro[0:5], xlab = "PC #", ylab = "Amount of explained variance", 
     main = "Cumulative variance plot")

# first two principle components plot (the results seems clustered..)
plot(pca$x[,1],pca$x[,2], xlab="PC1 (44.3%)", ylab = "PC2 (19%)", main = "PC1 / PC2 - plot")

# plot in 2 dimentions (seems not good results...)
# install.packages("factoextra")

# * Plot PCA ----

# Make the PCA plots in a loop

# Remove indicators

grouping_variables <- names(pca_data_5)[!(names(pca_data_5) %in% indicators)]

# Remove numerics
grouping_variables <- grouping_variables[!(grouping_variables %in% c("ecoregion_id",
                                                                     "predominant.threat.type",
                                                                     "headline.threat.type",
                                                                     "ecoregion.area.km.sq",
                                                                     "RLI_records",
                                                                     "LPI_records",
                                                                     "scenario",
                                                                     "predominant.threat.count",
                                                                     "mean.scientific.publications"))]

pca_plots <- list()

for (i in seq_along(grouping_variables)) {
  
  vargroup <- pca_data_5[,grouping_variables[i]]
  varname <- grouping_variables[i]
  
  pca_plot <- fviz_pca_biplot(pca, geom.ind = "point", pointshape = 21, 
                       pointsize = 2, 
                       col.ind = "black", 
                       palette = "jco", 
                       addEllipses = TRUE,
                       label = "var",
                       col.var = "black",
                       repel = TRUE,
                       fill.ind = as.factor(vargroup),
                       legend.title = varname) +
    ggtitle("2D PCA-plot from 30 feature dataset") +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(file.path(current_analysis_outputs, paste(varname, ".png", sep = "")),
         pca_plot, device = "png")
  
  pca_plots[[i]] <- pca_plot
  
}

pca_plots[[1]]
pca_plots[[2]]
pca_plots[[3]]
pca_plots[[4]]
pca_plots[[5]]
pca_plots[[6]]
pca_plots[[7]]
pca_plots[[8]]
pca_plots[[9]]
pca_plots[[10]]


# Try a model ----

# TODO: use tidymodels ----

library(Rmisc)
library(car)
library(lm.beta)



hfp_rli <- scatterplot(`RLI 2005` ~ `HFP 2005-` , data = model_input_data_1)

hfp_rli_lm <- lm(`RLI 2005` ~ `HFP 2005-`+ RLI_records + included.in.HFP, 
                 data = model_input_data_1)

plot(hfp_rli_lm)

influence.measures(hfp_rli_lm)

summary(hfp_rli_lm)
avPlots(hfp_rli_lm, ask = F)


predictions <- predict(hfp_rli_lm, data = analysis_input_data$included.in.HFP)

relationship <- cut(predictions, breaks = 3, 
                    labels = c("negative", "neutral", "positive"))

predictions <- data.frame(ecoregion_id = analysis_input_data$ecoregion_id, 
                          prediction = predictions,
                          threat = analysis_input_data$included.in.HFP)

prediction_map <- left_join(ecoregion_map_renamed, 
                            predictions,
                            by = "ecoregion_id") 


relationship_map_LU <-  ggplot(prediction_map) +
                     geom_sf(aes(fill = prediction), colour = "black", 
                              size = 0.05, show.legend = 'fill') +
                     scale_fill_viridis_c(alpha = .8,
                                           na.value = "grey70") +
                     theme(axis.line = element_line(),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            panel.background = element_blank()) +
                     labs(fill = "test") +
                     theme(legend.position = "right")

relationship_map_LU

# Old analysis code ----

#' # All indicators with 2005 data
#' 
#' indicator_matrix_2005 <- indicator_matrix[, c(1,2,10,16,23,27,28)]
#' 
#' indicator_correlations_2005 <- ggcorr(indicator_matrix_2005)
#' 
#' # All indicators with 2010 data
#' 
#' indicator_matrix_2010 <- indicator_matrix[, c(3,6,11,17,25)]
#' 
#' indicator_correlations_2010 <- ggcorr(indicator_matrix_2010)
#' 
#' # RLI all timepoints and 1990 HFP
#' 
#' indicator_matrix_rli <- indicator_matrix[, c(5,21,22,23,24,25)]
#' 
#' indicator_correlations_rli <- ggcorr(indicator_matrix_rli)
#' 
#' # BHI all timepoints and 1990 HFP
#' 
#' indicator_matrix_bhi <- indicator_matrix[, c(5,2, 3, 4)]
#' 
#' indicator_correlations_bhi <- ggcorr(indicator_matrix_bhi)
#' 
#' # Extinctions all timepoints and 1990 HFP
#' 
#' indicator_matrix_ex <- indicator_matrix[, c(5,13,14,15,16,17,18)]
#' 
#' indicator_correlations_ex <- ggcorr(indicator_matrix_ex)
#' 
#' # Indicators from different data sources 2005 - 2010
#' 
#' indicator_matrix_ind <- indicator_matrix[, c(1,2,6,10,20,24,27)]
#' 
#' indicator_correlations_ind <- ggcorr(indicator_matrix_ind, palette = "BuPu",
#'                                      label = TRUE)
#' 
#' 
#' # RLI Birds 1980 and following extinctions
#' 
#' indicator_matrix_er <- indicator_matrix[, c(21, 13, 14,15,16,17,18)]
#' 
#' indicator_correlations_er <- ggcorr(indicator_matrix_er)
#' 
#' 
#' # Transformed
#' 
#' indicator_matrix_t <- indicators_transformed %>%
#'   select(-ecoregion_id) %>%
#'   na.omit(.)
#' 
#' pearson_correlations_all_years_t <- as.data.frame(cor(indicator_matrix_t, 
#'                                                       method = "pearson"))
#' 
#' spearman_correlations_all_years_t <- as.data.frame(cor(indicator_matrix_t, 
#'                                                        method = "spearman"))
#' 
#' # mutate(HFP_adjusted_old = ifelse(raw_indicator_value > 33.29037,
#' #                                  33.29037, raw_indicator_value)) %>%
#' # mutate(RLI_adjusted_old = ifelse(raw_indicator_value == 0, NA,
#' #                ifelse(raw_indicator_value > 0 & 
#' #                         raw_indicator_value < 0.9538,
#' #                       0.9538, raw_indicator_value))) %>%
#' # mutate(year = as.numeric(year)) %>%
#' # mutate(HFP_adjusted_old = ifelse(indicator != 
#' #                           "mean human footprint index", NA,
#' #                           HFP_adjusted_old)) %>%
#' # mutate(RLI_adjusted_old = ifelse(grepl("red list index Aves",
#' #                                 indicator), 
#' #                                 RLI_adjusted_old, NA)) %>%
#' # mutate(HFP_adjusted =
#' #          pmin(pmax(raw_indicator_value,quantile(raw_indicator_value,
#' #                                          .005, na.rm = TRUE)),
#' #               quantile(raw_indicator_value, .995,
#' #                        na.rm = TRUE))) %>%
#' # mutate(HFP_scaled_adjusted = scale_to_1(HFP_adjusted)) %>%
#' # mutate(HFP_scaled_adjusted_inverted =
#' #          1 - HFP_scaled_adjusted) %>%
#' # mutate(scaled = scale_to_1(raw_indicator_value)) %>%
#' # mutate(inverted = 1 - raw_indicator_value) %>%
#' # mutate(scaled_inverted = 1 - scaled) %>%
#' # mutate(RLI_adjusted = ifelse(raw_indicator_value == 0, NA,
#' #                              pmin(pmax(raw_indicator_value,
#' #                   quantile(raw_indicator_value, .05, 
#' #                            na.rm = TRUE))))) %>%
#' # mutate(RLI_adjusted_inverted = 1 - RLI_adjusted)
#' 
#' # Scatterplots ----
#' 
#' indicator_values_master_1 <- indicator_values_master
#' 
#' indicator_values_master_2 <- indicator_values_master_1 %>%
#'   mutate(key2 = paste(ecoregion_id, indicator_year, sep = " ")) %>%
#'   mutate(key = make.names(key2)) %>%
#'   mutate(key = str_remove(key, "X")) %>%
#'   select(-key2) 
#' 
#' indicators_scaled <- indicators_scaled %>%
#'   mutate(variable = str_remove(variable, "-"))
#' 
#' indicator_values_scaled <- indicators_scaled %>%
#'   mutate(key = paste(ecoregion_id, variable, sep = ".")) %>%
#'   rename(centred_indicator_value = value)
#' 
#' unique(indicator_values_scaled$variable)
#' unique(indicator_values_master$indicator_year)
#' 
#' indicator_values_master <- indicator_values_master_2 %>%
#'   merge(indicator_values_scaled[c("key", 
#'                                   "centred_indicator_value")], by = "key") 
#' 
#' unique(indicator_values_master$indicator_year)
#' 
#' # * Subset by indicator ----
#' 
#' indicator_values_master_list <- split(indicator_values_master, 
#'                                       indicator_values_master$indicator)
#' 
#' indicator_values_master_wide <- list()
#' indicator_by_time_plots <- list()
#' 
#' for (i in seq_along(indicator_values_master_list)) {
#'   
#'   wide_indicator_data <- indicator_values_master_list[[i]] %>%
#'     select(ecoregion_id, indicator_year, 
#'            centred_indicator_value) %>%
#'     distinct(.) %>%
#'     spread(key = indicator_year, 
#'            value = centred_indicator_value) %>%
#'     select(-ecoregion_id)
#'   
#'   indicator_plot <- ggpairs(wide_indicator_data)
#'   
#'   indicator_name <- indicator_values_master_list[[i]][1,2]
#'   
#'   indicator_values_master_wide[[i]] <- wide_indicator_data
#'   indicator_by_time_plots[[i]] <- indicator_plot
#'   
#'   scatterplot_directory <- file.path(indicator_outputs, 
#'                                      "scatterplots")
#'   
#'   if( !dir.exists( scatterplot_directory ) ) {
#'     
#'     dir.create( scatterplot_directory, recursive = TRUE )
#'     
#'   }
#'   
#'   ggsave(file.path(scatterplot_directory, paste(location, eco_version,
#'                                                 indicator_name, 
#'                                                 "time_scatterplots.png", 
#'                                                 sep = "_")),
#'          indicator_plot,  device = "png")
#'   
#' }
#' 
#' 
#' # * Subset by time ----
#' 
#' 
#' # Non-transformed
#' 
#' indicators_by_year <- split(indicator_values_master, indicator_values_master$year)
#' indicators_by_year_names <- paste(location, names(indicators_by_year), sep = "_")
#' 
#' scatterplots_by_year <- list()
#' 
#' for (i in seq_along(indicators_by_year)){
#'   
#'   scatterplots_by_year[[i]] <- produce_scatterplots(indicators_by_year[[i]],
#'                                                     indicators_by_year_names[[i]],
#'                                                     save = TRUE, "centred")
#'   
#' }
#' 
#' scatterplots_by_year[[4]]
#' 
#' 
#' # Centred and inverted 
#' 
#' indicators_wic <- indicators_wic %>%
#'   rename(ecoregion_id = Ecoregion_id) %>%
#'   select(-geometry, -REALM)
#' 
#' indicators_wic_long <- melt(indicators_wic, 
#'                             id.vars = 'ecoregion_id')
#' 
#' 
#' indicator_values_wic <- indicators_wic_long  %>%
#'   rename(indicator_year = variable,
#'          raw_indicator_value = value) %>% # TODO: Change to transformed
#'   mutate(year = str_sub(indicator_year, start= -4)) %>%
#'   mutate(indicator_year = as.character(indicator_year)) %>%
#'   mutate(indicator = removeNumbers(indicator_year)) %>%
#'   mutate(indicator = str_replace_all(indicator,
#'                                      '[[:punct:]]',' '))
#' 
#' indicators_by_year <- split(indicator_values_wic, indicator_values_wic$year)
#' indicators_by_year_names <- paste(location, names(indicators_by_year), sep = "_")
#' 
#' scatterplots_by_year <- list()
#' 
#' for (i in seq_along(indicators_by_year)){
#'   
#'   scatterplots_by_year[[i]] <- produce_scatterplots(indicators_by_year[[i]],
#'                                                     indicators_by_year_names[[i]],
#'                                                     save = FALSE)
#'   
#' }
#' 
#' scatterplots_by_year[[1]]
#' 
#' 
#' # RLI and HFP 2005 only ---
#' 
#' indicator_values_rli_hfp <- indicator_values_wic %>%
#'   filter(indicator_year == "HFP.2005-" |
#'            indicator_year ==  "RLI.2005") 
#' 
#' rli_hfp_corr <- produce_scatterplots(indicator_values_rli_hfp, "test", save = FALSE,"raw")
#' 
#' indicators_by_year <- split(indicator_values_wic, indicator_values_wic$year)
#' indicators_by_year_names <- paste(location, names(indicators_by_year), sep = "_")
#' 
#' scatterplots_by_year <- list()
#' 
#' for (i in seq_along(indicators_by_year)){
#'   
#'   scatterplots_by_year[[i]] <- produce_scatterplots(indicators_by_year[[i]],
#'                                                     indicators_by_year_names[[i]],
#'                                                     save = FALSE)
#'   
#' }
#' 
#' scatterplots_by_year[[1]]
#' 
#' # Oceania ----
#' 
#' indicators_wic <- indicators_wic_oceania %>%
#'   rename(ecoregion_id = Ecoregion_id) %>%
#'   select(-geometry, -REALM)
#' 
#' indicators_wic_long <- melt(indicators_wic, 
#'                             id.vars = 'ecoregion_id')
#' 
#' 
#' indicator_values_wic <- indicators_wic_long  %>%
#'   rename(indicator_year = variable,
#'          raw_indicator_value = value) %>% # TODO: Change to transformed
#'   mutate(year = NA) %>%
#'   mutate(indicator = as.character(indicator_year)) 
#' 
#' oceania_comparisons <- indicator_values_wic %>%
#'   filter(indicator == "BII"| 
#'            indicator == "BHI"|
#'            indicator == "HFP"|
#'            indicator == "Endangered"|
#'            indicator == "RLI_Amph"|
#'            indicator == "RLI_Birds"|
#'            indicator == "RLI_Mammals")
#' 
#' length(unique(oceania_comparisons$ecoregion_id))
#' 
#' oceania_bhi_amph <- oceania_comparisons %>%
#'   filter(indicator == "RLI_Amph"| 
#'            indicator == "BHI")
#' 
#' produce_scatterplots(oceania_bhi_amph, "oceania_bhi_amph", save = TRUE)
#' 
#' oceania_scatterplots <- produce_scatterplots(oceania_comparisons, "oceania", save =FALSE)
#' 
#' plot(indicators_wic_oceania$HFP,indicators_wic_oceania$RLI_Amph)
#' 
#' oceania_master <- indicator_values_master %>%
#'   filter(realm == "Oceania")
#' 
#' oceania_master_bhi <- oceania_master %>%
#'   filter(indicator == "BHI plants") %>%
#'   select(-country) %>%
#'   distinct(.) %>%
#'   group_by(ecoregion_id, year) %>%
#'   arrange(.)
#' 
#' ggplot(oceania_master_bhi) +
#'   geom_line(aes(x = year, 
#'                 y = raw_indicator_value, 
#'                 color = ecoregion_id))
#' 
#' 
#' # Australasia ----
#' 
#' indicators_wic <- indicators_wic_australasia %>%
#'   rename(ecoregion_id = Ecoregion_id) %>%
#'   select(-geometry, -REALM)
#' 
#' indicators_wic_long <- melt(indicators_wic, 
#'                             id.vars = 'ecoregion_id')
#' 
#' 
#' indicator_values_wic <- indicators_wic_long  %>%
#'   rename(indicator_year = variable,
#'          raw_indicator_value = value) %>% # TODO: Change to transformed
#'   mutate(year = NA) %>%
#'   mutate(indicator = as.character(indicator_year)) 
#' 
#' australasia_comparisons <- indicator_values_wic %>%
#'   filter(indicator == "BII"| 
#'            indicator == "BHI"|
#'            indicator == "HFP"|
#'            indicator == "Endangered"|
#'            indicator == "RLI_Amph"|
#'            indicator == "RLI_Birds"|
#'            indicator == "RLI_Mammals")
#' 
#' length(unique(oceania_comparisons$ecoregion_id))
#' 
#' aus_bhi_amph <- australasia_comparisons %>%
#'   filter(indicator == "RLI_Amph"| 
#'            indicator == "BHI")
#' 
#' aus_bhi_amph$RLI_Amph_t <- 1/(max(aus_bhi_amph$RLI_Amph_t + 1, na.rm = TRUE) - 
#'                                 aus_bhi_amph$RLI_Amph_t)
#' 
#' 1/(max(indicators_wide_inv_centred$`proportion.at.risk.1980-`+1, 
#'        na.rm = TRUE) - 
#'      indicators_wide_inv_centred$`proportion.at.risk.1980-`) 
#' 
#' produce_scatterplots(aus_bhi_amph, "aus_bhi_amph", save = TRUE)
#' 
#' oceania_scatterplots <- produce_scatterplots(oceania_comparisons, "oceania", save =FALSE)
#' 
#' plot(indicators_wic_oceania$HFP,indicators_wic_oceania$RLI_Amph)
#' 
#' oceania_master <- indicator_values_master %>%
#'   filter(realm == "Oceania")
#' 
#' oceania_master_bhi <- oceania_master %>%
#'   filter(indicator == "BHI plants") %>%
#'   select(-country) %>%
#'   distinct(.) %>%
#'   group_by(ecoregion_id, year) %>%
#'   arrange(.)
#' 
#' ggplot(oceania_master_bhi) +
#'   geom_line(aes(x = year, 
#'                 y = raw_indicator_value, 
#'                 color = ecoregion_id))
#' 
#' # Transformed
#' 
#' indicators_by_year <- split(indicators_values_transformed, indicator_values_transformed$year)
#' indicators_by_year_names <- paste(location, "transformed", names(indicators_by_year), sep = "_")
#' 
#' scatterplots_by_year <- list()
#' 
#' for (i in seq_along(indicators_by_year)) {
#'   
#'   scatterplots_by_year[[i]] <- produce_scatterplots(indicators_by_year[[i]],
#'                                                     indicators_by_year_names[[i]],
#'                                                     save = FALSE)
#'   
#' }
#' 
#' scatterplots_by_year[[2]]
#' 
#' # * Subset by realm ----
#' 
#' # * Subset by time by realm ----
#' 
#' # Split data into nested list (level 1 = year, level 2 = realm)
#' 
#' indicators_by_year_by_realm <- list()
#' 
#' for (i in seq_along(indicators_by_year)) {
#'   
#'   indicators_timepoint <- indicators_by_year[[i]]
#'   
#'   indicators_by_year_by_realm[[i]] <- split(indicators_timepoint, 
#'                                             indicators_timepoint$realm)
#' }
#' 
#' indicators_by_year_by_realm <- flatten(indicators_by_year_by_realm)
#' realm_names <- names(indicators_by_year_by_realm)
#' timepoints <- names(indicators_by_year)
#' time_names <- rep(timepoints, each = length(unique(realm_names)))
#' indicators_by_year_by_realm_names <- paste(realm_names, time_names , sep = "_")
#' 
#' scatterplots_by_year_by_realm <- list()
#' 
#' for (i in seq_along(indicators_by_year_by_realm)){
#'   
#'   scatterplots_by_year_by_realm[[i]] <- produce_scatterplots(indicators_by_year_by_realm[[i]],
#'                                                              indicators_by_year_by_realm_names[[i]],
#'                                                              save = TRUE)
#'   
#' }
#' 
#' names(scatterplots_by_year_by_realm) <- indicators_by_year_by_realm_names
#' 
#' # Look at Australiasia 2005
#' 
#' scatterplots_by_year_by_realm[[27]]
#' 
#' # Look at Oceania 2005
#' 
#' scatterplots_by_year_by_realm[[31]]
#' 
#' # * Subset effect of land use on extinction risk ----
#' 
#' land_use_atrisk_values <- indicator_values_master %>%
#'   filter(indicator_year == "human footprint index 1990"|
#'            indicator_year == "proportion at risk 1990"|
#'            indicator_year == "proportion at risk 2000"|
#'            indicator_year == "proportion at risk 2005"|
#'            indicator_year == "proportion at risk 2010"|
#'            indicator_year == "proportion at risk 2015")
#' 
#' hfp_vs_atrisk_scatterplots <- produce_scatterplots(land_use_atrisk_values, 
#'                                                    paste(location, eco_version,
#'                                                          "HFP1990_vs_AtRisk_time_series",
#'                                                          sep = "_"), save = TRUE)
#' 
#' # * Subset effect of land use on extinction  ----
#' 
#' land_use_extinction_values <- indicator_values_master %>%
#'   filter(indicator_year == "human footprint index 1990"|
#'            indicator_year == "proportion extinct 1990"|
#'            indicator_year == "proportion extinct 2000"|
#'            indicator_year == "proportion extinct 2005"|
#'            indicator_year == "proportion extinct 2010"|
#'            indicator_year == "proportion extinct 2015")
#' 
#' hfp_vs_extinction_scatterplots <- produce_scatterplots(land_use_extinction_values, 
#'                                                        paste(location, eco_version,
#'                                                              "HFP1990_vs_Extinction_time_series",
#'                                                              sep = "_"), save = TRUE)
#' 
#' 
#' hfp_vs_extinction_scatterplots
#' 
#' # * Subset only unrelated indicators  ----
#' 
#' subset_values <- indicator_values_master %>%
#'   filter(indicator_year == "human footprint index 1990"|
#'            indicator_year == "human footprint index 2010"|
#'            indicator_year == "abundance biodiversity intactness index 2005"|
#'            indicator_year == "richness biodiversity intactness index 2005"|
#'            indicator_year == "proportion extinct 2005")
#' 
#' subset_scatterplots <- produce_scatterplots(subset_values, 
#'                                             paste(location, eco_version,
#'                                                   "subset_indicators",
#'                                                   sep = "_"), save = FALSE)
#' 
#' 
#' subset_scatterplots
#' 
#' 
#' ggplot(data = indicator_values_master_wide_scaled, 
#'        aes(x = abundance.biodiversity.intactness.index.2005, 
#'            y = human.footprint.index.1990)) +
#'   geom_point() + 
#'   coord_trans(x = "log10", y = "log10")
#' 
#' # Pairwise correlations and scatterplots ----
#' 
#' #' TODO: Can we use ggpair with other ggplot syntax? Including attributes but
#' #' also transformations? https://www.r-graph-gallery.com/199-correlation-matrix-with-ggally.html
#' 
#' correlations <- list()
#' scatterplots <- list()
#' years <- list()
#' 
#' for (i in seq_along(indicator_values_time_list)) {
#'   
#'   # Get the year name
#'   
#'   year <- range(indicator_values_time_list[[i]]$year)
#'   
#'   # Remove year column from dataframe
#'   
#'   step1 <- indicator_values_time_list[[i]] %>% 
#'     select(ecoregion_id, 
#'            raw_indicator_value, 
#'            indicator_year, 
#'            -year) %>%
#'     distinct(.)
#'   
#'   # Convert into wide format
#'   
#'   # step2 <- dcast(step1, ecoregion_code ~ indicator, 
#'   #                value.var = "raw_indicator_value", sum)
#'   
#'   step2 <- step1 %>%
#'     spread(key = indicator_year, 
#'            value = raw_indicator_value) 
#'   
#'   # Convert into matrix without id variables
#'   
#'   step3 <- as.matrix(step2[2:ncol(step2)])
#'   
#'   # Check if there's enough indicators to analyse correlation (ie more than 2)
#'   
#'   if(ncol(step3) == 1) {
#'     
#'     print(paste("Year", year, "does not have values for more than one indicator to test correlations between"))
#'     
#'     rm(step1, step2, step3)
#'     
#'   } else {
#'     
#'     # Analyse correlations between indicators 
#'     
#'     step4 <- step3[complete.cases(step3),]
#'     
#'     correlations[[i]] <- as.data.frame(cor(step4, method = "pearson"))
#'     
#'     scatterplots[[i]] <- ggpairs(as.data.frame(step4))
#'     
#'     years[[i]] <- year
#'     
#'     rm(step1, step2, step3, step4)
#'     
#'   }
#'   
#' }
#' 
#' # Remove the empty lists with no data in them
#' 
#' names(correlations) <- years
#' 
#' correlations <- list.clean(correlations, fun = is.null, recursive = FALSE)
#' 
#' names(scatterplots) <- years
#' 
#' if(save_outputs == "yes") {
#'   
#'   for (i in seq_along(scatterplots)) {
#'     
#'     ggsave(file.path(outputs, paste(as.character(years[[i]]),"_scatterplot",
#'                                     ".png", sep = "")), 
#'            scatterplots[[i]],  device = "png")
#'     
#'   }
#' }
