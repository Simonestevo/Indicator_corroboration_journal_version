
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
timepoint <- "2005"
load_map <- FALSE


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

indicator_properties <- read.csv(file.path(analysis_inputs,
                                          "indicator_properties.csv"))

raw_indicators_long_all <- readRDS(file.path(analysis_inputs,
                          "global_ecoregions_2017_indicator_values_master.rds"))

raw_indicators_long_all$indicator_year <- str_replace(raw_indicators_long_all$indicator_year, " ", "_")

# TEMPORARY: Remove LPI ----

raw_indicators_long_all <- raw_indicators_long_all %>%
                           filter(indicator != "LPI")

raw_indicators_long <- raw_indicators_long_all %>%
                       dplyr::select(ecoregion_id, indicator_year,
                                     raw_indicator_value) %>%
                       distinct(.)

raw_indicators_wide <- raw_indicators_long %>%
                       spread(key = indicator_year, 
                             value = raw_indicator_value) 

# indicators <- unique(raw_indicators_long$indicator_year)

# raw_indicators_wide <- readRDS(file.path(analysis_inputs,
#                           "global_ecoregions_2017_indicator_values_master_wide.rds"))

# if (!is.na(timepoint)) {
#   
#   raw_indicators_long <- raw_indicators_long_all %>%
#     filter(year == timepoint) %>%
#     dplyr::select(ecoregion_id, indicator_year,
#                   raw_indicator_value) %>%
#     distinct(.)
#   
#   raw_indicators_wide <- raw_indicators_long %>%
#     spread(key = indicator_year, 
#            value = raw_indicator_value) 
#   
#   indicators <- unique(raw_indicators_long$indicator_year)
#   
# }

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

if(load_map == TRUE) {

ecoregion_map_all <- readRDS(paste(
                     file.path("N:/Quantitative-Ecology/Simone/extinction_test/inputs",
                               "ecoregions_2017"),
                               "Ecoregions2017valid.rds"))

ecoregion_map <- ecoregion_map_all %>%
  dplyr::select(ECO_ID, ECO_NAME, OBJECTID, REALM, geometry)

ecoregion_map_renamed <- ecoregion_map %>% dplyr::rename(ecoregion_id = ECO_ID)

rm(ecoregion_map_all, ecoregion_map)
  
}

# Countries

ecoregion_countries <- readRDS(file.path(parent_outputs, "version_3",
                        "2020-08-10_database_output_files",
                        "ecoregion_country_data.rds"))



# Ecoregion data cleaning ----

# Convert numeric to factors

ecoregions_wide <- raw_ecoregions_wide %>%
                   mutate(LPI_records = as.numeric(ifelse(is.na(LPI_records), 
                                               0, LPI_records)))

ecoregions_wide$scientific.publications.factor <- cut(ecoregions_wide$mean.scientific.publications, breaks = 3, 
                            labels = c("Few_publications", 
                                       "Moderate_publications",
                                       "Many_publications"))

ecoregions_wide$area.factor <- cut(ecoregions_wide$ecoregion.area.km.sq, breaks = 3, 
                                  labels = c("Small_ecoregions", 
                                             "Medium_ecoregions",
                                             "Large_ecoregions"))

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
               scenario.numeric = as.factor(as.numeric(scenario))) %>%
        mutate(scenario = ifelse(scenario == "Few (< 590) RLI records & Threat external to HFP",
                                 "Few RLI HFP exclusive",
                                 ifelse(scenario == "Few (< 590) RLI records & Threat related to HFP",
                                 "Few RLI HFP inclusive",
                                 ifelse(scenario == "Moderate (between 590 & 1200) RLI records & Threat external to HFP",
                                 "Moderate RLI HFP exclusive",
                                 ifelse(scenario == "Moderate (between 590 & 1200) RLI records & Threat related to HFP",
                                 "Moderate RLI HFP inclusive",
                                 ifelse(scenario == "Many (> 1200) RLI records & Threat external to HFP",
                                 "Many RLI HFP exclusive",
                                 ifelse(scenario == "Many (> 1200) RLI records & Threat related to HFP",
                                 "Many RLI HFP inclusive",
                                 ifelse(scenario == "Few (< 590) RLI records & NA",
                                 "Few RLI HFP NA",
                                 NA))))))))


ecoregions_wide$endemics.factor <- cut(ecoregions_wide$number.of.endemics, 
                                       breaks = c(-1,0,100), 
                                   labels = c("no endemics", 
                                              "endemics"))

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

#PCA input data ----

indicators_wide <- as.data.frame(reverse.code(keys,raw_indicators_wide,
                                 mini = cols_min, maxi = cols_max))

# Remove the little negative thing at the end of reversed column names, otherwise
# they become difficult to use with dplyr

names(indicators_wide) <- str_replace(names(indicators_wide), "-", "")

# Change full stops to underscores

names(indicators_wide) <- str_replace(names(indicators_wide), " ", "_")

indicators <- names(indicators_wide)

indicators_all <- indicators[!indicators %in% "ecoregion_id"]

# Manage outliers ----

pca_input_data <- indicators_wide %>%
  # filter(LPI_2005 < quantile(LPI_2005,
  #                              0.99, na.rm = TRUE)) %>%
  # filter(`extinct 2005-` < quantile(`extinct 2005-`, 
  #                              0.99, na.rm = TRUE)) %>%
  filter(HFP_2005 < quantile(HFP_2005, 
                             0.99, na.rm = TRUE))

# test

# ECO <- east_australia # This ecoregion has a low HFP value (bad outcome)
# HFP <- raw_indicators_wide %>% filter(ecoregion_id == ECO) %>% select(`HFP 2005`)
# HFP
# 
# # Once inverted, it should have a high ecoregion value to represent the bad outcome
# HFP_inv <- indicators_wide %>% filter(ecoregion_id == ECO) %>% 
#   select(`HFP 2005-`)
# HFP_inv

# * Transform ----

# library(moments)
# skewness(indicators_wide$RLI_2005, na.rm = TRUE) # negatively skewed a lot (-2.48)
# skewness(indicators_wide$HFP_2005, na.rm = TRUE) # negatively skewed a bit (-0.89)
# 
# x <- indicators_wide$RLI_2005
# indicators_wide$RLI_2005_log <- 1/(max(x+1) - x) 
# indicators_wide$HFP_2005_log <- log10(indicators_wide$HFP_2005)
# 
# 
# skewness(indicators_wide$RLI_2005_log, na.rm = TRUE) # negatively skewed a lot (-2.48)
# hist(indicators_wide$RLI_2005_log)
# skewness(indicators_wide$HFP_2005_log, na.rm = TRUE)

# * Centre ----

#TODO: Do we need to transform any variables? bc probably need to do so before scaling
# https://www.datanovia.com/en/lessons/transform-data-to-normal-distribution-in-r/

indicators_wide_centred <- indicators_wide %>%
                           mutate_at(c(2:ncol(indicators_wide)), 
                                     funs(c(scale(.)))) 

summary(indicators_wide_centred)

# ** Centred boxplots ----

#' IMPORTANT - this boxplot shows with outliers already removed

indicator_boxplot_data <- reshape2::melt(indicators_wide_centred, 
                          id.vars = 'ecoregion_id')

boxplots <- ggplot(indicator_boxplot_data) +
            geom_boxplot(aes(x = variable, y = value)) +
            theme(axis.text.x = element_text(angle= 45,hjust=1))

boxplots

boxplot(indicators_wide_centred$HFP_2005)

# Manage outliers for the boxplots/correlation data

indicators_wide_centred_trunc <- indicators_wide_centred %>%
  # filter(LPI_2005 < quantile(LPI_2005,
  #                              0.99, na.rm = TRUE)) %>%
  # filter(LPI_2015 < quantile(LPI_2015,
  #                            0.99, na.rm = TRUE)) %>%
  # filter(`extinct 2005-` < quantile(`extinct 2005-`, 
  #                              0.99, na.rm = TRUE)) %>%
  filter(HFP_2005 < quantile(HFP_2005, 
                             0.99, na.rm = TRUE)) %>%
  filter(HFP_2013 < quantile(HFP_2013, 
                             0.99, na.rm = TRUE))

indicator_boxplot_data_2 <- reshape2::melt(indicators_wide_centred_trunc, 
                                           id.vars = 'ecoregion_id')

boxplots_2 <- ggplot(indicator_boxplot_data_2) +
              geom_boxplot(aes(x = variable, y = value)) +
              theme(axis.text.x = element_text(angle = 45,hjust = 1))

boxplots_2

ggsave(file.path(current_analysis_outputs, "indicator_boxplots_2.png"),
       boxplots_2, device = "png")

rm(raw_indicators_long, raw_indicators_wide)

# Finalise correlation data ----

correlation_input_data <- indicators_wide_centred_trunc %>%
  merge(ecoregions_wide, by = "ecoregion_id") 

# Time by time scatterplots ----

# Select the correct time point columns

BHI_data <- correlation_input_data %>%
            dplyr::select(all_of(c("ecoregion_id", "realm", 
                                   "BHI_plants_2005", 
                                   "BHI_plants_2015")))

# BIIri_data <- correlation_input_data %>%
#               dplyr::select(all_of(c("BIIri_2005", "BIIri_2015")))
#                 
# BIIab_all <- correlation_input_data %>%
#   dplyr::select(all_of(c("BIIab_2005", "BIIab_2015")))
# 
# RLI_all <- correlation_input_data %>%
#   dplyr::select(all_of(c("RLI_2005", "RLI_2015")))

extinct_data <- correlation_input_data %>%
  dplyr::select(all_of(c("ecoregion_id", "realm",
                         "extinct_2005", "extinct_2015")))

threatened_data <- correlation_input_data %>%
  dplyr::select(all_of(c("ecoregion_id", "realm",
                         "threatened_2005", "threatened_2015")))

LPI_data <- correlation_input_data %>%
  dplyr::select(all_of(c("ecoregion_id", "realm",
                         "LPI_2005", "LPI_2015")))

time_correlation_input_list <- list(BHI_data, extinct_data, threatened_data,
                                    LPI_data)


time_scatterplots <- list()

for (i in seq_along(time_correlation_input_list)) {
  
  scatterplot_data <- time_correlation_input_list[[i]]
  
  labx <- names(scatterplot_data)[3]
  laby <- names(scatterplot_data)[4]
  names(scatterplot_data) <- c("ecoregion_id", "realm","varx", "vary")
  vargroup <- "realm"
  
  scatterplot <- ggplot(scatterplot_data,aes(x = varx, 
                                             y = vary, 
                                             color = realm)) +
    geom_point() +
    #geom_text(label = scatterplot_data$ecoregion_id) + 
    #geom_smooth(method=lm) +
    labs(x= labx,
         y = laby)
  
  ggsave(file.path(current_analysis_outputs, paste(labx, "_x_", laby,"_NL", ".png", sep = "")),
         scatterplot, device = "png")
  
  time_scatterplots[[i]] <- scatterplot
  
}

time_scatterplots[[4]]


# Grouping variables ----

# Get the names of grouping variables

correlation_input_data_all <- correlation_input_data

correlation_input_data <- correlation_input_data_all %>%
                          merge(pca_data_6[c("ecoregion_id", "cluster")])

grouping_variables_all <- names(correlation_input_data)[!(names(correlation_input_data) %in% 
                                                            indicators_all)]

# Remove numerics
grouping_variables <- grouping_variables_all[!(grouping_variables_all %in% c("ecoregion_id",
                                            "predominant.threat.type",
                                            "headline.threat.type",
                                            "ecoregion.area.km.sq",
                                            "RLI_records",
                                            "LPI_records",
                                            "scenario.numeric",
                                            "predominant.threat.count",
                                            "mean.scientific.publications",
                                            "number.of.endemics",
                                            "mean.human.population.density"))]


# Subset to a single timepoint ----

if (!is.na(timepoint)) {
  
  indicators <- names(correlation_input_data)[str_detect(names(correlation_input_data),
                                                                   timepoint)]
  
  indicators_timepoint_names <- c(grouping_variables_all,
                                  indicators)
  
  correlation_input_data <- correlation_input_data %>%
    dplyr::select(all_of(indicators_timepoint_names))
  
}


# Scatterplots ----

# correlation_input_data <- correlation_input_data_all
# correlation_input_data_all <- correlation_input_data

vargroup <- "realm"

indicator_combinations <- combn(indicators, 2)

indicator_list <- as.list(as.data.frame(indicator_combinations))

# Loop through each combination of indicators and make a scatterplot.  Note
# these are really ugly with all the labels and stuff - they're just to identify
# the weird ecoregions, we can tidy up anything that would actually go into the paper

scatterplots <- list()

for (i in seq_along(indicator_list)) {
  
scatterplot_data <- correlation_input_data[, c("ecoregion_id", 
                                            vargroup,
                                            indicator_list[[i]][1],
                                            indicator_list[[i]][2])]
labx <- names(scatterplot_data)[3]
laby <- names(scatterplot_data)[4]
names(scatterplot_data) <- c("ecoregion_id", "vargroup", "varx", "vary")

scatterplot <- ggplot(scatterplot_data,aes(x = varx, 
                                           y = vary, 
                                           color = vargroup)) +
               geom_point() +
               # geom_text(label = correlation_input_data$ecoregion_id) + 
               # geom_smooth(method=lm) +
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
scatterplots[[28]]

# Get the ecoregions with the strangest results (eg really high HFP, low RLI)

weird_ecoregions <- correlation_input_data %>%
                    filter(ecoregion_id %in% c(325, 400, 637, 20, 632, 548))


# Correlation plots ----

all_grouping_variables <- paste("all", grouping_variables, sep = ".")
nice_grouping_variables <- str_to_title(gsub(".", " ", grouping_variables, fixed=TRUE))
nice_grouping_variables <- str_replace(nice_grouping_variables, "Hfp", "HFP")
nice_grouping_variables <- str_replace(nice_grouping_variables, "Lpi", "LPI")
nice_grouping_variables <- str_replace(nice_grouping_variables, "Rli", "RLI")
nice_grouping_variables <- str_remove(nice_grouping_variables, " Factor")


all_groups <- list()
all_heatmaps <- list()

for (i in seq_along(grouping_variables)) {
  
vargroup <- grouping_variables[i]
varall <- all_grouping_variables[i]

x <- c("ecoregion_id", indicators, vargroup)

group_indicators <- correlation_input_data %>%
                    dplyr::select(all_of(x))

all <- group_indicators
    
group_indicator_list <- split(group_indicators, group_indicators[,vargroup])

all[, ncol(all)] <- varall

group_indicator_list[[varall]] <- all

group_correlations <- list()
group_correlations_dataframes <- list()
group_matrices <- list()

for (j in seq_along(group_indicator_list)) {
  
   group_matrix <- group_indicator_list[[j]] %>%
                   dplyr::select(-ecoregion_id, -vargroup) %>%
                   na.omit(.)
   
   names(group_matrix) <- str_remove(names(group_matrix), "_2005")

   group_matrices[[j]] <- group_matrix
  
}

names(group_matrices) <- names(group_indicator_list)

for (j in seq_along(group_matrices)) {
  
  subset_matrix <- group_matrices[[j]]
  
  n <- nrow(subset_matrix)
  
  name <- paste(names(group_matrices)[j], "n =", n, sep = " ")
  name <- str_replace(name, "&", "and")
  name <- str_replace(name, "/", "_and_")
  name <- str_replace(name, "<", " less than ")
  name <- str_replace(name, ">", " more than ")
  name <- gsub("[()]", "", name)
  
  if (nrow(subset_matrix) < 3) {
    
    print(paste("insufficient data for", names(group_matrices)[j], sep = " "))
    
  } else {
    
    p.mat <- cor_pmat(subset_matrix)
    
    correlation_plot <- ggcorrplot(cor(subset_matrix, method = "spearman"),
                                   title = name, 
                                   p.mat = p.mat, 
                                   type = "lower", insig = "blank",
                                   outline.color = "white",
                                   colors = c("#453781FF", "white", "#287D8EFF"),
                                   lab = TRUE)
    
    correlation_df <- correlation_plot$data[,1:3] %>%
                      merge(indicator_properties[c("indicator","inputs")], 
                            by.x = "Var1", by.y = "indicator") %>%
                      rename(input1 = inputs) %>%
                      merge(indicator_properties[c("indicator","inputs")], 
                            by.x = "Var2", by.y = "indicator") %>%
                      rename(input2 = inputs) %>%
                            mutate(subgroup = names(group_matrices)[j],
                             combination = paste(Var1, "x", Var2, sep = " "),
                             inputs = paste(input1, "x", input2, sep = " ")) %>%
                      mutate(inputs = ifelse(inputs == "land use data x iucn red list",
                                             "iucn red list x land use data", inputs)) %>%
                      rename(coefficient = value) %>%
                      dplyr::select(-Var1, -Var2)
    
    group_directory <- file.path(current_analysis_outputs, paste(vargroup,
                                 "correlation matrices", sep = " "))
    
    if ( !dir.exists( group_directory ) ) {
      
      dir.create( group_directory, recursive = TRUE )
      
    }

    ggsave(file.path(group_directory,
                     paste(name, "indicator_correlation_matrix.png", sep = "_")),
           correlation_plot, device = "png")
  }
  
  group_correlations[[j]] <- correlation_plot
  
  group_correlations_dataframes[[j]] <- correlation_df

  }

all_groups[[i]] <- group_correlations

correlation_dataframe <- do.call(rbind, group_correlations_dataframes)

heatmap <- ggplot(correlation_dataframe, aes(x = subgroup,
                                             y = combination,
                                             fill = coefficient)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(limits = c(-1,1), low = "#453781FF", high = "#287D8EFF") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = paste(nice_grouping_variables[i],"ecoregion subgroups", sep = " "),
       y = "Pairwise indicator comparisons") +
  guides(fill = guide_legend(title = "Correlation\ncoefficient (r)")) +
  theme(panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.background = element_rect(fill = "grey97"),
        legend.position = "bottom",
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8),
        axis.text.x = element_text(size= 7),
        axis.text.y = element_text(size= 7)) +
  facet_wrap(~inputs) + 
  theme(strip.text.x = element_text(size = 7, hjust = 0))

heatmap <- heatmap +  scale_x_discrete(label = function(x) stringr::str_trunc(x, 28)) 

ggsave(file.path(group_directory,
                 paste(grouping_variables[i], "correlation_heatmap.png", sep = "_")),
       heatmap, device = "png")

all_heatmaps[[i]] <- heatmap

}

names(all_groups) <- grouping_variables
names(all_heatmaps) <- grouping_variables
all_heatmaps[[11]]


# PCA ----

# Subset to a single timepoint

if (!is.na(timepoint)) {
  
  indicators <- names(pca_input_data)[str_detect(names(pca_input_data),
                                                         timepoint)]
  
  pca_input_data <- pca_input_data %>%
    dplyr::select(all_of(c("ecoregion_id", indicators)))
  
}

# Prepare data

pca_data_1 <- gather(pca_input_data, indicator, indicator_value, 
                     "BHI_plants_2005":"threatened_2005", factor_key=TRUE)

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
                            endemics.factor,
                            scenario.numeric,
                            predominant.threat.count,
                            mean.human.population.density,
                            number.of.endemics,
                            everything()) 

pca_data_3 <- spread(pca_data_2,indicator,indicator_value) 

pca_data_4 <- pca_data_3[complete.cases(pca_data_3[,16:ncol(pca_data_3)]),]

pca_data_5 <- pca_data_4 

# * Conduct the PCA ----

pca <- prcomp(pca_data_5[,c(16:17,19:28)], center = TRUE, scale = TRUE)
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

# https://rstudio-pubs-static.s3.amazonaws.com/323416_ab58ad22d9e64ba2831569cf3d14a609.html

# Make the PCA plots in a loop

# Remove indicators

grouping_variables <- names(pca_data_5)[!(names(pca_data_5) %in% indicators)]

# Remove numerics
grouping_variables <- grouping_variables[!(grouping_variables %in% c("ecoregion_id",
                                                                     "predominant.threat.type",
                                                                     "headline.threat.type",
                                                                     # "ecoregion.area.km.sq",
                                                                     # "RLI_records",
                                                                     # "LPI_records",
                                                                     #"mean.scientific.publications",
                                                                     "scenario",
                                                                     "predominant.threat.count"))]

pca_plots <- list()

for (i in seq_along(grouping_variables)) {
  
  vargroup <- pca_data_5[,grouping_variables[i]]
  varname <- grouping_variables[i]
  
  pca_plot <- fviz_pca_biplot(pca, geom.ind = "point", pointshape = 21, 
                       pointsize = 2, 
                       col.ind = "black", 
                       palette = "viridis", 
                       addEllipses = TRUE,
                       label = "var",
                       col.var = "black",
                       repel = TRUE,
                       fill.ind = as.factor(vargroup),
                       legend.title = varname) +
    theme(plot.title = element_blank(),
          legend.position = "bottom") +
    theme(legend.text=element_text(size = 7)) +
    guides(fill=guide_legend(nrow=6,byrow=TRUE))
  
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

# Try combining cluster analysis on PCA ----

summary(pca)


pca_clusters <- kmeans(pca$x[,1:4], centers = 3, nstart = 100)

plot(kmeans_input_data[, c("BIIri_2005", "threatened_2005")],
     col = pca_clusters$cluster,
     pch = as.character(pca_data_5$island.status),
     main = paste("k-means clustering of indicator data with", number_clusters, "clusters"),
     xlab = "bii", ylab = "rli")


# Map the ecoregions in their clusters

cluster_map_data <- as.data.frame(cbind(pca_data_5$ecoregion_id, pca_clusters$cluster))

names(cluster_map_data) <- c("ecoregion_id","cluster")

cluster_map_data <- ecoregion_map_renamed %>%
                    merge(cluster_map_data, by = "ecoregion_id")

cluster_map <- plot(cluster_map_data["cluster"])


cluster_one <- cluster_map_data %>% 
  filter(cluster == 1) %>% 
  dplyr::select(ecoregion_id)

cluster_one_ecoregions <- correlation_input_data[correlation_input_data$ecoregion_id %in% 
                                                     cluster_one$ecoregion_id ,] 


cluster_two <- cluster_map_data %>% 
  filter(cluster == 2) %>% 
  dplyr::select(ecoregion_id)

cluster_two_ecoregions <- correlation_input_data[correlation_input_data$ecoregion_id %in% 
                                                   cluster_two$ecoregion_id ,] 


cluster_three <- cluster_map_data %>% 
  filter(cluster == 3) %>% 
  dplyr::select(ecoregion_id)

cluster_three_ecoregions <- correlation_input_data[correlation_input_data$ecoregion_id %in% 
                                                     cluster_three$ecoregion_id ,] 

cluster_three <- cluster_map_data %>% 
                 filter(cluster == 3) %>% 
                 dplyr::select(ecoregion_id)

cluster_three_ecoregions <- correlation_input_data[correlation_input_data$ecoregion_id %in% 
                                          cluster_three$ecoregion_id ,]


# Cluster biplot data

pca_data_6 <-  pca_data_5 %>%
               merge(cluster_map_data, by = "ecoregion_id")


pca_plot <- fviz_pca_biplot(pca, geom.ind = "point", pointshape = 21, 
                            pointsize = 2, 
                            col.ind = "black", 
                            palette = "viridis", 
                            addEllipses = TRUE,
                            label = "var",
                            col.var = "black",
                            repel = TRUE,
                            fill.ind = as.factor(pca_data_6$cluster),
                            legend.title = varname) +
  theme(plot.title = element_blank(),
        legend.position = "bottom") +
  theme(legend.text=element_text(size = 7)) +
  guides(fill=guide_legend(nrow=6,byrow=TRUE))

# Get the response variable  (or what you think it is)

#vargroup <- as.numeric(pca_data_5$island.status == "island")

vargroup <- as.numeric(as.factor(pca_data_5$island.status))

#vargroup <- pca_data_5$scenario.numeric

#vargroup <- as.numeric(as.factor(vargroup))

number_clusters <- length(unique(vargroup))

# Plot the different dimensions to see how clearly the response variable is delineated 

plot(pca$x[, c(1, 2)], col = as.numeric(as.factor(pca_data_5$island.status)), 
     xlab = "PC1", ylab = "PC2")

plot(pca$x[, c(1, 3)], col = as.numeric(as.factor(pca_data_5$island.status)), 
     xlab = "PC1", ylab = "PC2")

# Plot explained variability

par(mfrow = c(1, 2))

# Calculate variability of each component
pr.var <- pca$sdev ^2

# Variance explained by each principal component: pve
pve <- pr.var/sum(pr.var)

# Plot variance explained for each principal component
plot(pve, xlab = "Principal Component", 
     ylab = "Proportion of Variance Explained", 
     ylim = c(0, 1), type = "b")

# Plot cumulative proportion of variance explained
plot(cumsum(pve), xlab = "Principal Component", 
     ylab = "Cumulative Proportion of Variance Explained", 
     ylim = c(0, 1), type = "b")


# Try a model ----

# TODO: use tidymodels ----

library(Rmisc)
library(car)
library(lm.beta)
library(MASS)

# We want the non-inverted values so need to rebuild the input data

response <- raw_indicators_wide[,c("ecoregion_id", "threatened_2005")]

model_input_data <- raw_indicators_wide[,c("ecoregion_id",
                                                "BHI_plants 2005",
                                                "BIIab_2005",
                                                "BIIri_2005",
                                                "extinct_2005",
                                                "HFP_2005",
                                                "RLI_2005")] 

model_input_data <- model_input_data %>%
                    dplyr::filter(HFP_2005 < quantile(HFP_2005, 
                                               0.99, na.rm = TRUE)) %>%
                    mutate(ecoregion_id = as.character(ecoregion_id))

model_input_data <- model_input_data %>%
                    merge(ecoregions_wide, by = "ecoregion_id")

scaling_index <- unlist(lapply(model_input_data, is.numeric))

eco_id <- model_input_data$ecoregion_id
non_numeric <- model_input_data[,!scaling_index]

model_input_data <- as.data.frame(scale(model_input_data[,scaling_index], 
                                        scale = TRUE, center = TRUE)) 

model_input_data <- cbind(eco_id, model_input_data, non_numeric)

model_input_data <- model_input_data %>%
                    dplyr::select(-eco_id) %>%
                    dplyr::select(ecoregion_id, everything()) %>%
                    mutate(ecoregion_id = as.numeric(ecoregion_id))

model_input_data <- response %>%
                    merge(model_input_data, by = "ecoregion_id")


model_input_data <- model_input_data %>%
                    filter(RLI_records > 600)

scatt <- scatterplot(threatened_2005 ~ BIIri_2005 , data = model_input_data)

test_glm <- glm(BIIri_2005 ~ 1, family = gaussian, 
                 data = model_input_data)


step_glm <- stepAIC(test_glm, scope = list(upper = ~ HFP_2005 + RLI_records + 
                                           threatened_2005 +
                                           RLI_2005 +
                                           extinct_2005 +
                                           scenario +
                                           included.in.HFP + realm +
                                           LPI_records + headline.threat.type +
                                           Biome + mean.scientific.publications +
                                           island.status +
                                           number.of.endemics +
                                           mean.human.population.density,
                                           lower = ~ 1))

summary(step_glm)

# Map scenarios ----

# * Global ----

if(load_map == TRUE) {
  
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
  
  # * Realms ----
  
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
  
}

# K-means cluster analysis ----

# rownames(pca_input_data) <- pca_input_data$ecoregion_id
# 
# kmeans_input_data <- pca_input_data %>% 
#   dplyr::select(- ecoregion_id)
# 
# # Remove incomplete cases
# 
# kmeans_input_data <- kmeans_input_data[complete.cases(kmeans_input_data),]
# 
# # Scale 
# 
# kmeans_input_data <- scale(kmeans_input_data)
# 
# # Find best/truest number of clusters
# 
# wss <- 0
# 
# set.seed(1)
# 
# # Look over 1 to 15 possible clusters
# for (i in 1:15) {
#   # Fit the model: km.out
#   km.out <- kmeans(kmeans_input_data, centers = i, nstart = 20, iter.max = 200)
#   # Save the within cluster sum of squares
#   wss[i] <- km.out$tot.withinss
# }
# 
# # Produce a scree plot
# plot(1:15, wss, type = "b", 
#      xlab = "Number of Clusters", 
#      ylab = "Within groups sum of squares")
# 
# # Looks like 3 clusters is best
# 
# number_clusters <- 3
# 
# wisc.km <- kmeans(scale(kmeans_input_data), centers = number_clusters, nstart = 100)
# 
# # Compare k-means to actual groupings
# 
# x <- table(wisc.km$cluster, pca_data_5$lpi.records.factor)
# 
# true <- colSums(x)
# 
# model <- rowSums(x)
# 
# check_model <- cbind(true,model)
# 
# check_model
# 
# # View the resulting model
# wisc.km$tot.withinss
# 
# # Plot of BII and RLI by cluster membership
# 
# plot(kmeans_input_data[, c("BIIri_2005", "threatened_2005")],
#      col = wisc.km$cluster,
#      #pch = as.character(pca_data_5$Biome),
#      main = paste("k-means clustering of indicator data with", number_clusters, "clusters"),
#      xlab = "bii", ylab = "rli")
# 
# par(mfrow = c(2, 3))


# # Try hierarchical cluster ----
# 
# # Remove ID column
# 
# hclust_input_data <- kmeans_input_data
# 
# # Make cluster model
# 
# hclust.indicators <- hclust(dist(hclust_input_data), method = "complete")
# plot(hclust.indicators)
# abline(h = 9, col = "red")
# 
# cut_indicators <- cutree(hclust.indicators, k = number_clusters)
# 
# cut_indicators_df <- cbind(rownames(hclust_input_data), cut_indicators)
# 
# x <- table(cut_indicators,vargroup)
# 
# true <- colSums(x)
# 
# model <- rowSums(x)
# 
# check_model <- cbind(true,model)
# 
# check_model

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
