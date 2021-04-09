
# Using R version 4.0.3

## Clear the space
rm(list = ls()) # clear memory

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

# TODO ----

#' TODO: Update confidence intervals from correlation plots
#' TODO: Generate lists of ecoregions that are excluded
#' TODO: Prepare stats summary for indicators and variables (or can use 
#' boxplots? add jitter?)
#' TODO: Split input prep from data visualisation
#' TODO: Remove superfluous steps
#' TODO: Make another folder called something like MS figures? for stuff that's
#' definitely going in
#' TODO: Fix ecoregion lookup function

# PCA and Clustering
library(factoextra)
library(FactoMineR)
library(corrplot)
library(ape)
#library(MASS)
library(ggpubr)

# Data handling and table reshaping
library(tidyverse)
library(reshape2)
library(devtools)
library(data.table)
library(rlist)
library(e1071)
library(psych)
library(arules)

# Plotting
library(ggplot2)
library(RColorBrewer)
library(ggbiplot) # Won't install
library(plot3D)
library(plotly)
library(viridis)
library(png)
library(gridExtra)
library(GGally)
library(ggcorrplot)
library(visdat)
library(cowplot)

# Maps
library(sf)
library(leaflet)


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
timepoints <- c("2005", "2008")
load_map <- TRUE
indicators_to_remove <- c("RLIother", "RLILU")

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

make_heatmap <- function(data) {

heatmap <- ggplot(data, aes(x = subgroup,
                            y = combination,
                            fill = coefficient)) +
  geom_tile(color = "white") +
  geom_text(aes(label=coefficient), size = 2) +
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
  facet_grid(ind_group ~., scales = "free_y") +
  theme(strip.text.x = element_text(size = 7, hjust = 0))

heatmap <- heatmap +  scale_x_discrete(label = function(x) stringr::str_trunc(x, 28)) 

return(heatmap)

}

lookup_ecoregion <- function(eco_id_number, data) {
  
eco_name <- ecoregion_map_renamed %>% filter(ecoregion_id == eco_id_number) %>%
            select(ecoregion_id, ECO_NAME) %>%
            st_drop_geometry()

print(eco_name)

eco_vals <- data %>% 
  filter(ecoregion_id == eco_id_number) %>%
  select(all_of(indicators)) 

averages <- data %>%
            select(all_of(indicators)) %>%
            colMeans(na.rm = TRUE)

id_col <- c(paste("eco", eco_id_number, "scores", sep = " "), "global_mean")
eco_vals <- rbind(eco_vals, averages)

eco_vals <- cbind(id_col, eco_vals)

eco_vals

}

x <- group_dataframes[[1]][[1]]
df <- x


# matrix should be indicators as columns, rownames are eco_id
# group is string for group, subgroup is a string denoting the subgroup (eg mangroves)

make_subgroup_scatterplot <- function(df, group, subgroup, group_directory) {
  
  
  subgroup <- str_replace(subgroup, " - ", "_")
  subgroup <- str_replace(subgroup, " & ", "and")
  subgroup <- str_replace(subgroup, " , ", "_")
  subgroup <- str_replace(subgroup, " . ", "_")
  subgroup <- str_replace(subgroup, " ", "_")
  subgroup <- str_replace(subgroup, "/", "_")
  
  # Set rownames back to column
  # df <- tibble::rownames_to_column(df)
  # 
  # df <- df %>%
  #       dplyr::rename(ecoregion_id = rowname)
  
  ecoregion_id <- as.numeric(df$ecoregion_id)
  
  # Get all possible combinations of indicators from the df
  indicator_combos <- combn(names(df[-1]), 2)
  
  # Now loop through the combos
  i <- i + 1
  subgroup_scatterplots <- list()
  
  for (i in seq_along(1:ncol(indicator_combos))) {
   
  print(i)
  combo <- as.vector(indicator_combos[,i])
  
  print(combo)
  
  plotname <- paste(combo[1], combo[2], subgroup, sep = "_")
  
  plotname
  
  x_axis <- df[ , grepl( combo[1] , names( df ) ) ]
  range(x_axis)
  
  y_axis <- df[ , grepl( combo[2] , names( df ) ) ]
  range(y_axis)
  
  plot_df <- as.data.frame(cbind(ecoregion_id, x_axis, y_axis))

  scatterplot <- ggplot(plot_df, aes(x = x_axis, y = y_axis)) +
    geom_point() +
    labs(x = combo[1],
         y = combo[2]) + 
    stat_cor(method = "spearman") +
    ggtitle(paste(group, subgroup, sep = " - ")) +
    geom_text(label = plot_df$ecoregion_id, nudge_x = 0.1)
  
  ggsave(file.path(group_directory,
                   paste(plotname,
                         "scatterplot.png",
                         sep = "_")),
         scatterplot,  device = "png")
  
  subgroup_scatterplots[[i]] <- scatterplot
  
  }

  subgroup_scatterplots

}

# y <- make_subgroup_scatterplot(x, "cluster", "cluster one", scatter_directories[[1]])
# 
# y[[18]]
# y[[1]]

spearman_CI <- function(x, y, alpha = 0.05){
  
  rs <- cor(x, y, method = "spearman")
  ci <- sort(tanh(atanh(rs) + c(-1,1)*sqrt((1+rs^2/2)/(n-3))*qnorm(p = alpha/2)))
  
  z <- data.frame(rs = rs, lower_ci = ci[1], upper_ci = ci[2])
  return(z)
}

# Ecoregion data ----

raw_ecoregions_wide <- readRDS(file.path(analysis_inputs,
                               "global_ecoregions_2017_ecoregion_values_master_wide.rds"))

# Keep a copy of unaltered ecoregion data

ecoregions_wide <- raw_ecoregions_wide 

 
# ecoregions_wide$area.factor <- cut(ecoregions_wide$ecoregion.area.km.sq, breaks = 3, 
#                                    labels = c("Small_ecoregions", 
#                                               "Medium_ecoregions",
#                                               "Large_ecoregions"))

ecoregions_wide$area.factor <- discretize(ecoregions_wide$ecoregion.area.km.sq, 
                                   method = "frequency",
                                   breaks = 3, 
                                   labels = c("Small_ecoregions", 
                                              "Medium_ecoregions",
                                              "Large_ecoregions"),
                                   ordered_result = TRUE)

# Split the ecoregions by lpi data by fewer than 20 (minimum population sample
# required for SDMs) and more than 20

ecoregions_wide$lpi.records.factor <- cut(ecoregions_wide$LPI_records, 
                                          breaks = c(-Inf,19, Inf),
                                          labels = c("less than 20",
                                                     "more than 20"))

# Split ecoregions by RLI species fewer than 400 and more than 400 (based
# on min sample in Henriques et al 2020)

ecoregions_wide$rli.records.factor <- cut(ecoregions_wide$RLI_records, 
                                          breaks = c(-Inf,400, Inf),
                                          labels = c("less than 400",
                                                     "more than 400"))
table(ecoregions_wide$rli.records.factor)

ecoregions_wide <- ecoregions_wide %>% 
                   mutate(included.in.HFP = ifelse(included.in.HFP == 1, 
                          "Threat related to HFP",
                          ifelse(included.in.HFP == 0,
                          "Threat external to HFP",
                                         NA)))

ecoregions_wide <- ecoregions_wide %>%
  mutate(scenario = as.factor(paste(rli.records.factor, included.in.HFP, sep = " & ")),
         scenario.numeric = as.factor(as.numeric(scenario))) 

ecoregions_wide$endemics.factor <- cut(ecoregions_wide$number.of.endemics, 
                                       breaks = c(-Inf, 0 , Inf),
                                       labels = c("no endemics",
                                                  "endemics"))


table(ecoregions_wide$endemics.factor)

ecoregions_wide$High.beta.area.factor <- discretize(ecoregions_wide$High.beta.area, 
                                      method = "interval",
                                      breaks = 4,
                                      labels = c("very low beta",
                                                 "low beta",
                                                 "medium beta",
                                                 "high beta"))

table(ecoregions_wide$High.beta.area)

# Remove any unwanted grouping variables

names(ecoregions_wide)

ecoregions_wide <- ecoregions_wide %>%
                   dplyr::select(-scenario, -scenario.numeric, 
                                 -mean.scientific.publications,
                                 -headline.threat.type)

# Convert any characters into factors

ecoregions_wide <- ecoregions_wide %>% 
                   mutate(across(where(is.character), as.factor)) 


# Get the names of factors

grouping_variables <- names(dplyr::select_if(ecoregions_wide, is.factor))


# Get names of numeric variables 

numeric_variables <- names(dplyr::select_if(ecoregions_wide, is.numeric))
numeric_variables <- numeric_variables[!str_detect(numeric_variables, "ecoregion_id")]


# Indicator data ----

# Indicator data

indicator_properties <- read.csv(file.path(analysis_inputs,
                                          "indicator_properties.csv"))

#TODO: get this out of indicator_properties
indicator_relationships <- read.csv(file.path(analysis_inputs,
                                           "indicator_input_relationships.csv"))

raw_indicators_long_all <- readRDS(file.path(analysis_inputs,
                          "global_ecoregions_2017_indicator_values_master.rds"))

raw_indicators_long_all$indicator_year <- str_replace(raw_indicators_long_all$indicator_year, " ", "_")


raw_indicators_long <- raw_indicators_long_all %>%
                       dplyr::select(ecoregion_id, indicator_year,
                                     raw_indicator_value) %>%
                       distinct(.)

raw_indicators_wide <- raw_indicators_long %>%
                       spread(key = indicator_year, 
                             value = raw_indicator_value) 

dim(raw_indicators_wide)

raw_names <- names(raw_indicators_wide)
new_raw_names <- str_replace(raw_names, " ", "_")
names(raw_indicators_wide) <- new_raw_names

# Get single timepoint indicator names

# Get single timepoint indicators

indicators_05 <- names(raw_indicators_wide)[str_detect(names(raw_indicators_wide),
                                                  timepoints[[1]])]

indicators_08 <- names(raw_indicators_wide)[str_detect(names(raw_indicators_wide),
                                                  timepoints[[2]])]

# Remove LPI 2008 record which is the only indicator with timepoints for 05 and 08

indicators_08 <- indicators_08[!str_detect(indicators_08,
                                           "LPI_2008")]

indicators <- c(indicators_05, indicators_08)

# Threat scheme

threat_scheme <- read.csv(file.path("N:/Quantitative-Ecology/Simone/extinction_test/inputs/iucn_threats\\iucn_threat_classification_scheme.csv"))

headline_threats <- threat_scheme %>%
  dplyr::select(headline_name) %>%
  distinct(.) %>%
  pull()

headline_threats <- c(headline_threats, "All threats")

# Ecoregion map ----

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

ecoregions <- ecoregion_map_renamed %>% 
              st_drop_geometry()
# Countries

ecoregion_countries <- readRDS(file.path(parent_outputs, "version_3",
                        "2020-08-10_database_output_files",
                        "ecoregion_country_data.rds"))

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

dim(indicators_wide)

# * Fix column names ----

# Remove the little negative thing at the end of reversed column names, otherwise
# they become difficult to use with dplyr

names(indicators_wide) <- str_replace(names(indicators_wide), "-", "")

# Change full stops to underscores

names(indicators_wide) <- str_replace(names(indicators_wide), " ", "_")

indicators <- names(indicators_wide)

indicators_all <- indicators[!indicators %in% "ecoregion_id"]

# * Manage outliers ---- 

# A couple of big outliers in LPI values
hist(indicators_wide$LPI_2005, breaks = 20)
max(indicators_wide$LPI_2005, na.rm = TRUE)

# Get the 95th percentile

lpi_95 <- quantile(indicators_wide$LPI_2005, 
                   probs = 0.95, 
                   na.rm = TRUE)
lpi_95

lpi_mx_eco <- raw_indicators_wide %>%
  filter(LPI_2005 == max(LPI_2005, na.rm = TRUE)) %>%
  select(ecoregion_id) %>%
  pull(.)

# Looks much better
hist(indicators_wide$LPI_2005)
max(indicators_wide$LPI_2005, na.rm = TRUE)

# Note, max positive LPI score goes to ecoregion 675 Po Basin, which also gets the third most
# negative (ie contradictory) score for HFP
ecoregion_map_renamed %>% filter(ecoregion_id == lpi_mx_eco)


# Just remove LPI outliers, given we know it is a dodgier dataset and volatile
indicators_wide_2 <- indicators_wide %>%
  mutate(LPI_2005 = ifelse(LPI_2005 > lpi_95,
                           NA, LPI_2005)) 

dim(indicators_wide_2)

# HFP values - have a couple of exceptionally low values

hist(indicators_wide$HFP_2005, breaks = 20)
min(indicators_wide$HFP_2005, na.rm = TRUE)
hfp_mx <- max(raw_indicators_wide$HFP_2005, na.rm = TRUE)
hfp_mx_eco <- raw_indicators_wide %>%
              filter(HFP_2005 == max(HFP_2005, na.rm = TRUE)) %>%
              select(ecoregion_id) %>%
              pull(.)

# Looking at HFP map, bermuda cells do have very high scores (eg 48), which
# matches with the ecoregions WWF description of severe degradation, so
# don't remove outliers

lookup_ecoregion(hfp_mx_eco)

# Extinctions also have some very high values, however these can be verified and are
# correct, so don't remove

hist(indicators_wide$extinct_2008, breaks = 20)
min(indicators_wide$extinct_2008, na.rm = TRUE)

# * Remove unwanted indicators ----

indicators_wide_3 <- indicators_wide_2 %>%
                   dplyr::select(-number_extinct_2008, -number_extinct_2016,
                                 -AmphRLI_2008, -BirdRLI_2008, -BirdRLI_2016,
                                 -MammRLI_2008)

dim(indicators_wide_3)

# PCA ----

# * Prepare PCA input data ----

pca_input_data <- indicators_wide_3

head(pca_input_data)

# * Subset to a single timepoint ----

  indicators_05 <- names(pca_input_data)[str_detect(names(pca_input_data),
                                                 timepoints[[1]])]
  
  indicators_08 <- names(pca_input_data)[str_detect(names(pca_input_data),
                                                 timepoints[[2]])]
  
  # Remove LPI 2008 record which is the only indicator with timepoints for 05 and 08
  
  indicators_08 <- indicators_08[!str_detect(indicators_08,
                                            "LPI_2008")]
  
  indicators <- c(indicators_05, indicators_08)
  
  # Subset so we only have indicators from the years we need
  
  pca_input_data <- pca_input_data %>%
                    dplyr::select(all_of(c("ecoregion_id", indicators)))
  
  dim(pca_input_data)
  head(pca_input_data)
  
# Prepare data
  
# Add the grouping ecoregion variables back in

pca_data_1 <- pca_input_data %>%
        merge(ecoregions_wide[c("ecoregion_id",
                                "Biome",
                                "disturbance.year",
                                "High.beta.area.factor",
                                "island.status",
                                "predominant.threat.type",
                                "realm",
                                "lpi.records.factor",
                                "rli.records.factor")],
              by = "ecoregion_id", all.y = TRUE)

dim(pca_data_1)

head(pca_data_1)
dim(pca_data_1)

# Remove incomplete cases - NOTE - while including LPI, removes 3/4 of the data

pca_data_2 <- pca_data_1[complete.cases(pca_data_1[,2:9]),]
dim(pca_data_2)

pca_with_lpi_excluded_ecoregions <- pca_data_1[!complete.cases(pca_data_1[,2:9]),]
dim(pca_with_lpi_excluded_ecoregions)

# Check that together the datasets = number of ecoregions (should = TRUE)
(nrow(pca_with_lpi_excluded_ecoregions) + nrow(pca_data_2)) == nrow(ecoregion_map_renamed) - 1 # exclude rock and ice

head(pca_with_lpi_excluded_ecoregions)

pca_with_lpi_excluded_ecoregions <- pca_with_lpi_excluded_ecoregions %>% 
                                    select(ecoregion_id) %>%
                                    merge(ecoregions,
                                          by = "ecoregion_id")

names(pca_with_lpi_excluded_ecoregions) <- tolower(pca_with_lpi_excluded_ecoregions)

write.csv(pca_with_lpi_excluded_ecoregions, file.path(
          current_analysis_outputs, "ecoregions_excluded_from_lpi_pca.csv"))

# * Indicators only PCA ----

# Only conducts the PCA on the main 8 indicators (excludes RLI land use and RLI non-land use)

indicators_for_pca <- indicators[!str_detect(indicators, indicators_to_remove[1])]
indicators_for_pca <- indicators_for_pca[!str_detect(indicators_for_pca , 
                                                      indicators_to_remove[2])]

indicator_only_pca_data <- pca_data_2[,c("ecoregion_id", indicators_for_pca)]

names(indicator_only_pca_data)

# ** PCA ----

# Tidy up names

pl.data <- indicator_only_pca_data[,indicators_for_pca]
rownames(pl.data) <- indicator_only_pca_data[,1]
dimC <- dim(pl.data)

# Scale the data to a mean of 0 and sd of 1
pl.data <- scale(pl.data)
summary(pl.data)

names(pl.data) <- colnames(pl.data)
pl.data <- as.data.frame(pl.data)

sd(pl.data$BHI_plants_2005)

# Save a copy of the inputs
write.csv(pl.data, file.path(current_analysis_outputs, "indicator_only_pca_data.csv"))

# Get formula for pca

pc.f <- formula(paste("~", paste(names(pl.data), collapse = "+")))

# Run pca

pl.pca <- princomp(pc.f, cor=TRUE, data=pl.data)


# Print out PCA loadings
pl.pca$loadings

# Print out eigenvalues
pl.pca$sd^2

# Print PCA summary
summary(pl.pca)

# Plot results - look to see number of PCA axes to retain
plot(pl.pca, type = "lines")

### How many dimensions do we include?
# Plot points
text(pl.pca$scores, labels = as.character(row.names(pl.data)), pos=1, cex=0.7)

# Can see ecoregion 374 sticks out

lookup_ecoregion(374)

# Biplot of PCA - can see LPI gives little contribution

tiff(file = file.path(current_analysis_outputs, "LPI_inclusive_pca_biplot.tiff"), 
     units = "in", width=10, height=5, res = 300)

biplot(pl.pca, cex=0.8, col=c(1,8))

dev.off()

# * No LPI indicator PCA ----

# Repeat PCA but remove LPI because it reduces the data points and has little
# influence on any of the components

# Remove LPI before filtering for complete cases

pca_data_4 <- pca_data_1 %>%
              dplyr::select(-LPI_2005)

pca_data_5 <- pca_data_4[complete.cases(pca_data_4[,2:8]),]
dim(pca_data_5)

# Only conducts the PCA on the main 7 indicators (excludes LPI)

indicators_for_pca_2 <- indicators_for_pca[!str_detect(indicators_for_pca,"LPI_2005")]

indicator_only_pca_2_data <- pca_data_5[,c("ecoregion_id", indicators_for_pca_2)]

names(indicator_only_pca_2_data)
dim(indicator_only_pca_2_data)

# ** PCA ----

# Tidy up names

pl2.data <- indicator_only_pca_2_data[,indicators_for_pca_2]
rownames(pl2.data) <- indicator_only_pca_2_data[,1]
dimC <- dim(pl2.data)

# Scale the data to a mean of 0 and sd of 1
pl2.data <- scale(pl2.data)
summary(pl2.data)

names(pl2.data) <- colnames(pl2.data)
pl2.data <- as.data.frame(pl2.data)
dim(pl2.data)
sd(pl2.data$BHI_plants_2005)

# Save a copy of the inputs
write.csv(pl2.data, file.path(current_analysis_outputs, "indicator_only_pca_2_data.csv"))

# Get formula for pca

pc.f2 <- formula(paste("~", paste(names(pl2.data), collapse = "+")))

# Run pca

pl2.pca <- princomp(pc.f2, cor=TRUE, data=pl2.data)

# Print out PCA loadings
pca_loadings <- pl2.pca$loadings

pca_loadings <- as.table(pca_loadings)

write.csv(pca_loadings, file.path(current_analysis_outputs, 
                                  "pca_loadings.csv"))

# Print out eigenvalues
pl2.pca$sd^2

# Print PCA summary
summary(pl2.pca)

# Look at the amount of explained variance

cumpro <- cumsum(pl2.pca$sdev^2 / sum(pl2.pca$sdev^2))
plot(cumpro[0:5], xlab = "PC #", 
     ylab = "Amount of explained variance",
     main = "Cumulative variance plot")

# Plot results - look to see number of PCA axes to retain
plot(pl2.pca, type="lines")

### Looks like need to consider the first three axes still?

# Plot points
text(pl2.pca$scores, labels=as.character(row.names(pl2.data)), pos=1, cex=0.7)

# Biplot of PCA
biplot(pl2.pca, cex=0.8, col=c(1,8))

# Alt approach to plotting

PoV <- pl2.pca$sdev^2/sum(pl2.pca$sdev^2)
fviz_eig(pl2.pca)

# Put on row named
#Rename Col 1
View(pl2.data)

pcx <- pl2.pca$scores[,1]
pcy <- pl2.pca$scores[,2]
pcz <- pl2.pca$scores[,3]

pcxlab <- paste("PC1 (", round(PoV[1] * 100, 2), "%)")
pcylab <- paste("PC2 (", round(PoV[2] * 100, 2), "%)")
pczlab <- paste("PC3 (", round(PoV[3] * 100, 2), "%)")

View(pl2.pca$scores)

tiff(file = file.path(current_analysis_outputs, "pca_3D_plot.tiff"), 
     units = "in", width=10, height=5, res = 300)

# 3D as points
scatter3D(pcx, pcy, pcz, bty = "g", pch = 20, cex = 2, 
          col = gg.col(100), theta = 150, phi = 0, main = "PCA Scores", xlab = pcxlab,
          ylab =pcylab, zlab = pczlab)
text3D(pcx, pcy, pcz,  labels = rownames(pl2.pca$scores), add = TRUE, colkey = FALSE, cex = 0.7)

dev.off()

###Use prcomp() instead - this uses singular value decomposition 

res.pca2 <- prcomp(pl2.data, scale = TRUE)
res.pca2$x

# Visualize eigenvalues (scree plot). Show the percentage of variances explained by each principal component.
fviz_eig(res.pca2)

#Graph of individuals. Individuals with a similar profile are grouped together.
fviz_pca_ind(res.pca2,
             col.ind = "contrib", # Color by congtribution
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,     # Avoid text overlapping
             #label=SP_as_col$Year
) +
  labs(title ="PCA", x = "PC1", y = "PC2")

tiff(file = file.path(current_analysis_outputs, "pca_biplot.tiff"), 
     units = "in", width=10, height=5, res = 300)

# Graph of variables
fviz_pca_var(res.pca2,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#36648B", "#FFA500", "#8B2500"),
             repel = TRUE     # Avoid text overlapping
)

dev.off()

# * Get ecoregion clusters ----

# Compute hierarchical clustering on principal components
res.pca24 <- PCA(pl2.data, ncp = 5, graph = FALSE)
res.hcpc2 <- HCPC(res.pca24, graph = FALSE)

data.hcpc2 <- as.data.frame(setDT(res.hcpc2$data.clust, keep.rownames = TRUE)[])

data.hcpc2 <- data.hcpc2 %>%
              dplyr::rename(cluster = clust,
                            ecoregion_id = rn)

names(data.hcpc2)
unique(data.hcpc2$cluster)

# Vertical dendrogram

tiff(file = file.path(current_analysis_outputs, "pca_cluster_dendrogram.tiff"), 
     units = "in", width=10, height=5, res = 300)

fviz_dend(res.hcpc2, 
          cex = 0.7,                     # Label size
          palette = "jco",               # Color palette see ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
          rect_border = "jco",           # Rectangle color
          labels_track_height = 0.8      # Augment the room for labels
)

dev.off()

# Factor map

fviz_cluster(res.hcpc2,
             repel = TRUE,            # Avoid label overlapping
             show.clust.cent = TRUE, # Show cluster centers
             palette = "jco",         # Color palette see ?ggpubr::ggpar
             ggtheme = theme_bw(),
             main = "Factor map"
)


# Map & plot clusters ----

cluster_map_data <- ecoregion_map_renamed %>%
                    merge(data.hcpc2[c("ecoregion_id", "cluster")], 
                          by = "ecoregion_id")

cluster_map <-  ggplot(cluster_map_data) +
  geom_sf(aes(fill = cluster), colour = "black", 
          size = 0.05, show.legend = 'fill') +
  scale_fill_viridis_d(alpha = .8,
                       na.value = "grey70") +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  labs(fill = "Cluster") +
  theme(legend.position = "right")

cluster_map

ggsave(file.path(current_analysis_outputs,
                 paste(location, eco_version,
                       "cluster_map.png",
                       sep = "_")), 
       cluster_map,  device = "png")

# Add clusters to our data


# Look at the ecoregions in each cluster

cluster_one_ecoregions <- cluster_map_data %>% 
                          filter(cluster == 1) %>%
                          st_drop_geometry(.)

cluster_two_ecoregions <- cluster_map_data %>% 
                          filter(cluster == 2) %>%
                          st_drop_geometry(.)

cluster_three_ecoregions <- cluster_map_data %>% 
                          filter(cluster == 3) %>%
                          st_drop_geometry(.)

write.csv(cluster_one_ecoregions, file.path(current_analysis_outputs, 
                                            "cluster_one_ecoregions.csv"))
write.csv(cluster_two_ecoregions, file.path(current_analysis_outputs, 
                                            "cluster_two_ecoregions.csv"))
write.csv(cluster_three_ecoregions, file.path(current_analysis_outputs,
                                              "cluster_three_ecoregions.csv"))

# Cluster biplot data

cluster_biplot_data <-  pca_data_2 %>%
  merge(cluster_map_data, by = "ecoregion_id")


pca_cluster_plot <- fviz_pca_biplot(pl2.pca, geom.ind = "point", pointshape = 21, 
                                    pointsize = 2, 
                                    col.ind = "black", 
                                    palette = c("#453781FF","#287D8EFF", "#DCE319FF"),
                                    alpha = 0.4,
                                    addEllipses = TRUE,
                                    label = "var",
                                    col.var = "black",
                                    repel = TRUE,
                                    fill.ind = as.factor(cluster_biplot_data$cluster),
                                    legend.title = "Cluster") +
                                    theme(plot.title = element_blank(),
                                          legend.position = "right") +
                                    theme(legend.text=element_text(size = 7)) +
                                    guides(fill=guide_legend(nrow=6,byrow=TRUE))

pca_cluster_plot

ggsave(file.path(current_analysis_outputs,
                 paste(location, eco_version,
                       "cluster_biplot.png",
                       sep = "_")), 
       pca_cluster_plot,  device = "png") 

rm(cluster_map)


# Prepare correlation data ----

# * Centre ----

#TODO: Do we need to transform any variables? bc probably need to do so before scaling
# https://www.datanovia.com/en/lessons/transform-data-to-normal-distribution-in-r/

indicators_wide_centred <- indicators_wide %>%
                           mutate_at(c(2:ncol(indicators_wide)), 
                                     funs(c(scale(.)))) 

summary(indicators_wide_centred)

# ** Centred boxplots ----

#' IMPORTANT - this boxplot shows with outliers already removed

# Get 2005 only
indicator_boxplot_data_wide <- indicators_wide_centred %>%
  dplyr::select(all_of(c("ecoregion_id", indicators)))

# Convert back into long format
indicator_boxplot_data <- reshape2::melt(indicator_boxplot_data_wide, 
                          id.vars = 'ecoregion_id')

boxplots <- ggplot(indicator_boxplot_data) +
            geom_boxplot(aes(x = variable, y = value)) +
            theme(axis.text.x = element_text(angle= 45,hjust=1))

boxplots

ggsave(file.path(current_analysis_outputs, "indicator_boxplots.png"),
       boxplots, device = "png")

rm(raw_indicators_wide)

# Finalise correlation data ----

correlation_input_data_all <- indicators_wide_centred %>%
  merge(ecoregions_wide, by = "ecoregion_id", all.y = TRUE) %>%
                          merge(cluster_map_data[c("ecoregion_id", "cluster")],
                                by = "ecoregion_id") %>%
                          dplyr::select(-geometry)

# Check it looks right (should have all timepoints)

summary(correlation_input_data_all)

saveRDS(correlation_input_data_all,
        file.path(current_analysis_outputs, "correlation_input_data.RDS"))

# Time by time scatterplots ----

# Select the correct time point columns

BHI_data <- correlation_input_data_all %>%
            dplyr::select(all_of(c("ecoregion_id", "realm", 
                                   "BHI_plants_2005", 
                                   "BHI_plants_2015")))

 
RLI_birds_data <- correlation_input_data_all %>%
  dplyr::select(all_of(c("ecoregion_id", "realm", 
                         "BirdRLI_2008", 
                         "BirdRLI_2016")))


extinct_data <- correlation_input_data_all %>%
  dplyr::select(all_of(c("ecoregion_id", "realm",
                         "extinct_2008", "extinct_2016")))

threatened_data <- correlation_input_data_all %>%
  dplyr::select(all_of(c("ecoregion_id", "realm",
                         "threatened_2008", "threatened_2016")))

LPI_data <- correlation_input_data_all %>%
  dplyr::select(all_of(c("ecoregion_id", "realm",
                         "LPI_2005", "LPI_2015")))

time_correlation_input_list <- list(BHI_data, RLI_birds_data, extinct_data, 
                                    threatened_data, LPI_data)


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
         y = laby) + stat_cor(method = "spearman")
  
  ggsave(file.path(current_analysis_outputs, paste(labx, "_x_", laby,"_NL", ".png", sep = "")),
         scatterplot, device = "png")
  
  time_scatterplots[[i]] <- scatterplot
  
}

time_scatterplots[[4]]

# Grouping variables ----

# Get the names of grouping variables


# 
# grouping_variables <- names(correlation_input_data)[!(names(
#                           correlation_input_data) %in% 
#                                 indicators_all)]
# 
# grouping_variables <- grouping_variables[!grouping_variables %in% "ecoregion_id"]

# Get categorical grouping variables only 

# groups <- correlation_input_data_all %>%
#           dplyr::select(all_of(grouping_variables))
# 
# lapply(groups, class)

# Subset to a single timepoint ----

dim(correlation_input_data_all)

correlation_input_data <- correlation_input_data_all %>%
    dplyr::select(all_of(c("ecoregion_id","cluster",indicators,
                           grouping_variables)))

names(correlation_input_data)

# Cluster boxplots ----

# Plot data by cluster

library(forcats)
library(hrbrthemes)

cluster_boxplot_data_wide <- correlation_input_data_all %>%
  dplyr::select(all_of(c("ecoregion_id", "cluster", 
                         indicators, numeric_variables)))

summary(cluster_boxplot_data_wide)

# Make sure all the values are scaled

cluster_boxplot_data_wide <- scale(cluster_boxplot_data_wide[,c(indicators, 
                                                                numeric_variables)])

cluster_boxplot_data_wide <- cbind(correlation_input_data_all[,c("ecoregion_id",
                                                                 "cluster")], 
                                   cluster_boxplot_data_wide)
summary(cluster_boxplot_data_wide)

cluster_boxplot_data_wide <- cluster_boxplot_data_wide %>%
                             dplyr::select(-predominant.threat.count,
                                           -number.of.endemics) %>%
                             rename(human_pop_density = mean.human.population.density)

# Melt back into long format

cluster_boxplot_data <- reshape2::melt(cluster_boxplot_data_wide, 
                                       measure.vars = c(indicators, 
                                                        "ecoregion.area.km.sq",
                                                        "High.beta.area",
                                                        "LPI_records",
                                                        "human_pop_density",
                                                        "RLI_records"))

head(cluster_boxplot_data)

names(cluster_boxplot_data) <- c("ecoregion_id", "cluster", 
                                 "Indicator", "Indicator value")

cluster_boxplots <- ggplot(cluster_boxplot_data) +
                    geom_boxplot(aes(x = Indicator, y = `Indicator value`,
                                     fill = Indicator)) +
                    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                    facet_wrap( ~ cluster) +
                    geom_hline(yintercept = 0) +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) 
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(size = 6)
  )

cluster_boxplots

ggsave(file.path(current_analysis_outputs, "cluster_boxplots.png"),
       cluster_boxplots, device = "png")

violin_plots <- cluster_boxplot_data %>%
ggplot(aes(x = Indicator, y = `Indicator value`,
           fill = Indicator,
           color = Indicator)) +
  geom_violin(width=2.1, size=0.2) +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) +
  theme_ipsum() +
  theme(
    legend.position = "none"
  ) + 
  #coord_flip() +
  facet_wrap( ~ cluster, nrow = 3) +
  theme(axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 6, angle = 45, hjust = 1)) +
  geom_hline(yintercept = 0, color = "red")

violin_plots

# Cluster barplots ----

categorical_variables <- names(dplyr::select_if(ecoregions_wide, is.factor))

cluster_barplot_data_wide <- correlation_input_data_all %>%
                            dplyr::select(all_of(c("ecoregion_id", "cluster",
                                                   categorical_variables)))  

cluster_barplot_data <- reshape2::melt(cluster_barplot_data_wide, 
                                       measure.vars = categorical_variables)

cluster_barplot_data <- cluster_barplot_data %>%
                        rename(group = variable,
                               subgroup = value) %>%
                        group_by(cluster, group, subgroup) %>%
                        summarise(ecoregion_count = n_distinct(ecoregion_id)) %>%
                        ungroup()

cluster_barplot_data <- as.data.frame(cluster_barplot_data)

group_cluster_barplot_data <- split(cluster_barplot_data, cluster_barplot_data$group)
length(group_cluster_barplot_data)


barplots <- list()

for ( i in seq_along(group_cluster_barplot_data)) {

barplots[[i]] <- ggplot(group_cluster_barplot_data[[i]]) +
            geom_col(aes(x = cluster, 
                         y = ecoregion_count,
                         fill = subgroup), 
                     position = "fill") +
            facet_wrap(~ group)  +
            theme(strip.text.x = element_text(size = 10),
                  axis.text.x = element_text(size = 8),
                  legend.position = "bottom",
                  legend.text = element_text(size = 8),
                  axis.ticks = element_blank()) +
            xlab("Cluster") + 
            ylab("Ecoregion categories") + 
            labs(color ='Ecoregion categories') +
            scale_fill_viridis(discrete=TRUE) +
            scale_color_viridis(discrete=TRUE) 

}

i <- 0

i <- i +1
barplots[[10]]

png(paste(current_analysis_outputs,"barplots1.png", sep = "/"), 
    units="in", width=9, height=12, res=400)

gridone <- plot_grid(barplots[[1]], barplots[[2]], barplots[[3]],
           align = "v", 
           nrow = 3,
           ncol = 1)

gridone

dev.off()

png(paste(current_analysis_outputs,"barplots2.png", sep = "/"), 
    units="in", width=9, height=12, res=400)

gridtwo <- plot_grid(barplots[[4]], barplots[[5]], barplots[[6]],
                     align = "v", 
                     nrow = 3,
                     ncol = 1)

gridtwo

dev.off()

png(paste(current_analysis_outputs,"barplots3.png", sep = "/"), 
    units="in", width=9, height=12, res=400)

gridthree <- plot_grid(barplots[[7]], barplots[[8]], barplots[[9]],
                     align = "v", 
                     nrow = 3,
                     ncol = 1)

gridthree

dev.off()

png(paste(current_analysis_outputs,"barplots4.png", sep = "/"), 
    units="in", width=9, height=12, res=400)

gridfour <- plot_grid(barplots[[10]], barplots[[11]],
                     align = "v", 
                     nrow = 2,
                     ncol = 1)

gridfour
dev.off()

# Scatterplots ----

# Look at how groups are currently classed
lapply(correlation_input_data, class)

# Ensure ecoregion_id is numeric so it doesn't get converted

correlation_input_data$ecoregion_id <- as.numeric(correlation_input_data$ecoregion_id)

# Convert characters to factors (should only affect grouping vars)

correlation_input_data <- correlation_input_data %>% 
                          mutate(across(where(is.character), as.factor))

# Check it worked correctly
lapply(correlation_input_data, class)

vargroup <- "disturbance.year"

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
                                           y = vary)) +
               geom_point() +
               geom_text(label = correlation_input_data$ecoregion_id) + 
               # geom_smooth(method=lm) +
               labs(x= labx,
                    y = laby) +
               stat_cor(method = "spearman")

ggsave(file.path(current_analysis_outputs, paste(labx, "_x_", laby, ".png", sep = "")),
       scatterplot, device = "png")

scatterplots[[i]] <- scatterplot

}

i <- i + 1

scatterplots[[1]]


# Correlation plots ----



# * SAMPLE SIZE ISSUE ----
#' TODO:  https://www.personality-project.org/r/html/r.test.html - do we need
#' a better way to compare correlations, esp with all different sample size?
#' https://www.youtube.com/watch?v=8mj_DeHtSU4&ab_channel=RajeshDorbala
#' https://stats.stackexchange.com/questions/18887/how-to-calculate-a-confidence-interval-for-spearmans-rank-correlation

# Add clusters to the grouping variables
new_grouping_variables <- c("cluster", grouping_variables)

# Tidy up the names for plotting
all_grouping_variables <- paste("all", new_grouping_variables, sep = ".")
nice_grouping_variables <- str_to_title(gsub(".", " ", new_grouping_variables, fixed=TRUE))
nice_grouping_variables <- str_replace(nice_grouping_variables, "Hfp", "HFP")
nice_grouping_variables <- str_replace(nice_grouping_variables, "Lpi", "LPI")
nice_grouping_variables <- str_replace(nice_grouping_variables, "Rli", "RLI")
nice_grouping_variables <- str_remove(nice_grouping_variables, " Factor")

# * Individual subgroup plots ----

# Loop through the data and calculate correlation coefficients for each subgroup
# of every grouping variable


groups <- list()

group_correlation_plots <- list()

group_directories <- list()

group_coefficients <- list()

group_matrices <- list()

group_dataframes <- list()


for (i in seq_along(new_grouping_variables)) {

# Format the names of the subgroups (eg mangroves, tropical forests) 
# in the grouping variable (eg. Biome)
  
vargroup <- new_grouping_variables[i]
varall <- all_grouping_variables[i]

print(paste( "Producing correlation plots for", vargroup, "grouping variable", sep = " "))

# Create a directory for your outputs

group_directory <- file.path(current_analysis_outputs, paste(vargroup,
                             "correlation matrices", sep = " "))

if ( !dir.exists( group_directory ) ) {
  
  dir.create( group_directory, recursive = TRUE )
  
}

## TEMPORARY
# Remove LPI because it reduces the amount of data we have by about 3/4

indicators <- indicators[!indicators %in% "LPI_2005"]

# Get names of columns we want to subset

x <- c("ecoregion_id", indicators, vargroup)

# Subset correlation inputs to just the indicators, grouping variable and eco ids

group_indicators <- correlation_input_data %>%
                    dplyr::select(all_of(x))

all <- group_indicators

# Split the full dataset by the grouping variable subgroup values (eg mangroves), 
# each into a separate matrix

group_indicator_list <- split(group_indicators, group_indicators[,vargroup])

# Add the full dataset to the list too (so we can compare subgroups to global)

all[, ncol(all)] <- varall

group_indicator_list[[varall]] <- all

subgroups <- names(group_indicator_list)

# Now to format the data for each sub-group, ready to calculate correlations


subgroup_correlation_plots <- list()

subgroup_coefficients <- list()

subgroup_matrices <- list()

subgroup_dataframes <- list()

for (j in seq_along(group_indicator_list)) {
  
   group_dataframe <- group_indicator_list[[j]] %>%
                      dplyr::select(-vargroup) %>%
                      na.omit(.)
  
   group_matrix <- group_indicator_list[[j]] %>%
                   dplyr::select(-ecoregion_id, -vargroup) %>%
                   na.omit(.)
   
   # Remove the years from the column names to make it tidier
   
   names(group_matrix) <- str_remove(names(group_matrix), "_2005")
   
   names(group_matrix) <- str_remove(names(group_matrix), "_2008")
   
   names(group_dataframe) <- c("ecoregion_id",names(group_matrix))

   subgroup_matrices[[j]] <- group_matrix
   subgroup_dataframes[[j]] <- group_dataframe
   
   n <- nrow(group_matrix)
   
   # Create name for plot and output file
   
   name <- paste(names(group_indicator_list)[j], "n =", n, sep = " ")
   name <- str_replace(name, "&", "and")
   name <- str_replace(name, "/", "_and_")
   name <- str_replace(name, "<", " less than ")
   name <- str_replace(name, ">", " more than ")
   name <- gsub("[()]", "", name)
   
   # Check how many rows, don't calculate correlation if too few ecoregions in group
   if (nrow(group_matrix) < 3) {
     
     print(paste("insufficient data for", names(subgroup_matrices)[j], sep = " "))
     
   } else {
     
     # Get the correlation coefficients
     
     rs <- cor(group_matrix, method = "spearman")
     
     subgroup_coefficients[[j]] <- rs
     
     p.mat <- cor_pmat(group_matrix)
     
     # Make a correlation heatmap for the individual sub group
     
     correlation_plot <- ggcorrplot(cor(group_matrix, method = "spearman"),
                                    title = name, 
                                    p.mat = p.mat, 
                                    type = "lower", insig = "blank",
                                    outline.color = "white",
                                    colors = c("#453781FF", "white", "#287D8EFF"),
                                    lab = TRUE)
     
     subgroup_correlation_plots[[j]] <- correlation_plot
     
     ggsave(file.path(group_directory,
                      paste(name, "indicator_correlation_matrix.png", sep = "_")),
            correlation_plot, device = "png")
  
   }

groups[[i]] <- subgroups

group_correlation_plots[[i]] <- subgroup_correlation_plots

group_directories[[i]] <- group_directory

group_coefficients[[i]] <- subgroup_coefficients

group_matrices[[i]] <- subgroup_matrices

group_dataframes[[i]] <- subgroup_dataframes

  }
}

names(group_correlation_plots) <- new_grouping_variables
names(groups) <- new_grouping_variables
names(group_matrices) <- new_grouping_variables
names(group_coefficients) <- new_grouping_variables
names(group_dataframes) <- new_grouping_variables

# Check nested lists are correct lengths
length(group_correlation_plots)
lapply(group_correlation_plots, length)

# * Get coefficient dataframes ----

group_correlation_dataframes <- list()

#subgroup_correlations_dataframes <- list()

for (i in seq_along(group_correlation_plots)) {
  
  # Get the plots for one grouping variable (eg Biomes)
  
  subgroup_plots <- group_correlation_plots[[i]]
  
  names(subgroup_plots) <- groups[[i]]
  
  subgroup_plots <- list.clean(subgroup_plots)
  
  subgroups <- names(subgroup_plots)
  
  subgroup_correlations_dataframes <- list()
  
  for (j in seq_along(subgroup_plots)) {
    
    # Get the plot for one sub-group of the grouping variable (eg Mangroves)
    
    correlation_plot <- subgroup_plots[[j]]
    
    subgroup <- subgroups[j]
    
    # Get the coefficients, and add some additional grouping variables 
    # for the different types of indicator combinations (eg independent or related inputs)
    
    correlation_df1 <- correlation_plot$data[,1:5] %>%
      merge(indicator_properties[c("indicator","inputs")], 
            by.x = "Var1", by.y = "indicator") %>%
      dplyr::rename(input1 = inputs) %>%
      merge(indicator_properties[c("indicator","inputs")], 
            by.x = "Var2", by.y = "indicator") %>%
      dplyr::rename(input2 = inputs) %>%
      mutate(subgroup = subgroup,
             combination = paste(Var1, "x", Var2, sep = " "),
             inputs = paste(input1, "x", input2, sep = " ")) %>%
      mutate(inputs = ifelse(inputs == "land use data x iucn red list",
                             "iucn red list x land use data", inputs)) %>%
      dplyr::rename(coefficient = value) %>%
      dplyr::select(-Var2) %>%
      merge(indicator_relationships, by = "combination") %>%
      dplyr::rename(ind_group = Var1) 
    
    correlation_df2 <- correlation_plot$data[,1:5] %>%
      merge(indicator_properties[c("indicator","inputs")], 
            by.x = "Var1", by.y = "indicator") %>%
      dplyr::rename(input1 = inputs) %>%
      merge(indicator_properties[c("indicator","inputs")], 
            by.x = "Var2", by.y = "indicator") %>%
      dplyr::rename(input2 = inputs) %>%
      mutate(subgroup = subgroup,
             combination = paste(Var2, "x", Var1, sep = " "),
             inputs = paste(input1, "x", input2, sep = " ")) %>%
      mutate(inputs = ifelse(inputs == "land use data x iucn red list",
                             "iucn red list x land use data", inputs)) %>%
      dplyr::rename(coefficient = value) %>%
      dplyr::select(-Var1) %>%
      merge(indicator_relationships, by = "combination") %>%
      dplyr::rename(ind_group = Var2)
    
    correlation_df <- rbind(correlation_df1, correlation_df2)
    
    subgroup_correlations_dataframes[[j]] <- correlation_df
    
    all_combos <- correlation_df$combination
    
  }
  
  # Combine all possible subgroups and correlation coefficients for the group
  # into one dataframe
  
  subgroup_full_dataframe <- do.call(rbind, subgroup_correlations_dataframes)
  
  #group_correlation_dataframes[[i]] <- subgroup_correlations_dataframes
  group_correlation_dataframes[[i]] <- subgroup_full_dataframe
  
}

length(group_correlation_dataframes)
lapply(group_correlation_dataframes, length)

names(group_correlation_dataframes) <- new_grouping_variables

# * Get confidence intervals ----

group_confidence_intervals <- list()

for (i in seq_along(group_matrices)) {
  
  group_matrix <- group_matrices[[i]] 
  
  subgroups <- groups[[i]]
  
  group_names <- c("cluster", new_grouping_variables)
  
  subgroup_confidence_intervals <- list()
  
  for (j in seq_along(group_matrix)) {
    
    subgroup_matrix <- group_matrix[[j]]
    
    indicator_combos <- combn(names(subgroup_matrix), 2)
    
    # Now loop through the combos
    
    out <- list()
    
    for (k in seq_along(1:ncol(indicator_combos))) {
      
      combo <- as.vector(indicator_combos[,k])
      combo
      
      col1 <- subgroup_matrix[,combo[1]]
      col2 <- subgroup_matrix[,combo[2]]
      
      n <- nrow(subgroup_matrix)
      
      df <- spearman_CI(col1, col2, alpha = 0.05) # single indicator pair coefficient
      
      df <- cbind(group_names[i],subgroups[j], combo[1], combo[2], df, n)
      
      names(df) <- c("group", "subgroup", "ind_1", "ind_2", "rs", "lower_ci",
                     "upper_ci", "n")
      
      df <- df %>%
            mutate(pair = paste(combo[1], combo[2], sep = " x "))
      
      out[[k]] <- df
      
    }
    
    l2df <- do.call(rbind, out) # list of all subgroup pair coefficients
    
    subgroup_confidence_intervals[[j]] <- l2df
  }
  
  l1df <- do.call(rbind, subgroup_confidence_intervals)

  group_confidence_intervals[[i]] <- subgroup_confidence_intervals
  
}

# Turn into one dataframe per grouping variable

confidence_intervals <- list()

for (i in seq_along(group_confidence_intervals)) {
  
  confidence_intervals[[i]] <- do.call(rbind, group_confidence_intervals[[i]])
  
}

# Remove non-significant coefficients ----

x <- confidence_intervals[[1]] %>%
     mutate(id = paste(pair, subgroup, sep = "_"))

head(x)
class(x)

pval <- group_correlation_dataframes[[1]] %>%
        select(combination, coefficient, pvalue, signif, subgroup) %>%
        distinct(.) %>%
        mutate(id = paste(combination, subgroup, sep = "_"))
head(pval)
class(pval)

test <- x %>%
        merge(pval[c("id","coefficient", "pvalue", "signif")], by = "id") %>%
        dplyr::select(rs, coefficient, pvalue,signif, everything())

head(test)

# * Make LU * BD caterpillar plots ----

caterpillar_plots <- list()

for (i in seq_along(confidence_intervals)) {
  
y <- confidence_intervals[[i]]


names(y) <- c("grouping var", "subgroup", "ind1", "ind2", "rs", "lower_ci",
              "upper_ci", "n", "pair")

y <- y %>%
    mutate(subgroup_label = paste(subgroup, "n =", n, sep = " ")) %>%
    filter(pair == "BHI_plants x extinct"|
           pair == "BHI_plants x RLI"|
           pair == "BHI_plants x threatened"|
           pair == "BIIab x extinct"|
           pair == "BIIab x RLI"|
           pair == "BIIab x threatened"|
           pair == "BIIri x extinct"|
           pair == "BIIri x RLI"|
           pair == "BIIri x threatened"|
           pair == "HFP x extinct"|
           pair == "HFP x RLI"|
           pair == "HFP x threatened")


y <- y %>%
     mutate(id = paste(pair, subgroup, sep = "_"))

pval <- group_correlation_dataframes[[i]] %>%
        select(combination, coefficient, pvalue, signif, subgroup) %>%
        distinct(.) %>%
        mutate(id = paste(combination, subgroup, sep = "_"))


plotdata <- y %>%
  merge(pval[c("id","coefficient", "pvalue", "signif")], by = "id") %>%
  dplyr::select(rs, coefficient, pvalue,signif, everything()) %>%
  filter(signif != 0)

catplot <-  ggplot(plotdata, aes(rs, pair)) +
            geom_point(aes(col = pair)) +
            geom_linerange(aes(xmin = lower_ci, xmax = upper_ci,
                               col = pair)) +
            facet_wrap(~ subgroup_label) +
            geom_vline(xintercept = 0, col = "red") +
            ggtitle(new_grouping_variables[[i]]) +
            theme(axis.text.y = element_blank(),
                  strip.text.x = element_text(size = 5),
                  axis.text.x = element_text(size = 5),
                  legend.position = "none",
                  legend.text = element_text(size = 5),
                  axis.ticks = element_blank()) +
            xlab("Spearman's rank correlation coefficient") + 
            ylab("Ecoregion categories") + 
            labs(color ='Ecoregion categories') +
            scale_fill_viridis(discrete=TRUE) +
            scale_color_viridis(discrete=TRUE) + 
            geom_text(aes(label = pair), size= 2, vjust = 1.5)

ggsave(file.path(group_directories[[i]],
                 paste(new_grouping_variables[[i]],
                       "caterpillar_plot.png", sep = "_")),
       catplot, device = "png")

caterpillar_plots[[i]] <- catplot

}

i <- i +1
caterpillar_plots[[i]]

names(caterpillar_plots) <- new_grouping_variables

# * Make related indicator caterpillar plots ----

related_caterpillar_plots <- list()

for (i in seq_along(confidence_intervals)) {
  
  y <- confidence_intervals[[i]]
  
  names(y) <- c("grouping var", "subgroup", "ind1", "ind2", "rs", "lower_ci",
                "upper_ci", "n", "pair")
  
  y <- y %>%
    mutate(subgroup_label = paste(subgroup, "n =", n, sep = " ")) %>%
    filter(pair == "BHI_plants x BIIab"|
             pair == "BHI_plants x BIIri"|
             pair == "BHI_plants x HFP"|
             pair == "BIIab x BIIri"|
             pair == "BIIab x HFP"|
             pair == "BIIri x HFP"|
             pair == "extinct x RLI"|
             pair == "extinct x threatened"|
             pair == "RLI x threatened")
  
  y <- y %>%
    mutate(id = paste(pair, subgroup, sep = "_"))
  
  pval <- group_correlation_dataframes[[i]] %>%
    select(combination, coefficient, pvalue, signif, subgroup) %>%
    distinct(.) %>%
    mutate(id = paste(combination, subgroup, sep = "_"))
  
  
  plotdata <- y %>%
    merge(pval[c("id","coefficient", "pvalue", "signif")], by = "id") %>%
    dplyr::select(rs, coefficient, pvalue,signif, everything()) %>%
    filter(signif != 0)
  
  catplot <-  ggplot(plotdata, aes(rs, pair)) +
    geom_point(aes(col = pair)) +
    geom_linerange(aes(xmin = lower_ci, xmax = upper_ci,
                       col = pair)) +
    facet_wrap(~ subgroup_label) +
    geom_vline(xintercept = 0, col = "red") +
    ggtitle(new_grouping_variables[[i]]) +
    theme(axis.text.y = element_blank(),
          strip.text.x = element_text(size = 5),
          axis.text.x = element_text(size = 5),
          legend.position = "none",
          legend.text = element_text(size = 5),
          axis.ticks = element_blank()) +
    xlab("Spearman's rank correlation coefficient") + 
    ylab("Ecoregion categories") + 
    labs(color ='Ecoregion categories') +
    scale_fill_viridis(discrete=TRUE) +
    scale_color_viridis(discrete=TRUE) + 
    geom_text(aes(label = pair), size= 2, vjust = 1.5)
  
  ggsave(file.path(group_directories[[i]],
                   paste(new_grouping_variables[[i]],
                         "related_caterpillar_plot.png", sep = "_")),
         catplot, device = "png")
  
  related_caterpillar_plots[[i]] <- catplot
  
}

i <- i + 1
related_caterpillar_plots[[i]]

names(related_caterpillar_plots) <- new_grouping_variables





# * Make subgroup scatterplots ----

# Make directories for the scatterplots (this can probably go in the make plots
# loop itself now have resolved issues with naming)

scatter_directories <- list()

for (i in 1:length(names(group_matrices))) {

group <- names(group_matrices)[i]

print(group)

scatter_directory <- file.path(current_analysis_outputs, paste(group,
                      "scatterplots", sep = "_"))

if ( !dir.exists( scatter_directory ) ) {
  
  dir.create( scatter_directory, recursive = FALSE )
  
}

scatter_directories[[i]] <- scatter_directory

}

# For each grouping variable, for each subgroup, for each indicator pair,
# make a scatterplot

group_scatterplots <- list()

for (i in seq_along(group_dataframes)){
  
  subgroup_dataframes <- group_dataframes[[i]]
  
  scatter_directory <- scatter_directories[[i]]
  
  subgroup_scatterplots <- list()
  
  for (j in seq_along(subgroup_dataframes)) {
    
  subgroup <- groups[[i]][j]
  
  subgroup
  
  subgroup <- str_replace(subgroup, " - ", "_")
  subgroup <- str_replace(subgroup, " & ", "and")
  subgroup <- str_replace(subgroup, " , ", "_")
  subgroup <- str_replace(subgroup, " . ", "_")
  subgroup <- str_replace(subgroup, " ", "_")
  
  subgroup
  
  subgroup_scatterplots[[j]] <- make_subgroup_scatterplot(subgroup_dataframes[[j]], 
                                                     names(group_dataframes)[i], 
                                                     subgroup,
                                                     scatter_directory)
    
  }

  group_scatterplots[[i]] <- subgroup_scatterplots

}



# * Grouping variable heatmaps ----

# Make a heatmap that includes correlations for all subgroups in one plot 
# (eg all biome subgroups (magroves, tropical forests etc) in one plot)

for (i in seq_along(group_correlation_dataframes)) {

group_df <- group_correlation_dataframes[[i]]
group_directory <- group_directories[[i]]

# Split the independent and related combinations

independent <- group_df %>%
               dplyr::filter(input_relationship == "independent") %>%
               dplyr::filter(ind_group == "BHI_plants"|
                      ind_group == "BIIab"|
                      ind_group == "BIIri"|
                      ind_group == "HFP")

independent_heatmap <- make_heatmap(independent)

# Save in the group directory

ggsave(file.path(group_directory,
                 paste(new_grouping_variables[i], 
                       "independent_correlation_heatmap.png", sep = "_")),
       independent_heatmap, device = "png")

all_independent_heatmaps[[i]] <- independent_heatmap


related <- group_df %>%
  filter(input_relationship == "related")

related_heatmap <- make_heatmap(related)
related_heatmap

ggsave(file.path(group_directory,
                 paste(new_grouping_variables[i], "related_correlation_heatmap.png", 
                       sep = "_")),
       related_heatmap, device = "png")

all_related_heatmaps[[i]] <- related_heatmap

}

names(correlation_plots) <- new_grouping_variables
names(all_independent_heatmaps) <- new_grouping_variables
names(all_related_heatmaps) <- new_grouping_variables

i <- 1
all_independent_heatmaps[[i]]

i <- i + 1
all_independent_heatmaps[[i]]

i <- 1
all_related_heatmaps[[i]]

i <- i + 1
all_related_heatmaps[[i]]

# Make an ecoregion map

eco_pal <- colorNumeric("PuBu", domain = cluster_map_data$ecoregion_id)
hfp_pal <- colorNumeric("PuBu", domain = indicator_map_input_data$HFP_2005)

cluster_map_data %>% 
  leaflet() %>% 
  addTiles() %>% 
  addPolygons(weight = 1, 
              color = ~ eco_pal(ecoregion_id),
              fillOpacity = 0.8,
              group = "ecoregions",
              # add labels that display indicator value and ecoregion id
              label = ~paste(ecoregion_id, ECO_NAME, 
                             sep = " "),
              # highlight polygons on hover
              highlightOptions = highlightOptions(weight = 5, color = "white",
                                                  bringToFront = TRUE)) %>%

tmap_save(beta, file.path(indicator_outputs, paste(location,
                                                   "proportion_high_beta_map.png", sep = "_")),
          width=1920, height=1080, asp=0)



# * Plot PCA 

# https://rstudio-pubs-static.s3.amazonaws.com/323416_ab58ad22d9e64ba2831569cf3d14a609.html

# Make the PCA plots in a loop

# Remove indicators
#' 
#' grouping_variables <- names(pca_data_5)[!(names(pca_data_5) %in% indicators)]
#' 
#' # Remove numerics
#' grouping_variables <- grouping_variables[!(grouping_variables %in% c("ecoregion_id",
#'                                                                      "predominant.threat.type",
#'                                                                      "headline.threat.type",
#'                                                                      # "ecoregion.area.km.sq",
#'                                                                      # "RLI_records",
#'                                                                      # "LPI_records",
#'                                                                      #"mean.scientific.publications",
#'                                                                      "scenario",
#'                                                                      "predominant.threat.count"))]
#' 
#' pca_plots <- list()
#' 
#' for (i in seq_along(grouping_variables)) {
#'   
#'   vargroup <- pca_data_5[,grouping_variables[i]]
#'   varname <- grouping_variables[i]
#'   
#'   pca_plot <- fviz_pca_biplot(pca, geom.ind = "point", pointshape = 21, 
#'                        pointsize = 2, 
#'                        col.ind = "black", 
#'                        palette = "viridis", 
#'                        addEllipses = TRUE,
#'                        label = "var",
#'                        col.var = "black",
#'                        repel = TRUE,
#'                        fill.ind = as.factor(vargroup),
#'                        legend.title = varname) +
#'     theme(plot.title = element_blank(),
#'           legend.position = "bottom") +
#'     theme(legend.text=element_text(size = 7)) +
#'     guides(fill=guide_legend(nrow=6,byrow=TRUE))
#'   
#'   ggsave(file.path(current_analysis_outputs, paste(varname, ".png", sep = "")),
#'          pca_plot, device = "png")
#' 
#'   pca_plots[[i]] <- pca_plot
#'   
#' }
#' 
#' pca_plots[[1]]
#' pca_plots[[2]]
#' pca_plots[[3]]
#' pca_plots[[4]]
#' pca_plots[[5]]
#' pca_plots[[6]]
#' pca_plots[[7]]
#' pca_plots[[8]]
#' pca_plots[[9]]
#' pca_plots[[10]]
#' 
#' 
#' 
#' # Get the response variable  (or what you think it is)
#' 
#' #vargroup <- as.numeric(pca_data_5$island.status == "island")
#' 
#' vargroup <- as.numeric(as.factor(pca_data_5$island.status))
#' 
#' #vargroup <- pca_data_5$scenario.numeric
#' 
#' #vargroup <- as.numeric(as.factor(vargroup))
#' 
#' number_clusters <- length(unique(vargroup))
#' 
#' # Plot the different dimensions to see how clearly the response variable is delineated 
#' 
#' plot(pca$x[, c(1, 2)], col = as.numeric(as.factor(pca_data_5$island.status)), 
#'      xlab = "PC1", ylab = "PC2")
#' 
#' plot(pca$x[, c(1, 3)], col = as.numeric(as.factor(pca_data_5$island.status)), 
#'      xlab = "PC1", ylab = "PC2")
#' 
#' # Plot explained variability
#' 
#' par(mfrow = c(1, 2))
#' 
#' # Calculate variability of each component
#' pr.var <- pca$sdev ^2
#' 
#' # Variance explained by each principal component: pve
#' pve <- pr.var/sum(pr.var)
#' 
#' # Plot variance explained for each principal component
#' plot(pve, xlab = "Principal Component", 
#'      ylab = "Proportion of Variance Explained", 
#'      ylim = c(0, 1), type = "b")
#' 
#' # Plot cumulative proportion of variance explained
#' plot(cumsum(pve), xlab = "Principal Component", 
#'      ylab = "Cumulative Proportion of Variance Explained", 
#'      ylim = c(0, 1), type = "b")
#' 
#' 
#' # Try a model ----
#' 
#' # TODO: use tidymodels ----
#' 
#' library(Rmisc)
#' library(car)
#' library(lm.beta)
#' library(MASS)
#' 
#' # We want the non-inverted values so need to rebuild the input data
#' 
#' response <- raw_indicators_wide[,c("ecoregion_id", "threatened_2005")]
#' 
#' model_input_data <- raw_indicators_wide[,c("ecoregion_id",
#'                                                 "BHI_plants 2005",
#'                                                 "BIIab_2005",
#'                                                 "BIIri_2005",
#'                                                 "extinct_2005",
#'                                                 "HFP_2005",
#'                                                 "RLI_2005")] 
#' 
#' model_input_data <- model_input_data %>%
#'                     dplyr::filter(HFP_2005 < quantile(HFP_2005, 
#'                                                0.99, na.rm = TRUE)) %>%
#'                     mutate(ecoregion_id = as.character(ecoregion_id))
#' 
#' model_input_data <- model_input_data %>%
#'                     merge(ecoregions_wide, by = "ecoregion_id")
#' 
#' scaling_index <- unlist(lapply(model_input_data, is.numeric))
#' 
#' eco_id <- model_input_data$ecoregion_id
#' non_numeric <- model_input_data[,!scaling_index]
#' 
#' model_input_data <- as.data.frame(scale(model_input_data[,scaling_index], 
#'                                         scale = TRUE, center = TRUE)) 
#' 
#' model_input_data <- cbind(eco_id, model_input_data, non_numeric)
#' 
#' model_input_data <- model_input_data %>%
#'                     dplyr::select(-eco_id) %>%
#'                     dplyr::select(ecoregion_id, everything()) %>%
#'                     mutate(ecoregion_id = as.numeric(ecoregion_id))
#' 
#' model_input_data <- response %>%
#'                     merge(model_input_data, by = "ecoregion_id")
#' 
#' 
#' model_input_data <- model_input_data %>%
#'                     filter(RLI_records > 600)
#' 
#' scatt <- scatterplot(threatened_2005 ~ BIIri_2005 , data = model_input_data)
#' 
#' test_glm <- glm(BIIri_2005 ~ 1, family = gaussian, 
#'                  data = model_input_data)
#' 
#' 
#' step_glm <- stepAIC(test_glm, scope = list(upper = ~ HFP_2005 + RLI_records + 
#'                                            threatened_2005 +
#'                                            RLI_2005 +
#'                                            extinct_2005 +
#'                                            scenario +
#'                                            included.in.HFP + realm +
#'                                            LPI_records + headline.threat.type +
#'                                            Biome + mean.scientific.publications +
#'                                            island.status +
#'                                            number.of.endemics +
#'                                            mean.human.population.density,
#'                                            lower = ~ 1))
#' 
#' summary(step_glm)
#' 
#' # Create an interactive map ----
#' 
#' indicator_map_data_all <- left_join(ecoregion_map_renamed, 
#'                                     correlation_input_data[c("ecoregion_id", "realm",
#'                                                              indicators)],
#'                                     by = "ecoregion_id") 
#' 
#' ecoregion_subset <- ecoregion_countries %>%
#'   filter(CNTRY_NAME == "Australia") %>%
#'   unique(.)
#' 
#' 
#' indicator_map_input_data <- indicator_map_data_all[indicator_map_data_all$ecoregion_id %in% 
#'                                    ecoregion_subset$ECO_ID,] 
#' 
#' # Create a colour palette for proportion of extinctions
#' 
#' indicator_map_input_data <- indicator_map_input_data_sf
#' 
#' ext_pal <- colorNumeric("PuBu", domain = indicator_map_input_data$extinct_2005)
#' hfp_pal <- colorNumeric("PuBu", domain = indicator_map_input_data$HFP_2005)
#' 
#' #' TODO: Add real, not scaled indicator values
#' #' TODO: Add ecoregion labels?
#' #' TODO: Try adding all the indicator values into one label
#' 
#' indicator_map_input_data %>% 
#'   leaflet() %>% 
#'   addTiles() %>% 
#'   addPolygons(weight = 1, color = ~ext_pal(extinct_2005),
#'               fillOpacity = 0.8, group = "Proportion of species extinct",
#'               # add labels that display indicator value and ecoregion id
#'               label = ~paste(ecoregion_id, "proportion extinct =", extinct_2005, 
#'                               sep = " "),
#'               # highlight polygons on hover
#'               highlightOptions = highlightOptions(weight = 5, color = "white",
#'                                                   bringToFront = TRUE)) %>%
#'   addPolygons(weight = 1, color = ~hfp_pal(HFP_2005),
#'               fillOpacity = 0.8, group = "Human footprint index",
#'               # add labels that display indicator value and ecoregion id
#'               label = ~paste(ecoregion_id, "HFP =", HFP_2005, sep = " "),
#'               # highlight polygons on hover
#'               highlightOptions = highlightOptions(weight = 5, color = "white",
#'                                                   bringToFront = TRUE)) %>%
#'   addLayersControl(baseGroups = c("OSM", "Carto", "Esri"), 
#'                    overlayGroups = c("Proportion of species extinct", 
#'                                      "Human footprint index"))
#' 
#' 
#' # Map indicators ----
#' 
#' 
#' # indicator_map_data <- left_join(ecoregion_map_renamed,
#' #                                  data)
#' # 
#' # indicator_map_data <- indicator_map_data_all %>% filter(realm == "Australasia")
#' 
#' #indicator_map_data_all_all <- indicator_map_data_all
#' 
#' for (i in seq_along(indicators)) {
#' 
#' legend_title <- indicators[i]
#' 
#' indicator_map_data <- indicator_map_data_all[,c("ecoregion_id",indicators[i])]
#' 
#' 
#' if(legend_title == "extinct_2008") {
#'   
#' 
#' indicator_map_data <- indicator_map_data %>%
#'                       mutate(extinct_2008 = ifelse(extinct_2008 < -1,
#'                                                  -1, extinct_2008))
#' 
#' print(paste(legend_title, "truncated fill values for visualisation", sep = " "))
#' 
#' }
#' 
#' names(indicator_map_data) <- c("ecoregion_id", "indicator","geometry")
#' 
#' if(legend_title == "threatened_2008" | legend_title == "RLI_2008") {
#'   
#'   
#'   indicator_map_data <- indicator_map_data %>%
#'     mutate(indicator = ifelse(indicator < -1,
#'                                  -1, indicator))
#'   
#'   print(paste(legend_title, "truncated fill values for visualisation", sep = " "))
#'         
#' }
#' 
#' 
#' 
#' map <-  ggplot(indicator_map_data) +
#'         geom_sf(aes(fill = indicator), colour = "black", 
#'                 size = 0.05, show.legend = 'fill') +
#'         scale_fill_viridis_c(alpha = .8,
#'                              na.value = "grey70") +
#'         theme(axis.line = element_line(),
#'               panel.grid.major = element_blank(),
#'               panel.grid.minor = element_blank(),
#'               panel.background = element_blank()) +
#'         labs(fill = legend_title) +
#'         theme(legend.position = "right")
#' 
#' ggsave(file.path(current_analysis_outputs,
#'                  paste(location, eco_version, legend_title,
#'                        "map.png",
#'                        sep = "_")), 
#'        map,  device = "png")
#' }
#' 
#' 
#' # Map scenarios ----
#' 
#' # * Global ----
#' 
#' if(load_map == TRUE) {
#'   
#'   ecoregion_values_scenario <- ecoregions_wide %>%
#'     dplyr::select(ecoregion_id, scenario)
#'   
#'   ecoregion_values_map_data <- left_join(ecoregion_map_renamed, 
#'                                          ecoregion_values_scenario,
#'                                          by = "ecoregion_id")
#'   
#'   global_scenario_map <-  ggplot(ecoregion_values_map_data) +
#'     geom_sf(aes(fill = scenario), colour = "black", 
#'             size = 0.05, show.legend = 'fill') +
#'     scale_fill_viridis_d(alpha = .8,
#'                          na.value = "grey70") +
#'     theme(axis.line = element_line(),
#'           panel.grid.major = element_blank(),
#'           panel.grid.minor = element_blank(),
#'           panel.background = element_blank()) +
#'     labs(fill = "test") +
#'     theme(legend.position = "right")
#'   
#'   global_scenario_map
#'   
#'   ggsave(file.path(current_analysis_outputs,
#'                    paste(location, eco_version,
#'                          "data_threat_scenarios.png",
#'                          sep = "_")), 
#'          scenario_map,  device = "png")
#'   
#'   # * Realms ----
#'   
#'   ecoregion_values_scenario_oceania <- ecoregions_wide %>%
#'     dplyr::select(ecoregion_id, scenario, realm) 
#'   
#'   ecoregion_values_map_data <- left_join(ecoregion_map_renamed, 
#'                                          ecoregion_values_scenario_oceania,
#'                                          by = "ecoregion_id") %>%
#'     filter(realm == "Oceania")
#'   
#'   scenario_map <-  ggplot(ecoregion_values_map_data) +
#'     geom_sf(aes(fill = scenario), colour = "black", 
#'             size = 0.05, show.legend = 'fill') +
#'     scale_fill_viridis_d(alpha = .8,
#'                          na.value = "grey70") +
#'     theme(axis.line = element_line(),
#'           panel.grid.major = element_blank(),
#'           panel.grid.minor = element_blank(),
#'           panel.background = element_blank()) +
#'     labs(fill = "test") +
#'     theme(legend.position = "right") +
#'     scale_x_continuous(limits = c(110, 300)) + 
#'     scale_y_continuous(limits = c(-50, 70)) 
#'   
#'   scenario_map
#'   
#'   ggsave(file.path(current_analysis_outputs,
#'                    paste(location, eco_version,
#'                          "data_threat_scenarios.png",
#'                          sep = "_")), 
#'          scenario_map,  device = "png")
#'   
#' }
#' 
#' # K-means cluster analysis ----
#' 
#' rownames(pca_input_data) <- pca_input_data$ecoregion_id
#' 
#' kmeans_input_data <- pca_input_data %>%
#'   dplyr::select(- ecoregion_id)
#' 
#' # Remove incomplete cases
#' 
#' kmeans_input_data <- kmeans_input_data[complete.cases(kmeans_input_data),]
#' 
#' # Scale
#' 
#' kmeans_input_data <- scale(kmeans_input_data)
#' 
#' # Find best/truest number of clusters
#' 
#' wss <- 0
#' 
#' set.seed(1)
#' 
#' # Look over 1 to 15 possible clusters
#' for (i in 1:15) {
#'   # Fit the model: km.out
#'   km.out <- kmeans(kmeans_input_data, centers = i, nstart = 20, iter.max = 200)
#'   # Save the within cluster sum of squares
#'   wss[i] <- km.out$tot.withinss
#' }
#' 
#' # Produce a scree plot
#' plot(1:15, wss, type = "b",
#'      xlab = "Number of Clusters",
#'      ylab = "Within groups sum of squares")
#' 
#' # Looks like 3 clusters is best
#' 
#' number_clusters <- 3
#' 
#' wisc.km <- kmeans(scale(kmeans_input_data), centers = number_clusters, nstart = 100)
#' 
#' # Compare k-means to actual groupings
#' 
#' x <- table(wisc.km$cluster, pca_data_5$lpi.records.factor)
#' 
#' true <- colSums(x)
#' 
#' model <- rowSums(x)
#' 
#' check_model <- cbind(true,model)
#' 
#' check_model
#' 
#' # View the resulting model
#' wisc.km$tot.withinss
#' 
#' # Plot of BII and RLI by cluster membership
#' 
#' plot(kmeans_input_data[, c("BIIri_2005", "threatened_2005")],
#'      col = wisc.km$cluster,
#'      #pch = as.character(pca_data_5$Biome),
#'      main = paste("k-means clustering of indicator data with", number_clusters, "clusters"),
#'      xlab = "bii", ylab = "rli")
#' 
#' par(mfrow = c(2, 3))


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

# 3D plotting code ----

# Plot3D with plotly
#manual_palette <- c("#551A8B", "#8968CD", "#AB82FF", "#FF8C00", "#CD6600", "#8B4500")

# df3D <- data.frame(comp1=pl.pca$scores[,1],
#                    comp2=pl.pca$scores[,2],
#                    comp3=pl.pca$scores[,3])
# fig <- plot_ly(df3D, x = ~comp1, y = ~comp2, z = ~comp3, color = ~comp3,  mode = 'lines+markers',
#                # Hover text:
#                text = ~rownames(pl.pca$scores))
# fig <- fig %>% add_markers()
# fig <- fig %>% add_text(textposition = "top right")
# fig <- fig %>% layout(scene = list(xaxis = list(title = pcxlab),
#                                    yaxis = list(title = pcylab),
#                                    zaxis = list(title = pczlab)),
#                       annotations = list(
#                         x = 1.13,
#                         y = 1.05,
#                         text = 'PC3 Score',
#                         showarrow = FALSE
#                       ))
# fig

# Old analysis code ----

#PCA code

pca <- prcomp(pca_data_5[,c(22:28)], center = TRUE, scale = TRUE)
summary(pca)
print(pca)

pca_loadings <- as.data.frame(print(pca$rotation))

write.csv(pca_loadings, file.path(current_analysis_outputs, "pca_loadings.csv"))

var <- get_pca_var(pca)
variable_contributions <- var$contrib

write.csv(variable_contributions, file.path(current_analysis_outputs,
                                            "pca_variable_contributions.csv"))

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
