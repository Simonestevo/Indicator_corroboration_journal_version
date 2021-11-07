
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

# PCA, Clustering, Models

library(corrplot)
library(ape)
library(MASS)
library(ggpubr)
library(betareg)
library(performance)

# Data handling and table reshaping
library(tidyverse)
library(tidylog)
library(reshape2)
library(devtools)
library(data.table)
library(rlist)
library(e1071)
library(psych)
library(arules)
library(ppsr)
library(broom)

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
# library(forcats)

# Maps
library(sf)
library(leaflet)


# Set input and output locations ----
analysis_date <- "2021-05-15"
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

# Folder for final manuscript figures

manuscript_outputs <- "N:\\Quantitative-Ecology\\Simone\\extinction_test\\outputs\\version_3\\2021-05-04_manuscript_outputs"

dir.create(manuscript_outputs, recursive = TRUE ) # create a new directory for today's outputs

dirs <- list.dirs(file.path(manuscript_outputs, analysis_date), full.names = TRUE)

main_outputs <- dirs[str_detect(dirs, "main")]

# Folder for supporting info figures

supp_outputs <- dirs[str_detect(dirs, "supporting")]

# Folder with indicator data

data_outputs <- dirs[str_detect(dirs, "data")]

current_analysis_outputs <- file.path(analysis_outputs,paste(date,
                                                             "_analysis_output_files",sep=""))

dir.create(current_analysis_outputs)

files_paths <- list.files(data_outputs, full.names = TRUE)

indicators_raw_lpi <- readRDS(files_paths[str_detect(files_paths, "1_indicators_cleaned_lpi.rds")])
indicators_raw_no_lpi <- readRDS(files_paths[str_detect(files_paths, "2_indicators_cleaned_no_lpi.rds")])
indicators_std_no_lpi <- readRDS(files_paths[str_detect(files_paths, "3_indicators_cleaned_standardised_no_lpi.rds")])
indicators_std_lpi <- readRDS(files_paths[str_detect(files_paths, "4_indicators_cleaned_standardised_lpi.rds")])
ecoregions_numeric_raw <- readRDS(files_paths[str_detect(files_paths, "5_ecoregions_numeric.rds")])
ecoregions_numeric_std <- readRDS(files_paths[str_detect(files_paths, "6_ecoregions_numeric_standardized.rds")])
ecoregions_categorical <- readRDS(files_paths[str_detect(files_paths, "7_ecoregions_categorical.rds")])
ecoregion_clusters <- read.csv(files_paths[str_detect(files_paths, 
                                                      "8_ecoregion_clusters.csv")])
# Get the names of factors

grouping_variables <- names(ecoregions_categorical)
grouping_variables <- grouping_variables[!str_detect(grouping_variables, "ecoregion_id")]

# Get names of numeric variables 

numeric_variables <- names(ecoregions_numeric_raw)
numeric_variables <- numeric_variables[!str_detect(numeric_variables, "ecoregion_id")]

# Make correlation matrix of numeric predictors

ecoregions_numeric_raw <- ecoregions_numeric_raw[complete.cases(
                          ecoregions_numeric_raw),]

numeric_cor_inputs <- as.matrix(ecoregions_numeric_raw %>% 
                      dplyr::select(all_of(numeric_variables)))

numeric_cor <- cor(numeric_cor_inputs, method = "spearman")

View(numeric_cor)

# Prepare model data ----

# https://www.zeileis.org/papers/ERCIM-2010.pdf

model_input_data <- indicators_raw_no_lpi %>%
  merge(ecoregion_clusters[c("ecoregion_id",
                             "cluster")],
        by = "ecoregion_id") %>% 
  merge(ecoregions_categorical,
        by = "ecoregion_id") %>%
  merge(ecoregions_numeric_raw,
        by = "ecoregion_id") %>%
  dplyr::select(-included_in_hfp)

head(model_input_data)

write.csv(model_input_data, "9_model_input_data.csv")


model_input_data$BHI_plants_2005 <- scale(model_input_data$BHI_plants_2005,
                                          center = TRUE,
                                          scale = TRUE)

model_input_data$ecoregion_area_km_sq <- scale(model_input_data$ecoregion_area_km_sq,
                                               center = TRUE,
                                               scale = TRUE)

model_input_data$rli_records <- scale(model_input_data$rli_records,
                                      center = TRUE,
                                      scale = TRUE)

model_input_data$scaled_rli_records <- scale(model_input_data$scaled_rli_records,
                                             center = TRUE,
                                             scale = TRUE)

model_input_data$BIIab_2005 <- scale(model_input_data$BIIab_2005,
                                     center = TRUE,
                                     scale = TRUE)

model_input_data$scaled_threat_count <- scale(model_input_data$scaled_threat_count,
                                              center = TRUE,
                                              scale = TRUE)

boxplot(RLI_2008 ~ area_factor, data = model_input_data) 


# Beta regression ----
# https://stats.stackexchange.com/questions/233366/how-to-fit-a-mixed-model-with-response-variable-between-0-and-1

test <- glmmadmb(RLI_2008 ~ BHI_plants_2005 + (1|ecoregion_area_km_sq), 
         model_input_data, family="beta")

bhi_mod <- betareg(RLI_2008 ~ BHI_plants_2005, data = model_input_data)
summary(bhi_mod)
qqnorm(resid(bhi_mod))
qqline(resid(bhi_mod))

bhi_rli_mod <- betareg(RLI_2008 ~ BHI_plants_2005 * rli_records, data = model_input_data)
summary(bhi_rli_mod)
qqnorm(resid(bhi_rli_mod))
qqline(resid(bhi_rli_mod))

bii_mod <- betareg(RLI_2008 ~ BIIab_2005 | scaled_rli_records, data = model_input_data)
summary(bii_mod)
qqnorm(resid(bii_mod))
qqline(resid(bii_mod))

bhi_clus_mod <- betareg(RLI_2008 ~ BHI_plants_2005 + cluster, data = model_input_data)
summary(bhi_clus_mod)
qqnorm(resid(bii_mod))
qqline(resid(bii_mod))

bii_clus_mod <- betareg(RLI_2008 ~ BIIab_2005 + cluster, data = model_input_data)
summary(bii_clus_mod)
qqnorm(resid(bii_clus_mod))
qqline(resid(bii_clus_mod))

bii_threat_mod <- betareg(RLI_2008 ~ BIIab_2005 + scaled_threat_count, data = model_input_data)
summary(bii_threat_mod)
qqnorm(resid(bii_threat_mod))
qqline(resid(bii_threat_mod))

bhi_threat_mod <- betareg(RLI_2008 ~ BHI_plants_2005 * scaled_threat_count, data = model_input_data)
summary(bii_threat_mod)
qqnorm(resid(bii_threat_mod))
qqline(resid(bii_threat_mod))

bhi_island_mod <- betareg(RLI_2008 ~ BHI_plants_2005 * island_status, data = model_input_data)
summary(bii_threat_mod)
qqnorm(resid(bii_threat_mod))
qqline(resid(bii_threat_mod))

## Don't think we can include disturbance yr as it is also related to BHI through the LUH data
# bhi_island_time_mod <- betareg(RLI_2008 ~ BHI_plants_2005 * island_status * disturbance_year, data = model_input_data)
# summary(bii_threat_mod)
# qqnorm(resid(bii_threat_mod))
# qqline(resid(bii_threat_mod))

bhi_island_threats_mod <- betareg(RLI_2008 ~ BHI_plants_2005 * island_status * scaled_threat_count, data = model_input_data)
summary(bii_threat_mod)
qqnorm(resid(bii_threat_mod))
qqline(resid(bii_threat_mod))

bhi_threat_type_mod <- betareg(RLI_2008 ~ BHI_plants_2005 + predominant_threat_type, 
                               data = model_input_data)
summary(bhi_threat_type_mod)
qqnorm(resid(bhi_threat_type_mod))
qqline(resid(bhi_threat_type_mod))

bii_threat_area_mod <- betareg(RLI_2008 ~ BIIab_2005 + scaled_threat_count + 
                                 ecoregion_area_km_sq, data = model_input_data)
summary(bii_threat_area_mod)
qqnorm(resid(bii_threat_area_mod))
qqline(resid(bii_threat_area_mod))

bii_threat_area_data_mod <- betareg(RLI_2008 ~ BIIab_2005 + scaled_threat_count + 
                                      ecoregion_area_km_sq + scaled_rli_records, data = model_input_data)
summary(bii_threat_area_mod)
qqnorm(resid(bii_threat_area_mod))
qqline(resid(bii_threat_area_mod))

bii_data_mod <- betareg(RLI_2008 ~ BIIab_2005 + rli_records, data = model_input_data)

summary(bii_threat_area_mod)
qqnorm(resid(bii_threat_area_mod))
qqline(resid(bii_threat_area_mod))

hfp_mod <- betareg(RLI_2008 ~ HFP_2005, 
                   data = model_input_data)

summary(hfp_mod)
qqnorm(resid(hfp_mod))
qqline(resid(hfp_mod))

null_mod <- betareg(RLI_2008 ~ 1 | 1, data = model_input_data)
summary(null_mod)
qqnorm(resid(null_mod))
qqline(resid(null_mod))


test_mod <- betareg(RLI_2008 ~ BHI_plants_2005 * island_status * scaled_threat_count|scaled_threat_count, link = "cloglog", data = model_input_data)

test_mod2 <- betareg(RLI_2008 ~ BHI_plants_2005 * island_status * scaled_threat_count|scaled_threat_count + island_status, link = "cloglog", data = model_input_data)


test_mod3 <- betareg(RLI_2008 ~ BHI_plants_2005 * island_status * scaled_threat_count|scaled_threat_count + island_status, link = "loglog", data = model_input_data)



beta_AIC <- as.data.frame(AIC(bii_threat_mod, bhi_threat_mod, bii_threat_area_mod, null_mod, 
                              bhi_mod, bhi_rli_mod, bii_mod, bii_threat_area_data_mod, bii_data_mod,
                              hfp_mod, bhi_threat_type_mod, bhi_island_threats_mod, test_mod, test_mod2, test_mod3))

beta_AIC <- beta_AIC %>% 
  arrange(AIC)

beta_AIC

best_model <- test_mod3

comp <- compare_performance(bii_threat_mod, bhi_threat_mod, bii_threat_area_mod, null_mod, 
                            bhi_mod, bhi_rli_mod, bii_mod, bii_threat_area_data_mod, bii_data_mod,
                            hfp_mod, bhi_threat_type_mod, bhi_island_threats_mod, test_mod, test_mod3, test_mod2)

comp

write.csv(beta_AIC, file.path(main_outputs, "draft_model_AIC.csv"))

# Try a dredge

install.packages("MuMIn")
library(MuMIn)

options(na.action = "na.fail")

dredge_input <- model_input_data %>% 
                dplyr::select(BHI_plants_2005, RLI_2008, BIIab_2005, HFP_2005, biome, disturbance_year, island_status, predominant_threat_type, land_use, realm,
                              ecoregion_area_km_sq, high_beta_area, rli_records, lpi_records, scaled_threat_count)

gm <- betareg(RLI_2008 ~ ., link = "loglog", data = dredge_input)

dredge_models <- dredge(gm)

subset(dredge_models, delta < 4)
# Visualize the model selection table:
par(mar = c(3,5,6,4))
plot(dredge_models, labAsExpr = TRUE)

par(mar = c(3,5,6,4))
plot(dd, labAsExpr = TRUE)
# Model average models with delta AICc < 4
model.avg(dd, subset = delta < 4)

# https://www.stata.com/manuals/rbetareg.pdf

modsum <- tidy(best_model)

binned_residuals(best_model)
plot(best_model, which = 1:4, type = "pearson")
plot(best_model,  which = 5, type = "deviance", sub.caption = "")
plot(best_model,  which = 1, type = "deviance", sub.caption = "")
coef(best_model, model = "precision")
library(lmtest)
coeftest(best_model)
lrtest(best_model)

unique(fitted(best_model))

plot(fitted(best_model),
     residuals(best_model))

x <- cbind(model_input_data$RLI_2008,fitted(best_model))

par(mar=c(5,5,1,1))
cex.lab <- cex.axis <- 2
plot(model_input_data$RLI_2008 ~ model_input_data$BHI_plants_2005, xlab="BHI", ylab="RLI")
axis(1, at=c(1,2,3), cex.axis=cex.axis)
points(fitted(best_model), model_input_data$BHI_plants_2005, col="red", pch=0, cex=2)

write.csv(modsum, file.path(main_outputs, "draft_model_summary.csv"))

# Model type comparisons

# GLM guassian

glmtestg <- glm(RLI_2008 ~ BHI_plants_2005 * island_status * scaled_threat_count, 
               family = "gaussian", data = model_input_data)

# GLM bernoulli

test <- model_input_data %>% 
        mutate(success = round(RLI_2008 * 100),
               fail = 100 - success)

glmtestb <- glm(cbind(success, fail) ~ BHI_plants_2005 * island_status * scaled_threat_count, 
                family = "binomial", data = test)

# GAM betareg

library(mgvc)

gambr <- gam(RLI_2008 ~ BHI_plants_2005 * island_status * scaled_threat_count, 
                        family = betar(), data = model_input_data)

compare_performance(best_model, glmtestb, glmtestg)
