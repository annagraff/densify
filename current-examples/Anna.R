rm(list=ls())

source("../package_densify/R/flat_taxonomy_matrix.R")
source("../package_densify/R/densify_steps.R")
source("../package_densify/R/densify_score.R")
source("../package_densify/R/densify_prune.R")

library(readr)
library(plyr)
library(phytools)
library(tidyverse)
library(testthat)
library(GetoptLong)
library(vegan)

# paper tests

# Read language and value information from the ZIP file:
wals_languages <- read.csv(unz("../old-data/wals_dataset.cldf.zip", "languages.csv"),
                           header=TRUE, sep=",")
wals_values <- read.csv(unz("../old-data/wals_dataset.cldf.zip", "values.csv"),header=TRUE, sep=",")
wals_parameters <- read.csv(unz("../old-data/wals_dataset.cldf.zip", "parameters.csv"),header=TRUE, sep=",")

wals_parameters$Feature <- apply(wals_parameters, 1, function(x) paste(x[1]," ",x[2],sep=""))
wals_parameters <- wals_parameters %>% select(c("ID","Feature"))

names(wals_languages)[1]<-"Language_ID"
names(wals_parameters)[1]<-"Parameter_ID"

wals <- wals_languages %>% merge(wals_values, by = "Language_ID", all = TRUE) %>%
  merge(wals_parameters, by = "Parameter_ID", all = TRUE) %>%
  filter(Glottocode!="") %>%
  select(c("Glottocode","Feature","Value")) %>% na.omit()

wals <- as.data.frame(pivot_wider(wals,
                       names_from = Feature,
                       values_from = Value,
                       values_fill = "?",
                       values_fn = function(x){if (length(unique(x))==1){unique(x)} else ("?")}))

lgs <- wals$Glottocode
rownames(wals) <- lgs
wals <- select(wals,-Glottocode)

saveRDS(wals,file="../package_densify/data/wals.rda")

glottolog_languoids <- readr::read_csv("../old-data/glottolog_languoid_v4.8/languoid.csv")
saveRDS(glottolog_languoids,file="../package_densify/data/glottolog_languoids.rda")


wals[wals=="?"] <- NA
wals[wals=="NA"] <- NA
head(wals)


# test F1
taxonomy_matrix <- build_flat_taxonomy_matrix(id = glottolog_languoids$id, parent_id = glottolog_languoids$parent_id)
head(taxonomy_matrix)

tic()
set.seed(2023)
documentation <- densify_steps(original_data = wals, max_steps = nrow(wals)+ncol(wals), mean_type = "log_odds", taxonomy = TRUE, taxonomy_matrix = taxonomy_matrix, tax_weight_factor = 0.99, coding_weight_factor = 0.99)
head(documentation)
toc()

exponent_prop_coded_data <- 1
exponent_available_data_points <- 1
exponent_lowest_taxon_coding_score <- 1
exponent_lowest_variable_coding_score <- 0
exponent_taxonomic_diversity <- 1

optimum <- densify_score(documentation = documentation, 
                         exponent_prop_coded_data = exponent_prop_coded_data, 
                         exponent_available_data_points = exponent_available_data_points, 
                         exponent_lowest_taxon_coding_score = exponent_lowest_taxon_coding_score,
                         exponent_lowest_variable_coding_score = exponent_lowest_variable_coding_score,
                         exponent_taxonomic_diversity = exponent_taxonomic_diversity)





####### Testing the functions with the Grambank data (Grambank statistical)
original_data <- read.csv("../data/Original_Grambank_6_august_2023.csv",row.names = "glottocode") %>% select(-X)

# test F2
# original_data must be a data frame with the glottocodes as row names and variable names as column names. any question marks, empty entries, "NA"s must be NAs
original_data[original_data=="?"] <- NA
original_data[original_data=="NA"] <- NA

# make sure that all variables are/remain sufficiently variable --> the second-most-frequent variable state must contain at least 3 languages
# if any variables are not/no longer sufficiently variable for our areal analysis, remove them
nrlevels <- data.frame(variable=colnames(original_data),
                         number_of_variable_states=apply(original_data,2,function(x)length(table(as.factor(x)))),
                         count_second_largest_variable_state=apply(original_data,2,function(x)sort(table(as.factor(x)),decreasing=T)[2]))

  if (sum(nrlevels$count_second_largest_variable_state%in%c(NA,1,2))!=0){ # only act if there is a variable that needs removal
    uninformative_variables <- rownames(filter(nrlevels,count_second_largest_variable_state%in%c(NA,1,2)))
    original_data <- original_data[,-which(colnames(original_data)%in%uninformative_variables)] # update matrix by pruning away uninformative variables
  }

max_steps = ncol(original_data)+nrow(original_data)
mean_type = "log_odds"
taxonomy = T
tax_weight_factor <- 0.999 # we must multiply all taxonomic weights by a number below (but close to) 1, because if tax_weight ever reaches 1, the odds mean will be undefined!
coding_weight_factor <- 0.99 # we should also multiply all non-taxonomic weights by a factor (identical to the tax_weight_factor) to not privilege coding over taxonomy!
taxonomy_matrix <- taxonomy_matrix

set.seed(2023)
d <- densify_steps(original_data, max_steps, mean_type, taxonomy, taxonomy_matrix, tax_weight_factor, coding_weight_factor)
write.csv(documentation_GBO2208_L_T_2023_0999_099,"documentation files/documentation_GBO2208_L_T_2023_0999_099.csv")


dL <- read.csv("documentation files/documentation_GBL1808_L_T_2023_0999_099.csv") %>% select(-X)
dS <- read.csv("documentation files/documentation_GBS1808_L_T_2023_0999_099.csv") %>% select(-X)

names(dL)<-names(documentation_GBO2208_L_T_2023_0999_099)
names(dS)<-names(documentation_GBO2208_L_T_2023_0999_099)


# test F3
exponent_prop_coded_data <- 1
exponent_available_data_points <- 1
exponent_lowest_taxon_coding_score <-1
exponent_lowest_variable_coding_score <-1
exponent_taxonomic_diversity <- 1

optimum <- densify_score(d,
                         exponent_prop_coded_data=exponent_prop_coded_data,
                         exponent_available_data_points=exponent_available_data_points,
                         exponent_lowest_taxon_coding_score=exponent_lowest_taxon_coding_score,
                         exponent_lowest_variable_coding_score=exponent_lowest_variable_coding_score,
                         exponent_taxonomic_diversity=exponent_taxonomic_diversity)

# test F4
# these two are equivalent
pruned_matrix <- densify_prune(original_data, d, optimum)
#pruned_matrix <- original_data
  # ######### some plots
  # hist(apply(pruned_matrix,1,function(x)(length(na.omit(x))))/ncol(pruned_matrix),col = "cadetblue2",xlab="coding density per language", ylab="frequency")
  # hist(apply(pruned_matrix,2,function(x)(length(na.omit(x))))/nrow(pruned_matrix),col = "cadetblue2",xlab="coding density per variable", ylab="frequency")
  # 
  # proportion coded languages
  sum(!is.na(pruned_matrix))/(ncol(pruned_matrix)*nrow(pruned_matrix))
  nrow(pruned_matrix)
  register_pruned <- filter(taxonomy_matrix, id %in% rownames(pruned_matrix))
  length(unique(register_pruned$level1))
  ncol(pruned_matrix)
  
# write.csv(pruned_matrix,"original-220823-pruned-GBMO-11111.csv")
  # write.csv(pruned_matrix,"logical-220823-pruned-GBML-11111.csv")
  # write.csv(pruned_matrix,"statistical-220823-pruned-GBMS-11111.csv")
  