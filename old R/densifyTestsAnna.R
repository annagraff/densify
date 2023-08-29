rm(list=ls())
library(tidyverse) 

source("F1_glottocode_taxonomy.R")
source("F2_densify_steps.R")
source("F3_densify_score.R")
source("F4_densify_prune.R")


####### Testing the functions with the Grambank data (Grambank statistical)
gb_original_data <- read.csv("../data/Original_Grambank_6_august_2023.csv",row.names = "glottocode") %>% select(-X)

# 
original_register <- glottocode_taxonomy(path_to_file="../data/glottolog_4.8_languoid/languoid.csv")

  # test F2
  # original_data must be a data frame with the glottocodes as row names and variable names as column names. any question marks, empty entries, "NA"s must be NAs
  original_data <- gb_original_data
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

  max_steps = nrow(original_data)+ncol(original_data)
  mean_type = "log_odds"
  taxonomy = T
  tax_weight_factor <- 0.999 # we must multiply all taxonomic weights by a number below (but close to) 1, because if tax_weight ever reaches 1, the odds mean will be undefined!
  coding_weight_factor <- 0.99 # we should also multiply all non-taxonomic weights by a factor (identical to the tax_weight_factor) to not privilege coding over taxonomy!
  original_register <- original_register

  set.seed(2023)
  documentation_L_T_2023_0999_099 <- densify_steps(original_data, max_steps, mean_type, taxonomy, original_register, tax_weight_factor, coding_weight_factor)
  write.csv(documentation_L_T_2023_0999_099,"documentation files/documentation_L_T_2023_0999_099.csv")
  
  documentation_L_T_2023_0999_099<-read.csv("documentation files/documentation_L_T_2023_0999_099.csv")
  
  
  # test F3
  exponent_prop_coded_data <- 1
  exponent_available_data_points <-1
  exponent_lowest_language_score <- 1
  exponent_taxonomic_diversity <- 1
  
  optimum <- densify_score(documentation_L_T_2023_0999_099,exponent_prop_coded_data=exponent_prop_coded_data, exponent_available_data_points=exponent_available_data_points, exponent_lowest_language_score=exponent_lowest_language_score, exponent_taxonomic_diversity=exponent_taxonomic_diversity)

# test F4
# these two are equivalent
pruned_matrix <- densify_prune(original_data, documentation_L_T_2023_0999_099, optimum)

######### some plots
hist(apply(pruned_matrix,1,function(x)(length(na.omit(x))))/ncol(pruned_matrix),col = "cadetblue2",xlab="coding density per language", ylab="frequency")
hist(apply(pruned_matrix,2,function(x)(length(na.omit(x))))/nrow(pruned_matrix),col = "cadetblue2",xlab="coding density per variable", ylab="frequency")
  
# proportion coded languages
sum(!is.na(pruned_matrix))/(ncol(pruned_matrix)*nrow(pruned_matrix))
nrow(pruned_matrix)
register_pruned <- filter(original_register, glottocode %in% rownames(pruned_matrix))
length(unique(register_pruned$glottolog.node1))
ncol(pruned_matrix)


  