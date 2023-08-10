########
# WALS data (dummy)
########

library(dplyr)  # needed for %>% e.g.
library(tidyr)  # needed for pivot_wider e.g.

### create WALS dataframe (wals_data)
download.file("https://cdstar.shh.mpg.de/bitstreams/EAEA0-7269-77E5-3E10-0/wals_dataset.cldf.zip",
              destfile = "wals_dataset.cldf.zip")

# Read language and value information from the ZIP file:
wals_languages <- read.csv(unz("wals_dataset.cldf.zip", "languages.csv"),
                           header=TRUE, sep=",")
wals_values <- read.csv(unz("wals_dataset.cldf.zip", "values.csv"),header=TRUE, sep=",")
wals_parameters <- read.csv(unz("wals_dataset.cldf.zip", "parameters.csv"),header=TRUE, sep=",")

wals_parameters$Feature <- apply(wals_parameters, 1, function(x) paste("WALS_",x[1]," ",x[2],sep=""))
wals_parameters <- wals_parameters %>% select(c("ID","Feature"))

names(wals_languages)[1]<-"Language_ID"
names(wals_parameters)[1]<-"Parameter_ID"

wals_merged <- wals_languages %>% merge(wals_values, by = "Language_ID", all = TRUE) %>%
  merge(wals_parameters, by = "Parameter_ID", all = TRUE) %>%
  filter(Glottocode!="") %>%
  select(c("Glottocode","Feature","Value")) %>% na.omit()

wals_data<-pivot_wider(wals_merged,
                       names_from = Feature,
                       values_from = Value,
                       values_fill = "?",
                       values_fn = function(x){if (length(unique(x))==1){unique(x)} else ("?")})

lgs <- wals_data$Glottocode
wals_data <- as.data.frame(select(wals_data, -Glottocode))
rownames(wals_data) <- lgs

############################################################################################################

####### Testing the functions with the dummy data (WALS)

# test F1
original_register <- glottocode_taxonomy(path_to_file="data/glottolog_languoid_v4.8/languoid.csv")

# test F2
# original_data must be a data frame with the glottocodes as row names and variable names as column names. any question marks, empty entries, "NA"s must be NAs
original_data <- wals_data
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
mean_type = "arithmetic"
taxonomy = T
tax_weight_factor <- 0.99 # we must multiply all taxonomic weights by a number below (but close to) 1, because if tax_weight ever reaches 1, the odds mean will be undefined!
coding_weight_factor <- 0.999 # we should also multiply all non-taxonomic weights by a factor (identical to the tax_weight_factor) to not privilege coding over taxonomy!
original_register <- original_register

set.seed(2023)
documentation_A_T_2023 <- densify_steps(original_data, max_steps, mean_type, taxonomy, original_register, tax_weight_factor, coding_weight_factor)

write.csv(documentation_A_T_2023,"R/documentation files/documentation_A_T_2023b.csv")

documentation_A_T_2023 <- read.csv("R/documentation files/documentation_A_T_2023.csv")


# test F3
exponent_prop_coded_data <- 1
exponent_available_data_points <- 1
exponent_lowest_language_score <- 0
exponent_taxonomic_diversity <- 0

documentation = documentation_A_T_2023

optimum <- densify_score(documentation,exponent_prop_coded_data=exponent_prop_coded_data, exponent_available_data_points=exponent_available_data_points, exponent_lowest_language_score=exponent_lowest_language_score)

# test F4
pruned_df <- densify_prune(original_data, documentation_A_T_2023, optimum)

######### some plots
hist(apply(pruned_df,1,function(x)(length(na.omit(x))))/ncol(pruned_df),col = "cadetblue2",xlab="coding density per language", ylab="frequency")

# proportion coded languages
sum(!is.na(pruned_df))/(ncol(pruned_df)*nrow(pruned_df))
nrow(pruned_df)
register_pruned <- filter(original_register, glottocode %in% rownames(pruned_df))
length(unique(register_pruned$glottolog.node1))
ncol(pruned_df)
# rm(pruned_df)

# Also prep the pruned matrix:
densify_prep <- function(original_data) {
  # save row names for later
  languages <- rownames(original_data)
  # replace NAs by 0 (no data available) and non-NA entries by 1 (data available)
  full_matrix <- as.matrix(original_data)
  full_matrix[!is.na(full_matrix)] <- 1
  full_matrix[is.na(full_matrix)] <- 0
  # convert dataframe entries to numeric
  full_matrix <- apply(full_matrix, 2, as.numeric)
  # rename row names
  rownames(full_matrix) <- languages
  return(full_matrix)
}

pruned_matrix = densify_prep(pruned_df)

####

# histogram original matrix
hist(apply(original_data,1,function(x)(length(na.omit(x))))/ncol(original_data),col = "cadetblue2",xlab="coding density per language", ylab="frequency")
sum(!is.na(original_data))/(ncol(original_data)*nrow(original_data))

nrow(original_data)
register_original <- filter(original_register, glottocode %in% rownames(original_data))
length(unique(register_original$glottolog.node1))
ncol(original_data)

# number of families
register <- filter(original_register, glottocode %in% rownames(original_data))
register <- register[match(rownames(original_data),register$glottocode),]
length(unique(register$glottolog.node1))
