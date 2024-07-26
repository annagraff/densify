# densify

densify ist an R package for densifying (sparse) matrices.

# Installation (for users)

To install densify, run the following code in R:

~~~~
install.packages("devtools")
library("devtools")
install_github('annagraff/densify', build_vignettes = T)
library(densify)
~~~~

# Performing matrix densification with densify

## Preparing input

The data frame that requires subsetting must have rows representing taxa or observations (with taxon names provided in a dedicated column) and columns representing variables (and variable names as column names). Any cells with empty entries, not applicable or question marks must be coded as `NA`. If taxonomic structure is relevant to the pruning process, a taxonomy must be provided as a `phylo` object or as an adjacency table (i.e. a data frame containing columns `id` and `parent_id`, with each row encoding one parent-child relationship). The `glottolog_languoids` dataframe provided by the package can be used directly for this purpose.

~~~~
# prepare example data: WALS and Glottolog
data(WALS)
data(glottolog_languoids)

# any question marks, empty entries, "NA"s must be coded as NA
WALS[WALS=="?"] <- NA
WALS[WALS=="NA"] <- NA
head(WALS)

# all taxa must be present in the taxonomy used for pruning
WALS <- WALS[which(WALS$Glottocode %in% glottolog_languoids$id), ]
~~~~

## Densifying, visualizing, ranking and pruning the input
The `densify()` function iteratively prunes the input matrix. It can be modulated by several parameters to, among others: specify to what extent taxonomic diversity of the sample should weigh in the pruning process; choose among various methods to calculate importance weights (e.g. via arithmetic or logit-transformed means of row-wise coding density); and choose the minimum variation required for a feature to be retained (e.g. constant variables might be of no interest). For a detailed discussion of the parameters, refer to the function documentation or to the vignette hosted in the software repository.

~~~~
set.seed(2024)
example_result <- densify(data = WALS, 
                          cols = !Glottocode,
                          taxonomy = glottolog_languoids, 
                          taxon_id = "Glottocode", 
                          density_mean = "log_odds", 
                          min_variability = 3,
                          limits = list(min_coding_density = 1),
                          density_mean_weights = list(coding = 1, taxonomy = 1))
~~~~

The output of the function is a `densify_result` object, documenting several summary statistics of all resulting sub-matrices. These summary statistics can be used to define a scoring function, which is used by `rank_results()`, `visualize()` and `prune()` to identify the optimum: `rank_results()` returns the relative ranks of all generated sub-matrices given the scoring function, `visualize()` visualizes their relative ranking, and `prune()` extracts the optimal sub-matrix (ranked first).

~~~~
head(example_result)

# use rank_results() to obtain a vector indicating the rank of each sub-matrix
# with the default scoring function
example_ranks_1 <- rank_results(example_result, 
                   scoring_function = n_data_points*coding_density)
                                
# with a scoring function that gives high weight to taxonomic diversity:
example_ranks_2 <- rank_results(example_result, 
                   scoring_function = n_data_points*coding_density*taxonomic_index^3)

# use visualize() to illustrate the quality scores and optimum given each scoring function
visualize(example_result, 
          scoring_function = n_data_points*coding_density)
visualize(example_result, 
          scoring_function = n_data_points*coding_density*taxonomic_index^3)

# use prune() to obtain the optimum sub-matrix given each scoring function
example_optimum_1 <- prune(example_result, 
                     scoring_function = n_data_points*coding_density)
                           
example_optimum_2 <- prune(example_result, 
                     scoring_function = n_data_points*coding_density*taxonomic_index^3)
~~~~

For more details on each function, refer to the help pages (see below) and/or resort to the publication paper and vignette. 
~~~~
?densify
?rank_results
?visualize
?prune
~~~~

## Contributing

To report bugs, seek support or suggest improvements (e.g. additional functionalities or changes to functionality or API arguments), please open an issue on GitHub. Similarly, to make direct package contributions, please open an issue on GitHub for discussion before submitting a merge request.

### Code of Conduct

Please note that `densify` is released with a Contributor Code of Conduct. By contributing to this project, you agree to abide by its terms.
