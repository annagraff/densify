---
title: 'densify: An R package to prune sparse data frames of typological linguistic data'
tags:
- sparse matrices
- "sub-sampling"
- linguistic data
- diversity samples
- R
date: "22 January 2024"
output: pdf_document
citation_author: Graff et al.
authors:
- name: Anna Graff
  orcid: "0000-0002-7703-3471"
  equal-contrib: yes
  corresponding: yes
  affiliation: 1, 2, 3
- name: Marc Lischka
  orcid: XXXXXXXX
  equal-contrib: yes
  affiliation: 4
- name: Taras Zakharko
  orcid: XXXXX
  affiliation: 1, 3
- name: Reinhard Furrer
  orcid: "0000-0002-6319-2332"
  affiliation: 3, 4
- name: Balthasar Bickel
  orcid: "0000-0002-9087-0565"
  affiliation: 1, 3
year: 2024
bibliography: sources.bib
affiliations:
- name: University of Zurich, Department of Comparative Language Science
  index: 1
- name: University of Zurich, Department of Evolutionary Biology and Environmental
    Studies
  index: 2
- name: University of Zurich, Center for the Interdisciplinary Study of Language Evolution
  index: 3
- name: University of Zurich, Department of Mathematics
  index: 4
---

<!--General comment: I am not familiar with the journal but I wonder whether we shouldn't explain a bit more what linguistic typological data is? Perhaps given an example and illustrate the problem? 

RF: yes. probably better. at least much better than to include detailed list of arguments. Even if we are slightly above the 1000 word limit. See also comments below.
-->

# Summary

The R package `densify` provides a procedure to prune input data frames containing empty cells (or cells with values {?} or {NA}) to denser sub-matrices with fewer rows, columns and empty cells. The pruning process, split into three functions, trades off a series of variably weighted concerns, including data retention, coding density (proportion of non-empty cells) and taxonomic diversity of rows (if the entities are represented in a taxonomic structure). Users can adapt the relative weights given to these concerns through various parameters for the densification process to best fit their needs. As such, the software is useful for several purposes, including the densification of sparse input matrices and the subsampling of large input matrices according to a procedure that is sensitive to taxonomic structure.

# Statement of Need

While the software will run on any data frame (with rows representing any entities with or without taxonomic structure and columns representing variables), it was primarily designed to prune data frames of linguistic typological data, where the rows are languages (taxa) and the columns typological variables (sometimes also called characters, parameters or features).

Linguistic typological data is increasingly available in large-scale databases, and many analyses that aim at exploring diversity or testing hypotheses rely on such databases. Some of these resources have information for nearly all variables across nearly all languages in the database (i.e. they have complete or near-complete coding density), but the dataframe may be too large for certain computationally intensive analyses (e.g., PHOIBLE [@phoible], Grambank [@grambank]). Other databases (e.g., WALS [@wals], AUTOTYP [@AUTOTYP], Lexibank [@Lexibank]) exhibit variables that are coded for different sets of languages, resulting in sparse language-variable matrices. Combining data from various databases via language identifiers like glottocodes usually increases sparsity because the language samples do not match. 

When datasets are too large or too sparse for computational applications, it can thus be necessary to generate and subsequently operate on a subset of the data represented in a smaller and denser matrix. Thereby, researchers might be particularly interested in maintaining taxonomic diversity in the languages represented, preferentially removing languages belonging to clades represented by many other languages in the sample and penalizing the removal of language isolates or languages which represent small language families.

While certain packages exist to generate sub-matrices from varying input matrices according to principled criteria (e.g., `admmDensestSubmatrix` [@R-admmDensestSubmatrix], which identifies the densest sub-matrix of an input graph of a specified size; or `FSelector` [@R-FSelector] and `varrank` [@R-varrank], which perform attribute subset selection based on various tests and entropy measures to identify the most relevant attributes of a data input), `densify` provides a flexible, explicit, and taxonomy-sensitive pruning algorithm that focuses both on the removal of rows and columns and does not require the size of the sub-matrix to be specified a priori.

# Usage

The package `densify` provides the data from The World Atlas of Language Structures (WALS) [@wals] and the language taxonomy provided by Glottolog v. 4.8 [@glottolog] as example data. The accompanying package vignette features a detailed demonstration of the utility and flexibility of `densify` to subsample an input matrix according to varying needs, using this data.

## Preparing input

The data frame that requires subsetting must have rows representing taxa or observations (with taxon names provided in a dedicated column) and columns representing variables (and variable names as column names). Any cells with empty entries, not applicable or question marks must be coded as NA. If matrix densification should be sensitive to taxonomic structure, a taxonomy must be provided as a `phylo` object [cf. @R-ape] or as an adjacency table (i.e. a data frame containing columns `id` and `parent_id`, with each row encoding one parent-child relationship). Every taxon in the input data frame must be included in the taxonomy (as a tip or node). If a language taxonomy is to be provided, glottolog (`glottolog_languoids`) can be used directly.

``` r
install.packages("densify")
library(densify)

# prepare data: WALS and Glottolog
data(WALS)
data(glottolog_languoids)

# any question marks, empty entries, "NA"s must be coded as NA
WALS[WALS=="?"] <- NA
WALS[WALS=="NA"] <- NA
head(WALS)

# all taxa must be present in the taxonomy used for pruning
WALS <- WALS[which(WALS$Glottocode %in% glottolog_languoids$id), ]
```

## Pruning

Iterative pruning of the input matrix is performed by the `densify()` function, which requires the following information:

<!-- not sure all the arguments should be included here. When I read about the package, I am interested in the main functionality and the possible output. With cool interpretations thereof.
-->

-   The original data frame with observations in rows and variables in columns (`data`). No default.
-	A specification of which columns should be densified (`cols`). Default: all columns densified.
-   A taxonomy tree as a phylo object or a data frame with columns `id` and `parent_id` (`taxonomy`). No default.
-	The name fo the column identifying taxa (`taxon_id`).
-	A string specifying the scoring type used for calculating row-wise importance weights (`scoring`). Possible values are `"arithmetic"`, `"geometric"`, `"log_odds"`. Default: `"log_odds"`.
-   An optional integer specifying the threshold for variability in variables (`min_variability`). Default: `1`.
-	A logical argument denoting whether taxonomic diversity of taxa should be considered in the pruning process (`consider_taxonomic_diversity`). Default: `TRUE`, if taxonomy provided.

<!-- -   A taxonomy weight factor, required if `mean_type = "log_odds"` and `consider_taxonomic_diversity = TRUE`. Default: 0.99
-   A coding weight factor, required if `mean_type = "log_odds"`. Default: 0.99-->

For a more detailed discussion of the parameters, refer to the vignette hosted in the software repository.

The output of the function is a densify\_result object, documenting several summary statistics of the sub-matrix resulting from each iteration. These include the number of available data points, the overall proportion of coded data in the matrix, the coding densities of the least well-coded taxon and variable, the coding densities of the most well-coded taxon and variable, the average coding densities of the all taxa and variables, and a taxonomic diversity index (Shannon diversity of the highest taxonomic level).

``` r
set.seed(2023)
example_result <- densify(data = WALS, 
                          taxonomy = glottolog_languoids, 
                          taxon_id = "Glottocode", 
                          scoring = "log_odds", 
                          min_variability = 3, 
                          consider_taxonomic_diversity = TRUE)
head(example_result)
```
## Finding the optimal number of iterations and producing the sub-matrix

The functions `prune()` and `rank_results()` can identify the optimal sub-matrix via a quality score, computed via a user-defined function composed of any combination of available statistics from the densify\_result object. The default score thereby maximizes the product of the number of available data points and the overall coding density.

The function `prune()` thereby retrieves the optimal sub-matrix, while `rank_results()` returns the relative ranks of all sub-matrices given the score. 

``` r
# use prune() to obtain the optimum submatrix, given different scores
example_optimum_1 <- prune(example_result, 
                           score = n_data_points*coding_density)
example_optimum_2 <- prune(example_result, 
                           score = n_data_points*coding_density*taxonomic_index^3)

# use rank_results() to obtain the same optima alongside 
example_ranks_1 <- rank_results(example_result, 
                                score = n_data_points*coding_density)
example_ranks_2 <- rank_results(example_result, 
                                score = n_data_points*coding_density*taxonomic_index^3)
```

The relative ranking of sub-matrices given the specified quality score can also be visualized using `visualize()` or `plot()`.

``` r
# use visualize() to illustrate quality scores and optimum (first score)
visualize(example_result, score = n_data_points*coding_density)

# use plot() to illustrate quality scores and optimum (second score)
plot(example_result, score = n_data_points*coding_density*taxonomic_index^3)
```

# Conclusions

The R package `densify` provides users with a flexible and explicit method to generate sub-matrices from an input matrix in a mathematically principled way. This paper shows case examples using a standard sparse linguistic dataset (WALS) and the standard linguistic taxonomy provided by Glottolog.
Additional examples and usage details are found in the vignette and hosted in the software repository on GitHub.

# Acknowledgements

The authors declare that there are no conflicts of interest.

# References
