---
title: 'densify: An R package to prune sparse data frames of typological linguistic
  data'
tags:
- sparse matrices
- "sub-sampling"
- linguistic data
- diversity samples
- R
date: "29 August 2023"
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
  affiliation: 3, 4
- name: Taras Zakharko
  orcid: XXXXX
  affiliation: 1, 3
- name: Reinhard Furrer
  orcid: XXXXX
  affiliation: 3, 4
- name: Balthasar Bickel
  orcid: "0000-0002-9087-0565"
  affiliation: 1, 3
year: 2023
bibliography: sources.bib
affiliations:
- name: University of Zurich, Department of Comparative Language Science
  index: 1
- name: University of Zurich, Department of Evolutionary Biology and Environmental
    Studies
  index: 2
- name: University of Zurich, Center for the Interdisciplinary Study of Language Evolution
  index: 3
- name: University of Zurich, Institute of Mathematics
  index: 4
---

<!--General comment: I am not familiar with the journal but I wonder whether we shouldn't explain a bit more what linguistic typological data is? Perhaps given an example and illustrate the problem? -->

# Summary

The R package `densify` provides a procedure to prune input data frames containing empty cells (or cells with value "?" or NA) to denser sub-matrices with less rows, columns and empty cells. The pruning process, split into three functions, trades off a series of variably weighted concerns, including data retention, coding density (proportion of non-empty cells) and taxonomic diversity of rows (if the entities are represented in a taxonomic structure). Users can adapt the relative weight given to these concerns through various parameters for the densification process to best fit their needs. As such, the software is useful for several purposes, including the densification of sparse input matrices and the subsampling of large input matrices according to a procedure that is sensitive to taxonomic structure.

# Statement of Need

While the software will run on any data frame (with rows representing any entities with or without taxonomic structure and columns representing variables), it was designed to prune data frames of linguistic typological data, where the rows are languages (taxa) and the columns typological variables (sometimes also called characters, parameters or features).

Linguistic typological data is increasingly available in large-scale databases, and many analyses that aim at testing hypotheses at global scales rely on such databases. Some of these resources have information for nearly all variables across nearly all languages in the database (i.e. they have complete or near-complete coding density), but the dataframne may be too large for certain computationally intensive analyses (e.g. PHOIBLE [@phoible], Grambank [@grambank]). Other databases (e.g. WALS [@wals], AUTOTYP [@AUTOTYP], Lexibank [@Lexibank]) exhibit variables that are coded for very different sets of languages, resulting in sparse language-variable matrices. Combining data from various databases via language identifiers like glottocodes usually increases sparsity because the language samples do not match.

Although there are methods available to operate on such matrices (e.g. `Matrix` [@R-Matrix], `MatrixExtra` [@R-MatrixExtra]), sparse matrices may lack power for certain research questions, such that researchers may prefer to operate on a subset of the data represented in a denser matrix. <!-- BB: Why are the matrix packages relevent here. There is a presupposition I am missing.:: ANNA: This is a requirement of the journal for the statement of need section ("A Statement of need section that clearly illustrates the research purpose of the software and places it in the context of related work. A list of key references, including to other software addressing related needs.": https://joss.readthedocs.io/en/latest/submitting.html) -->

When generating subsets of a language-variable matrix, researchers might be particularly interested in maintaining taxonomic diversity in the languages represented, preferentially removing languages belonging to clades represented by many other languages in the sample and penalizing the removal of so-called language isolates or languages which represent language families with few members.

While certain packages exist to generate sub-matrices from varying input matrices according to principled criteria (e.g. `admmDensestSubmatrix` [@R-admmDensestSubmatrix], which identifies the densest sub-matrix of an input graph of a specified size; or `FSelector` [@R-FSelector], which performs attribute subset selection based on various tests and entropy measures to identify the most relevant attributes of a data input), `densify` provides a flexible, explicit, and taxonomy-sensitive pruning algorithm that focuses both on the removal of rows and columns and does not require specification of the size of the sub-matrix.

# Usage

`densify` provides the data from The World Atlas of Language Structures (WALS [@wals]) and the language taxonomy provided by Glottolog v. 4.8 [@glottolog] as example data. The vignette features a detailed demonstration of the utility and flexibility of `densify` to subsample an input matrix according to varying needs, using this data.

## Preparing input

`densify` requires the input data frame to be a data frame with rows as taxa or observations (and taxon names as row names) and columns as variables (and variable names as column names). Any cells with empty entries, not applicable or question marks must be coded as NA. If matrix densification should consider taxonomic structure, a flat taxonomy must be provided, listing every taxon present in the initial data frame along with all the nodes connecting it to the root. Such a taxonomy can be generated with the `build_flat_taxonomy_matrix()` function, if all nodes and tips are provided alongside each of their parent nodes. For generating a language taxonomy, glottolog can be used directly.

<!-- I think you need to explain how a non-flat taxonomy looks like, i.e. what the input format is. Readers might assume newick or something like that. Also, to maximize usability, I am wondering whether you could point to tools that convert a newick or similar representation into the ID-parentID format of glottolog? ::-->

``` r
install.packages("densify")
library(densify)

# prepare data: WALS and glottolog
data(wals)
data(glottolog_languoids)

# the input data frame must be a data frame with the taxon names as row names and variable names as column names; any question marks, empty entries, "NA"s must be coded as NAs
wals[wals=="?"] <- NA
wals[wals=="NA"] <- NA
head(wals)

# generate flat taxonomy using glottolog
taxonomy_matrix <- build_flat_taxonomy_matrix(id = glottolog_languoids$id, parent_id = glottolog_languoids$parent_id)
head(taxonomy_matrix)
```

## Pruning

Iterative pruning of the input matrix is performed by the `densify_steps` function, which requires the following information:

-   The original data frame. Default: `wals`
-   An integer specifying the number of iterations performed. Default: 1
-   An integer specifying the threshold for variability in variables (columns). Default: 1
-   The mean type, which can be: "arithmetic", "geometric", "log_odds". Default: "log_odds".
-   Should taxonomy be accounted for in the pruning process? Default: `FALSE`
-   The taxonomy matrix, required if taxonomy = `TRUE`, optional otherwise. Default: NULL.
-   A taxonomy weight factor, required if mean_type = "log_odds" and taxonomy = `TRUE`. Default: 0.99
-   A coding weight factor, required if mean_type = "log_odds". Default: 0.99

For a more detailed discussion of the parameters, refer to the vignette hosted in the software repository.

The output of the function is an iteration log, documenting key characteristics of the sub-matrix resulting from each iteration. 

``` r
set.seed(2023)
example_iteration_log <- densify_steps(original_data = wals, max_steps = nrow(wals)+ncol(wals)-2, variability_threshold = 3, mean_type = "log_odds", taxonomy = TRUE, taxonomy_matrix = taxonomy_matrix, tax_weight_factor = 0.99, coding_weight_factor = 0.99)
head(iteration_log)
```
## Retrieving the optimal number of iterations and producing the sub-matrix

Given the pruning iteration log output, the second function `densify_score` will identify the optimal sub-matrix via a quality score, computed via user-defined exponents relating to the overall proportion of coded data in the matrix, the number of available data points, the coding density of the least well-coded taxon, the coding density of the least well-coded variable, and a taxonomic diversity index (Shannon diversity of the highest taxonomic level).

Finally, the function `densify_prune` will prune the original matrix to the optimal submatrix, given the outputs of `densify_steps` and `densify_score`.

<!--# I think it would be good to explain what the purpose of the documentation is, how its contents look like, and what it means.-->

``` r
example_optimum <- densify_score(iteration_log = example_iteration_log, 
                         exponent_prop_coded_data = 1, 
                         exponent_available_data_points = 1,
                         exponent_lowest_taxon_coding_score = 1,
                         exponent_lowest_variable_coding_score = 1,
                         exponent_taxonomic_diversity = 1,
						 plot = TRUE)

pruned_wals <- densify_prune(original_data = wals, 
							iteration_log = example_iteration_log, 
							optimum = example_optimum)
```

# Conclusions

The R package `densify` provides users with a flexible and explicit method to generate sub-matrices from an input matrix in a mathematically principled way. This paper shows a case example using a standard sparse linguistic dataset (WALS) and the standard linguistic taxonomy provided by Glottolog. Additional examples and usage details are hosted in the vignette and found in the software repository on GitHub.

# Acknowledgements

The authors declare that there are no conflicts of interest.

# References
