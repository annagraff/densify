---
title: 'densify: An R package to prune sparse data frames of typological linguistic data'
tags:
- sparse matrices
- "sub-sampling"
- linguistic data
- diversity samples
- R
date: "5 June 2024"
output: pdf_document
citation_author: Graff et al.
authors:
- name: Anna Graff
  orcid: "0000-0002-7703-3471"
  corresponding: yes
  affiliation: 1, 2, 3
- name: Marc Lischka
  orcid: 0009-0007-9493-2392
  affiliation: 4
- name: Taras Zakharko
  orcid: 0000-0001-7601-8424
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
- name: University of Zurich, Department of Mathematical Modeling and Machine Learning
  index: 4
---

# Summary

The R package `densify` provides a procedure to prune input data frames containing empty cells (or cells with values {?} or {NA}) to denser sub-matrices with fewer empty cells. The pruning process trades off a series of variably weighted concerns, including data retention, coding density (proportion of non-empty cells) and taxonomic diversity of rows (representing for example phylogenetic relations). Users can adapt the relative weights given to these concerns through various parameters so that the densification process best fits their needs. As such, the software is useful for several purposes, including the densification of sparse input matrices and the subsampling of large input matrices according to a procedure that is sensitive to taxonomic structure. 

# Statement of Need
Linguistic typological data is increasingly available in large-scale databases, and many analyses that aim at exploring diversity or testing hypotheses rely on such databases. Some of these resources have information for nearly all features (sometimes also called characters, parameters or variables) across nearly all languages in the database (i.e. they have complete or near-complete coding density), but the dataframe may be too large for certain computationally intensive analyses (e.g., PHOIBLE [@phoible], Grambank [@grambank]). Other databases (e.g., WALS [@wals], AUTOTYP [@AUTOTYP], Lexibank [@Lexibank]) exhibit features that are coded for different sets of languages, resulting in sparse language-feature matrices. Combining data from various databases via language identifiers like glottocodes [@glottolog] usually increases sparsity because the language samples do not match. 

When datasets are too large or too sparse for computational applications, it can thus be necessary to generate and subsequently operate on a subset of the data represented in a smaller and denser matrix. Thereby, researchers might be particularly interested in maintaining taxonomic diversity in the languages represented, preferentially removing languages belonging to clades represented by many other languages in the sample and penalizing the removal of language isolates or languages which represent small language families.

While certain packages exist to generate sub-matrices from varying input matrices according to principled criteria (e.g., `admmDensestSubmatrix` [@R-admmDensestSubmatrix], which identifies the densest sub-matrix of an input graph of a specified size; or `FSelector` [@R-FSelector] and `varrank` [@R-varrank], which perform attribute subset selection based on various tests and entropy measures to identify the most relevant attributes of a data input), `densify` adds sensitivity to taxonomic structure and generally more flexibility in parameter settings for the user. The algorithm focuses both on the removal of rows and columns and does not require the size of the sub-matrix to be specified a priori. 

The software therefore addresses recurring problems researchers face when working with typological linguistic data, and it was primarily designed to handle such data frames (with rows representing languages and columns representing typological features). However, it will run on any data frame with rows representing any entities with or without taxonomic structure and columns representing variables. It may thus be of use for other applications as well.

# Usage

The package `densify` provides the data from The World Atlas of Language Structures (WALS) [@wals] and the language taxonomy provided by Glottolog v. 5.0 [@glottolog] as example data. The accompanying package vignette features a detailed demonstration of the utility and flexibility of `densify` to subsample an input matrix according to varying needs, using this data.

Input data must be prepared in a dataframe with rows representing taxa or observations with taxon names specified in a dedicated column, and columns representing variables with variable names as column names. Cells that are empty or contain non-applicable or question mark entries must be coded as `NA`. If densification should be sensitive to taxonomic structure, a taxonomy must be provided as (i) a `phylo` object [cf. @R-ape], (ii) as an adjacency table (i.e. a data frame containing columns `id` and `parent_id`, with each row encoding one parent-child relationship), or (iii) by the `glottolog_languoids` dataframe provided by the package. Every taxon in the input data frame must be included in the taxonomy (as a tip or node). 

Iterative pruning of the input matrix is performed by the `densify()` function, which can be modulated by several parameters (described in detail in the function documentation and the vignette). It returns a specially formatted tibble (a `densify_result` object), which describes the result of each densification step alongside some summary statistics. The function `rank_results()` ranks the densify results accoring to a specifiable, subjectively useful scoring function that recruits these summary statistics. The optimal sub-matrix according to the scoring function receives rank 1 and can be directly retrieved by the function `prune()`. `visualize()`, an alias of `plot()`, visually compares the quality scores between different pruning steps. The default scoring function used by `rank_results()`, `prune()` and `visualize()` maximizes the product of the number of available data points and the overall coding density, but it can be adjusted to include other measures and trade off their relative weight.

# Conclusions

The R package `densify` provides users with a flexible and explicit method to generate sub-matrices from an input matrix in a mathematically principled way. The package documents case examples using a standard sparse linguistic dataset (WALS) and the standard linguistic taxonomy provided by Glottolog.

Examples and further usage details for this software are found in the vignette hosted in the software repository on GitHub.

# Acknowledgements

The authors declare that there are no conflicts of interest.

# References
