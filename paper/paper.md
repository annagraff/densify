---
title: "densifieR: An R package to prune sparse data frames of typological linguistic data"
tags:
  - sparse matrices
  - sub-sampling
  - linguistic data
  - diversity samples
  - R
authors:
  - name: Anna Graff
    orcid: 0000-0002-7703-3471
    equal-contrib: true
    corresponding: true 
    affiliation: "1, 2, 3"
  - name: Marc Lischka
    orcid: XXXXXXXX
    equal-contrib: true
    affiliation: "3, 4"
  - name: Taras Zakharko
    orcid: XXXXX
    affiliation: "1, 3"
  - name: Reinhard Furrer
    orcid: XXXXX
    affiliation: "3, 4"
  - name: Balthasar Bickel
    orcid: XXXXX
    affiliation: "1, 3"

affiliations:
 - name: 'University of Zurich, Department of Comparative Language Science'
   index: 1
 - name: 'University of Zurich, Department of Evolutionary Biology and Environmental Studies'
   index: 2
 - name: 'University of Zurich, Center for the Interdisciplinary Study of Language Evolution'
   index: 3
 - name: 'University of Zurich, Institute of Mathematics'
   index: 4
citation_author: Graff et al.
date: 29 August 2023
year: 2023
bibliography: sources.bib
output: rticles::joss_article
csl: apa.csl
journal: JOSS

---

# Summary

The R package ``densify`` takes a data frame as an input, and generates a pruned data frame consisting of a subset of rows and columns. The pruning process follows criteria that the user can adapt to best fit her needs, including data retention, codedness and taxonomic diversity of rows (if the entities are represented in a taxonomic structure). As such, the software is useful for several purposes, including the densification of sparse input matrices and the subsampling of large input matrices according to a procedure that is sensitive to taxonomic structure.

# Statement of Need

While the software will run on any data frame (with rows representing any entities with or without taxonomic structure), it was designed to prune data frames of typological linguistic data. Linguistic data is increasingly available in large-scale databases, and many analyses that aim at testing hypotheses at global scales, for which project-specific data collection may not be feasible, resort to such databases. Some of these resources have complete or near-complete variable coding density for all languages data is available for, but may be too large for certain computationally intensive analyses (e.g. PHOIBLE [@phoible], Grambank [@grambank]). Other databases (e.g. WALS [@wals], AUTOTYP [@AUTOTYP], Lexibank [@Lexibank]) exhibit variables that are coded for very different sets of languages, resulting in sparse language-variable matrices. Combining data from various databases via language identifiers like glottocodes usually results in further sparsity. Although there are methods available to operate on such matrices (e.g. ``Matrix`` [@R-Matrix], ``MatrixExtra`` [@R-MatrixExtra]), sparse matrices may lack power for certain research questions, such that researchers may prefer to operate on a subset of the data represented in a denser matrix. When generating subsets of a language-variable matrix, researchers might be particularly interested in maintaining taxonomic diversity in the languages represented, penalizing the removal of language isolates or members of small language families more strongly than languages belonging to clades represented by many other languages in the sample.

While certain packages exist to generate sub-matrices from varying input matrices according to principled criteria (e.g. ``admmDensestSubmatrix`` [@R-admmDensestSubmatrix], which identifies the densest sub-matrix of an input graph of a specified size; or ``FSelector`` [@R-FSelector], which performs attribute subset selection based on various tests and entropy measures to identify the most relevant attributes of a data input), ``densify`` provides a flexible, explicit and taxonomy-sensitive pruning algorithm that focuses both on the removal of rows and columns and does not require specification of the submatrix.

# Usage

``densify`` provides the data from The World Atlas of Language Strucutres (WALS [@wals]) and the language taxonomy provided by Glottolog v. 4.8 [@glottolog] as example data. The vignette features a detailed demonstration of the utility and flexibility of ``densify`` to subsample an input matrix according to varying needs using this data.

## Preparing input

``densify`` requires the input data frame to be a data frame with taxon names as row names and variable names as column names. Any cells with empty entries, not applicable or question marks must be coded as NA. If matrix densification should consider taxonomic structure, a flat taxonomy must be provided, listing every taxon present in the initial data frame along with all the nodes connecting it to the root. Such a taxonomy can be generated with the `build_flat_taxonomy_matrix` function, if all nodes and tips are provided alongside each of their parent nodes. For generating a language taxonomy, glottolog can be used directly.

```r
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
taxonomy\_matrix <- build\_flat\_taxonomy\_matrix(id = glottolog\_languoids$id, parent\_id = glottolog\_languoids$parent\_id)
head(taxonomy\_matrix)
```
## Pruning
Iterative pruning of the input matrix is performed by the `densify_steps` function, which requires the following information:

  *	The original data frame. Default: wals
  *	An integer specifying the number of iterations performed. Default: 1
  *	An integer specifying the threshold for variability in columns. Default: 1
  *	The mean type, which can be: "arithmetic", "geometric", "log\_odds". Default: "log\_odds".
  *	Is taxonomy accounted for in the pruning process? Default: `FALSE`
  *	The taxonomy matrix, required if taxonomy = `TRUE`, optional otherwise. Default: NULL.
  *	The tax\_weight\_factor, required if mean\_type = "log_odds". Default: 0.99
  *	The coding\_weight\_factor, required if mean\_type = "log_odds". Default: 0.99

For a more detailed discussion of the parameters, refer to the vignette hosted in the software repository.

```r
set.seed(2023)
documentation <- densify\_steps(original\_data = wals, max\_steps = nrow(wals)+ncol(wals)-2, variability_threshold=3, mean\_type = "log\_odds", taxonomy = TRUE, taxonomy\_matrix = taxonomy\_matrix, tax\_weight\_factor = 0.99, coding\_weight\_factor = 0.99)
head(documentation)
```

## Retrieving the optimal number of iterations and producing the sub-matrix
Given the pruning documentation output, `densify_score` will identify the optimal sub-matrix via a quality score, computed via user-defined exponents relating to the overall proportion of coded data in the matrix, the number of available data points, the coding density of the least well-coded taxon, the coding density of the least well-coded variable, and a taxonomic diversity index (Shannon diversity of the highest taxonomic level).

```r
exponent\_prop\_coded\_data <- 1
exponent\_available\_data\_points <- 1
exponent\_lowest\_taxon\_coding\_score <-1
exponent\_lowest\_variable\_coding\_score <-1
exponent\_taxonomic\_diversity <- 1

optimum <- densify_score(documentation = documentation, 
                         exponent\_prop\_coded\_data = exponent\_prop\_coded\_data, 
                         exponent\_available\_data\_points = exponent\_available\_data\_points, 
                         exponent\_lowest\_taxon\_coding\_score = exponent\_lowest\_taxon\_coding\_score,
                         exponent\_lowest\_variable\_coding\_score = exponent\_lowest\_variable\_coding\_score,
                         exponent\_taxonomic\_diversity = exponent\_taxonomic\_diversity)
                         
pruned_wals <- densify_prune(wals, documentation, optimum)
```

# Conclusions
The R package ``densify`` provides users a flexible and explicit method to generate sub-matrices from an input matrix in a mathematically principled way. This paper shows a case example using a standard sparse linguistic dataset (WALS) and the standard linguistic taxonomy provided by Glottolog. More examples and usage details are hosted in the software repository on GitHub.

# Acknowledgements
The authors declare that there are no conflicts of interest.

# References

