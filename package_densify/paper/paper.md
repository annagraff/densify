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
citation_author: Graff et al et. al.
date: 25 August 2023
year: 2023
bibliography: sources.bib
output: rticles::joss_article
csl: apa.csl
journal: JOSS

---

# Summary

The R package ``densifieR`` takes a data frame as an input, and generates a pruned data frame consisting of a subset of rows and columns. The pruning process follows criteria that the user can adapt to best fit her needs, including data retention, codedness and taxonomic diversity of rows (if the entities are represented in a taxonomic strucutre). As such, the software is useful for several purposes, including the densification of sparse input matrices and the subsampling of large input matrices according to a procedure that is sensitive to taxonomic structure.

# Statement of Need

While the software will run on any data frame (with rows representing any entities with or without taxonomic structure), we believe it will be particularly useful for handling data frames of typological linguistic data. Linguistic data is increasingly available in large-scale databases, and many analyses that aim at testing hypotheses at global scales, for which data collection may not be feasible, resort to such databases. Some of these resources have complete or near-complete variable coding density for all languages data is available for, but may be too large for certain computationally intensive analyses (e.g. PHOIBLE, Grambank). Other databases (e.g. WALS, AUTOTYP, Lexibank) exhibit variables that are coded for very different sets of languages, resulting in sparse or even extremely sparse language-variable matrices. Combining data from various databases via language identifiers like glottocodes usually results in further sparsity. Such sparse matrices may lack power for statistical analyses. In either case, researchers might be particularly interested in maintaining taxonomic diversity in the languages represented in a sub-matrix, penalizing the removal of language isolates or members of small language families more strongly than languages belonging to clades represented by many other languages in the matrix.

To the best of our knowledge, there do not yet exist any principled approaches to generate sub-matrices from sparse input matrices.


Many contemporary applications and software packages are optimized for large-scale networks. For example, ``igraph`` [@R-igraph], ``network`` [@R-networkpkg], and ``sna`` [@R-sna] were developed to analyze social media networks [@Jones_2017], epidemiological networks [@Christakis_2011], and political networks [@Hobbs_2016], respectively. In contrast, discourse networks in educational contexts are substantially smaller, typically with only 3-8 students [@Wagner_2018]. Consequently, certain parameters that are relevant for these larger networks are not necessarily applicable, and analysis of discourse networks demands additional parameters beyond what is available in graph theory [@Lou_2001; @Chai_2019].

# Usage

``densifieR`` provides the data from The World Atlas of Language Strucutres (WALS; CITATION) and the language taxonomy provided by Glottolog v. 4.8 (CITATION) as example data. This data is used to demonstrate the utility and flexibility of ``densifieR`` to subsample an input matrix according to varying needs.

## Preparing input

An ``igraph`` object is the core input to many of the modular analytical functions offered in ``discourseGT``. Prior to generating an ``igraph`` object, a weighted edge list needs to be generated from imported raw data, structured as two columns containing sequential nodes or individual students who start or continue a discussion episode (@Chai_2019). This is addressed by the `tabulate_edges` function. By default, the weight of an edge is defined as the number of times an edge has occurred between two nodes. Weights can be redefined based on other available criteria, but this must be done manually. 

```r
install.packages("densifieR")
library(densifier)

# Prepare data: WALS and glottolog
data(wals)
data(glottolog_languoids)

# the input data frame must be a data frame with the taxon names as row names and variable names as column names; any question marks, empty entries, "NA"s must be coded as NAs
wals[wals=="?"] <- NA
wals[wals=="NA"] <- NA
head(wals)

# Generate flat taxonomy using glottolog
taxonomy_matrix <- build_flat_taxonomy_matrix(glottolog_languoids$id, glottolog_languoids$parent_id)
head(taxonomy_matrix)
```
## Pruning
Iterative pruning of the input matrix is performed by the `densify_steps` function, which requires the following information:

  *	The original data frame. Default: wals
  *	An integer specifying the number of iterations performed. Default: 1
  *	The mean type, which can be: "arithmetic", "geometric", "log_odds". Default: "log_odds".
  *	Is taxonomy accounted for in the pruning process? Default: `FALSE`
  *	The taxonomy matrix, required if taxonomy = `TRUE`, optional otherwise? Default: NULL.
  *	The tax_weight_factor, required if mean_type = "log_odds". Default: 0.99
  *	The coding_weight_factor, required if mean_type = "log_odds". Default: 0.99

For a more detailed discussion of the parameters, refer to the vignette hosted in the software repository.
```r
set.seed(2023)
documentation <- densify_steps(original_data = wals, max_steps = (nrow(wals)+ncol(wals), mean_type = "log_odds", taxonomy = TRUE, taxonomy_matrix = taxonomy_matrix, tax_weight_factor = 0.99, coding_weight_factor = 0.99)
head(documentation)
```

The graph settings specified by `prepareGraphs` will influence the analytical output of downstream functions.

## Retrieving the optimal number of iterations and producing the sub-matrix
``discourseGT`` offers graph theory-based analytics via two separate functions: `coreNetAnalysis()` and `subgroupsNetAnalysis()`. 

`coreNetAnalysis()` will perform core graph theory operations, such as the counting number of nodes and edges and calculating edge weights, average graph degree, centrality, and other graph theory parameters [@Chai_2019].

```r
coreNet <- coreNetAnalysis(prepNet)
```

# Conclusions
The R package ``densifieR`` provides users a flexible and explicit method to generate sub-matrices from an input matrix in a mathematically principled way. This paper shows a case example using a standard sparse linguistic dataset (WALS) and the standard linguistic taxonomy provided by Glottolog. More examples and usage details are hosted in the software repository on GitHub.

# Acknowledgements
We thank Reinhard Furrer for his input on several package functions and usability.
FINANCIAL SUPPORT??
The authors declare that there are no conflicts of interest.

# References
A list of key references, including to other software addressing related needs. Note that the references should include full names of venues, e.g., journals and conferences, not abbreviations only understood in the context of a specific discipline.

