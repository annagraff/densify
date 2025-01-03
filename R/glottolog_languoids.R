#' Glottolog Languoid Taxonomy, Version 5.0
#'
#' Glottolog Languoid Taxonomy, Version 5.0
#'
#' A data frame containing linguistic taxonomic information from Glottolog version 5.0.
#' You can use this dataset as the `taxonomy` argument for the [densify()] function if your
#' data relies on glottocodes and you wish to perform taxonomy-aware data pruning.
#'
#' The dataset was extracted from the official Glottolog CLDF release and converted into a
#' wide tabular representation. Each row describes a single variety (Glottolog languoid).
#' The columns `id`, `parent_id`, and `family_id` correspond to the Glottolog code of the
#' languoid, its parent, and its root node, respectively. The `level` column denotes the
#' languoid classification (e.g., dialect, language, or family). Singletons and top-level
#' families have `parent_id` set to `NA`.
#'
#'
#' @name glottolog_languoids
#' @docType data
#' @source Glottolog, Version 5.0
#' @references
#' - Hammarström, Harald & Forkel, Robert & Haspelmath, Martin & Bank, Sebastian. 2023.
#' Glottolog 5.0.
#' Leipzig: Max Planck Institute for Evolutionary Anthropology.
#' https://doi.org/10.5281/zenodo.10804357
#' (Available online at http://glottolog.org, Accessed on 2023-03-13.)
#'
NULL
