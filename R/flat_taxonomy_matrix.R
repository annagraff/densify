#' Flatten a taxonomy and encode it as a matrix with one column per taxonomic level
#'
#' This function converts a taxonomy tree (a phylo object or an adjacency table) into a flat
#' representation. This is used internally by [densify] to traverse the taxonomy tree. The
#' package user does not need to invoke this function directly, but we chose to expose it since it
#' can be useful for other applications.
#'
#' The input taxonomy can be provided either as a `phylo` object (e.g. via [ape::read.nexus]), as a
#' data frame that contains columns `id` and `parent_id` (each row encodes one parent-child
#' relationship), or by passing in `id` and `parent_id` vectors explicitly. The vectors `id` and
#' `parent_id` must be of type character.
#'
#' A taxonomy tree is said to be "flattened" if the depth of each of it's leaf nodes is equal to the
#' tree height. Such tree can be constructed by injecting dummy nodes where appropriate. Example:
#'
#' \preformatted{
#'    A                              A      F''     level 1
#'   / \         dummy nodes        / \     |
#'   B  \       ------------->      B  E'   F'      level 2
#'  / \  \                         / \  \   |
#' C   D  E  F                    C   D  E  F       leaf level
#' }
#'
#' A flat taxonomy can be trivially encoded in a tabular format as it has the full set of nodes at every
#' level. Each column encodes a taxonomic level. The main purpose of this representation is the
#' ability to find the parent of every leaf at every level in O(1).
#'
#' Note that there can be multiple flattened representations for a given tree, since it is not
#' guaranteed that nodes at different levels can be aligned in the unique way. In our
#' implementation we always align at the top level, which means that the dummy nodes are always
#' inserted at the bottom, between the leaf and it's closest parent in the original tree.
#'
#' @docType package
#'
#' @examples
#' as_flat_taxonomy_matrix(glottolog_languoids)
#'
#' @param x An object of class `phylo` (e.g. result of [ape::read.tree]) or a data frame with
#'   columns `id` and `parent_id`.
#'
#' @param id Character vector of node identifiers. This cannot be supplied at the same time as `x`.
#'
#' @param parent_id Character vector of parent identifiers. This cannot be supplied at the same
#'   time as `x`.
#'
#' @param .x Optional label for the argument `x` to use in error messages
#'
#' @return A data frame with one row per node id and as many columns as there are levels in the tree.
#'
#' @export
as_flat_taxonomy_matrix <- function(x, id, parent_id, .x = rlang::caller_arg(x)) {
  # process the arguments
  x <- check_as_flat_taxonomy_matrix_args(x, id = id, parent_id = parent_id, .x = .x)

  id <- x$id
  parent_id <- x$parent_id

  # we build the matrix by recursively joining last ids by parent
  mat <- list(id)

  next_level0 <- NULL

  repeat {
    # get the next level of parents
    next_level <- parent_id[vctrs::vec_match(mat[[1L]], id)]

    # check for recursion
    !identical(next_level0, next_level) || rlang::abort(c(
      "recursive taxonomy detected!",
      "!" = {
        nodes <- unique(next_level[!is.na(next_level)])
        nodes <- paste0("'", nodes, "'")
        if (length(nodes) > 5L) nodes <- c(nodes[1:4L], "...")
        sprintf("possible recursion in %s", paste0(nodes, collapse = ", "))
      }
    ))
    next_level0 <- next_level

    missing <- is.na(next_level)
    if (all(missing)) break

    # add the next parent level
    mat <- c(list(next_level), mat)

    # shift missing values by column to align the roots
    missing <- which(missing)
    for (i in seq_len(length(mat) - 1L)) {
       mat[[i]][missing] <- mat[[i + 1]][missing]
       mat[[i + 1]][missing] <- NA
    }
  }

  # and make this into a data frame
  out <- c(list(id), mat)
  names(out) <- c("id", paste0("level", seq_along(mat)))
  out <- vctrs::new_data_frame(out)

  # replace first level by nodes, if first level is empty
  if(nrow(out) > 0L && all(out$level1 == "")) {
    nms <- names(out)[-length(names(out))]
    out <- out[,-2]
    names(out) <- nms
  }

  out
}

# helper fucntion for parsing and checking the arguments to as_flat_taxonomy_matrix
check_as_flat_taxonomy_matrix_args <- function(x, id, parent_id, .x = rlang::caller_arg(x)) {

  # check arguments
  (missing(x) != any(missing(id), missing(parent_id))) &&
  (missing(id) == missing(parent_id)) || rlang::abort(c(
    "either `x` or both `id` and `parent_id` must be supplied",
    "x" = if(!missing(x)) {
      "`x` cannot be suplied at the same time as `id` or `parent_id`"
    } else if (missing(parent_id)) {
      "missing required argument `parent_id` (must be supplied with `id`)"
    } else if (missing(id)) {
      "missing required argument `id` (must be supplied with `parent_id`)"
    },
    "",
    "i" = "see ?as_flat_typology_matrix for details"
  ))

  # explicit arguments
  if (!missing(id) && !missing(parent_id)) {
    data <- vctrs::vec_recycle_common(id = id, parent_id = parent_id)

    is.factor(data$id) || is.character(data$id) || rlang::abort(c(
      "taxonomy data `id` must be a vector of character values",
      "!" = sprintf("<%s> found", class(data$id)[[1L]]),
      "",
      "i" = "see ?as_flat_taxonomy_matrix for details"
    ))

    is.factor(data$parent_id) || is.character(data$parent_id) || rlang::abort(c(
      "taxonomy data `parent_id` must be a vector of character values",
      "!" = sprintf("<%s> found", class(data$parent_id)[[1L]]),
      "",
      "i" = "see ?as_flat_taxonomy_matrix for details"
    ))

    data
  } else
  # a phylo object
  if (inherits(x, "phylo")) {
    labels <- c(x$tip.label, if(is.null(x$node.label)) paste0("node", seq_len(x$Nnode)), x$node.label)

    # add roots to the edge list
    edge <- x$edge
    roots <- setdiff(edge[, 1L], edge[, 2L])
    if(length(roots) > 0L) edge <- unique(rbind(cbind(NA_integer_, roots), edge))

    list(id = labels[edge[, 2L]], parent_id = labels[edge[, 1L]])
  } else
  if (rlang::is_bare_list(x) || is.data.frame(x)) {
    rlang::has_name(x, "id") &&
    rlang::has_name(x, "parent_id") || rlang::abort(c(
      sprintf("taxonomy data `%s` must contain named components `id` and `parent_id`", .x),
      "",
      "i" = "see ?as_flat_taxonomy_matrix for details"
    ))

    data <- vctrs::vec_recycle_common(id = x$id, parent_id = x$parent_id, .arg = .x)

    is.factor(data$id) || is.character(data$id) || rlang::abort(c(
      sprintf("taxonomy data `%s$id` must be a vector of character values", .x),
      "!" = sprintf("<%s> found", class(data$id)[[1L]]),
      "",
      "i" = "see ?as_flat_taxonomy_matrix for details"
    ))

    is.factor(data$parent_id) || is.character(data$parent_id) || rlang::abort(c(
      sprintf("taxonomy data `%s$parent_id` must be a vector of character values", .x),
      "!" = sprintf("<%s> found", class(data$parent_id)[[1L]]),
      "",
      "i" = "see ?as_flat_taxonomy_matrix for details"
    ))

    data
  } else {
    rlang::abort(c(
      sprintf("`%s` is expected to be a phylo object or data frame with columns `id` and `parent_id`", .x),
      "",
      "i" = "see ?as_flat_taxonomy_matrix for details"
    ))
  }
}


