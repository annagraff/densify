#' Flattens a taxonomy and encodes it in a table (with one column per taxonomic level)
#'
#' A taxonomy tree is said to be "flattened" if the depth of each of it's leaf nodes is equal to the
#' tree height. Such tree can be constructed by injecting dummy nodes where appropriate. Example:
#'
#'    A                              A      F''     level 1
#'   / \         dummy nodes        / \     |
#'   B  \       ------------->      B  E'   F'      level 2
#'  / \  \                         / \  \   |
#'  C  D  E  F                     C  D  E  F       leaf level
#'
#' A flat taxonomy can be trivially encoded in a tabular format as it has nodes at every level. Each
#' column encodes a taxonomic level. The main purpose of this representation is to query a leaf
#' parent at level k in O(1).
#'
#' Note that there can be multiple flattened representations for a given tree, since there it is not
#' guaranteed that nodes at different levels can be aligned in the unique way. In our implementation
#' we always align at the top level, which means that the dummy nodes are always inserted at the
#' bottom, between the leaf and it's closest parent in the original tree.
#'
#' @param id Vector of node ids
#' @param parent_id Vector of node parent ids (must be of same size and type as `id`)
#'
#' @return a data frame with one row per node id and as many columns as there are levels in the tree
#' @export

build_flat_taxonomy_matrix <- function(id, parent_id) {
  # TODO: argument checks
  length(id) == length(parent_id) || rlang::abort("`id` and `parent_id` must have same size")

  # we build the matrix by recursively joining last ids by parent
  mat <- list(id)

  repeat {
    # get the next level of parents
    next_level <- parent_id[vctrs::vec_match(mat[[1L]], id)]

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

  out
}





