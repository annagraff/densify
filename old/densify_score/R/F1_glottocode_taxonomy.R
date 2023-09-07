######################################################
# F1
#
#' glottocode_taxonomy
#'
#' Prepare glottolog info for all languages.
#'
#' Our goal is to flatten out the glottolog tree for every languoid (glottocode) in the
#' database and represent the node as linear levels. Here, the top-level items
#' (roots of the glottolog tree) get mapped to node level 1 and the bottom-level items
#' (leaves) are mapped to the last node. Intermediate dummy nodes are created in between
#' and filled with NA if necessary.
#'
#' For true singletons (leaves which are also roots), a virtual root is created during
#' the flattening process, e.g. Basque -> Basque
#'
#' @param path_to_file The path to the input CSV file containing glottolog information.
#'
#' @return A data frame with the flattened glottolog tree for each languoid, including node levels.
#'
#' @examples
#' # Assuming `path_to_file` is the path to the glottolog CSV file
#' glottocode_taxonomy(path_to_file)
#'
#' @import phytools
#' @import tidyverse
#' @import testthat
#' @import GetoptLong
#' @import readr
#' @importFrom dplyr select mutate filter
#'
#' @export
#'
#' @seealso \code{\link{phytools}}, \code{\link{tidyverse}}, \code{\link{testthat}},
#' \code{\link{GetoptLong}}, \code{\link{readr}}
#'
#' @references Add references here (if available).
#'
#' @keywords glottolog, taxonomy, flattening, nodes, leaves, data frame
#'
#' @family Linguistic Data Processing
######################################################

glottocode_taxonomy <- function(path_to_file){
  local({
    library(phytools)
    library(tidyverse)
    library(testthat)
    library(GetoptLong)
    library(readr)
    # load the list of glottolog languoids
    languoids <-  read.csv(path_to_file, na.strings=c("NA", "", "na")) %>%
      select(id=id, parent_id=parent_id, glottolog.isocode = iso639P3code, glottocode=id, glottolog.name =name, longitude=longitude, latitude=latitude)

    expect_false(any(duplicated(languoids$id)), info="Glottolog contains duplicate nodes")
    expect_false(any(duplicated(languoids$glottocode)), info="Glottolog contains duplicate glottocodes")

    # identify roots and leaves
    languoids <- mutate(languoids, is_root = is.na(parent_id), is_leaf = id %in% setdiff(id, parent_id))

    # the tree is described in the table as an adjacency list
    # we want to record the path from the root to every node
    # this is done recustively, in the following visitor function
    languoids$path <- list(integer())

    visit_node_constructing_paths <- function(node_id, path) {
      # identify the table row
      node_row <- match(node_id, languoids$id)

      expect(!is.na(node_row), qq("Node @{node_id} is absent from the table"))
      expect(identical(languoids$path[[node_row]], integer()), qq("Node @{node_id} already visited!"))

      # construct and record the path
      path <- c(path, node_id)
      languoids$path[[node_row]] <<- path

      # exit if leaf
      if(languoids$is_leaf[node_row]) return()

      # visit all the children
      child_ids <- filter(languoids, parent_id ==node_id)$id

      expect_true(length(child_ids)>0, info=qq("Intermediate node @{node_id} must have at least one child"))

      # and apply tree construction recursively
      for(child_id in child_ids) visit_node_constructing_paths(child_id, path)
    }

    # visit all roots
    for(root_id in languoids$id[languoids$is_root]) visit_node_constructing_paths(root_id, integer())

    expect_false(any(sapply(languoids$path, length)==0), info="Some languoids not visited during graph traversal")

    # now we want to normalize the paths — each path must be of the same length, with roots in front and leaves in back
    # we insert dummy empty nodes if necessary
    longest_path_length <- max(sapply(languoids$path, length))

    languoids <- mutate(languoids,
                        # path normalization
                        normalized_path = mapply(path, is_leaf, SIMPLIFY=FALSE, FUN=function(path, is_leaf) {
                          # if the path is already longest, there is nothing to do
                          if(length(path)==longest_path_length) return(path)

                          # if there is only one element, we need to insert a dummy root (if it's leaf) and a dummy leaf (if it's a root)
                          if( (length(path) == 1)) path <- if(is_leaf) c(path, path) else c(path, 0)

                          # now insert dummy nodes 0 between the leaf and the upper part of the path
                          path <- c(path[-length(path)], integer(longest_path_length-length(path)), path[length(path)])

                          path
                        }),
                        # get the node level of the respective languoid
                        glottolog.node.level = mapply(normalized_path, id, FUN=function(path, id)  length(path) - match(id, rev(path)) + 1)
    )

    # map the paths to glottocodes and create the flattened out hierarchy table

    flat_tree <- languoids$normalized_path %>%
      sapply(function(path)  languoids$glottocode[match(path, languoids$id)]) %>%
      t() %>%
      as.data.frame(stringsAsFactors=FALSE)

    names(flat_tree) <- paste0("glottolog.node", 1:ncol(flat_tree))

    # and combine this with the languoid infos
    languoids <-
      select(languoids, glottocode, glottolog.name, longitude, latitude) %>%
      cbind(flat_tree) %>%
      filter(!is.na(glottocode))

    languoids
  })
}

