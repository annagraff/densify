rm(list=ls())
library(tidyverse)

########
# WALS data (dummy)
########

### create WALS dataframe (wals_data)
download.file("https://cdstar.shh.mpg.de/bitstreams/EAEA0-7269-77E5-3E10-0/wals_dataset.cldf.zip",
              destfile = "input/wals_dataset.cldf.zip")

# Read language and value information from the ZIP file:
wals_languages <- read.csv(unz("input/wals_dataset.cldf.zip", "languages.csv"),
                           header=TRUE, sep=",")
wals_values <- read.csv(unz("input/wals_dataset.cldf.zip", "values.csv"),header=TRUE, sep=",")
wals_parameters <- read.csv(unz("input/wals_dataset.cldf.zip", "parameters.csv"),header=TRUE, sep=",")

wals_parameters$Feature <- apply(wals_parameters, 1, function(x) paste("WALS_",x[1]," ",x[2],sep=""))
wals_parameters <- wals_parameters %>% select(c("ID","Feature"))

names(wals_languages)[1]<-"Language_ID"
names(wals_parameters)[1]<-"Parameter_ID"

wals_merged <- wals_languages %>% merge(wals_values, by = "Language_ID", all = TRUE) %>%
  merge(wals_parameters, by = "Parameter_ID", all = TRUE) %>%
  filter(Glottocode!="") %>%
  select(c("Glottocode","Feature","Value"))

wals_data<-pivot_wider(wals_merged,
                       names_from = Feature,
                       values_from = Value,
                       values_fill = "?",
                       values_fn = function(x){if (length(unique(x))==1){unique(x)} else ("?")})
##################################################################################################################################
##################################################################################################################################

######################################################
# F1
#
# Prepare glottolog info for all languages
#
# Our goal is to flatten out the glottolog tree for every languoid (glottocode) in the
# database and represent the node as linear levels. Here, the top-level items
# (roots of the glottolog tree) get mapped to node level 1 and the bottom-level items
# (leafs) are mapped to last node. Intermediate dummy nodes are created in between
# and filled with NA if necessary
#
# For true singletons (leaves which are also roots), a virtual root is created during
# the flattening process, e.g. Basque -> Basque
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

    # identify roots and leafs
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

    # now we want to normalize the paths — each path must be of the same lenght, with roots in front and
    # leafs in back
    # we insert dummy empty nodes if nessesary
    longest_path_length <- max(sapply(languoids$path, length))

    languoids <- mutate(languoids,
                        # path normalization
                        normalized_path = mapply(path, is_leaf, SIMPLIFY=FALSE, FUN=function(path, is_leaf) {
                          # if the path is already longes, there is nothing to do
                          if(length(path)==longest_path_length) return(path)

                          # if there is only one element, we need to insert a dummy root (if its  leaf) and a dummy leaf (if its a root)
                          if( (length(path) == 1)) path <- if(is_leaf) c(path, path) else c(path, 0)

                          # now insert dummy nodes 0 between the leaf and the upper part of the path
                          path <- c(path[-length(path)], integer(longest_path_length-length(path)), path[length(path)])

                          path
                        }),
                        # get the node level of the respective languoid
                        glottolog.node.level = mapply(normalized_path, id, FUN=function(path, id)  length(path) - match(id, rev(path)) + 1)
    )

    # map the paths to glottocode sand create the flattened out hierarchy table
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

### Test for F1
register <- glottocode_taxonomy(path_to_file="input/glottolog_4.7_languoid/languoid.csv")

##################################################################################################################################
##################################################################################################################################

######################################################
# F2
#
# Iterative matrix densification according to specified criteria.
#
# The output of this densification is a log-file, which specifies details about the matrix after
# each iteration.
######################################################

original_data <- wals_data
pruning_prep <- function(original_data) {
  # Replace NAs by 0 (no data available) and non-NA entries by 1 (data available)
  full_matrix <- full_matrix[full_matrix=="?"] <- NA
  full_matrix <- as.matrix(original_data)
  full_matrix <- full_matrix[full_matrix=="?"] <- NA  # ML: Why do we do this twice?
  full_matrix[!is.na(full_matrix)] <- 1
  full_matrix[is.na(full_matrix)] <- 0

  # This is currently emptying the entire matrix, so it's commented out.
  # # Remove the invalid glottocode row
  # full_matrix <- full_matrix[-which(rownames(full_matrix) == "NA."), ]

  # Convert dataframe entries to numeric
  full_matrix <- apply(full_matrix, 2, as.numeric)

  # Remove dimnames attribute from resulting matrix
  dimnames(full_matrix) <- NULL

  return(full_matrix)
}

# Call the function
result <- pruning_prep(original_data)

# Validate the output
expectedOutput <- matrix(c(0, 1, 1, 1, 1, 0), nrow = 3)
identical(result, expectedOutput)


