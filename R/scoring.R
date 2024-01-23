
# Select the rows or columns with the lowest index
identify_lowest_scores <- function(weights) {
  row_min <- min(weights$rows$scores)
  col_min <- min(weights$cols$scores)


  # rows and columns with smallest scores (return both if they are equal)
  rows <- if(row_min <= col_min) which(row_min == weights$rows$scores) else integer(0)
  cols <- if(col_min <= row_min) which(col_min == weights$cols$scores) else integer(0)
  nrows <- length(rows)
  ncols <- length(cols)

  vctrs::new_data_frame(list(
    axis  = c(rep(1L, length(rows)), rep(2L, length(cols))),
    index = c(rows, cols),
    score = c(weights$rows$scores[rows], weights$cols$scores[cols])
  ))
}




# Calculates the taxonomic diversity index
#
# The index calculation is per group (usually the top-level clade, i.e. family), with the parameter
# `group_indices` mapping groups to leaf taxa. This parameter can also be used to limit
# the calculation to certain groups only (we use it to speed up updates)
calculate_taxonomic_diversity <- function(taxonomy, group_indices) {
  # no groups, nothing to compute
  if (length(group_indices) == 0L) return(vctrs::vec_init(double(), nrow(taxonomy)))

  # replace each node by how often it occurs in the sample
  # these counts form the basis for the diversity measure (we prioritize retaining data from smaller branches)
  #
  # uneven taxonomies are equalized by counting NA nodes in the flat taxonomy as unique
  counts <- calculate_counts(taxonomy[unlist(group_indices), , drop = FALSE])

  # balance the taxonomic structure by only considering the unique levels
  prods <- rowapply(counts, function(x) prod(unique(x)), .ptype = double())

  # the taxonomic weights for each taxon are reciprocal geometric means of the (balanced) node level counts
  # note: fractional exponent effectively pads each product to the full taxonomic depth
  weights <- 1/(prods^(1/ncol(taxonomy)))

  # scatter the weights back into the correct rows using group indices
  weights <- vctrs::vec_assign(vctrs::vec_init(double(), nrow(taxonomy)), unlist(group_indices), weights)

  # normalize the weights within individual groups
  for (group in group_indices) weights[group] <- weights[group]/sum(weights[group])

  weights
}


# Row-wise score using arithmetic mean
row_scores_arithmetic <- function(matrix) {
  rowMeans(matrix)
}

# Row-wise score using geometric mean
row_scores_geometric <- function(matrix) {
  ncol(matrix) > 1L || return(matrix)

  # there is no rowProds
  prod <- matrix[, 1L, drop = TRUE]
  for(i in seq_len(ncol(matrix) - 1L)) {
    prod <- prod*matrix[, i+1, drop = TRUE]
  }

  # geometric mean
  prod^(1/ncol(matrix))
}


# Row-wise score using log odds
row_scores_log_odds <- function(matrix) {
  ncol(matrix) > 1L || return(matrix)

  # clamp the value range
  matrix[matrix > 0.99] <- 0.99
  matrix[matrix < 0.01] <- 0.01


  # logistic transformation
  matrix <- matrix(qlogis(matrix), ncol = ncol(matrix), nrow = nrow(matrix))
  means <- rowMeans(matrix)
  scores <- exp(means)/(1+exp(means))

  scores
}
