# Set up ininial pruning state
#
# The pruning state tracks currently pruned data, variable levels, taxonomy, and importance weights
#
# - data and taxonomy are matched by rows (i.e. i-th row in the data describes the
# i-th taxon in the taxonomy)
init_pruning_state <- function(data, vars, ids, taxonomy, ..., scoring_fn, density_mean_weights) {
  # build a data frame that only contains variables
  # row names track taxa ids
  var_frame <- new_data_frame(unclass(data)[vars])
  rownames(var_frame) <- ids

  # initial state
  state <- list(
    # the indicator matrix (row names are taxa ids, col names are variables)
    matrix = binary_indicator_matrix(var_frame),
    # the flat taxonomy for the active data (matrix of taxonomic levels)
    taxonomy = taxonomy,
    # data variable levels (for variability computation)
    data_levels = encode_unique_levels(var_frame),
    # row and column mapping to the original data (vars tracks mapping of variable columns to data columns)
    indices = list(rows = seq_len(nrow(data)), cols = seq_len(ncol(data)), vars = vars),
    # parameters
    params = list(
      scoring_fn = scoring_fn,
      density_mean_weights = density_mean_weights
    )
  )

  # calculate initial importance weights
  state$weights <- update_importance_weights_for_pruning_state(NULL, state)

  state
}

# Prune the rows and/or columns and return the new pruning state
#
#  - column indices refer to variable indices, not original data indices
prune_indices <- function(state, indices) {
  updated_families <- character(0)

  # elements to prune
  rows <- indices$index[indices$axis == 1L]
  cols <- indices$index[indices$axis == 2L]

  # prune rows
  if (length(rows) > 0L) {
    updated_families <- state$taxonomy[rows, 1L, drop = TRUE]

    # update the data state
    state$matrix <- vec_slice(state$matrix, -rows)
    if (!is.null(state$taxonomy)) state$taxonomy <- vec_slice(state$taxonomy, -rows)
    if (!is.null(state$data_levels)) state$data_levels <- vec_slice(state$data_levels, -rows)

    # update the weights (we do not update the scores, since they will be recalculated anyway)
    state$weights$rows$weights <- vec_slice(state$weights$rows$weights, -rows)

    # update the data row indices
    state$indices$rows <- state$indices$rows[-rows]
  }

  # prune columns
  if (length(cols) > 0L) {
    # update the data state
    state$matrix <- state$matrix[, -cols, drop = FALSE]
    if (!is.null(state$data_levels)) state$data_levels <- state$data_levels[-cols]

    # update the weights
    state$weights$cols$weights <- vec_slice(state$weights$cols$weights, -cols)

    # update the data column indices (cols need to be translated to original data columns via vars mapping)
    state$indices$cols <- state$indices$cols[-state$indices$vars[cols]]

    # update the vars column to data column mapping (prefix sum is used to update the indices)
    state$indices$vars <- (state$indices$vars - cumsum(which_mask(cols, length(state$indices$vars))))[-cols]
  }

  # update the weights and scores
  state$weights <- update_importance_weights_for_pruning_state(state$weights, state)

  state
}


# Prunes non-informative taxa and variables below the threshold
prune_non_informative_data <- function(state, check_taxa = TRUE, check_vars = TRUE, min_variability = NA_integer_, ..., .changes) {
  changes <- NULL

  # repeat until we prune everything we can
  while ((check_taxa || check_vars) && prod(dim(state$matrix)) > 0L) {
    # check rows (taxa) with no data coverage
    if (check_taxa  && length(indices <- which(rowSums(state$matrix) == 0)) > 0L) {
      # update the change list
      if (!missing(.changes)) {
        changes <- vec_rbind(changes, data_frame(
          axis = 1L,
          index = indices,
          type = "taxon",
          id = dimnames(state$matrix)[[1L]][indices],
          score = state$weights$rows$scores[indices],
          reason = "no data"
        ))
      }

      state <- prune_indices(state, data_frame(axis = 1L, index = indices))
      check_vars <- TRUE
    }
    check_taxa <- FALSE

    # check columns (vars) with no data
    if (check_vars  && length(indices <- which(colSums(state$matrix) == 0)) > 0L) {
      # update the change list
      if (!missing(.changes)) {
        changes <- vec_rbind(changes, data_frame(
          axis = 2L,
          index = indices,
          type = "variable",
          id = dimnames(state$matrix)[[2L]][indices],
          score = state$weights$cols$scores[indices],
          reason = "no data"
        ))
      }

      state <- prune_indices(state, data_frame(axis = 2L, index = indices))
    }

    if (check_vars && !is.na(min_variability) && length(indices <- which(calculate_variability(state$data_levels, 2L) < min_variability)) > 0L ) {
      # update the change list
      if (!missing(.changes)) {
        changes <- vec_rbind(changes, data_frame(
          axis = 2L,
          index = indices,
          type = "variable",
          id = dimnames(state$matrix)[[2L]][indices],
          score = state$weights$cols$scores[indices],
          reason = "low variability"
        ))
      }

      state <- prune_indices(state, data_frame(axis = 2L, index = indices))
      check_taxa <- TRUE
    }
    check_vars <- FALSE
  }

  if (!missing(.changes) && !is.null(changes)) .changes(changes)

  state
}

# Recalculates the importance weights for a pruning state
#
# - `updated_families` is used to optimize taxonomic diversity index calculation
#
# The importance weights consists of the following factors
#
# - absolute coding density (percentage of variables with data) for individual taxa (rows)
# - relative coding density (matrix-weighted) for individual taxa (rows)
# - [if requested] taxonomic diversity index  for individual taxa (rows)
# - relative coding density (matrix-weighted) for individual variables (columns)
update_importance_weights_for_pruning_state <- function(weights, state, updated_families = character(0)) {
  m <- state$matrix

  # data is fully pruned
  if (nrow(m) == 0 || ncol(m) == 0) {
    row_weights <- matrix(double(), nrow = nrow(m))
    col_weights <- matrix(double(), nrow = ncol(m))

    return(list(
      rows = list(weights = row_weights, scores = row_weights),
      cols = list(weights = col_weights, scores = col_weights)
    ))
  }

  # coding densities must be recalculated at every step
  row_absolute <- rowSums(m)/ncol(m)

  # we use the eigenvalues of the Gramian to calculate the degree of linear independence in the indicator matrix
  # this gives us the relative coding weights of rows
  gram <- crossprod(m)
  ev   <- eigen(gram)$vectors[, 1L]

  row_weighted <- (m %*% ev)/sum(ev)
  col_weighted <- (t(row_weighted) %*% m)/sum(row_weighted)

  # taxonomic weights are only computed if needed
  row_taxonomic <- if (state$params$density_mean_weights$taxonomy > 0) {
    # previous weights
    prev_weights <- if(!is.null(weights)) weights$rows$weights[, 3L, drop = TRUE]

    # only recalculate if something changed
    if (is.null(prev_weights) || length(updated_families) > 0L) {
      # top-level clade groups (families) for taxonomic diversity calculation
      # only recalculate updated families
      groups <- vec_group_loc(state$taxonomy[, 1L])
      group_indices <- if (is.null(weights)) groups$loc else groups$loc[groups$key %in% unique(updated_families)]

      # calculate new weights
      updated_weights <- calculate_taxonomic_diversity(state$taxonomy, group_indices)

      # update the taxonomic weights
      if(is.null(prev_weights)) {
        updated_weights
      } else {
        update_at <- which(!is.na(updated_weights))
        prev_weights[update_at] <- updated_weights[update_at]

        prev_weights
      }
    } else {
      prev_weights
    }
  }

  # scale the weights
  row_absolute  <- row_absolute*state$params$density_mean_weights$coding
  row_weighted  <- row_weighted*state$params$density_mean_weights$coding
  row_taxonomic <- row_taxonomic*state$params$density_mean_weights$taxonomy

  # new weights
  row_weights <- cbind(as.vector(row_absolute), as.vector(row_weighted), as.vector(row_taxonomic), deparse.level = 0L)
  col_weights <- cbind(as.vector(col_weighted), deparse.level = 0L)

  list(
    rows = list(weights = row_weights, scores = state$params$scoring_fn(row_weights)),
    cols = list(weights = col_weights, scores = as.vector(col_weights))
  )
}


# Pruning state statistics
get_pruning_state_stats <- function(state) {
  m <- state$matrix

  row_sums <- rowSums(m)
  col_sums <- colSums(m)
  total    <- sum(row_sums)
  row_prop <- row_sums/length(row_sums)
  col_prop <- col_sums/length(col_sums)


  stats <- if ((nrows <- nrow(m)) > 0L &&  (ncols <- ncol(m)) > 0L) {
    row_sums <- rowSums(m)
    col_sums <- colSums(m)
    total    <- sum(row_sums)
    row_prop <- row_sums/ncols
    col_prop <- col_sums/nrows

    list(
      n_data_points  = total,
      n_rows = nrows,
      n_cols = ncols,
      coding_density = total/(nrow(m)*ncol(m)),
      row_coding_density_min = min(row_prop),
      row_coding_density_median = stats::median(row_prop),
      row_coding_density_max = max(row_prop),
      col_coding_density_min = min(col_prop),
      col_coding_density_median = stats::median(col_prop),
      col_coding_density_max = max(col_prop)
    )
  } else {
    list(
      n_data_points  = 0,
      n_rows = 0,
      n_cols = 0,
      coding_density = NA_real_,
      row_coding_density_min = NA_real_,
      row_coding_density_avg = NA_real_,
      row_coding_density_max = NA_real_,
      col_coding_density_min = NA_real_,
      col_coding_density_avg = NA_real_,
      col_coding_density_max = NA_real_
    )
  }

  # Shannon diversity index for the top-level taxon
  stats$taxonomic_index <- if(!is.null(taxonomy <- state$taxonomy)) {
    counts <- if (nrow(taxonomy) > 0L) vec_count(taxonomy[,1])$count else integer()
    props <- counts/sum(counts)

    -sum(props * log(props))
  } else {
    NA_real_
  }

  new_data_frame(stats)
}


pruning_state_stat_vars <- names(get_pruning_state_stats(list(matrix = matrix())))

