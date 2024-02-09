# Calculate the frequency of the nth most frequent value for all variables
calculate_variability <- function(vars, n) {
  vapply(unname(vars), function(var) {
    if (!is.integer(var)) {
      var <- vec_slice(var, !vec_detect_missing(var))
      var <- vec_group_id(var)
    }
    counts <-  sort(tabulate(var), decreasing = TRUE)

    if (n > length(counts)) 0L else counts[[n]]
  }, 0L)
}

# Replace values in the table by their unique level ids
encode_unique_levels <- function(data) {
  out <- list()
  ids0 <- vec_init(integer(), nrow(data))

  # replace all columns by vectors of their unique levels
  for (coldata in unclass(data)) {
    ids <- ids0

    slice <- which(!vec_detect_missing(coldata))
    ids[slice] <- vec_group_id(vec_slice(coldata, slice))

    out <- c(out, list(ids))
  }

  names(out) <- paste0(rep(".", length(out)), seq_along(out))
  new_data_frame(out)
}


# A streamlined apply(mat, 1L, fn)
rowapply <- function(matrix, fn, ..., .ptype = list()) {
  out <- vec_init(.ptype, nrow(matrix))
  for (i in seq_len(nrow(matrix))) out[[i]] <- fn(matrix[i, ])
  out
}

# A which.min() that returns all indices
# which_min_all <- function(x) {
#   which(x == min(x))
# }


# Computes a binary indicator matrix (1 for presence, 0 for absence of data) for the dataset
binary_indicator_matrix <- function(data, ...) {
  # save the data frame dimension
  dim <- c(nrow(data), ncol(data))

  # detect all missing entries
  is_missing <- unlist(lapply(unname(data), vec_detect_missing))

  # construct a binary indicator matrix
  structure(1 - is_missing, dim = dim, dimnames = dimnames(data))
}

# Replaces each value with how often it occurs in the data
calculate_counts <- function(x) {
  dims <- dim(x)
  dim(x) <- NULL

  groups <- vec_group_loc(x)
  groups$count <- lengths(groups$loc)
  groups$count[is.na(groups$key)] <- 1L

  out <- list_unchop(as.list(groups$count), indices = groups$loc, ptype = integer())
  dim(out) <- dims

  out
}

# A variable setter
setter <- function(var, ..., .compose = function(old, new) new) {
  env <- parent.frame()
  var <- as.character(substitute(var))

  function(value) env[[var]] <- .compose(env[[var]], value)
}

type_label <- function(x) {
  if (is.data.frame(x)) class(x)[[L]] else vec_ptype_full(x)
}


# returns a boolean mask corresponding to the indices (reverse of which)
which_mask <- function(idx, n = max(idx)) {
  mask <- logical(n)
  mask[idx] <- TRUE

  mask
}

# check that the value is not NULL or NA
is_value <- function(x) !is_null(x) && !is.na(x)

# x is between a and b
is_between <- function(x, a, b) {
  is_bare_numeric(x) && length(x) == 1L && is_true(x >= a && x <= b)
}

