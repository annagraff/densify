######################################################
# F4
#
# Function to actually densify the original input matrix, given the optimum quality score.
#
######################################################

densify_prune <- function(original_data, documentation, optimum = 1) {
  # Validate inputs
  if (!is.data.frame(original_data)) {
    stop("The 'original_data' parameter must be a data frame.")
  }
  if (!is.data.frame(documentation)) {
    stop("The 'documentation' parameter must be a data frame.")
  }
  if (!is.integer(optimum) || optimum < 1) {
    stop("The 'optimum' parameter must be a strictly positive integer.")
  }

  # Ensure 'optimum' is not greater than the number of iterations in 'documentation'
  num_iterations <- nrow(documentation)
  if (optimum > num_iterations) {
    warning("The 'optimum' value is larger than the number of iterations in documentation.")
    optimum <- num_iterations
  }

  # Extract relevant documentation up to the optimum iteration
  documentation <- slice(documentation, 1:optimum)

  # Process removed languages and variables
  prune_lgs <- unique(unlist(strsplit(documentation$removed_lg, ";")))[unique(unlist(strsplit(documentation$removed_lg, ";"))) != "NA"]
  prune_vars <- unique(unlist(strsplit(documentation$removed_var, ";")))[unique(unlist(strsplit(documentation$removed_var, ";"))) != "NA"]

  # Create the pruned matrix by removing specified languages and variables
  pruned_matrix <- original_data[which(rownames(original_data)%in%prune_lgs==F),which(colnames(original_data)%in%prune_vars==F)]

  return(pruned_matrix)
}
