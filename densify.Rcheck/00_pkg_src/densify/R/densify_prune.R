######################################################
# F4
#
# Function to actually densify the original input matrix, given the optimum quality score.
#
######################################################

#' Function to densify and prune the original input matrix
#'
#' This function takes the original data matrix and the documentation generated from the densification process. It densifies the matrix based on the optimum quality score and removes languages and variables that were pruned during the iterative process.
#'
#' @param original_data A data frame with the glottocodes as row names and variable names as column names.
#' @param documentation A data frame containing the documentation generated during the densification process.
#' @param optimum The iteration number representing the optimum quality score. The function will densify the matrix and prune languages/variables up to this iteration. Default is set to 1 (output is original data frame).
#'
#' @return A pruned and densified data frame.
#'
#' @seealso \code{\link{densify_steps}}
#'
#' @examples
#' # Assuming 'original_data' and 'documentation' are prepared
#' # pruned_densified <- densify_prune(original_data, documentation, optimum = 2)
#' @import dplyr
#' @export

densify_prune <- function(original_data, documentation, optimum = 1) {
  # Validate inputs:
  if (!is.data.frame(original_data)) {
    stop("The 'original_data' parameter must be a data frame.")
  }
  if (!is.data.frame(documentation)) {
    stop("The 'documentation' parameter must be a data frame.")
  }
  if (!is.integer(optimum) || optimum < 1) {
    stop("The 'optimum' parameter must be a strictly positive integer.")
  }

  # Ensure 'optimum' is not greater than the number of iterations in 'documentation':
  num_iterations <- nrow(documentation)
  if (optimum > num_iterations) {
    warning("The 'optimum' value was larger than the number of iterations in documentation.")
    optimum <- num_iterations
  }

  # Extract relevant documentation up to the optimum iteration:
  documentation <- dplyr::slice(documentation, 1:optimum)

  # Process removed languages and variables:
  prune_taxa <- unique(unlist(strsplit(documentation$removed_tax, ";")))[unique(unlist(strsplit(documentation$removed_tax, ";"))) != "NA"]
  prune_vars <- unique(unlist(strsplit(documentation$removed_var, ";")))[unique(unlist(strsplit(documentation$removed_var, ";"))) != "NA"]

  # Create the pruned matrix by removing specified languages and variables:
  pruned_matrix <- original_data[which(rownames(original_data)%in%prune_taxa==F),which(colnames(original_data)%in%prune_vars==F)]

  return(pruned_matrix)
}
