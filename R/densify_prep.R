#' pre-processing linguistic input data frame
#'
#' This function pre-processes a linguistic input data frame by replacing NAs with 0
#' (indicating no data available) and non-NA entries with 1 (indicating data available).
#'
#' @param original_data input data frame
#' @return prepped matrix
#' @export
densify_prep <- function(original_data) {
  # Replace NAs by 0 (no data available) and non-NA entries by 1 (data available)
  full_matrix <- as.matrix(original_data)
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

