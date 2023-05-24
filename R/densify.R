#' pre-processing linguistic input data frame
#'
#' This function pre-processes a linguistic input data frame by replacing NAs with 0
#' (indicating no data available) and non-NA entries with 1 (indicating data available).
#'
#' @param OriginalData input data frame
#' @return prepped matrix
#' @export
densify_prep <- function(OriginalData) {
  # Replace NAs by 0 (no data available) and non-NA entries by 1 (data available)
  FullMatrix <- as.matrix(OriginalData)
  FullMatrix[!is.na(FullMatrix)] <- 1
  FullMatrix[is.na(FullMatrix)] <- 0

  # This is currently emptying the entire matrix, so it's commented out.
  # # Remove the invalid glottocode row
  # FullMatrix <- FullMatrix[-which(rownames(FullMatrix) == "NA."), ]

  # Convert dataframe entries to numeric
  FullMatrix <- apply(FullMatrix, 2, as.numeric)

  # Remove dimnames attribute from resulting matrix
  dimnames(FullMatrix) <- NULL

  return(FullMatrix)
}

