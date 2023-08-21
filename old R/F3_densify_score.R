######################################################
# F3
#
# Function to define a quality score for each iteration as documented in the output created in F2 and identify the optimal number of iterations post-hoc. Exponents can be modified.
# The function also yields a plot of the quality score at each iteration.
#
######################################################

#' densify_score
#'
#' Function to define a quality score for each iteration as documented in the output created in F2 and to identify the optimal number of iterations post-hoc. Exponents can be modified.
#'
#' @param documentation A data frame containing the output documentation from the `densify_steps` function.
#'
#' @param exponent_prop_coded_data Exponent for the proportion of coded data. Default is 1.
#' @param exponent_available_data_points Exponent for the number of available data points. Default is 1.
#' @param exponent_lowest_taxon_coding_score Exponent for the lowest language score. Default is 1.
#' @param exponent_taxonomic_diversity Exponent for the taxonomic diversity index. Default is 1.
#'
#' @return The optimal number of iterations with the highest quality score.
#'
#' @examples
#' # Assuming `documentation` is a data frame generated from `densify_steps` function
#' densify_score(documentation)
#'
#' @export

densify_score <- function(documentation, exponent_prop_coded_data = 1, exponent_available_data_points = 1, exponent_lowest_taxon_coding_score = 1, exponent_taxonomic_diversity = 1) {
  # Validate input parameters:
  if (!is.numeric(exponent_prop_coded_data) ||
      !is.numeric(exponent_available_data_points) ||
      !is.numeric(exponent_lowest_taxon_coding_score) ||
      !is.numeric(exponent_taxonomic_diversity)) {
    stop("Exponents must be numeric.")
  }

  # Calculate the quality score for each iteration:
  quality <- as.numeric(documentation$prop_coded_data)^exponent_prop_coded_data *
    as.numeric(documentation$available_data_points)^exponent_available_data_points *
    as.numeric(documentation$worst_lg_abs_coding_density)^exponent_lowest_taxon_coding_score *
    as.numeric(documentation$taxonomic_index)^exponent_taxonomic_diversity

  # Plot the quality scores:
  plot(quality, xlab = "Iteration", ylab = "Quality Score", main = "Quality Score at Each Iteration")
  legend("topright", legend = "Quality Score", col = "black", lty = 1, cex = 0.8)

  # Identify the first iteration number with the highest quality score:
  optimum <- min(which(quality == max(na.omit(quality))))
  return(optimum)
}
