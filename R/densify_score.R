#' Quality score computation for each iteration
#'
#' Defines a quality score for each iteration as documented in the output created in the `densify_steps` function and identifies the optimal number of iterations post-hoc according to modifiable exponents.
#'
#' @param iteration_log A data frame containing the output iteration_log from the `densify_steps` function.
#'
#' @param exponent_prop_coded_data Exponent for the proportion of coded data. Default is 1.
#'
#' @param exponent_available_data_points Exponent for the number of available data points. Default is 1.
#'
#' @param exponent_lowest_taxon_coding_score Exponent for the lowest taxon score. This exponent is not mandatory, so there is no default.
#'
#' @param exponent_lowest_variable_coding_score Exponent for the lowest variable score. This exponent is not mandatory, so there is no default.
#'
#' @param exponent_taxonomic_diversity Exponent for the taxonomic diversity index. This exponent is not mandatory, so there is no default.
#'
#' @param plot Logical denoting whether the quality scores should be plotted. Default is `FALSE`.
#'
#' @return The optimal number of iterations with the highest quality score.
#'
#' @examples
#' # Assuming `iteration_log` is a data frame generated from `densify_steps` function
#' # densify_score(iteration_log)
#'
#' @seealso \code{\link{densify_steps}}
#'
#' @importFrom stats na.omit
#' @import ggplot2
#' @export

densify_score <- function(iteration_log, exponent_prop_coded_data = 1, exponent_available_data_points = 1, exponent_lowest_taxon_coding_score = NULL, exponent_lowest_variable_coding_score = NULL, exponent_taxonomic_diversity = NULL, plot = FALSE) {
  # Validate input parameters:
  if (!is.numeric(exponent_prop_coded_data) ||
      !is.numeric(exponent_available_data_points)) {
    stop("Exponents must be numeric.")
  }

  if (!is.null(exponent_lowest_taxon_coding_score)){
    if(!is.numeric(exponent_lowest_taxon_coding_score)){
      stop("Exponents must be numeric.")
    }
  }

  if (!is.null(exponent_lowest_variable_coding_score)){
    if(!is.numeric(exponent_lowest_variable_coding_score)){
      stop("Exponents must be numeric.")
    }
  }

  if(!is.null(exponent_taxonomic_diversity) & "taxonomic_index" %in% names(iteration_log == F)){
    stop("densify_steps() was carried out without providing a taxonomy matrix, so there cannot be an exponent for taxonomic diversity.")
  }

  if (is.null(exponent_taxonomic_diversity) & "taxonomic_index" %in% names(iteration_log)){
    warning("Attention, an exponent for taxonomic diversity would be applicable but was not set, and no default exponent applies.\n")
  }

  if (is.null(exponent_lowest_taxon_coding_score)){
    warning("Attention, an exponent for lowest taxon coding score was not set, and no default exponent applies.\n")
  }

  if (is.null(exponent_lowest_variable_coding_score)){
    warning("Attention, an exponent for lowest variable coding score was not set, and no default exponent applies.\n")
  }


  # Calculate the quality score for each iteration:
  quality <- as.numeric(iteration_log$prop_coded_data)^exponent_prop_coded_data *
    as.numeric(iteration_log$available_data_points)^exponent_available_data_points

  # if exponent_taxonomic_diversity is provided, include this in the quality score
  if(!is.null(exponent_lowest_taxon_coding_score)){
    quality <- (quality * as.numeric(iteration_log$worst_tax_abs_coding_density)^exponent_lowest_taxon_coding_score)
  }

  # if exponent_lowest_variable_coding_score is provided, include this in the quality score
  if(!is.null(exponent_lowest_variable_coding_score)){
    quality <- (quality * as.numeric(iteration_log$worst_var_abs_coding_density)^exponent_lowest_variable_coding_score)
  }

  # if exponent_taxonomic_diversity is provided, include this in the quality score
  if(!is.null(exponent_taxonomic_diversity)){
    quality <- (quality * as.numeric(iteration_log$taxonomic_index)^exponent_taxonomic_diversity)
  }

  # Plot the quality scores:
  if(plot == TRUE){
    gg<-ggplot(data.frame(iteration = 1:length(quality), quality_score = quality), aes(x = iteration, y = quality_score)) +
      geom_point(pch=1) +
      theme_bw() +
      labs(x = "Iteration", y = "Quality Score", title = "Quality Score at Each Iteration")
    print(gg)
  }

  # Identify the first iteration number with the highest quality score:
  optimum <- min(which(quality == max(na.omit(quality))))
  return(optimum)
}
