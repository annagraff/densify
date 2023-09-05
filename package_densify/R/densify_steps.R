######################################################
# F2
#
#' Iterative matrix densification according to specified criteria.
#'
#' The output of this densification is a log-file, which specifies details about the matrix after
#' each iteration.
#'
#' @param original_data A data frame with the taxa (e.g. languages, as glottocodes) as row names and variable names as column names.
#'   Any question marks, empty entries, or "NA"s must be represented as NAs. Default is the data from WALS.
#'
#' @param max_steps An integer specifying the maximum number of iterations attempted during densification.
#'   Default is 1.
#'
#' @param variability_threshold specifying how many taxa the second-most-frequent variable state of any variable must contain for the variable to be maintained in the matrix.
#'   Default is 1.
#'
#' @param mean_type A character string specifying the type of mean to be used for calculating the final weights.
#'   Possible values are "arithmetic", "geometric", or "log_odds". Default is "log_odds".
#'
#' @param taxonomy A logical indicating whether taxonomic diversity should be factored in for the pruning/retention of rows.
#'   If `TRUE`, taxonomic diversity is considered; if `FALSE`, it is not. Default is `FALSE`.
#'
#' @param taxonomy_matrix A data frame containing a taxonomy matrix, which can be retrieved from the `build_flat_taxonomy_matrix` function or provided directly.
#'   This parameter must be specified if `taxonomy = TRUE` and is optional otherwise (if specified, taxonomic index is computed for each iteration).
#'
#' @param tax_weight_factor A numeric value between 0 and 1 that determines the relative weight given to taxonomy in taxon pruning.
#'   This parameter must be specified if `taxonomy = TRUE` and `mean_type = "log_odds"`. Default is 0.99.
#'
#' @param coding_weight_factor A numeric value between 0 and 1 that determines the relative weight given to coding quality
#'   (absolute coding score and weighted coding score) in taxon pruning. This parameter must be specified if
#'  `mean_type = "log_odds"`. Default is 0.99.
#'
#' @return A data frame with details about the matrix after each iteration of densification.
#'
#' @examples
#' # Assuming `original_data` and `taxonomy_matrix` are appropriate data frames
#' # densify_steps(original_data = wals,
#' # max_steps = 3,
#' # variability_threshold = 3,
#' # mean_type = "log_odds",
#' # taxonomy = T,
#' # taxonomy_matrix,
#' # tax_weight_factor = 0.99,
#' # coding_weight_factor = 0.99)
#'
#' @import vegan
#' @import tidyverse
#' @importFrom stats qlogis var

#' @export
######################################################

densify_steps <- function(original_data = wals, max_steps = 1, variability_threshold = 1, mean_type = "log_odds", taxonomy = F, taxonomy_matrix, tax_weight_factor = 0.99, coding_weight_factor = 0.99){

  if (taxonomy == F){
    warning("Attention, taxonomy disregarded for pruning.\n")
  }

  if (is.data.frame(taxonomy_matrix) == F){
    warning("Attention, no taxonomy matrix has been provided.\n")
  }

  # prepare original_data and taxonomy_matrix, if applicable:
  densify_prep <- function(original_data) {
    # save row names for later
    taxon_identity <- rownames(original_data)
    # replace NAs by 0 (no data available) and non-NA entries by 1 (data available)
    full_matrix <- as.matrix(original_data)
    full_matrix[!is.na(full_matrix)] <- 1
    full_matrix[is.na(full_matrix)] <- 0
    # convert dataframe entries to numeric
    full_matrix <- apply(full_matrix, 2, as.numeric)
    # rename row names
    rownames(full_matrix) <- taxon_identity
    return(full_matrix)
  }

  full_matrix <- densify_prep(original_data)

  if(is.data.frame(taxonomy_matrix)){ # if a taxonomy matrix is provided:
    # reduce both the matrix and the taxonomy to the glottocodes present in both files
    full_matrix <- full_matrix[which(rownames(full_matrix)%in%taxonomy_matrix$id),]

    # reorganise taxonomy_matrix: trim it to the taxa in the initial matrix in the corresponding order; reorganise it so that each tip is listed not just in "id" but also in the last listed level
    taxonomy_reorganised <- taxonomy_matrix %>% dplyr::filter(id %in% rownames(full_matrix))
    taxonomy_reorganised <- taxonomy_reorganised[match(rownames(full_matrix),taxonomy_reorganised$id),]
    max_node_nr <- ncol(taxonomy_reorganised)-1
  }

  # prepare iterative coding scheme:

  # prepare documentation file with information on initial (full) matrix
  documentation <- data.frame(iterations = "0 (full)",
                              available_data_points = sum(full_matrix),
                              prop_coded_data = sum(full_matrix)/(nrow(full_matrix)*ncol(full_matrix)),
                              n_tax = nrow(full_matrix),
                              n_var = ncol(full_matrix),
                              worst_tax_abs = paste(names(rowSums(full_matrix))[which(rowSums(full_matrix)==min(rowSums(full_matrix)))],collapse=";"),
                              worst_tax_abs_coding_density = min(rowSums(full_matrix)/ncol(full_matrix)),
                              worst_var_abs = paste(names(colSums(full_matrix))[which(colSums(full_matrix)==min(colSums(full_matrix)))],collapse=";"),
                              worst_var_abs_coding_density = min(colSums(full_matrix)/nrow(full_matrix)),
                              removed_tax = "NA",
                              removed_var = "NA",
                              na_distribution_index = sqrt(var(rowSums(full_matrix)/ncol(full_matrix))+var(colSums(full_matrix)/nrow(full_matrix))))

  if(is.data.frame(taxonomy_matrix)){ # if a taxonomy matrix is provided, score the taxonomic index for the input data frame
    documentation <- cbind(documentation, taxonomic_index = vegan::diversity(table(taxonomy_reorganised$level1)))
  }


  # determine weighted row and column scores (r_weights and c_weights)
  cross_prod = crossprod(full_matrix)
  ev_cross_prod = eigen(cross_prod)$vectors[,1]
  r_weights <- (full_matrix %*% (ev_cross_prod))/sum(ev_cross_prod) # the weights are weighted averages
  c_weights <- (t(r_weights) %*% full_matrix)/sum(r_weights) # the weights are weighted averages

  updated_matrix <- full_matrix
  r_weights_updated <- data.frame(weighted_coding_score = r_weights, abs_prop_non_NA = rowSums(full_matrix)/ncol(full_matrix))
  c_weights_updated <- c_weights

  # create a data frame collecting different kinds of weights per taxon
  weight_collection<-r_weights_updated

  # initialise iterative densifying procedure
  for (iterations in 1:max_steps){

    cat("Iteration number",iterations,"\n")

    # We stopped checking the following, because our code can handle cases with zero
    # and it handles them in a sensible way.

    # if (min(c_weights)==0|min(weight_collection$weighted_coding_score)==0){
    #   warning("The lowest weighted row or column mean is 0: min(r_weights) = ", min(weight_collection$weighted_coding_score), "; min(c_weights) = ", min(c_weights_updated))
    # }
    if (mean_type %in% c("arithmetic","geometric","log_odds") == F){
      stop("Mean type must be arithmetic, geometric or log_odds.")
      break
    }
    # fill up current matrix with taxonomic_weight measures entirely only in the first iteration and if taxonomy is considered for pruning
    if (iterations==1 & taxonomy==T){
      cat("Computing initial taxonomic weights.\n")
      families <- unique(taxonomy_reorganised$level1)
      for (f in 1:length(families)){
        taxa <- dplyr::filter(taxonomy_reorganised, .data$level1%in%families[f])$id
        if(length(taxa)==1){ # if taxon is sole representative of its family, --> weight 1
          weight_collection[taxa,"taxonomic_weight"] <- 1
        }
        else{
          # rather than listing the identity of each node, list how often this node appears within the sample. these counts are needed to derive a taxonomic diversity measure
          # where the taxonomy is exhausted (nodes start appearing as NA) --> record 1 rather than NA (this is to make taxa with different taxonomic depth comparable by the measure we employ)
          node_freq_count <- taxonomy_reorganised %>% dplyr::filter(.data$level1==families[f])
          for (node_level in 2:(max_node_nr+1)){ # transform node columns as described above
            node_freq_count[is.na(node_freq_count[,node_level]),node_level] <- "none"
            count_table <- table(node_freq_count[,node_level])
            node_freq_count[,node_level]<-apply(node_freq_count[,node_level],1,function(node_id){if(node_id!="none"){as.integer(count_table[node_id])}else{1}})
          }

          # sometimes taxa have nodes that can be collapsed, given the taxon sample. where this is the case, collapse the nodes accordingly
          # (e.g. if we have 2 low-level sub-taxa as the sole representatives of a larger clade within a family, they will have several (identical) nodes that are unnecessary and which distort their relative taxonomic position to taxa in another clade of this family)
          for (tax in 1:nrow(node_freq_count)){
            relevant_levels<-(unique(as.numeric(node_freq_count[tax,2:(max_node_nr+1)])))
            new_levels<-(c(relevant_levels,rep(1,(max_node_nr-length(relevant_levels))))) # this collapses any unnecessary nodes
            node_freq_count[tax,2:(max_node_nr+1)]<-as.list(new_levels) # update node_freq_count with new node count
            node_freq_count[tax,"prov_tax_weight"]<-1/(prod(new_levels)^(1/length(new_levels))) # this is 1/(geometric mean of (modified) node level counts)
          }
          weight_collection[node_freq_count$id,"taxonomic_weight"]<-prop.table(node_freq_count$prov_tax_weight)
        }
      }
    }

    cat("Computing final weights for all taxa\n")

    # compute a final weight for each taxon via the odds-mean of the absolute proportion of non-NA, the weighted coding score and the taxonomic weight
    if (mean_type == "arithmetic"){ # arithmetic mean
      weight_collection$mean_score <- apply(weight_collection,1,mean)
    } else if (mean_type == "geometric"){ # geometric mean
      weight_collection$mean_score <- apply(weight_collection,1,function(x) prod(x) ** (1/length(x)))
    } else if (mean_type == "log_odds" & taxonomy == T) { # log odds mean if taxonomy considered
      weight_collection$abs_prop_non_NA <- weight_collection$abs_prop_non_NA*coding_weight_factor # multiply all weights by coding_weight_factor
      weight_collection$weighted_coding_score <- weight_collection$weighted_coding_score*coding_weight_factor # multiply all weights by coding_weight_factor
      weight_collection$taxonomic_weight <- weight_collection$taxonomic_weight*tax_weight_factor # multiply all weights by tax_weight_factor so that no values are equal to 1 (which would make computing the log-odds impossible)
      mn<-apply(weight_collection["weighted_coding_score"],1,qlogis)+apply(weight_collection["abs_prop_non_NA"],1,qlogis)+apply(weight_collection["taxonomic_weight"],1,qlogis)/3
      weight_collection$mean_score<-exp(mn)/(1+exp(mn))
    } else if (mean_type == "log_odds" & taxonomy == F) { # log odds mean if taxonomy not considered
      weight_collection$abs_prop_non_NA <- weight_collection$abs_prop_non_NA*coding_weight_factor # multiply all weights by coding_weight_factor
      weight_collection$weighted_coding_score <- weight_collection$weighted_coding_score*coding_weight_factor # multiply all weights by coding_weight_factor
      mn<-apply(weight_collection["weighted_coding_score"],1,qlogis)+apply(weight_collection["abs_prop_non_NA"],1,qlogis)/2
      weight_collection$mean_score<-exp(mn)/(1+exp(mn))
    }

    # identify the worst taxon, family and variable; there may be ties among the worst taxa/variables, thus randomly sample which to remove at such an iteration
    worsttax <- dplyr::filter(weight_collection,.data$mean_score==min(weight_collection$mean_score))
    worsttax <- worsttax[sample(nrow(worsttax),1),]
    worstvar <- sample(colnames(updated_matrix)[which(c_weights_updated==min(c_weights_updated))],1)

    cat("Identifying taxa and variables with lowest score:", rownames(worsttax), "and", worstvar ,"\n")

    # remove taxon if the worst taxon is currently worse than or equally bad as the worst variable(s)
    if(worsttax$mean_score <= min(c_weights_updated)){
      cat("Remove the taxon with lowest coding score:",rownames(worsttax),"\n")

      # track this taxon will be removed; track no variable will be removed
      removed_taxa <- rownames(worsttax)
      removed_vars <- NA

      updated_matrix <- updated_matrix[-which(rownames(updated_matrix)==removed_taxa),] # update matrix by pruning away worst taxon
      weight_collection <- weight_collection[-which(rownames(weight_collection)==removed_taxa),] # update weight collection by pruning away worst taxon
      if(is.data.frame(taxonomy_matrix)){taxonomy_reorganised <- taxonomy_reorganised[taxonomy_reorganised$id%in%removed_taxa==F,]} # update taxonomy_reorganised by pruning away worst taxon if taxonomy provided

    }

    if ((is.data.frame(updated_matrix))){
      if (nrow(updated_matrix)==0) { cat("Trimming aborted - there are no more taxa left.")
        break }
    }

    # remove variable if the worst variable is currently worse than the worst taxon/taxa
    if(worsttax$mean_score > min(c_weights_updated)){
      cat("Remove the variable with the lowest coding score:", worstvar,"\n")

      # track this variable will be removed; track no taxon will be removed
      removed_taxa <- NA
      removed_vars <- worstvar

      updated_matrix <- updated_matrix[,-(which(colnames(updated_matrix)==worstvar))] # update matrix by pruning away worst variable
      c_weights_updated <- c_weights_updated[-(which(colnames(c_weights_updated)==worstvar))] # update weights by pruning away worst variable
    }

    if (!is.matrix(updated_matrix)){
      cat("Trimming aborted - there are no more variables left.")
      break
    }

    # make sure that all variables are/remain sufficiently variable --> the second-most-frequent variable state must contain at least N taxa (N=variability_threshold)
    # if any variables are not/no longer sufficiently variable, remove them
    cat("Ensuring variable variablity.\n")

    pruned_matrix <- original_data[which(rownames(original_data)%in%rownames(updated_matrix)),which(colnames(original_data)%in%colnames(updated_matrix))]
    nrlevels <- data.frame(variable=colnames(pruned_matrix),
                           number_of_variable_states=apply(pruned_matrix,2,function(x)length(table(as.factor(x)))),
                           count_second_largest_variable_state=apply(pruned_matrix,2,function(x)sort(table(as.factor(x)),decreasing=T)[2]))

    if (sum(nrlevels$count_second_largest_variable_state%in%c(NA,1:variability_threshold))!=0){ # only act if there is a variable that needs removal
      uninformative_variables <- rownames(dplyr::filter(nrlevels,.data$count_second_largest_variable_state%in%c(NA,1:variability_threshold)))
      removed_vars <- c(removed_vars,uninformative_variables)[-which(is.na(c(removed_vars,uninformative_variables)))]
      updated_matrix <- updated_matrix[,-which(colnames(updated_matrix)%in%uninformative_variables)] # update matrix by pruning away uninformative variables
      c_weights_updated <- c_weights_updated[-which(colnames(updated_matrix)%in%uninformative_variables)] # update weights by pruning away uninformative variables
      cat("; remove the following uninformative variables", uninformative_variables,"\n")
    }

    if (!is.matrix(updated_matrix)){
      cat("Trimming aborted - there are no variables left.")
      break
    }

    # if removing variables results in uninformative taxa, these must also be removed
    if(min(rowSums(updated_matrix))==0){
      uninformative_taxa <- names(which(rowSums(updated_matrix)==0))
      removed_taxa <- c(removed_taxa,uninformative_taxa)[-which(is.na(c(removed_taxa,uninformative_taxa)))]
      updated_matrix <- updated_matrix[-which(rownames(updated_matrix)%in%uninformative_taxa),] # update matrix by pruning away worst taxon
      weight_collection <- weight_collection[-which(rownames(weight_collection)%in%uninformative_taxa),] # update weight collection by pruning away worst taxon
      if (is.data.frame(taxonomy_matrix)){taxonomy_reorganised <- taxonomy_reorganised[taxonomy_reorganised$id%in%uninformative_taxa==F,]} # update taxonomy_reorganised by pruning away worst taxon
      cat("Remove the following uninformative taxa:", uninformative_taxa,"\n")
    }

    if (nrow(updated_matrix)==0){
      cat("Trimming aborted - there are no taxa left.")
      break
    }

    # if any taxa were removed, the phylogenetic weights for remaining taxa of their families must be updated!
    if ((length(removed_taxa))!=0 & taxonomy == T){
      cat("There are phylogenetic weights to be updated.\n")
      fams_for_tax_weight_update <- unique(dplyr::filter(taxonomy_matrix,id%in%removed_taxa)$level1)
      for (fam in fams_for_tax_weight_update){ # proceed family-wise
        node_freq_count <- dplyr::filter(taxonomy_reorganised,.data$level1%in%fam) # filter out remaining taxa of this family
        if(nrow(node_freq_count)==0){} else if(nrow(node_freq_count)==1){ # if taxon is sole representative of its family, --> weight 1
          weight_collection[node_freq_count$id,"taxonomic_weight"] <- 1
        } else { # if taxon is one of several taxa from the same family in the current sample, the other taxa from that family need to be updated
          for (node_level in 2:(max_node_nr+1)){ # transform node columns as described above
            node_freq_count[is.na(node_freq_count[,node_level]),node_level] <- "none"
            count_table <- table(node_freq_count[,node_level])
            node_freq_count[,node_level] <- apply(node_freq_count[,node_level],1,function(node_id){if(node_id!="none"){as.integer(count_table[node_id])}else{1}})
          }
          # sometimes taxa have nodes that can be collapsed, given the taxon sample. where this is the case, collapse the nodes accordingly
          # (e.g. if we have 2 low-level dialects as the sole representatives of a clade within a family, they will have several (identical) nodes that are unnecessary and which distort their relative taxonomic position to taxa in another clade of this family)
          for (tax in 1:nrow(node_freq_count)){
            relevant_levels <- (unique(as.numeric(node_freq_count[tax,2:(max_node_nr+1)])))
            new_levels <- (c(relevant_levels,rep(1,(max_node_nr-length(relevant_levels))))) # this collapses any unnecessary nodes
            node_freq_count[tax,2:(max_node_nr+1)] <- as.list(new_levels) # update node_freq_count with new node count
            node_freq_count[tax,"prov_tax_weight"] <- 1/(prod(new_levels)^(1/length(new_levels))) # this is 1/(geometric mean of (modified) node level counts)
          }
          ## update the taxonomic weight in weight_collection (multiplied by tax_weight_factor to ensure log-odds can be computed if applicable)
          weight_collection[which(rownames(weight_collection)%in%node_freq_count$id),"taxonomic_weight"]<-prop.table(node_freq_count$prov_tax_weight)
        }
      }
    }

    # update documentation
    if(is.data.frame(taxonomy_matrix)){
      documentation<-rbind(documentation,
                           c(iterations,
                             sum(updated_matrix),
                             sum(updated_matrix)/(nrow(updated_matrix)*ncol(updated_matrix)),
                             nrow(updated_matrix),
                             ncol(updated_matrix),
                             paste(names(rowSums(updated_matrix))[which(rowSums(updated_matrix)==min(rowSums(updated_matrix)))],collapse=";"),
                             min(rowSums(updated_matrix)/ncol(updated_matrix)),
                             paste(names(colSums(updated_matrix))[which(colSums(updated_matrix)==min(colSums(updated_matrix)))],collapse=";"),
                             min(colSums(updated_matrix)/nrow(updated_matrix)),
                             paste(removed_taxa,collapse=";"),
                             paste(removed_vars,collapse=";"),
                             sqrt(var(rowSums(updated_matrix)/ncol(updated_matrix))+var(colSums(updated_matrix)/nrow(updated_matrix))),
                             vegan::diversity(table(taxonomy_reorganised$level1))))
    } else{
      documentation<-rbind(documentation,
                           c(iterations,
                             sum(updated_matrix),
                             sum(updated_matrix)/(nrow(updated_matrix)*ncol(updated_matrix)),
                             nrow(updated_matrix),
                             ncol(updated_matrix),
                             paste(names(rowSums(updated_matrix))[which(rowSums(updated_matrix)==min(rowSums(updated_matrix)))],collapse=";"),
                             min(rowSums(updated_matrix)/ncol(updated_matrix)),
                             paste(names(colSums(updated_matrix))[which(colSums(updated_matrix)==min(colSums(updated_matrix)))],collapse=";"),
                             min(colSums(updated_matrix)/nrow(updated_matrix)),
                             paste(removed_taxa,collapse=";"),
                             paste(removed_vars,collapse=";"),
                             sqrt(var(rowSums(updated_matrix)/ncol(updated_matrix))+var(colSums(updated_matrix)/nrow(updated_matrix)))))
    }

    cat("Recomputing weighted means.\n\n")

    # update the absolute coding proportion of each taxon
    abs_prop_non_NA<-(rowSums(updated_matrix)/ncol(updated_matrix))

    # update c_weights and r_weights
    cross_prod = crossprod(updated_matrix)
    ev_cross_prod = eigen(cross_prod)$vectors[,1]
    r_weights <- (updated_matrix %*% (ev_cross_prod))/sum(ev_cross_prod) # the weights are weighted averages
    r_weights_updated <- data.frame(weighted_coding_score=r_weights,abs_prop_non_NA=abs_prop_non_NA)
    c_weights_updated <- (t(r_weights_updated$weighted_coding_score) %*% updated_matrix)/sum(r_weights_updated$weighted_coding_score) # the weights are weighted averages

    # update the taxa's weight_collection to contain the correct prop_non_NA and weighted_coding_score
    weight_collection$abs_prop_non_NA <- abs_prop_non_NA
    weight_collection$weighted_coding_score <- signif(r_weights_updated$weighted_coding_score,digits=7)
  }
  return(documentation)
}
