######################################################
# F2
#
#' Iterative matrix densification according to specified criteria.
#'
#' The output of this densification is a log-file, which specifies details about the matrix after
#' each iteration.
#'
#' @param original_data A data frame with the glottocodes as row names and variable names as column names.
#'   Any question marks, empty entries, or "NA"s must be represented as NAs.
#'
#' @param max_steps An optional integer specifying the maximum number of iterations attempted during densification.
#'   Default is 1.
#'
#' @param mean_type A character string specifying the type of mean to be used for calculating the final weights.
#'   Possible values are "arithmetic", "geometric", or "log_odds". Default is "log_odds".
#'
#' @param taxonomy A logical indicating whether taxonomic diversity should be factored in for language pruning/retention.
#'   If `TRUE`, taxonomic diversity is considered; if `FALSE`, it is not. Default is `FALSE`.
#'
#' @param original_register A data frame containing the register retrieved from the `glottocode_taxonomy` function.
#'   This parameter must be specified if `taxonomy = TRUE`.
#'
#' @param tax_weight_factor A numeric value between 0 and 1 that determines the relative weight given to taxonomy in language pruning.
#'   This parameter must be specified if `taxonomy = TRUE` and `mean_type = "log_odds"`.
#'
#' @param coding_weight_factor A numeric value between 0 and 1 that determines the relative weight given to coding quality
#'   (absolute coding score and weighted coding score) in language pruning. This parameter must be specified if `taxonomy = TRUE`
#'   and `mean_type = "log_odds"`.
#'
#' @return A data frame with details about the matrix after each iteration of densification.
#'
#' @examples
#' # Assuming `original_data` and `original_register` are appropriate data frames
#' densify_steps(original_data, max_steps = 3, mean_type = "log_odds", taxonomy = TRUE, original_register, tax_weight_factor = 0.99, coding_weight_factor = 0.99)
#'
#' @import vegan
#' @export
######################################################

# original_data: must be a data frame with the glottocodes as row names and variable names as column names. any question marks, empty entries, "NA"s must be NAs
# max_steps: can be specified at will, defines how many iterations are attempted
# mean_type: "arithmetic", "geometric" or "log_odds"
# taxonomy: T or F; if T, taxonomic diversity is factored in for language pruning/retention, if F it is not
# original_register: register retrieved from F1; must be specified if taxonomy = T
# tax_weight_factor: must be between (not including) 0 and 1 and determines the relative weight given to taxonomy in language pruning; must be specified if taxonomy = T and if mean_type = log_odds
# coding_weight_factor: must be between (not including) 0 and 1 and determines the relative weight given to coding quality (absolute coding score and weighted coding score) in language pruning; must be specified if taxonomy = T and if mean_type = log_odds

densify_steps <- function(original_data, max_steps = 1, mean_type = "log_odds", taxonomy = F, original_register, tax_weight_factor = 0.99, coding_weight_factor = 0.99){
  library(vegan)

  # prepare original_data and original_register, if applicable:

  densify_prep <- function(original_data) {
    # save row names for later
    languages <- rownames(original_data)
    # replace NAs by 0 (no data available) and non-NA entries by 1 (data available)
    full_matrix <- as.matrix(original_data)
    full_matrix[!is.na(full_matrix)] <- 1
    full_matrix[is.na(full_matrix)] <- 0
    # convert dataframe entries to numeric
    full_matrix <- apply(full_matrix, 2, as.numeric)
    # rename row names
    rownames(full_matrix) <- languages
    return(full_matrix)
  }

  full_matrix <- densify_prep(original_data)

  if(taxonomy == T){
    # reduce both the matrix and the register to the glottocodes present in both files
    full_matrix <- full_matrix[which(rownames(full_matrix)%in%original_register$glottocode),]

    # reorganise register: trim it to the languages in the initial matrix in the corresponding order; reorganise it so that each tip is not listed in the last column but at whatever node appropriate for that tip
    register <- original_register %>% filter(glottocode %in% rownames(full_matrix)) %>% select(-c(glottolog.name,longitude,latitude))
    register <- register[match(rownames(full_matrix),register$glottocode),]
    for(i in 1:nrow(register)){
      register[i,sum(!is.na(register[i,2:ncol(register)]))+1]<-register$glottocode[i] # the 2nd column lists the first node (hence 2; 1)
    }
    register <- select(register,-ncol(register)) # delete the last column of register
    max_node_nr <- ncol(register)-1
  }

  # prepare iterative coding scheme:

  # prepare documentation file with information on initial (full) matrix
  documentation <- data.frame(iterations = "0 (full)",
                              available_data_points = sum(full_matrix),
                              prop_coded_data = sum(full_matrix)/(nrow(full_matrix)*ncol(full_matrix)),
                              n_lg = nrow(full_matrix),
                              n_var = ncol(full_matrix),
                              worst_lg_abs = paste(names(rowSums(full_matrix))[which(rowSums(full_matrix)==min(rowSums(full_matrix)))],collapse=";"),
                              worst_lg_abs_coding_density = min(rowSums(full_matrix)/ncol(full_matrix)),
                              worst_var_abs = paste(names(colSums(full_matrix))[which(colSums(full_matrix)==min(colSums(full_matrix)))],collapse=";"),
                              worst_var_abs_coding_density = min(colSums(full_matrix)/nrow(full_matrix)),
                              removed_lg = "NA",
                              removed_var = "NA",
                              taxonomic_index = vegan::diversity(table(register$glottolog.node1)),
                              na_distribution_index = sqrt(var(rowSums(full_matrix)/ncol(full_matrix))+var(colSums(full_matrix)/nrow(full_matrix))))

  # determine weighted row and column scores (r_weights and c_weights)
  cross_prod = crossprod(full_matrix)
  ev_cross_prod = eigen(cross_prod)$vectors[,1]
  r_weights <- (full_matrix %*% (ev_cross_prod))/sum(ev_cross_prod) # the weights are weighted averages
  c_weights <- (t(r_weights) %*% full_matrix)/sum(r_weights) # the weights are weighted averages

  updated_matrix <- full_matrix
  r_weights_updated <- data.frame(weighted_coding_score = r_weights, abs_prop_non_NA = rowSums(full_matrix)/ncol(full_matrix))
  c_weights_updated <- c_weights

  # create a data frame collecting different kinds of weights per language
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
    if (taxonomy == F){
      warning("Attention, taxonomy disregarded.\n")
    }

    # fill up current matrix with taxonomic_weight measures entirely only in the first iteration!
    if (iterations==1 & taxonomy==T){
      cat("Computing initial taxonomic weights.\n")
      families <- unique(register$glottolog.node1)
      for (f in 1:length(families)){
        lgs <- filter(register, glottolog.node1%in%families[f])$glottocode
        if(length(lgs)==1){ # if lg is sole representative of its family, --> weight 1
          weight_collection[lgs,"taxonomic_weight"] <- 1
        }
        else{
          # rather than listing the identity of each node, list how often this node appears within the sample. these counts are needed to derive a taxonomic diversity measure
          # where the taxonomy is exhausted (nodes start appearing as NA) --> record 1 rather than NA (this is to make languages with different taxonomic depth comparable by the measure we employ)
          node_freq_count <- register %>% filter(glottolog.node1==families[f])
          for (node_level in 2:(max_node_nr+1)){ # transform node columns as described above
            node_freq_count[is.na(node_freq_count[,node_level]),node_level] <- "none"
            count_table <- table(node_freq_count[,node_level])
            node_freq_count[,node_level]<-(unlist(lapply(node_freq_count[,node_level],function(node_id){if(node_id!="none"){as.integer(count_table[node_id])}else{1}})))
          }

          # sometimes languages have nodes that can be collapsed, given the language sample. where this is the case, collapse the nodes accordingly
          # (e.g. if we have 2 low-level dialects as the sole representatives of a clade within a family, they will have several (identical) nodes that are unnecessary and which distort their relative taxonomic position to languages in another clade of this family)
          for (lg in 1:nrow(node_freq_count)){
            relevant_levels<-(unique(as.numeric(node_freq_count[lg,2:(max_node_nr+1)])))
            new_levels<-(c(relevant_levels,rep(1,(max_node_nr-length(relevant_levels))))) # this collapses any unnecessary nodes
            node_freq_count[lg,2:(max_node_nr+1)]<-new_levels # update node_freq_count with new node count
            node_freq_count[lg,"prov_tax_weight"]<-1/(prod(new_levels)^(1/length(new_levels))) # this is 1/(geometric mean of (modified) node level counts)
          }
          weight_collection[node_freq_count$glottocode,"taxonomic_weight"]<-prop.table(node_freq_count$prov_tax_weight)
        }
      }
    }

    cat("Computing final weights for all languages\n")

    # compute a final weight for each language via the odds-mean of the absolute proportion of non-NA, the weighted coding score and the taxonomic weight
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

    # identify the worst language, family and variable; there may be ties among the worst languages/variables, thus randomly sample which to remove at such an iteration
    worstlg <- filter(weight_collection,mean_score==min(weight_collection$mean_score))
    worstlg <- worstlg[sample(nrow(worstlg),1),]
    worstvar <- sample(colnames(updated_matrix)[which(c_weights_updated==min(c_weights_updated))],1)

    cat("Identifying languages and variables with lowest score:", rownames(worstlg), "and", worstvar ,"\n")

    # remove language if the worst language is currently worse than or equally bad as the worst variable(s)
    if(worstlg$mean_score <= min(c_weights_updated)){
      cat("Remove the language with lowest coding score:",rownames(worstlg),"\n")

      # track this language will be removed; track no variable will be removed
      removed_lgs <- rownames(worstlg)
      removed_vars <- NA

      updated_matrix <- updated_matrix[-which(rownames(updated_matrix)==removed_lgs),] # update matrix by pruning away worst lg
      weight_collection <- weight_collection[-which(rownames(weight_collection)==removed_lgs),] # update weight collection by pruning away worst lg
      if(taxonomy==T){register <- register[register$glottocode%in%removed_lgs==F,]} # update register by pruning away worst lg

    }

    if ((is.data.frame(updated_matrix))){
      if (nrow(updated_matrix)==0) { cat("Trimming aborted - there are no more languages left.")
        break }
    }

    # remove variable if the worst variable is currently worse than the worst language(s)
    if(worstlg$mean_score > min(c_weights_updated)){
      cat("Remove the variable with the lowest coding score:", worstvar,"\n")

      # track this variable will be removed; track no language will be removed
      removed_lgs <- NA
      removed_vars <- worstvar

      updated_matrix <- updated_matrix[,-(which(colnames(updated_matrix)==worstvar))] # update matrix by pruning away worst variable
      c_weights_updated <- c_weights_updated[-(which(colnames(c_weights_updated)==worstvar))] # update weights by pruning away worst variable
    }

    if (!is.matrix(updated_matrix)){
      cat("Trimming aborted - there are no more variables left.")
      break
    }

    # make sure that all variables are/remain sufficiently variable --> the second-most-frequent variable state must contain at least 3 languages
    # if any variables are not/no longer sufficiently variable for our areal analysis, remove them
    cat("Ensuring variable variablity.\n")

    pruned_matrix <- original_data[which(rownames(original_data)%in%rownames(updated_matrix)),which(colnames(original_data)%in%colnames(updated_matrix))]
    nrlevels <- data.frame(variable=colnames(pruned_matrix),
                           number_of_variable_states=apply(pruned_matrix,2,function(x)length(table(as.factor(x)))),
                           count_second_largest_variable_state=apply(pruned_matrix,2,function(x)sort(table(as.factor(x)),decreasing=T)[2]))

    if (sum(nrlevels$count_second_largest_variable_state%in%c(NA,1,2))!=0){ # only act if there is a variable that needs removal
      uninformative_variables <- rownames(filter(nrlevels,count_second_largest_variable_state%in%c(NA,1,2)))
      removed_vars <- c(removed_vars,uninformative_variables)[-which(is.na(c(removed_vars,uninformative_variables)))]
      updated_matrix <- updated_matrix[,-which(colnames(updated_matrix)%in%uninformative_variables)] # update matrix by pruning away uninformative variables
      c_weights_updated <- c_weights_updated[-which(colnames(updated_matrix)%in%uninformative_variables)] # update weights by pruning away uninformative variables
      cat("; remove the following uninformative variables", uninformative_variables,"\n")
    }

    if (!is.matrix(updated_matrix)){
      cat("Trimming aborted - there are no variables left.")
      break
    }

    # if removing variables results in uninformative languages, these must also be removed
    if(min(rowSums(updated_matrix))==0){
      uninformative_languages <- names(which(rowSums(updated_matrix)==0))
      removed_lgs <- c(removed_lgs,uninformative_languages)[-which(is.na(c(removed_lgs,uninformative_languages)))]
      updated_matrix <- updated_matrix[-which(rownames(updated_matrix)%in%uninformative_languages),] # update matrix by pruning away worst lg
      weight_collection <- weight_collection[-which(rownames(weight_collection)%in%uninformative_languages),] # update weight collection by pruning away worst lg
      if (taxonomy==T){register <- register[register$glottocode%in%uninformative_languages==F,]} # update register by pruning away worst lg
      cat("Remove the following uninformative languages:", uninformative_languages,"\n")
    }

    if (nrow(updated_matrix)==0){
      cat("Trimming aborted - there are no languages left.")
      break
    }

    # if any languages were removed, the phylogenetic weights for remaining languages of their families must be updated!
    if ((length(removed_lgs))!=0 & taxonomy==T){
      cat("There are phylogenetic weights to be updated.\n")
      fams_for_tax_weight_update <- unique(filter(original_register,glottocode%in%removed_lgs)$glottolog.node1)
      for (fam in fams_for_tax_weight_update){ # proceed family-wise
        node_freq_count <- filter(register,glottolog.node1%in%fam) # filter out remaining languages of this family
        if(nrow(node_freq_count)==0){} else if(nrow(node_freq_count)==1){ # if lg is sole representative of its family, --> weight 1
          weight_collection[node_freq_count$glottocode,"taxonomic_weight"] <- 1
        } else { # if lg is one of several lgs from the same family in the current sample, the other langauages from that family need to be updated
          for (node_level in 2:(max_node_nr+1)){ # transform node columns as described above
            node_freq_count[is.na(node_freq_count[,node_level]),node_level] <- "none"
            count_table <- table(node_freq_count[,node_level])
            node_freq_count[,node_level] <- (unlist(lapply(node_freq_count[,node_level], function(node_id){if(node_id!="none"){as.integer(count_table[node_id])}else{1}})))
          }
          # sometimes languages have nodes that can be collapsed, given the language sample. where this is the case, collapse the nodes accordingly
          # (e.g. if we have 2 low-level dialects as the sole representatives of a clade within a family, they will have several (identical) nodes that are unnecessary and which distort their relative taxonomic position to languages in another clade of this family)
          for (lg in 1:nrow(node_freq_count)){
            relevant_levels <- (unique(as.numeric(node_freq_count[lg,2:(max_node_nr+1)])))
            new_levels <- (c(relevant_levels,rep(1,(max_node_nr-length(relevant_levels))))) # this collapses any unnecessary nodes
            node_freq_count[lg,2:(max_node_nr+1)] <- new_levels # update node_freq_count with new node count
            node_freq_count[lg,"prov_tax_weight"] <- 1/(prod(new_levels)^(1/length(new_levels))) # this is 1/(geometric mean of (modified) node level counts)
          }
          ## update the taxonomic weight in weight_collection (multiplied by tax_weight_factor to ensure log-odds can be computed if applicable)
          weight_collection[which(rownames(weight_collection)%in%node_freq_count$glottocode),"taxonomic_weight"]<-prop.table(node_freq_count$prov_tax_weight)
        }
      }
    }

    # update documentation
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
                           paste(removed_lgs,collapse=";"),
                           paste(removed_vars,collapse=";"),
                           vegan::diversity(table(register$glottolog.node1)),
                           sqrt(var(rowSums(updated_matrix)/ncol(updated_matrix))+var(colSums(updated_matrix)/nrow(updated_matrix)))))

    cat("Recomupting weighted means.\n\n")

    # update the absolute coding proportion of each langauge
    abs_prop_non_NA<-(rowSums(updated_matrix)/ncol(updated_matrix))

    # update c_weights and r_weights
    cross_prod = crossprod(updated_matrix)
    ev_cross_prod = eigen(cross_prod)$vectors[,1]
    r_weights <- (updated_matrix %*% (ev_cross_prod))/sum(ev_cross_prod) # the weights are weighted averages
    r_weights_updated <- data.frame(weighted_coding_score=r_weights,abs_prop_non_NA=abs_prop_non_NA)
    c_weights_updated <- (t(r_weights_updated$weighted_coding_score) %*% updated_matrix)/sum(r_weights_updated$weighted_coding_score) # the weights are weighted averages

    # update the languages' weight_collection to contain the correct prop_non_NA and weighted_coding_score
    weight_collection$abs_prop_non_NA <- abs_prop_non_NA
    weight_collection$weighted_coding_score <- signif(r_weights_updated$weighted_coding_score,digits=7)
  }
  return(documentation)
}
