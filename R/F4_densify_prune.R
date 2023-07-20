######################################################
# F4
#
# Function to actually densify the original input matrix, given the optimum quality score.
#
######################################################

densify_prune <- function(original_data, documentation, optimum){
  documentation <- slice(documentation, 1:optimum)
  prune_lgs <- unique(unlist(strsplit(documentation$removed_lg,";")))[unique(unlist(strsplit(documentation$removed_lg,";")))!="NA"]
  prune_vars <- unique(unlist(strsplit(documentation$removed_var,";")))[unique(unlist(strsplit(documentation$removed_var,";")))!="NA"]
  pruned_matrix <- original_data[which(rownames(original_data)%in%prune_lgs==F),which(colnames(original_data)%in%prune_vars==F)]
  return(pruned_matrix)
}
