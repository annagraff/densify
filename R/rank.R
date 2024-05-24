#' Densify results ranking based on a quality scoring function
#'
#' Ranks densify results by a quality scoring function from lowest to highest. Ties receive the same rank. The user can provide their own
#' quality scoring formula in `scoring_function`, which can use any statistics provided in a [densify_result] tibble.
#'
#' @param x a [densify_result] object
#' @param ... other arguments passed to methods
#' @param scoring_function an expression that computes a densification quality score. The default maximizes the product of
#'   the number of data points and the coding density. See [densify_result] for available statistics that can
#'   be used to compute a suitable scoring_function
#'
#' @return an integer vector of ranks corresponding to the densify iterations
#'
#' @importFrom generics rank_results
#' @export
rank_results.densify_result <- function(x, ..., scoring_function = n_data_points * coding_density) {
  check_dots_unnamed()

  n_data_points <- coding_density <- NULL

  score_quo <- enquo(scoring_function)

  # turn `list(...)` into list(`..1`, `.2`, ...)
  score_quo <- if (quo_is_call(score_quo, "list")) {
    lapply(as.list(quo_get_expr(score_quo))[-1], new_quosure,  quo_get_env(score_quo))
  } else {
    list(score_quo)
  }

  scores <- numeric()

  # evaluate the ranking scores
  for (quo in score_quo) {
    scoring_function <- eval_densify_results_quality_score(x, quo, .caller = "rank_results")
    scores <- cbind(scores, scoring_function)
  }

  vec_rank(scores, direction = "desc", incomplete = "na", ties = "min")
}

#' @importFrom generics rank_results
#' @export
generics::rank_results

#' Obtain densified data frame that maximizes a quality scoring function
#'
#' Returns the densified data frame from the densify iteration step that maximizes a quality scoring function. In case of ties,
#' first result is returned. The user can provide their own quality score formula in `scoring_function`, which can use any statistics
#' provided in a [densify_result] tibble.
#'
#' @param tree a [densify_result] object
#' @param ... other arguments passed to methods
#' @param scoring_function an expression that computes a densification quality score. The default maximizes the product of
#'   the number of data points and the coding density. See [densify_result] for available statistics that can
#'   be used to compute a suitable scoring_function
#'
#' @return the densified data frame
#'
#' @importFrom generics prune
#' @export
prune.densify_result <- function(tree, ..., scoring_function = n_data_points * coding_density) {
  check_dots_unnamed()

  n_data_points <- coding_density <- NULL

  ranks <- rank_results(tree, scoring_function = {{scoring_function}})

  # chose the highest ranked result
  best <- which(ranks == 1L)
  length(best) > 0L || cli::cli_abort("no ranking can be established, aborting")

  if (length(best) > 0) {
    best <- best[[1L]]
    cli::cli_warn("in {.fn prune}: multiple best matches, returning the first one at index {.field {best}}")
  }

  # return the pruned data set
  data <- as.data.frame(tree$data[[best]])
  if (nrow(data) == 0L) {
    cli::cli_warn("in {.fn prune}: returning empty data frame")
  }
  data
}

#' @importFrom generics prune
#' @export
generics::prune




eval_densify_results_quality_score <- function(densify_results, quo, ..., .caller = call_name(caller_call())) {
  local_error_call(caller_env())

  stats <- densify_results[names(densify_results) %in% pruning_state_stat_vars]
  expr <- quo_get_expr(quo)

  # evaluate the measure
  scoring_function <- eval_tidy(quo, stats)

  is_bare_numeric(scoring_function, nrow(stats)) || cli::cli_abort(c(
    "invalid expression {.code {as_label(score_expr)}}",
    i = "expression result must be a numeric vector of length {nrow(stats)}, got {pillar::obj_sum(scoring_function)}"
  ))

  # check that it uses variables correctly
  if(!any(extract_symbols(expr) %in% names(stats))) {
    known_stats <- sprintf("{.field %s}", names(stats))
    names(known_stats) <- rep_along(known_stats, "*")

    cli::cli_warn(c(
      "in {.fn {(.caller)}}: expression {.code {as_label(score_expr)}} does not appear to use any known densify stats",
      i = "available stats are:",
      known_stats
    ))
  }

  scoring_function
}


extract_symbols <- function(expr) {
  if (is_symbol(expr)) {
    as.character(expr)
  }
  else if (is_call(expr)) {
    unlist(lapply(expr[-1L], extract_symbols))
  } else {
    character()
  }

}
