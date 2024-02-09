#' Iterative matrix densification
#'
#' Iteratively densifies an input data frame, returning a specially formatted tible ([densify_result]) describing the result of each densification step.
#' Densified data can be retrieved using function [prune][prune.densify_result], or manually by invoking [as_tibble] on the values in the column `data`. The densification results
#' can also be inspected using [visualize][visualize.densify_result]
#'
#' @param data A data frame with observations/taxa (such as languages) in rows and variables in columns.
#'
#' @param cols  <[`tidy-select`][tidyselect::language]> specification of variable columns to densify (default is all columns are treated as variable
#'   columns). Columns not specified here will be preserved during densification.
#'
#' @param taxonomy A taxonomy tree, which can be a `phylo` object (e.g. result of [ape::read.tree] or a data frame with columns `id` and `parent_id`
#'   (such as [glottolog_languoids]). A taxonomy must be provided as a named argument if you want to consider the taxonomic diversity for densification.
#'
#' @param taxon_id  The name of the column with taxa identifiers. The `data` must contain such a column if a taxonomy is supplied. If not supplied,
#'   the function will try to make an educated guess based on the column contents.
#'
#' @param scoring A character string specifying the type of scoring used for calculating the importance weights. Possible values are are "arithmetic",
#'   "geometric", or "log_odds". Default is "log_odds".
#'
#' @param min_variability An integer specifying the minimal threshold of the second-most-frequent state of any variable. Variables below this threshold
#'   will be discarded. Supply `NA` to disable variability pruning. Default is `1`.
#'
#' @param limits A named list specifying the limit conditions for the densification process termination. Available limit conditions are
#'   `min_coding_density` (a number between 0 and 1 specifying the target coding density), `min_prop_rows` (a number between 0 and 1 specifying the
#'   minimal proportion of rows that have to be retained in the data), and `min_prop_cols` (a number between 0 and 1 specifying the
#'   minimal proportion of variable columns that have to be retained in the data). More then one condition can be specified simultaneously, `densify()`
#'   will stop if any condition is reached.
#'
#' @param scoring_weights A named list specifying additional weighting factors to be applied during importance score calculation. Supported parameters
#'   are `coding` and `taxonomy`, which apply to coding and taxonomic importance weights, respectively. The values for these parameters must be a
#'   number between 0 and 1. Specify 0 (or `NULL`, or `NA`) to disable the respective weight calculation. Use these parameters to tweak the relative
#'   importance between the available coding and taxonomic diversity. Note: the weighting factor is applied before the scoring function, experimentation
#'   will be required to find a good balance of parameters for `scoring = "log_odds"`. See the vignette for more details.
#'
#' @return A specially formatted tibble of iterative densification results, [densify_result]
#'
#' @examples
# # densify first 100 taxa and 100 variables from WALS
#' densify(WALS[1:100, 1:100], taxonomy = glottolog_languoids)
#'
#' @export
densify <- function(
  data,
  cols,
  ...,
  taxonomy,
  taxon_id,
  scoring = c("log_odds", "arithmetic", "geometric"),
  min_variability = 1,
  limits = list(min_coding_density = 1, min_prop_rows = NA, min_prop_cols = NA),
  scoring_weights = list(coding = 1, taxonomy = 1)
) {
  # argument validation and processing
  check_dots_empty()

  data <- check_data(data)
  taxonomy <- check_taxonomy(taxonomy)
  taxon_id <- check_taxon_id(taxon_id, data, taxonomy)
  vars <- check_variables(data, enquo(cols), taxon_id, .arg = "cols")
  scoring_fn <- get_scoring_function(arg_match(scoring))
  limits <- check_limits(limits, nrow(data), length(vars))
  weights <- check_scoring_weights(scoring_weights, !is.null(taxonomy))

  # match taxonomy and data by taxon
  if (!is_null(taxonomy)) taxonomy <- match_taxa(data, taxonomy, taxon_id)

  # init the pruning state
  state <- init_pruning_state(
    data,
    vars,
    ids = if(!is.null(taxon_id)) data[[taxon_id]],
    taxonomy = taxonomy,
    scoring_fn = scoring_fn,
    scoring_weights = weights
  )

  # factory for pruned data promise creation
  pruned_data_factory <- make_pruned_data_factory(data)

  # pruning steps will be recorded here
  densify_log <- make_log_entry_for_pruning_state(state, pruned_data_factory)

  # utility that displays current coding density of the pruning process
  current_coding_density <- function() {
    if ((n <- nrow(densify_log)) == 0 || !is.finite(densify_log$coding_density[[n]])) return("")

    sprintf("| current coding density %s%%", round(densify_log$coding_density[[n]], 2)*100)
  }

  # initial pruning of uninformative taxa
  state <- prune_non_informative_data(state, min_variability = min_variability, .changes = function(changes) {
    densify_log <<- vec_rbind(densify_log, make_log_entry_for_pruning_state(state, pruned_data_factory, changes))
  })

  # start the progress bar
  cli::cli_progress_bar(
    type = "custom",
    format = "{cli::pb_spin} Pruned {cli::pb_current} dimensions ({cli::pb_rate}) {current_coding_density()} | {cli::pb_elapsed}"
  )

  # prune data until limit conditions are satisfied
  while (check_densify_limit_condition(densify_log, limits)) {
    # select the datum dimension with the lowest score
    lowest_score <- identify_lowest_scores(state$weights)

    # resolve at random if we have multiple candidates
    if (nrow(lowest_score)) lowest_score <- lowest_score[sample.int(nrow(lowest_score), 1L), ]

    # record the changes
    changes <- data_frame(
      type = c("taxon", "variable")[[lowest_score$axis]],
      id = dimnames(state$matrix)[[lowest_score$axis]][[lowest_score$index]],
      score = lowest_score$score,
      reason = "low score"
    )

    # prune the datum
    state <- prune_indices(state, lowest_score)

    # prune uninformative taxa that might have resulted from removal
    state <- prune_non_informative_data(state,
      check_taxa = lowest_score$axis == 2L,
      check_vars = lowest_score$axis == 1L,
      min_variability = min_variability,
      .changes = setter(changes, .compose = vec_rbind)
    )

    # save the state
    densify_log <- vec_rbind(densify_log, make_log_entry_for_pruning_state(state, pruned_data_factory, changes))
    cli::cli_progress_update()
  }

  densify_result(densify_log)
}


#' Result of densify function
#'
#' A tibble that contains the information about each pruning step performed by [densify] along with some useful
#' statistics. Use [rank_results][rank_results.densify_result] to rank the densify results according to a quality score, [prune][prune.densify_result] to retrieve a
#' subjectively best result, and [visualize][visualize.densify_result] to visually compare the quality scores between different pruning steps.
#' You can also do custom evaluation on the results tibble, e.g. by using `dplyr` and retrieve the pruned datasets from the
#' column `data`.
#'
#' @format A tibble with the following columns:
#' \describe{
#'   \item{.step}{densify process iteration step}
#'   \item{n_data_points}{number of total non-missing data points in the pruned data}
#'   \item{n_rows}{number of rows (taxa) in the pruned data}
#'   \item{n_cols}{numbr of variables in the pruned data (note: number refers specifically to variables, extra columns are always preserved)}
#'   \item{coding_density}{proportion of non-missing data points relative to the data size}
#'   \item{row_coding_density_min}{minimum proportion of non-missing data points, considering each row}
#'   \item{row_coding_density_median}{median proportion of non-missing data points, considering each row}
#'   \item{row_coding_density_max}{maximum proportion of non-missing data points, considering each row}
#'   \item{col_coding_density_min}{minimum proportion of non-missing data points, considering each column}
#'   \item{col_coding_density_median}{median proportion of non-missing data points, considering each column}
#'   \item{col_coding_density_max}{maximum proportion of non-missing data points, considering each column}
#'   \item{taxonomic_index}{Shannon diversity index for top-level taxa (if taxonomy has been supplied)}
#'   \item{data}{a list containing promise objects that evaluate to pruned datasets (use [as.data.frame] to evaluate a promise)}
#'   \item{changes}{a tibble containing information about the pruned data dimensions}
#' }
#'
#'
#' @name densify_result
#' @rdname densify_result
NULL

densify_result <- function(densify_log) {
  # densify_log$data <- make_pruned_data_frame_promise_list(densify_log$data)

  densify_log <- tibble::new_tibble(df_list(.step = seq_len(nrow(densify_log)), densify_log))

  class(densify_log) <- c("densify_result", class(densify_log))
  densify_log
}


select.densify_result <- function(.data, ...) {
  cli::cli_warn(c(
    "executing {.fn select} on densify results drops the {.attr densify_result} class",
    i = "manually convert to tible with {.fn as_tibble} to silence this warning"
  ), call = caller_call())

  .data <- tibble::new_tibble(unclass(.data))
  NextMethod()
}

rename.densify_result <- function(.data, ...) {
  cli::cli_warn(c(
    "executing {.fn rename} on densify results drops the {.attr densify_result} class",
    i = "manually convert to tible with {.fn as_tibble} to silence this warning"
  ), call = caller_call())

  .data <- tibble::new_tibble(unclass(.data))
  NextMethod()
}

#' @importFrom pillar tbl_format_header
#' @export
tbl_format_header.densify_result <- function(x, setup, ...) {
  cli::cli_fmt({


    cli::cli_text("# Densify result with {nrow(x)} pruning steps")

    best_idx <- if (has_name(x, "coding_density")) which.max(x$coding_density)
    if (length(best_idx) == 1L) {
      best <- x[best_idx, ]

      cli::cli_text("# Highest achieved coding density is {round(best$coding_density, 2)*100}%")
    }
  })
}

#' @importFrom pillar tbl_format_footer
#' @export
tbl_format_footer.densify_result <- function(x, setup, ...) {
  arg <- caller_arg(x)

  out <- cli::cli_fmt({
    cli::cli_text("# ")
    cli::cli_text("# Use {.fn plot} or {.fn visualize} to visualize density statistics and pick the best result")
    cli::cli_text("# Use {.fn rank_results} or {.fn prune} to pick the best densify result based on your criteria")
    cli::cli_text("# Use {.fn eval} to manually retrieve any data frame from column {.code data}")
  })

  c(NextMethod(), out)
}


# Build a new entry in the densify log
#
# - changes is a data frame describing the changes in the state
# - data is the original data frame
make_log_entry_for_pruning_state <- function(state, pruned_data_factory, changes = data_frame()) {
  new_data_frame(list2(
    !!!get_pruning_state_stats(state),
    data = list(pruned_data_factory(state$indices$rows, state$indices$cols)),
    changes = list(tibble::new_tibble(changes))
  ), class = c("tbl_df", "tbl"))
}



#  █████  ██████   ██████       ██████ ██   ██ ███████  ██████ ██   ██ ██ ███    ██  ██████
# ██   ██ ██   ██ ██           ██      ██   ██ ██      ██      ██  ██  ██ ████   ██ ██
# ███████ ██████  ██   ███     ██      ███████ █████   ██      █████   ██ ██ ██  ██ ██   ███
# ██   ██ ██   ██ ██    ██     ██      ██   ██ ██      ██      ██  ██  ██ ██  ██ ██ ██    ██
# ██   ██ ██   ██  ██████       ██████ ██   ██ ███████  ██████ ██   ██ ██ ██   ████  ██████
#
#
# Functions for validating and processing densify() arguments

check_data <- function(data, ..., .arg = caller_arg(data)) {
  local_error_call(caller_env())
  .arg <- .arg

  is.data.frame(data) || cli::cli_abort(c(
    "{.arg {(.arg)}} must be a data frame",
    i = "got {.code {as_label(data)}}")
  )

  if (nrow(data) == 0L || ncol(data) == 0L) {
     cli::cli_warn(c(
      "!" = "in {.fun densify}: empty dataset supplied"
    ), call = caller_call())
  }

  data
}


check_taxonomy <- function(taxonomy, ..., .arg = caller_arg(taxonomy)) {
  local_error_call(caller_env())
  .arg <- .arg

  # no taxonomy specified, which is ok
  if(missing(taxonomy)) return(NULL)

  # flatten the taxonomy
  taxonomy <- as_flat_taxonomy_matrix(taxonomy, .x = .arg)

  # convert the taxonomy to matrix with id column removed
  identical(names(taxonomy)[1L], "id") && is.character(taxonomy$id) || abort("invalid taxonomy structure", .internal = TRUE)

  `rownames<-`(as.matrix(taxonomy[-1L]), taxonomy$id)
}


check_variables <- function(data, quoted_cols, taxon_id, ..., .arg) {
  local_error_call(caller_env())
  .arg <- .arg

  vars <- if(!quo_is_missing(quoted_cols)) {
    tidyselect::eval_select(quoted_cols, data, allow_rename = FALSE, error_call = caller_env())
  } else {
    cli::cli_warn(c(
      "!" = "in {.fun densify}: no {.arg {(.arg)}} argument specified, using all columns as variables",
      i = "use {.code {(.arg)} = <tidy column spec>} to silence this warning"
    ), call = caller_call())
    set_names(seq_len(ncol(data)), names(data))
  }

  if (any(taxon_id %in% names(vars))) {
    if (!quo_is_missing(quoted_cols)) {
      cli::cli_warn(c(
        "!" = "in {.fun densify}: removing taxon id {.var {taxon_id}} from the variable selection",
        i = "please adjust {.code {(.arg)} = <tidy column spec>} to silence this warning"
      ), call = caller_call())
    }

    vars <- vars[!names(vars) %in% taxon_id]
  }

  if (length(vars) == 0L) {
     cli::cli_warn(c(
      "!" = "in {.fun densify}: no variables selected, returning empty data frame",
      i = "please adjust {.code {(.arg)} = <tidy column spec>}"
    ), call = caller_call())
  }


  # build the densify object


  vars
}



check_taxon_id <- function(taxon_id, data, taxonomy, ..., .arg = caller_arg(taxon_id)) {
  local_error_call(caller_env())
  .arg <- .arg

  # taxon id is provided
  if (!missing(taxon_id)) {
    is_string(taxon_id) || cli::cli_abort(c(
      "{.arg {(.arg)}} must be be a name of a column",
      i = "got {.code {as_label(data)}}"
    ))

    has_name(data, taxon_id) || cli::cli_abort(c(
      "no column named {.var {taxon_id}} in the data",
      i = "{.code {(.arg)} = {as_label(taxon_id)} }"
    ))
  } else {
    # try to guess the taxon id

    # not relevant if no taxonomy is specified
    if (is.null(taxonomy)) return(NULL)

    # try to guess the taxon id
    overlaps <- lapply(as.list(data), function(values) {
      tryCatch(sum(vec_in(values, rownames(taxonomy), na_equal = FALSE)), error = function(cnd) 0L)
    })

    best <- which.max(overlaps)
    overlaps[[best]] > 0 || cli::cli_abort(c(
      "unable to automatically detect the column with the taxon id",
      i = "please specify {.code {(.arg)} = <column name>}"
    ))

    taxon_id <- names(data)[[best]]
    cli::cli_warn(c(
      "!" = "in {.fun densify}: using column {.arg {taxon_id}} as taxon id",
      i = "specify {.code {(.arg)} = <column name>} to silence this warning"
    ), call = caller_call())
  }

  # check that taxon id does not contain any NAs
  !vec_any_missing(data[[taxon_id]]) || {
    missing <- as.character(which(vec_detect_missing(data[[taxon_id]])))
    missing <- cli::cli_vec(missing, list("vec-trunc" = 5))

    cli::cli_abort(c(
      "taxon id {.arg {taxon_id}} contains missing values",
      i = "at location{?s} {.field {missing}}"
    ))
  }


  taxon_id
}

check_limits <- function(limits, nrows, ncols, ..., .arg = caller_arg(limits)) {
  local_error_call(caller_env())
  .arg <- .arg

  is.null(limits) || is.list(limits) || cli::cli_abort(c(
    "{.arg {(.arg)}} must be a list",
    i = "got {.code {as_label(limits)}}"
  ))

  known_limits <- names(eval(formals(densify)$limits))

  all(names(limits) %in% known_limits) || cli::cli_abort(c(
    "{.arg {(.arg)}} contains invalid parameter{?s} {.field {setdiff(names(limits), known_limits)}}",
    i = "valid parameter names are  {.field {known_limits}}"
  ))

  is.null(limits$min_coding_density) ||
  is.na(limits$min_coding_density) ||
  is_between(limits$min_coding_density, 0, 1) || cli::cli_abort(c(
    "{.arg {(.arg)}$min_coding_density} must be a numeric scalar in the range [0, 1]",
    i = "got {.code {as_label(limits$min_coding_density)}}"
  ))
  if(!is_true(is.finite(limits$min_coding_density))) limits$min_coding_density <- 0


  is.null(limits$min_prop_rows) ||
  is.na(limits$min_prop_rows) ||
  is_between(limits$min_prop_rows, 0, 1) || cli::cli_abort(c(
    "{.arg {(.arg)}$min_prop_rows} must be a numeric scalar in the range [0, 1]",
    i = "got {.code {as_label(limits$min_prop_rows)}}"
  ))
  limits$min_prop_rows <- limits$min_prop_rows*nrows
  if(!is_true(is.finite(limits$min_prop_rows))) limits$min_prop_rows <- 0

  is.null(limits$min_prop_cols) ||
  is.na(limits$min_prop_cols) ||
  is_between(limits$min_prop_cols, 0, 1) || cli::cli_abort(c(
    "{.arg {(.arg)}$min_prop_cols} must be a numeric scalar in the range [0, 1]",
    i = "got {.code {as_label(limits$min_prop_cols)}}"
  ))
  limits$min_prop_cols <- limits$min_prop_cols*ncols
  if(!is_true(is.finite(limits$min_prop_cols))) limits$min_prop_cols <- 0

  limits
}


check_scoring_weights <- function(weights, has_taxonomy, ..., .arg = caller_arg(weights)) {
  local_error_call(caller_env())
  .arg <- .arg

  # check if the weights has been supplied by the user or whether it's a default value
  # we need this to work around default taxonomic weight dependent on taxonomy argument
  if (!.arg %in% call_args_names(sys.call(-1L))) {
    weights <- list(coding = 1, taxonomy = ifelse(has_taxonomy, 1, 0))
  }

  is.null(weights) || is.list(weights) || is_named(weights) || cli::cli_abort(c(
    "{.arg {(.arg)}} must be a named list",
    i = "got {.code {as_label(weights)}}"
  ))

  known_weights <- c("coding", "taxonomy")

  all(names(weights) %in% known_weights) || cli::cli_abort(c(
    "{.arg {(.arg)}} contains invalid parameter{?s} {.field {setdiff(names(weights), known_weights)}}",
    i = "valid parameter names are  {.field {known_weights}}"
  ))

  is_null(weights$coding) ||
  is_na(weights$coding) ||
  is_between(weights$coding, 0, 1) || cli::cli_abort(c(
    "{.arg {(.arg)}$coding} must be a numeric scalar in the range [0, 1]",
    i = "got {.code {as_label(weights$coding)}}"
  ))
  if (!has_name(weights, "coding")) weights$coding <- 1
  if(!is_true(is.finite(weights$coding))) weights$coding <- 0

  is_null(weights$taxonomy) ||
  is_na(weights$taxonomy) ||
  is_between(weights$taxonomy, 0, 1) || cli::cli_abort(c(
    "{.arg {(.arg)}$taxonomy} must be a numeric scalar in the range [0, 1]",
    i = "got {.code {as_label(weights$taxonomy)}}"
  ))

  if (!has_name(weights, "taxonomy")) weights$taxonomy <- ifelse(has_taxonomy, 1, 0)
  if(!is_true(is.finite(weights$taxonomy))) weights$taxonomy <- 0

  if (!has_taxonomy && weights$taxonomy != 0) cli::cli_abort(c(
    "{.arg {(.arg)}$taxonomy} weight is non-zero, but no taxonomy supplied",
    i = "taxonomy must be supplied to perform taxonomic scoring"
  ))

  weights
}

get_scoring_function <- function(method) {
  switch(method,
    arithmetic = row_scores_arithmetic,
    geometric = row_scores_geometric,
    log_odds = row_scores_log_odds,
    cli::cli_abort("unknown scoring method {.val {as_label(method)}}", .internal = TRUE)
  )
}

match_taxa <- function(
  data,
  taxonomy,
  taxon_id,
  ...,
  .data_arg = caller_arg(data),
  .taxonomy_arg = caller_arg(taxonomy)
) {
  local_error_call(caller_env())
  .data_arg  <- .data_arg
  .taxonomy_arg <- .taxonomy_arg

  # locate taxa in the taxonomy
  matches <- tryCatch(
    vec_match(data[[taxon_id]], rownames(taxonomy)),
    vctrs_error_ptype2 = function(.) NULL
  )
  is.integer(matches) || cli::cli_abort(c(
    "data uses incompatible taxa identifiers",
    i = "{.arg {taxon_id}} (of type {.cls {type_label(data[[taxon_id]])}}) is not comparable to character taxa ids"
  ))

  !anyNA(matches) || {
    missing <- as.character(data[[taxon_id]][is.na(matches)])
    missing <- cli::cli_vec(missing, list("vec-trunc" = 5))

    cli::cli_abort(c(
      "unknown taxa id{?s} {.val {missing}}",
      i = "{.arg {taxon_id}} contains taxa ids not present in the taxonomy"
    ))
  }

  # return reordered and trimmed taxonomy
  vec_slice(taxonomy, matches)
}

# check that the limit conditions don't apply
check_densify_limit_condition <- function(result, limits) {
  last <- nrow(result)

  !any(
    result$coding_density[[last]] >= limits$min_coding_density,
    result$n_rows[[last]] <= limits$min_prop_rows,
    result$n_cols[[last]] <= limits$min_prop_cols,
    na.rm = TRUE
  )
}

