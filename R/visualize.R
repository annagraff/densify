#' Visualize the distribution of quality scores across densify results
#'
#' Produces a `ggplot2` visualization showing the distribution of a quality score across the iteration steps
#' of the densification process and shows the result with the maximum score. The user can provide their own
#' quality score formula in `scoring_function`, which can use any statistics provided in a [densify_result] tibble.
#'
#' @param x a [densify_result] object
#' @param ... other arguments passed to methods
#' @param scoring_function an expression that computes a densification quality score. The default maximizes the product of
#'   the number of data points and the coding density. See [densify_result] for available statistics that can
#'   be used to compute a suitable scoring_function
#'
#' @return a `ggplot2` plot
#'
#' @importFrom generics visualize
#' @export
visualize.densify_result <- function(x, ..., scoring_function = n_data_points * coding_density) {
  requireNamespace("ggplot2", quietly = TRUE) || cli::cli_abort(c(
      "{.fn visualize} requires {.pkg ggplot2} to be installed",
      i = "please install {.pkg ggplot2} using {.code install.packages('ggplot2')}"
    ))

  n_data_points <- coding_density <- step <- NULL

  scoring_function <- eval_densify_results_quality_score(x, (score_quo <- enquo(scoring_function)), .caller = "visualize")

  if (all(is.na(scoring_function))) {
    cli::cli_warn(c(
      "in {.fn visualize}: no scoring_function to visualize",
      if(anyNA(scoring_function)) i = "{.code {as_label(score_expr)}} produced NAs"
    ))

    return(NULL)
  }

  # visualization frame
  frame <- data_frame(idx = seq_len(nrow(x)), step = x$.step, scoring_function = scoring_function)

  # maximal value
  best <- which.max(frame$scoring_function)
  best_score <- frame$scoring_function[best]
  best_step <- frame$step[best]

  score_label <- paste("step ", best_step, "  \ndata frame [",x$n_rows[best_step],",", x$n_cols[best_step],"]", sep="")

  cli::cli_alert_info("use {.code tibble::as_tibble({caller_arg(x)}$data[[{frame$idx[[best]]}]])} to obtain the pruned data frame")

  ggplot2::ggplot(vec_slice(frame, vec_detect_complete(frame)), ggplot2::aes(x = step, y = scoring_function)) +
    ggplot2::geom_path(col="steelblue") +
    ggplot2::geom_vline(ggplot2::aes(xintercept = best_step), linetype = "dashed", col = "goldenrod1") +
    ggplot2::ylab(as_label(score_quo)) +
    ggplot2::scale_x_continuous(limits=range(frame$step), expand = c(0,0)) +
    ggplot2::annotate("text", x = best_step, y = best_score, label = score_label, hjust = "inward", size = 3) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45)) +
    ggplot2::theme_bw()
}


#' @rdname visualize.densify_result
#' @export
plot.densify_result <- visualize.densify_result

#' @importFrom generics visualize
#' @export
generics::visualize

