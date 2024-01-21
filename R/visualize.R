#' Visualise the distrubution of quality scores across densify results
#'
#' Produces a `ggplot2` visualisation showing the distribution of a quality score across the iteration steps
#' of the densification process and shows the result with the maximal score. The user can provide their own 
#' quality score formula in `score`, which can use any statistics provided in a [densify_result] tibble.
#'
#' @param x a [densify_result] object
#' @param score an expression that computes a densification quality score. Default score aims to maximise both
#'   the number of data points and the coding density. See [densify_result] for available statistics that can 
#'   be used to compute a suitable score
#'
#' @return a `ggplot2` plot
#'
#' @importFrom generics visualize
#' @export
visualize.densify_result <- function(x, ..., score = n_data_points * coding_density) {
  requireNamespace("ggplot2", quietly = TRUE) || cli::cli_abort(c(
      "{.fn visualize} requires {.pkg ggplot2} to be installed",
      i = "please install {.pkg ggplot2} using {.code install.packages('ggplot2')}"
    ))

  score <- eval_densify_results_quality_score(x, (score_quo <- enquo(score)), .caller = "visualize")

  if (all(is.na(score))) {
    cli::cli_warn(c(
      "in {.fn visualize}: no score to visualize",
      if(anyNA(score)) i = "{.code {as_label(score_expr)}} produced NAs"
    ))

    return(NULL)
  }


  # visualization frame
  frame <- data_frame(idx = seq_len(nrow(x)), step = x$.step, score = score)
  
  # maximal value
  best <- which.max(frame$score)
  best_score <- frame$score[best]
  best_step <- frame$step[best]

  
  score_label <- cli::format_inline(
    "  max score is {prettyNum(best_score)}\n",
    "  step {best_step}\n",
    "  data frame [{x$n_rows[best_step]}, {x$n_cols[best_step]}]"
  )

  cli::cli_alert_info("use {.code tibble::as_tibble({caller_arg(x)}$data[[{frame$idx[[best]]}]])} to obtain the pruned data frame")


  ggplot2::ggplot(vec_slice(frame, vec_detect_complete(frame)), ggplot2::aes(x = step, y = score, color = score)) + 
    ggplot2::geom_path() +
    ggplot2::geom_vline(ggplot2::aes(xintercept = best_step), linetype = "dashed", col = "goldenrod1") +
    ggplot2::ylab(as_label(score_quo)) +
    ggplot2::geom_text(ggplot2::aes(x = best_step, y = best_score, label = score_label), hjust = "inward", family = "sans", size = 3, fontface = "plain") +
    ggplot2::guides(color = ggplot2::guide_colourbar(title = "")) +
    ggplot2::scale_x_continuous(n.breaks = 20L, limits = range(frame$step)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45))
}

#' @rdname visualize.densify_result
#' @export
plot.densify_result <- visualize.densify_result
