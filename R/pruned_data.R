# Return a function that constructes pruned data frame promises
#
# We use this design to work around serialization issues of lazy expressions,
# R serializes references to the same data multiple times, which results in bad
# performance and very large files when writing densify_result frames to disk
#
# Using a nested environment with lazy bindings solves this
make_pruned_data_factory <- function(data) {
  parent_env <- new_environment(list(.data = data), parent = baseenv())

  function(rows, cols) {
    promise <- new_environment(list(.rows = rows, .cols = cols), parent = parent_env)
    env_bind_lazy(promise, data = vctrs::vec_slice(.data[.cols], .rows), .eval_env = promise)

    class(promise) <- "pruned_data_promise"
    promise
  }
}



#' @export
print.pruned_data_promise <- function(x, ...) {
  arg <- caller_arg(x)

  out <- cli::cli_fmt({
    cli::cli_bullets(c(
      "<promise: pruned data frame[{length(x$.rows)}, {length(x$.cols)}]>",
      "i" = "use {.code as.data.frame({arg})} or {.code as_tibble({arg})} to obtain the data frame"
    ))
  }, collapse = TRUE)
  cat(out)
}


#' @export
as.data.frame.pruned_data_promise <- function(x, ...) {
  x$data
}


#' @export
vec_ptype_abbr.pruned_data_promise <- function(x, ...) {
  sprintf("pruned df[%s, %s]", length(x$.rows), length(x$.cols))
}

#' @export
vec_ptype_full.pruned_data_promise <- function(x, ...) {
  sprintf("pruned data frame[%s, %s]", length(x$.rows), length(x$.cols))
}



