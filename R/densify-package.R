#' @keywords internal
#' @import rlang
#' @import vctrs
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL


.onLoad <- function(...) {
  s3_register("dplyr::select", "densify_result")
  s3_register("dplyr::rename", "densify_result")
}