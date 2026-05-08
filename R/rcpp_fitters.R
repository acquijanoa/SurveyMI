#' Check if native weighted fitters are available
#'
#' @return Logical scalar indicating whether compiled Rcpp fitters are loaded.
#' @export
has_rcpp_fitters <- function() {
  all(vapply(
    c("weighted_gaussian_cpp", "weighted_binomial_cpp", "weighted_poisson_cpp"),
    exists,
    logical(1),
    mode = "function"
  ))
}
