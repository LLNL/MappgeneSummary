#' Calculate Hill numbers of order 1, 2, and Inf from a matrix
#'
#' For both Hill 1 and 2, a maximally variant site will have a value of \eqn{n}
#' (which is 4 for the 4 possible nucleotides), and maximally invariant site is
#' 1.
#'
#' @param x a matrix with nrows samples and ncol counts
#'
#' @return a data.frame with Hill's N_1, N_2, N_Inf for each sample
#' @details log(N_1) is entropy, log(N_2) is inverse Simpson's
#'
#' @export

hill_12I <- function(x) {
  vegan::renyi(x, scales = c(1, 2, Inf), hill = TRUE)
}




#' Calculate the Simpson index for a set of frequencies
#'
#' A maximally variant site will have a value of 1, and a maximally invariant
#' site will have a value of \eqn{1/n} (or 0.25 for the 4 possible nucleotides).
#'
#' @param x a numeric vector
#'
#' @return the Simpson index
#'
#' @details The Simpson index is the probability that two bases are the
#'          same. Aka the reciprocal of N_2. Aka π for 1 nt.
#' @seealso vegan::diversity
#' @export

simpson <- function(x) {
  1 / hill_12I(x)["2"]
}




#' Calculate the Gini-Simpson index for a set of frequencies
#'
#' @param x a numeric vector
#'
#' @return the Gini-Simpson index
#'
#' @details The Gini-Simpson is 1 minus the probability that two bases are the
#'          same. Aka 1 minus Simpson's index. Aka 1 minus 1/N_2
#' @seealso vegan::diversity
#' @export

gini_simpson <- function(x) {
  1 - simpson(x)
}
