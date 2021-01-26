#' @title Hypergeometric.test
#'
#' @description Calculates the significance of an overlap of two sets using an hypergeometric test.
#'    It is a wrapper of the `phyper` function.
#' @param success Is the number of elements in the overlap of the sets.
#' @param universe_success Is the number of elements of the set of interest.
#' @param universe_failure Is the number of elements of the set of the other set.
#' @param size_collected The total of elements in the universe
#' @param lower.tail Should the test be calculated on the lower tail? (Hypothesis test is lower than)
#'
#' @return the p-value for the hypergeometric test.
#' @importFrom stats phyper
#' @export
#'
#' @examples
#' require(magrittr)
#' s = 10; S = 15; f = 10; T = 30
#' Hypergeometric.test(success = s,
#' universe_success = S,
#' universe_failure = f,
#' size_collected = T
#' )



Hypergeometric.test = function(success, universe_success, universe_failure, size_collected, lower.tail = FALSE){
  stats::phyper(q = success, m= universe_success, n = universe_failure, k = size_collected, lower.tail = TRUE, log.p = FALSE)
}
