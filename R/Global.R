#' Global Definition
#'

#' @docType package
#' @name NetSci
#' @importFrom  utils install.packages
#' @description Basic global variables to make sure the package runs.

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1") utils::globalVariables(c(".", "bin", "ID", "Target"))




