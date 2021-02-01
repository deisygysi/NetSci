#' @title Histogram_LCC
#' @description Plots the histogram to evaluate the significance of the Largest Connected Component (LCC).
#' @param LCC_L an output from the function LCC_Significance
#' @param Name title of the plot
#' @importFrom graphics abline hist title
#' @importFrom stats ecdf sd
#' @return An Histogram of the simulated LCC, and a red line of the actual LCC.
#' @export
#'
#' @examples
#' set.seed(666)
#' net  = data.frame(
#' Node.1 = sample(LETTERS[1:15], 15, replace = TRUE),
#' Node.2 = sample(LETTERS[1:10], 15, replace = TRUE))
#' net$value = 1
#' net =  CoDiNA::OrderNames(net)
#' net = unique(net)
#'
#' g <- igraph::graph_from_data_frame(net, directed = FALSE )
#' targets = c("N", "A", "I", "F")
#' LCC_Out = LCC_Significance(N = 1000,
#'                  Targets = targets,
#'                                   G = g,
#'                                   bins = 5,
#'                                   min_per_bin = 2)
#'                                   # in a real interactome, please use the default
#'
#' Histogram_LCC(LCC_Out, "Example")


Histogram_LCC = function(LCC_L, Name = NULL){
  Fn = ecdf(LCC_L$LCCZ)
  p = 1 - Fn(LCC_L$LCC)

  lim = c(LCC_L$LCC, LCC_L$LCCZ)
  LCC_L[[1]]%>% hist(., las = 1,
                     main = "",
                     xlim = c(min(lim -10), max(lim + 10)),
                     col = 'gray75', ylab = "")
  abline(v = LCC_L$LCC, col = "red")
  title(main = Name, sub = paste0("LCC: ",
                                  round(LCC_L$LCC,0),
                                  " (",
                                  round(LCC_L$mean,2),
                                  " +- ",
                                  round(LCC_L$sd,2),"; ",
                                  "p: " ,round(p,2),
                                  ")"))


}
