#' LCC Significance
#' @description Calculates the Largest Connected Component (LCC) from a given graph, and calculates its significance using a degree preserving approach.
#' Menche, J., et al (2015) <doi.org:10.1126/science.1065103>
#'
#' @param N Number of randomizations.
#' @param Targets Name of the nodes that the subgraph will focus on - Those are the nodes you want to know whether if forms an LCC.
#' @param G The  graph of interest (often, in NetMed it is an interactome - PPI).
#' @param bins the number os bins for the degree preserving randomization. When bins = 1, assumes a uniform distribution for nodes.
#' @param min_per_bin the minimum size of each bin.
#' @param hypothesis are you expecting an LCC greater or smaller than the average?
#'
#' @importFrom graphics abline hist title
#' @importFrom stats ecdf sd
#' @importFrom magrittr `%>%` `%<>%`
#' @importFrom  igraph graph_from_data_frame bipartite_mapping V E as_incidence_matrix
#' @return a list with the LCC
#' - $LCCZ all values from the randomizations
#' - $mean the average LCC of the randomizations
#' - $sd the sd LCC of the randomizations
#' - $Z The score
#' - $LCC the LCC of the given targets
#' - $emp_p the empirical p-value for the LCC
#' - $rLCC the relative LCC
#'
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
#' plot(g)
#' targets = c("I", "H", "F", "E")
#'  LCC_Significance(N = 100,
#'                    Targets = targets,
#'                                     G = g,
#'                                   bins = 1,
#'                                   min_per_bin = 2)
#'
#'

LCC_Significance = function(N = N, Targets = Targets, G,
                            bins =100, hypothesis = "greater",
                            min_per_bin = 20){
  LCC_1 = LCC_number(G, Targets)
  Vn = Targets
  LCC_total = LCC_1
  LCCZ = LCC_p(PPI_g = G, Vn = Vn, n = N, bins = bins, min_per_bin)
  muC = mean(LCCZ)
  sdC = sd(LCCZ)
  Z = (LCC_1-muC)/sdC
  p = pvals(x = LCCZ, val = LCC_1 )
  # Fn =stats::ecdf(LCCZ)
  # d = density(LCCZ)
  if(hypothesis == "greater"){
    p = p$p_gt
    # p = 1 - Fn(LCC_1)
    # p2 = sum(d$y[d$x>LCC_1])
  } else  if(hypothesis == "lower"){
    p = p$p_lt
    # p = Fn(LCC_1)
    # p2 = sum(d$y[d$x<LCC_1])
  }

  return(list(LCCZ = LCCZ, mean = muC,
              sd = sdC,
              Z= Z,
              LCC = LCC_1,
              emp_p = p,
              rLCC = LCC_1/length(Targets)))
}
