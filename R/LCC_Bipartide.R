#' LCC_Bipartide
#' @description LCC size for Rewire Bipartide Network
#' @param g Bipartide Graph to be rewired
#' @param N Number of resamples
#' @importFrom magrittr `%>%` `%<>%`
#' @importFrom  igraph components graph_from_data_frame bipartite_mapping degree V E as_incidence_matrix induced_subgraph

#'
#' @return a list with the LCC
#' - $LCCZ all values from the randomizations
#' - $mean the average LCC of the randomizations
#' - $sd the sd LCC of the randomizations
#' - $Z The Z-score
#' - $LCC the LCC of the original network
#' - $emp\_p the empirical p-value for the LCC
#' @export
#'
#' @examples
#' nodes = data.frame(c("D1", "D2", "D3", "D4", "D5",
#'                      "G1", "G2", "G3", "G4", "G6"),
#'                    type = c(TRUE, TRUE, TRUE, TRUE, TRUE,
#'                             FALSE, FALSE, FALSE, FALSE, FALSE))
#'
#' g = data.frame(from = c("D1", "D2", "D2", "D3", "D4", "D5", "D5", "D3"),
#'                to =   c("G1", "G1", "G2", "G3", "G4", "G4", "G1", "G6"))
#'
#' g = igraph::graph_from_data_frame(g, directed = FALSE, vertices = nodes)
#'
#' plot(g, layout = igraph::layout.bipartite)
#'
#' igraph::components(g)
#' LCC_BIP = LCC_Component(g)
#' Histogram_LCC(LCC_BIP)

LCC_Component = function(g, N = 1000){
  original =  igraph::components(g)
  cmp_o = original$csize %>% max
  LCC = list()

  for(i in 1:N){
    cmp = BiRewire::birewire.rewire.bipartite(g, verbose = FALSE) %>%
      igraph::components()
    LCC[[i]] = cmp$csize %>% max()
  }
  LCC %<>% unlist()

  muC = mean(LCC)
  sdC = sd(LCC)
  Z = (cmp_o-muC)/sdC
  p = pvals(x = LCC, val = cmp_o )

  out = list(LCC = cmp_o,
    LCCZ = LCC,
             mean = mean(LCC),
             sd = sd(LCC),
             Z = Z,
    emp_p = p$p_gt
               )
  return(out)
}
