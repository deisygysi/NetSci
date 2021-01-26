#' Proximity from target to source
#' @description Calculates the proximity (average or closest) from source to targets.
#' @param G The original graph (often an interactome).
#' @param source nodes from the network (in a drug repurpusing set-up those are the disease genes)
#' @param targets targets in the network (in a drug repurpusing set-up those are the drug-targets)
#'
#' @importFrom magrittr `%>%` `%<>%`
#' @importFrom  igraph graph_from_data_frame bipartite_mapping degree V E as_incidence_matrix induced_subgraph


#' @return the proximity value for the source-targets
#' @export
#'
#' @examples
#' #' set.seed(666)
#' net  = data.frame(
#' Node.1 = sample(LETTERS[1:15], 15, replace = TRUE),
#' Node.2 = sample(LETTERS[1:10], 15, replace = TRUE))
#' net$value = 1
#' net =  CoDiNA::OrderNames(net)
#' net = unique(net)
#'
#' g <- igraph::graph_from_data_frame(net, directed = FALSE )
#' T = c("G", "A", "D")
#' S = c("C", "M")
#' proximity_average(g, source = S, targets = T)


proximity_average <- function(G, source, targets){
  `%<>%` <- magrittr::`%<>%`
  `%>%` <- magrittr::`%>%`

  source = source [source %in% V(G)$name] %>% unique
  targets = targets [targets %in% V(G)$name]%>% unique
  out = NA
  if(length(source)> 0 & length(targets)> 0){
    d <-  igraph::distances(G, source, targets)
    out = apply(d, 2, min) %>% mean()
  }
  return(out)
}


