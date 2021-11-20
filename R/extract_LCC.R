#' Extract LCC from a graph
#'
#' @param g is the graph you want to extract the largest connected component
#'
#' @return a graph (from igraph) with only the largest connected component
#' @export
#' @importFrom magrittr `%>%` `%<>%`
#' @importFrom igraph shortest.paths distances graph_from_data_frame bipartite_mapping degree V E as_incidence_matrix induced_subgraph
#' @importFrom igraph components induced.subgraph
#' @importFrom dplyr  group_by mutate  ungroup filter pull
#'

#'
#' @examples
#' set.seed(12)
#' x = data.frame(n1 = sample(LETTERS[1:5]),
#'                n2 =  sample(LETTERS[1:20]))
#'
#' g = igraph::graph_from_data_frame(x, directed = FALSE)
#' g = igraph::simplify(g)
#' LCC = extract_LCC(g)
#'
#'
extract_LCC = function(g){
  mem = g %>% igraph::components()
  mem = mem$membership %>%
    as.data.frame()

  names(mem) = c("cluster")
  mem$nodes = row.names(mem)

  keep = mem %>%
    dplyr::group_by(cluster) %>%
    dplyr::mutate(n = n()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(n == max(n)) %>%
    dplyr::pull(nodes)

  g %<>% igraph::induced_subgraph(., keep)
  return(g)
}
