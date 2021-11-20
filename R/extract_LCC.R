#' Extract LCC from a graph
#'
#' @param g is the graph you want to extract the largest connected component
#'
#' @return a graph (from igraph) with only the largest connected component
#' @export
#'
#' @examples
#' set.seed(12)
#' x = data.frame(n1 = sample(LETTERS[1:5]),
#'                n2 =  sample(LETTERS[1:20]))
#'
#' g = igraph::graph_from_data_frame(x, directed = F)
#' g = simplify(g) 
#' LCC = extract_LCC(g)
#' 
#' 
extract_LCC = function(g){
  mem = g %>% components()
  mem = mem$membership %>% 
    as.data.frame() 
  
  names(mem) = c("cluster")
  mem$nodes = row.names(mem)
  
  keep = mem %>% 
    group_by(cluster) %>%
    mutate(n = n()) %>%
    ungroup() %>% 
    filter(n == max(n)) %>%
    pull(nodes)
  
  g %<>% induced_subgraph(., keep)
  return(g)
}