#' @title Separation
#' @description Calculates the separation of two set of targets on a network.
#' Often used to measure separation of disease modules in a interactome.
#' Separation is calculated as in Menche, J. et al (2015) <doi:10.1126/science.1257601>.
#' @param G The original graph (often an interactome).
#' @param ST Set-Target data. It is a data.frame with two columns. ID and Target.
#'
#' @importFrom magrittr `%>%` `%<>%`
#' @importFrom igraph shortest.paths graph_from_data_frame bipartite_mapping degree V E as_incidence_matrix induced_subgraph
#' @importFrom dplyr filter

#' @return the separation and distance of modules.
#' @export


#' @example
#' set.seed(12)
#' x = data.frame(n1 = sample(LETTERS[1:5]),
#'                n2 =  sample(LETTERS[1:20]))
#'
#' D1 = data.frame(gene = c("H", "I", "S", "N", "A"), disease = "D1")
#' D2 = data.frame(gene = c("E", "C",  "R" , "J", "Q", "O"), disease = "D2")
#' D3 = data.frame(gene = c("E", "G", "T", "P"), disease = "D3")
#' D4 = data.frame(gene = c("A", "B", "E"), disease = "D4")
#' D5 = data.frame(gene = c("D", "F", "L"), disease = "D5")
#' D6 = data.frame(gene = c("D", "F", "K"), disease = "D6")
#' D7 = data.frame(gene = c("A", "B" ,"F", "K"), disease = "D7")
#'
#' Diseases = rbind(D1, D2, D3, D4, D5, D6, D7)
#' Diseases %<>% dplyr::select(disease, gene)
#' g = igraph::graph_from_data_frame(x, directed = F)
#' g = simplify(g)
#'
#' separation(G = g, ST = Diseases)


separation = function(G, ST){
  names(ST)[1:2] = c("ID", "Target")
  ST$Target %<>% as.character()
  ST %<>% dplyr::filter(Target %in% V(G)$name)
  d = ST$ID %>% unique()
  asp = igraph::shortest.paths(G, v=unique(ST$Target),
                               to = unique(ST$Target))

  # S_aa = list()
  S_ab = S_ab.tmp = matrix(NA, ncol = length(d), nrow = length(d))
  pb <- txtProgressBar(min = 0, max = (length(d)), style = 3)
  for(i in 1:length(d)){
    setTxtProgressBar(pb, i)
    g_a = ST$Target[ST$ID == d[i]]
    # S_aa[[i]] = asp[g_a, g_a][upper.tri(asp[g_a, g_a])] %>% mean()
    tmp = asp[g_a, g_a]#[upper.tri(asp[g_a, g_a], diag = F)]
    diag(tmp) <- Inf
    tmp = apply(tmp,1,min)
    tmp = tmp[!is.infinite(tmp)] %>% mean
    S_ab.tmp[i,i] = tmp
    # S_ab.tmp[i,i] = tmp[!is.infinite(tmp)] %>% mean
    j = i + 1
    # for(j in (i + 1):length(d)){
    #   if(j <= length(d)){
    while(j <= length(d)){
      g_b = ST$Target[ST$ID == d[j]]
      tmp = asp[g_a, g_b]#[upper.tri(asp[g_a, g_a], diag = F)]
      # tmp = tmp[!is.infinite(tmp)]
      # diag(tmp) <- Inf
      t1 = apply(tmp,1,min)
      t2 = apply(tmp,2,min)
      # tmp = sum(t1[!is.infinite(t1)]) + sum(t2[!is.infinite(t2)])
      # tmp = tmp/(length(g_a)+length(g_b))
      tmp = c(t1, t2)
      tmp = mean(tmp[!is.infinite(tmp)])
      S_ab.tmp[i,j] = tmp
      j = j + 1
    }

    # S_ab.tmp[i,j] =  tmp[!is.infinite(tmp)] %>% mean
    # }
    # }
  }
  close(pb)
  message("Calculating S_ab...\n")
  ### calculate the sabs
  pb <- txtProgressBar(min = 0, max = (length(d)), style = 3)
  for(i in 1:length(d)){
    setTxtProgressBar(pb, i)
    for(j in 1:length(d)){
      S_ab[i,j] = S_ab.tmp[i,j] - mean(c(S_ab.tmp[i,i], S_ab.tmp[j,j]))
    }
  }
  close(pb)
  message("Done...\n")
  colnames(S_ab)= rownames(S_ab) = d
  colnames(S_ab.tmp)= rownames(S_ab.tmp) = d
  return(list(Sab = S_ab, Dab = S_ab.tmp))
}

