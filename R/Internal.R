#' Internal
#' @description  Internal functions. Not exported.
#'
#' @param PPIg
#' @param bins
#' @param nodes
#' @param n
#' @param min_per_bin
#' @importFrom magrittr `%>%` `%<>%`
#' @importFrom stats quantile ecdf
#' @importFrom dplyr summarise group_by n
#' @importFrom  igraph graph_from_data_frame bipartite_mapping degree V E as_incidence_matrix induced_subgraph
#' @keywords internal
#'
#'



LCC_Calc = function(PPIg,bins = 100, nodes, n, min_per_bin = 20){
  `%<>%` <- magrittr::`%<>%`
  `%>%` <- magrittr::`%>%`
  # bin = NULL
  DG = igraph::degree(PPIg) %>% as.data.frame()
  BINS = binr::bins(DG$., target.bins = bins,
                    exact.groups = FALSE,
                    minpts = min_per_bin)
  numbers_factors = cut(DG$., binr::bins.getvals(BINS),
                        labels = names(BINS$binct)) %>% as.numeric()
  DG$bin = numbers_factors
  DG$names = row.names(DG)
  nodes_ppi = subset(DG, DG$names %in% nodes)
  bins_get = nodes_ppi %>%
    dplyr::group_by(bin) %>%
    dplyr::summarise(., n = dplyr::n()) %>% as.data.frame()

  OUT = list()
  for ( j in 1:n){
    resampled = list()
    for ( i in 1:nrow(bins_get)){
      in_the_bin = subset(DG$names, DG$bin %in% bins_get$bin[i])
      test= sample(in_the_bin, size = bins_get$n[i])
      resampled[[i]] = test
    }
    new_nodes = unlist(resampled)

    OUT[[j]] = LCC_number(PPIg, new_nodes)
  }
  return(OUT)
}

LCC_p = function(n = 1000, PPI_g, Vn, bins = NULL, min_per_bin){
  `%<>%` <- magrittr::`%<>%`
  `%>%` <- magrittr::`%>%`

  Out = vector()
  if(is.null(bins)){
    Vn = length(Vn)
    for (i in 1:n) {
      Out[i]=LCC_aux(PPI_g, Vn)
    }
  } else{
    Out =  LCC_Calc(PPIg = PPI_g, bins = bins, nodes = Vn, n, min_per_bin)
    Out = unlist(Out)
  }
  return(Out)
}


LCC_number = function(PPIg, v){
  `%<>%` <- magrittr::`%<>%`
  `%>%` <- magrittr::`%>%`

  C_LCC = igraph::induced_subgraph(PPIg, v = v)
  LCC_total = igraph::components(C_LCC)$csize %>% max
  return(LCC_total)
}


LCC_aux = function(PPIg, Vn){
  `%<>%` <- magrittr::`%<>%`
  `%>%` <- magrittr::`%>%`

  v = sample(x = igraph::V(PPIg)$name, size = Vn)
  OUT = LCC_number(PPIg, v)
  return(OUT)
}


LCC = function(g){
  `%<>%` <- magrittr::`%<>%`
  `%>%` <- magrittr::`%>%`

  x = igraph::components(g)$csize %>% max()
  return(x)
}


### Functions for proximity
###

bins <- function(g, bins, min_bin , nodes){
  DG <- igraph::degree(g) %>% as.data.frame()
  names(DG) <- "degree"
  BINS <- binr::bins(DG$degree, target.bins = bins,
                     exact.groups = FALSE,
                     minpts = min_bin)
  numbers_factors <- cut(DG$degree, binr::bins.getvals(BINS),
                         labels = names(BINS$binct)) %>% as.numeric()
  DG$bin <- numbers_factors
  DG$names <- row.names(DG)
  nodes_ppi <- subset(DG, DG$names %in% nodes)
  bins_get <- nodes_ppi %>%
    dplyr::group_by(bin) %>%
    dplyr::summarise(., n = dplyr::n()) %>% as.data.frame()
  return(list(bins_get, DG))
}


boot_distance_aux <- function(n, bins_get, g, source, targets, DG, proximity){
  OUT <- list()
  for ( j in 1:n){
    resampled <- list()
    for ( i in 1:nrow(bins_get)){
      in_the_bin <- subset(DG$names, DG$bin %in% bins_get$bin[i])
      resampled[[i]] <- sample(in_the_bin, size = bins_get$n[i])
    }
    new_nodes <- unlist(resampled)
    if(proximity == "close"){
      OUT[[j]] <- proximity_close(g, new_nodes, targets)
    }
    else if(proximity == "average"){
      OUT[[j]] <- proximity_average(g, new_nodes, targets)
    }
  }
  return(OUT)
}

boot_distance = function(g, source, targets, bins, min_bin, n, proximity){
  BIN <- bins(g, bins, min_bin, nodes = targets)
  bins_get <- BIN[[1]]
  DG <- BIN[[2]]
  rm(BIN)

  x <- boot_distance_aux(n, bins_get, g, source, targets, DG, proximity) %>% unlist()
  return(x)
}

boot_distance_IC= function(x, alpha, d){
  IC = stats::quantile(x, c(1-alpha/2, alpha/2))

  d_fun <- stats::ecdf (x)
  p_gt = 1 - d_fun(d)
  p_lt = d_fun(d)
  Z = (d - mean(x))/sd(x)
  return(list(IC = IC, p_gt = p_gt, p_lt = p_lt, Z = Z))
}

boot_distance_p = function(g, source, targets, bins, min_bin, n, d, alpha, proximity){
  x = boot_distance(g, source, targets, bins, min_bin, n, proximity)
  y = boot_distance_IC(x, alpha, d)

  return(y)
}


