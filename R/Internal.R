#' Internal
#' @description  Internal functions. Not exported.
#'
#' @param PPIg
#' @param bins
#' @param nodes
#' @param n
#' @param min_per_bin
#' @importFrom magrittr `%>%` `%<>%`
#' @importFrom stats quantile ecdf approxfun integrate
#' @importFrom dplyr summarise group_by n
#' @importFrom  igraph graph_from_data_frame bipartite_mapping degree V E as_incidence_matrix induced_subgraph
#' @keywords internal
#'
#'



LCC_Calc = function(PPIg,
                    bins = 100,
                    nodes,
                    n,
                    min_per_bin = 20){
  `%<>%` <- magrittr::`%<>%`
  `%>%` <- magrittr::`%>%`
  # bin = NULL
  DG = igraph::degree(PPIg) %>% as.data.frame()
  BINS = binr::bins(DG$., target.bins = bins,
                    exact.groups = FALSE,
                    minpts = min_per_bin)
  numbers_factors = cut(DG$.,
                        binr::bins.getvals(BINS),
                        labels = names(BINS$binct)) %>%
    as.numeric()
  DG$bin = numbers_factors
  DG$names = row.names(DG)
  nodes_ppi = subset(DG, DG$names %in% nodes)
  bins_get = nodes_ppi %>%
    dplyr::group_by(bin) %>%
    dplyr::summarise(., n = dplyr::n()) %>%
    as.data.frame()

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

LCC_p = function(n = 1000,
                 PPI_g,
                 Vn,
                 bins = NULL,
                 min_per_bin){
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


boot_distance_aux <- function(n, bins_get, g,
                              source, targets, DG, proximity){
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
    }    else if(proximity == "weighted"){
      OUT[[j]] <- proximity_average_weighted(g, new_nodes, targets)
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

  # d_fun <- stats::ecdf (x)
  # p_gt = 1 - d_fun(d)
  # p_lt = d_fun(d)
  #
  # dens = stats::density(x)
  # p2_gt = sum(dens$y[dens$x>d])
  # p2_lt = sum(dens$y[dens$x<d])

  p = pvals(x = x, val = d)

  Z = (d - mean(x))/sd(x)
  return(list(IC = IC,
              # p_gt = p_gt,
              # p_lt = p_lt,
              p_gt = p$p_gt,
              p_lt = p$p_lt,
              Z = Z))
}

boot_distance_p = function(g, source, targets, bins, min_bin, n, d, alpha, proximity){
  x = boot_distance(g, source, targets, bins, min_bin, n, proximity)
  y = boot_distance_IC(x, alpha, d)

  return(y)
}

#############################

pvals = function(x, val){

  d = stats::density(x)

  xx <- d$x
  dx <- xx[2L] - xx[1L]
  yy <- d$y

  f <- stats::approxfun(xx, yy, rule = 2)
  C <- cubature::cubintegrate(f, min(xx, val), max(xx, val), method = "pcubature")$integral
  p.unscaled <- cubature::cubintegrate(f, val, max(xx), method = "pcubature")$integral
  p.scaled_gt <- p.unscaled / C

  p.unscaled <- cubature::cubintegrate(f,  min(xx), val, method = "pcubature")$integral
  p.scaled_lt <- p.unscaled / C

  p.scaled_gt = ifelse(p.scaled_gt > 1, 1, p.scaled_gt)
  p.scaled_gt = ifelse(p.scaled_gt < 0, 0, p.scaled_gt)

  p.scaled_lt = ifelse(p.scaled_lt > 1, 1, p.scaled_lt)
  p.scaled_lt = ifelse(p.scaled_lt < 0, 0, p.scaled_lt)


  return(list(p_gt = p.scaled_gt,
              p_lt = p.scaled_lt))
}


##### Fct for separation

saa = function(g1, g2, sps){

  if (identical(g1,g2)){
    tmp = sps[g1, g2]
    diag(tmp) <- Inf
    tmp = apply(tmp,1,min)
    tmp = tmp[!is.infinite(tmp)] %>%
      mean
  } else{
    tmp = sps[g1, g2]
    t1 = apply(tmp,1,min)
    t2 = apply(tmp,2,min)
    tmp = c(t1, t2)
    tmp = tmp[!is.infinite(tmp)] %>%
      mean
  }
  return(tmp)
}



resample = function(total,
                    n){
  samples = sample(1:total,
                   size = n,
                   replace = F)
  return(samples)
}

resample_saa = function(i){
  # require(magrittr)
  # require(igraph)
  `%>%`<- magrittr::`%>%`()
  `%>%`<- magrittr::`%<>%`()
  tmp = list()
  for(n in 1:N){
    tmp[[n]] = resample(n = nodes_ID$n[i],
                        total = nnodes)
  }


  saa_star_tmp = list()
  for(runs in 1:N){
    saa_star_tmp[[runs]] = saa(tmp[[runs]],
                               tmp[[runs]],
                               sps = all_sps)
  }

  saa_original = ST$Target[ST$ID == d[i]] %>%
    saa(.,., sps = all_sps)
  saa_star_tmp %<>% unlist()
  saa_stars = saa_star_tmp %>%
    t %>%
    as.data.frame() %>%
    dplyr::mutate(Disease = d[i],
                  Saa_Dis = saa_original)

  SAMPLES = tmp %>%
    unlist %>%
    matrix(., nrow = N, byrow = F)
  return(list(saa_stars = saa_stars, SAMPLES = SAMPLES))
}

# requires = function(){
#   require(magrittr)
#   require(igraph)
#   require(dplyr)
#   require(parallel)
# }

SAB_complete = function(i){
  # require(magrittr)
  `%>%`<- magrittr::`%>%`()
  `%>%`<- magrittr::`%<>%`()
  tmp =
    Sab_tmp[i,1:N] %>%
    as.numeric %>%
    pvals(., Sab_tmp$Saa_Dis[i])

  pval = tmp$p_lt %>% as.numeric()

  SAB = Sab_tmp[i,] %>%
    dplyr::select(x,
                  y,
                  Sab = Saa_Dis) %>%
    dplyr::mutate(pvalue_lt = pval)

  return(SAB)
}


sab_aux = function(j){
  X = 0; Sab_tmp = list()
  k = j
  while(k < length(d)){
    k = k + 1
    X = X + 1
    sab_star = list()
    for(resample_id in 1:N){
      sab_star[[resample_id]] = saa(SAMPLES[[k]][resample_id,],
                                    SAMPLES[[j]][resample_id,],
                                    sps = all_sps)
    }
    sab_original =
      saa(ST$Target[ST$ID == d[j]],
          ST$Target[ST$ID == d[k]],
          sps = all_sps)

    tmp2 = sab_star %>%
      unlist() %>%
      t %>%
      as.data.frame() %>%
      dplyr::mutate(Saa_Dis = sab_original)

    tmp3 = saa_stars %>%
      dplyr::filter(Disease %in% c(d[j], d[k])) %>%
      dplyr::select(-Disease) %>%
      colMeans()

    Sab_tmp[[X]]  = (tmp2 - tmp3) %>%
      dplyr:: mutate(x = d[j],
                     y = d[k])
  }
  return(Sab_tmp = Sab_tmp %>%
           dplyr::bind_rows())
}

