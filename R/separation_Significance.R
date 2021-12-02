#' Separation Significance
#'
#' @description Calculates the separation of two set of targets on a network and assigns a p-value to it.
#' Often used to measure separation of disease modules in a interactome.
#' Separation is calculated as in Menche, J. et al (2015) <doi:10.1126/science.1257601>.
#' p-values are calculates based on the permutation of nodes, you can set the full network to be
#' in the set for permutation or can select the ones you include as input.
#' @param G The original graph (often an interactome / PPI).
#' @param ST Set-Target data. It is a data.frame with two columns. ID and Target.
#' @param Threads How many threads you'd like to use (for parallel computation).
#' @param N = 1000 The number of permutations
#' @param correct_by_target TRUE by default. If you want to use the set of targets for the permutation or the full network.
#' @importFrom magrittr `%>%` `%<>%`
#' @importFrom igraph shortest.paths distances graph_from_data_frame bipartite_mapping degree V E as_incidence_matrix induced_subgraph
#' @importFrom dplyr filter select  group_by summarise bind_rows
#' @importFrom parallel makeCluster clusterExport clusterApplyLB stopCluster
#' @export
#'
#' set.seed(12)
#' require(magrittr)
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
#' separation_Significance(G = g, ST = Diseases)


separation_Significance =  function(G,
                                    ST,
                                    Threads = 10,
                                    N = 1000,
                                    correct_by_target = TRUE){
  # requires()
  NetSci.Sep <- new.env()
  G %<>% extract_LCC()
  names(ST)[1:2] = c("ID", "Target")
  ST$Target %<>% as.character()
  ST %<>% dplyr::filter(Target %in% V(G)$name)

  if(correct_by_target){
    ts = unique(ST$Target)
  } else{
    ts = V(G)$name
  }

  d = ST$ID %>% unique()
  message("Starting now. It might take some time, please be patient.\n")
  all_sps = igraph::distances(G, v = ts, to = ts)


  nnodes = nrow(all_sps)
  nodes_ID = ST %>%
    dplyr::group_by(ID) %>%
    dplyr::summarise(n = dplyr::n())

  SAMPLES = list()


  NetSci.Sep$nodes_ID = nodes_ID
  NetSci.Sep$N = N
  NetSci.Sep$nnodes = nnodes
  NetSci.Sep$all_sps = all_sps

  NetSci.Sep$ST = ST
  NetSci.Sep$SAMPLES = SAMPLES
  NetSci.Sep$d = d

  rm(all_sps)

  cl = parallel::makeCluster(Threads)
  parallel::clusterExport(cl, "resample_saa")
  parallel::clusterExport(cl , "saa")
  parallel::clusterExport(cl , "resample")
  parallel::clusterExport(cl , "pvals")

  parallel::clusterExport(cl , "nodes_ID", envir = NetSci.Sep)
  parallel::clusterExport(cl , "N", envir = NetSci.Sep)
  parallel::clusterExport(cl , "nnodes", envir = NetSci.Sep)
  parallel::clusterExport(cl , "all_sps", envir = NetSci.Sep)
  parallel::clusterExport(cl , "ST", envir = NetSci.Sep)
  parallel::clusterExport(cl , "d", envir = NetSci.Sep)
  parallel::clusterExport(cl , "SAMPLES", envir = NetSci.Sep)



  tmporary = parallel::clusterApplyLB(cl, 1:nrow(nodes_ID), resample_saa)

  message("Phew. The first part is done. Not ready yet.\n")
  SAMPLES = list(); saa_stars = list()
  for(diseases_all in 1:length(tmporary)){
    SAMPLES[[diseases_all]] = tmporary[[diseases_all]]$SAMPLES
    saa_stars[[diseases_all]] = tmporary[[diseases_all]]$saa_stars
  }
  saa_stars %<>% dplyr::bind_rows()

  NetSci.Sep$SAMPLES = SAMPLES
  NetSci.Sep$saa_stars = saa_stars

  # parallel::stopCluster(cl)
  #
  # cl = parallel::makeCluster(Threads)
  # parallel::clusterExport(cl, "resample_saa")
  # parallel::clusterExport(cl , "saa")
  # parallel::clusterExport(cl , "resample")
  # parallel::clusterExport(cl , "pvals")
  #
  # parallel::clusterExport(cl , "nodes_ID", envir = NetSci.Sep)
  # parallel::clusterExport(cl , "N", envir = NetSci.Sep)
  # parallel::clusterExport(cl , "nnodes", envir = NetSci.Sep)
  # parallel::clusterExport(cl , "all_sps", envir = NetSci.Sep)
  # parallel::clusterExport(cl , "ST", envir = NetSci.Sep)
  # parallel::clusterExport(cl , "d", envir = NetSci.Sep)
  parallel::clusterExport(cl , "saa_stars", envir = NetSci.Sep)
  parallel::clusterExport(cl , "SAMPLES", envir = NetSci.Sep)

  message("\n Starting the second part...\n")

  Sab_tmp  = parallel::clusterApplyLB(cl, 1:nrow(nodes_ID), sab_aux) %>%
    dplyr::bind_rows()

  Sab_tmp$Saa_Dis = ifelse(is.nan(Sab_tmp$Saa_Dis), Inf, Sab_tmp$Saa_Dis)
  message("Now we are almost there. Hold on :)\n")

  Sab_tmp[is.na(Sab_tmp)] <- Inf
  NetSci.Sep$Sab_tmp = Sab_tmp
  rm(Sab_tmp)

  parallel::clusterExport(cl , "Sab_tmp", envir = NetSci.Sep)


  SAB = parallel::clusterApplyLB(cl,
                                 1:nrow(Sab_tmp),
                                 SAB_complete) %>%
    dplyr::bind_rows()

  parallel::stopCluster(cl)
  return(SAB)
}
