#' @title Jaccard
#' @description Calculates the Jaccard index between different sets.
#'
#' @param Data A data.frame with 2 columns. The first refers to the set and the second the elements
#'
#' @return a data.frame with the set names and their Jaccard index
#' @importFrom magrittr `%>%` `%<>%`
#' @importFrom  igraph graph_from_data_frame bipartite_mapping V E as_incidence_matrix
#' @importFrom  Rfast transpose
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
#'
#' @examples
#' set.seed(123)
#' Data = data.frame(Class = sample(c("X", "Y", "Z"), replace = TRUE, size = 50),
#'                   Element = sample(LETTERS[1:15], replace = TRUE, size = 50))
#' Data = unique(Data)
#' Jaccard(Data)


Jaccard = function(Data){
  `%<>%` <- magrittr::`%<>%`
  `%>%` <- magrittr::`%>%`

  Data %<>% unique()

  if(nrow(Data) == 0){
    stop("Data must have more than 0 rows.")
  }
  if(is.data.frame(Data) == F){
    stop("Data must be a data.frame.")
  }

  levels = Data[1,] %>% unique() %>% length()
  if(levels <= 1){
    stop("You need more than one type.")
  }

  g =  Data %>%
    igraph::graph_from_data_frame(, directed = F)

  igraph::V(g)$type <- igraph::bipartite_mapping(g)$type
  #get the incidence matrix (or adjacency, depending on how the data was structured)
  A = igraph::as_incidence_matrix(g) %>% as.matrix()
  gg =   (A)  %*% Rfast::transpose(A)

  names(gg) = colnames(gg) = rownames(A)
  NORM = matrix(NA, ncol = ncol(gg), nrow = nrow(gg))
  #Normalize the values
  ADJ_for_DIS2DIS = gg
  pb <- utils::txtProgressBar(min = 0, max = (ncol(NORM)), style = 3)
  for( i in 1:ncol(NORM)){
    utils::setTxtProgressBar(pb, i)
    for(j in i:(nrow(NORM))){
      NORM[i,j] = NORM[j,i] = ADJ_for_DIS2DIS[i,j]/(ADJ_for_DIS2DIS[i,i]+ADJ_for_DIS2DIS[j,j]-ADJ_for_DIS2DIS[i,j])
    }
  }
  close(pb)

  # Genes = diag(gg) %>% as.data.frame()
  # Genes$ID = row.names(Genes)
  # Genes$prop = (Genes$./ sum(Genes$.) )%>% CoDiNA::normalize()
  # Genes$Count = Genes$.
  # Genes = Genes[,-1]
  #
  # Transform into a edge list
  rownames(NORM) = colnames(gg)
  colnames(NORM) = colnames(gg)
  G = NORM %>% wTO::wTO.in.line()
  names(G)[3]="Jaccard.Index"

  return(G)
}
