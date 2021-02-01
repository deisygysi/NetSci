#' avr_proximity_multiple_target_sets
#'
#' @param set Name of the sets you have targets for. (In a drug-target setup, those would be the drugs of interest).
#' @param G The original graph (often an interactome).
#' @param ST Set-Target data. It is a data.frame with two columns. ID and Target.
#' @param source The source nodes (disease genes).
#' @param N Number of randomizations.
#' @param bins the number os bins for the degree preserving randomization.
#' @param min_per_bin the minimum size of each bin.
#' @importFrom dplyr pull filter bind_rows
#' @description Calculates the average proximity from a set of targets to a set of source nodes.
#' It is calculate using a degree preserving randomization. It is calculated as described in
#' Guney, E. et al (2016) <doi.org:10.1038/ncomms10331>
#' @return proximity and its significance based on the degree preserving randomization.
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
#' S = c("N", "A", "F", "I")
#' T1 = data.frame(ID = "T1", Target = c("H", "M"))
#' T2 = data.frame(ID = "T2", Target = c("G", "O"))
#'
#' avr_proximity_multiple_target_sets(set = c('T1', 'T2'),
#' G = g,
#'  source = S,
#'  ST = rbind(T1,T2),
#'  bins = 5,
#'  min_per_bin = 2)
#'
#'
avr_proximity_multiple_target_sets = function(set, G, ST, source,
                                              N = 1000,
                                              bins = 100,
                                              min_per_bin  = 20){
  Out_NR = list()
  pb <- txtProgressBar(min = 0, max = (length(set)), style = 3)
  for (D in 1:length(set)){
    setTxtProgressBar(pb, D)
    DT = ST %>% dplyr::filter(ID == set[D]) %>% dplyr::pull(Target)

    s = source [source %in% V(G)$name] %>% unique
    t = DT [DT %in% V(G)$name]%>% unique

    if(length(s)> 0 & length(t)> 0){
      d = suppressMessages(proximity_average(G, targets = t, source = s))
      # targets are the DT -> They are resampled
      p = suppressMessages(boot_distance_p(G, targets = t, source = s,
                                           bins = bins,
                                           min_bin  = min_per_bin,
                                           n = N,
                                           d = d,
                                           alpha = 0.05,
                                           proximity  = "average"))
      Out_NR[[D]] = data.frame(Drug = set[D],
                               targets = length(DT),
                               targets_G = length(t),

                               proximity = d,
                               IC_97.5 = p$IC[1],
                               IC_2.5 = p$IC[2],
                               p_gt = p$p_gt,
                               p_lt = p$p_lt,
                               Z = p$Z, row.names =set[D] )
    }
  }
  close(pb)
  Out_NR = dplyr::bind_rows(Out_NR)
  return(Out_NR)
}
