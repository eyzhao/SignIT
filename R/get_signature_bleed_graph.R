#' Graph Representing Signature Bleed Relationships
#'
#' Creates an iGraph object representing bleed between mutation signatures.
#'
#' @param exposure_output       Output from get_exposures
#' @param min_bleed             Minimum bleed threshold for an interaction to be included in the graph
#'
#' @return An iGraph object where vertices are signatures and edges are bleed interactions.
#'
#' @import igraph
#' @import dplyr
#' @import tidyr

get_signature_bleed_graph <- function(exposure_output, min_bleed) {
  correl <- get_exposure_pairwise_correlations(exposure_output)
  
  g <- correl %>%
    unite(sig_pair, signature_1, signature_2, sep='||') %>%
    group_by(sig_pair) %>%
    summarise(bleed = -spearman[1]) %>%
    ungroup() %>%
    separate(sig_pair, c('from', 'to'), sep='\\|\\|') %>%
    graph_from_data_frame
  
  order_vector <- sapply(
    V(g)$name, 
    function(z) {
      which(exposure_output$signature_names == z)
    })
  
  g %>%
    delete.edges(which(E(g)$bleed < min_bleed)) %>%
    permute(order_vector)
}
