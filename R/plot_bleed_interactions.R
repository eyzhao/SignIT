#' Creates a circular plot with links between signatures that bleed into each other.
#'
#' @param exposures_mcmc_output           Output from get_exposures.
#' @param min_bleed                       Minimum bleed threshold, below which edges are removed.
#'
#' @import ggraph
#' @import ggplot2
#'
#' @return A ggplot object with bleed interactions.
#'
#' @export

plot_bleed_interactions <- function(exposures_mcmc_output, min_bleed=0.2) {
  bleed_graph <- get_signature_bleed_graph(exposures_mcmc_output, min_bleed)

  p <- bleed_graph %>%
    ggraph(layout = 'linear', circular=TRUE) +
    geom_edge_arc(aes(
      colour = bleed,
      edge_width = bleed
    )) +
    geom_node_text(aes(
      label = name,
      angle = -((-node_angle(x, y) + 90) %% 180) + 90,
      x = x * 1.08,
      y = y * 1.08
    ), hjust = 'outward') +
    scale_x_continuous(expand = c(0.4, 0.4)) +
    scale_y_continuous(expand = c(0.4, 0.4)) +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank()
    ) +
    theme_void() +
    labs(edge_colour = 'Signature\nBleed') +
    scale_edge_color_distiller(palette = 'Spectral') +
    guides(edge_width = FALSE)
  
  return(p)
}
