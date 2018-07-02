#' Plot the signature exposure posterior distribution along with a visualization of
#' signature bleed interactions.
#' 
#' @param exposures_mcmc_output     Output from get_exposures.
#' @param view                      Can be either "violin" or "boxplot" (default: "violin")
#' @param min_bleed                 Bleed threshold, between 0 and 1. Default: 0.2.
#' @param signature_trim            Text to be trimmed off of the x-axis signature
#'                                      labels.
#' @param legend                    Boolean - whether or not to include the legend.
#' @param plot_heights              Vector of 2 numbers denoting relative heights
#'                                      of exposure plot and bleed graph.
#' @param units                     Units to present exposures in. Can be
#'                                      'mutations', 'megabase', or 'fraction'.
#'
#' @return A cowplot merge of two ggplots showing the posterior distributions and signature bleed.
#'
#' @import ggplot2
#' @import dplyr
#' @import cowplot
#' @import ggraph
#'
#' @export


plot_exposure_posteriors_bleed <- function(
    exposures_mcmc_output,
    view='violin',
    min_bleed = 0.2,
    signature_trim='Signature',
    legend = TRUE,
    plot_heights = c(1,1),
    units = 'mutations'
) {
    bleed_graph <- get_signature_bleed_graph(exposures_mcmc_output, min_bleed)

    bleed_plot <- bleed_graph %>%
        ggraph(layout = 'linear') +
        geom_edge_arc(aes(
          colour = bleed,
          edge_width = bleed
        )) +
        geom_node_text(aes(
          label = '',
          angle = 90
        ), hjust = 1) +
        theme_void() +
        labs(edge_colour = 'Signature\nBleed') +
        scale_edge_color_distiller(palette = 'Spectral') +
        guides(edge_width = FALSE) +
        scale_y_continuous(limits = c(-15, 0), expand = c(-0.005, 0)) +
        scale_x_discrete(labels = exposures_mcmc_output$signature_names)

    posterior_plot <- plot_exposure_posteriors(
        exposures_mcmc_output,
        view=view,
        signature_trim=signature_trim,
        units = units
    ) +
        theme(
            axis.title.x = element_blank()
        )

    final_plot <- plot_grid(
        posterior_plot,
        bleed_plot + theme(legend.position = 'none'),
        ncol = 1,
        align = 'v',
        rel_heights = plot_heights
    )

    if (legend) {
        final_plot <- final_plot %>%
            plot_grid(
                get_legend(bleed_plot),
                rel_widths = c(5,1)
            )  
    }

    return(final_plot)
}
