#' Plot the signature exposure posterior distributions
#' 
#' Visualizes the posterior distributions of the MCMC solution to signature exposures.
#' By default, they are visualized as violin plots.
#'
#' @param exposures_mcmc_output     Output from get_exposures.
#' @param view                      Can be either "violin" or "boxplot".
#'
#' @return A ggplot plot of the posterior distributions.
#'
#' @import ggplot2
#' @import dplyr
#'
#' @export

plot_exposure_posteriors <- function(exposures_mcmc_output, view = 'violin', signature_trim='Signature') {
  plot <- exposures_mcmc_output$exposure_chain %>%
    mutate(
        signature = trim_signature_names(signature, signature_trim)
    ) %>%
    ggplot(aes(
      x = signature %>% as.factor,
      y = exposure
    )) +
    labs(x = 'Signature', y = 'Exposure') +
    rotate_x_axis_labels()
    
  if (view == 'boxplot') {
    plot <- plot + geom_boxplot()
  } else {
    plot <- plot + geom_violin()
  }

  return(plot)
}


