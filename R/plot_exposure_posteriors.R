#' Plot the signature exposure posterior distributions
#' 
#' Visualizes the posterior distributions of the MCMC solution to signature exposures.
#' By default, they are visualized as violin plots.
#'
#' @param exposures_mcmc_output     Output from get_exposures.
#' @param view                      Can be either "violin" or "boxplot".
#' @param units                     Units to present exposures in. Can be
#'                                      'mutations', 'megabase', or 'fraction'.
#'
#' @return A ggplot plot of the posterior distributions.
#'
#' @import ggplot2
#' @import dplyr
#'
#' @export

plot_exposure_posteriors <- function(exposures_mcmc_output, view = 'violin', signature_trim='Signature', units='mutations') {
  if (units == 'megabase') {
      unit_label = 'Mutations / Mb'
  } else if (units == 'fraction') {
      unit_label = 'Exposure Fraction'
  } else {
      unit_label = 'Exposure (Mutations)'
  }

  plot <- exposures_mcmc_output$exposure_chain %>%
    mutate(
        signature = trim_signature_names(signature, signature_trim),
        mutation_burden = sum(exposure),
        exposure = case_when(
            units == 'megabase' ~ exposure / 3234.83,
            units == 'fraction' ~ exposure / exposures_mcmc_output$n_mutations,
            TRUE ~ exposure
        )
    ) %>%
    ggplot(aes(
      x = signature %>% as.factor,
      y = exposure
    )) +
    labs(x = 'Signature', y = unit_label) +
    rotate_x_axis_labels()
    
  if (view == 'boxplot') {
    plot <- plot + geom_boxplot()
  } else {
    plot <- plot + geom_violin()
  }

  return(plot)
}


