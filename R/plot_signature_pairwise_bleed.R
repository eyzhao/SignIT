#' Plots pairwise correlations between signatures in MCMC samples
#'
#' This plot allows the user to explore their MCMC solutions for the presence
#' of mutation signature bleed effects.
#'
#' The MCMC sample space is likely to contain some degree of multicollinearity.
#' That is to say that some signatures are likely to correlate or anti-correlate
#' with others. This can result from a phenomenon known as "Signature Bleed,"
#' wherein signatures may feed into one anothers' signals if they share similar
#' mutational spectra.
#'
#' This function computes correlation metrics between each pair of signatures.
#' It then plots them as a grid. Signature bleed is associated with anticorrelation
#' between a signature pair, which suggests that good solutions can substitute
#' the contribution of one signature for the other.
#'
#' @param exposures_mcmc_output  Output from get_exposures.
#'
#' @return A ggplot object showing a heatmap of pairwise correlation coefficients.
#'
#' @import dplyr
#' @import ggplot
#'
#' @export

plot_signature_pairwise_bleed <- function(exposures_mcmc_output) {
  get_exposure_pairwise_correlations(exposures_mcmc_output) %>%
    mutate(
      signature_1 = signature_1 %>% `levels<-`(gsub('Signature ', '', levels(signature_1))),
      signature_2 = signature_2 %>% `levels<-`(gsub('Signature ', '', levels(signature_2)))
    ) %>%
    ggplot(aes(x = signature_1, y = signature_2, fill = spearman)) +
    geom_tile(width = 0.8, height = 0.8) +
    scale_fill_distiller(palette = 'Spectral') +
    rotate_x_axis_labels() +
    labs(x = 'First Signature', y = 'Second Signature', fill = 'Spearman Rho\n')
}
