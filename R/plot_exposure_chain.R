#' Plot the MCMC chains for the signature exposure solutions
#'
#' This function allows one to visualize the solution chains for the exposure
#' of each mutation signature. It plots each chain as a line graph, allowing
#' one to inspect for convergence.
#'
#' @param exposures_mcmc_output Output from running plot_exposures
#'
#' @return A ggplot object
#'
#' @import ggplot
#' @import magrittr
#'
#' @export

plot_exposure_chain <- function(exposures_mcmc_output) {
  exposures_mcmc_output$exposure_chain %>%
    mutate(
      chain = factor(chain),
      signature = trim_signature_names(signature)
    ) %>%
    ggplot(aes(
      y = exposure,
      x = iteration,
      group = chain,
      colour = chain
    )) +
    facet_grid(signature ~ .) +
    geom_line() +
    scale_colour_brewer(palette = 'Set1')
}
