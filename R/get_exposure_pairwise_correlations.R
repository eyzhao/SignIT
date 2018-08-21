#' Pairwise correlation coefficients between signatures
#'
#' As a measure of signature bleed, this method computes the Spearman correlation coefficient
#' between each pair of mutation signatures across the MCMC sample space.
#'
#' @param exposures_mcmc_output     Output from get_exposures()
#'
#' @return A dataframe of pairwise correlation values. One row per signature pair.
#'
#' @import dplyr
#' @import tibble
#'
#' @export

get_exposure_pairwise_correlations <- function(exposures_mcmc_output) {
  exposures_mcmc_output$exposure_chain %>%
    rename(signature_1 = signature) %>%
    plyr::ddply('signature_1', function(signature_1_chain) {
      exposures_mcmc_output$exposure_chain %>%
        rename(signature_2 = signature) %>%
        plyr::ddply('signature_2', function(signature_2_chain) {
          exp_table = tibble(sig_a = signature_1_chain$exposure, sig_b = signature_2_chain$exposure)
          tibble(
            spearman = with(exp_table, cor(sig_a, sig_b, method = 'spearman')),
            slope = lm(sig_b ~ sig_a, data = exp_table)$coefficients['sig_a']
          )
        })
    })
}


