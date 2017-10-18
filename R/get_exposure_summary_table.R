#' Creates a summary table of exposure statistics
#'
#' After running get_exposures.R, one receives an object containing the points sampled
#' by the MCMC chains. It is often necessary to review summary statistics on these
#' sample chains to examine the exposure ranges of individual signatures.
#' This function provides an interface for computing summary statistics including
#' mean, median, range, and "confidence bands".
#'
#' By default, the function reports mean and median, as well as range, 95% confidence bands,
#' and quartiles.
#'
#' @param exposures_mcmc_output     Output from get_exposures()
#'
#' @param alpha                     Vector defining confidence bands to report. Alpha is between 0
#'                                  and 1. The two confidence limits reported are
#'                                  at rank alpha/2 and 1 - (alpha/2). For example, 
#'                                  setting alpha = 0.05 results in 95% an interval capturing
#'                                  the middle 95% of the data. Alpha = 0.5 results in
#'                                  reporting of the 1st and 3rd quartiles. Multiple
#'                                  alpha values can be provided as vector elements. These
#'                                  values will each be computed and reported in separate columns.
#'
#' @return A dataframe with one row per signature and summary stats in the columns
#'
#' @importFrom plyr ddply
#' @import dplyr
#' @import tidyr
#'
#' @export

get_exposure_summary_table <- function(exposures_mcmc_output, alpha = c(0, 0.05, 0.5), fraction = FALSE) {
    quantile_breaks = c((alpha/2), 1-(alpha/2)) %>% sort

    if (fraction) {
        exposures_mcmc_output$exposure_chain$exposure <- exposures_mcmc_output$exposure_chain$exposure / exposures_mcmc_output$n_mutations
    }

    centre_table <- exposures_mcmc_output$exposure_chain %>%
        group_by(signature) %>%
        summarise(
            mean_exposure = mean(exposure),
            median_exposure = median(exposure)
        )

    quantile_table <- exposures_mcmc_output$exposure_chain %>%
        plyr::ddply(
            'signature', 
            function(z) { 
                q = quantile(z$exposure, quantile_breaks)
                data.frame(quantile = names(q), value = q) 
            }) %>%
        spread(quantile, value)

    return(centre_table %>% inner_join(quantile_table, by = 'signature'))
}
