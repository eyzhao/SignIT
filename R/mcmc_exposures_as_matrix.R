#' Convert MCMC sampled exposure values to matrix format
#'
#' MCMC samples of exposure solutions are initially in long format, with one row per
#' signature per iteration. This method converts this data frame to wide format,
#' with one row per iteration and one column per signature. It them converts this
#' to a matrix and returns it.
#'
#' @param exposures_mcmc_output     Output from get_exposures()
#'
#' @return Matrix of MCMC iterations. One row per iteration, one column per signature.
#'
#' @import dplyr
#' @import tidyr

mcmc_exposures_as_matrix <- function(exposures_mcmc_output) {
  exposures_mcmc_output$exposure_chain %>%
    spread(signature, exposure) %>%
    as.matrix()
}


