#' Convenience function which merges MCMC exposure iterations with density clusters
#'
#' @param exposures_mcmc_output     Output from get_exposures.
#' @param exposure_mcmc_clusters    Output from density_clustering.
#' @return A tibble of MCMC iterations in wide format (one column per signature) and a column indicating assigned cluster.
#'
#' @import dplyr
#' @import tibble

merge_exposures_with_clustering <- function(exposures_mcmc_output, exposure_mcmc_clusters) {
  exposures_mcmc_output %>% 
    mcmc_exposures_as_matrix %>%
    as_tibble %>%
    mutate(
      cluster = exposure_mcmc_clusters$cluster
    )
}
