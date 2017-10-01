#' Computes the residual sum of squares MCMC solutions
#'
#' Iterates through every MCMC exposure solution and computes the
#' residual sum of squares between that solution and the mutation catalog.
#'
#' @param exposures_mcmc_output     Output from get_exposures
#'
#' @return A vector of RSS values, one per MCMC iteration
#'
#' @import dplyr
#' @import tibble

get_rss <- function(exposures_mcmc_output) {
  catalog <- exposures_mcmc_output$mutation_catalog %>% .$count
  exposures <- exposures_mcmc_output %>% mcmc_exposures_as_matrix() %>% as_tibble %>% .[, exposures_mcmc_output$signature_names]
  reference <- exposures_mcmc_output$reference_signatures %>% select(-mutation_type) %>% as.matrix
  
  apply(exposures, 1, function(e) {
    error_rss <- (catalog - (reference %*% e))^2 %>% sum
  })
}

