#' Get the JAGS Model for Mutation Signature Analysis
#'
#' Creates the JAGS model for MCMC analysis of mutation signatures.
#' An NNLS solution is used as the initial guess.
#'
#' @param mutation_catalog     Data frame with columns mutation_type and count.
#' 
#' @param reference_signatures Reference signature data with mutation_type column and
#'                             one named column per signature. Can be obtained from
#'                             get_reference_signatures().
#' 
#' @param n_chains             The number of chains to run MCMC analysis with.
#'
#' @param n_adapt              The number of burn-in iterations to run per chain.
#'
#' @return JAGS model object with initialized values, ready to sample from
#'
#' @import rjags
#' @import nnls
#' @import dplyr

get_signature_model <- function(mutation_catalog, reference_signatures = NULL, n_chains = 4, n_adapt = 10000) {
  if (! mutation_types_match(reference_signatures, mutation_catalog, order_is_important = FALSE)) {
    stop('mutation_type does not match reference mutation types')
  }
  
  reference_signature_matrix <-  reference_signatures_as_matrix(reference_signatures, mutation_catalog)
  
  print(reference_signature_matrix)
  print(mutation_catalog[['count']])
  naive_nnls_solution <- nnls(reference_signature_matrix, mutation_catalog[['count']])$x
  naive_nnls_solution <- if_else(naive_nnls_solution == 0, 1, naive_nnls_solution) %>% ceiling()
  print(naive_nnls_solution)

  n_mutations <- sum(mutation_catalog[['count']])
  
  signature_data = list(
    mutation_catalog = mutation_catalog[['count']],
    exposure_dirichlet_prior = rep(1, length(naive_nnls_solution)),
    mutation_signatures = reference_signature_matrix,
    n_mutations = n_mutations
  )
  signature_init = list(
    exposures = naive_nnls_solution / sum(naive_nnls_solution)
  )
  
  signature_model = '
    model {
      mutation_catalog ~ dmulti((mutation_signatures %*% exposures) / n_mutations, n_mutations);
      exposures ~ ddirch(exposure_dirichlet_prior)
    }
  '
  
  jags.model(
    file=textConnection(signature_model), 
    data=signature_data,
    inits=signature_init,
    n.chains=n_chains, 
    n.adapt=n_adapt
  )
}


