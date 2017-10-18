#' Compute the MCMC sample space of exposure vectors
#'
#' At the heart of SignIT is the MCMC sampling of signature exposure solutions.
#' This function provides a convenient wrapper around the MCMC sampling steps.
#' It constructs the JAGS model, runs it, and then returns the result in a well formatted list.
#'
#' @param mutation_catalog      Data frame with columns mutation_type, count. The mutation types
#'                              must be equivalent to those in reference_signatures.
#'
#' @param file                  The mutation catalog can be stored as a TSV file, and its path can
#'                              can be provided in this parameter. Either mutation_catalog or file MUST
#'                              be provided.
#'
#' @param reference_signatures  Reference mutation signatures, such as that output from get_reference_signatures.
#'
#' @param n_chains              Number of MCMC chains
#'
#' @param n_iter                Number of iterations per chain. Total iterations will be n_chains times n_iter.
#'
#' @param n_adapt               Number of burn-in iterations
#'
#' @return List containing the MCMC samples, as well as other data such as reference signatures and mutation catalog.
#'
#' @import dplyr
#' @import readr
#'
#' @export

get_exposures <- function(mutation_catalog = NULL, file = NULL, reference_signatures = NULL, n_chains = 4, n_iter = 10000, n_adapt = 10000) {
  if (is.null(mutation_catalog) && is_null(file)) {
    stop('Must either provide mutation_catalog as data frame or path to TSV file')
  } else if (is.null(mutation_catalog)) {
    mutation_catalog <- read_tsv(
      file, 
      col_types = cols(
        mutation_type = col_character(),
        count = col_number()
      )
    ) %>%
    mutate(
        count = as.integer(round(count))
    )
  }
  
  if (! 'mutation_type' %in% names(mutation_catalog) || ! 'count' %in% names(mutation_catalog)) {
    stop('Mutation catalog is not properly formatted. Must have two columns: mutation_type (character) and count (integer)')
  }
  
  if (is.null(reference_signatures)) {
    reference_signatures <- get_reference_signatures()
  }
  
  signature_model <- get_signature_model(
    mutation_catalog, 
    reference_signatures = reference_signatures,
    n_chains = n_chains,
    n_adapt = n_adapt
  )
  
  signature_names <- get_signature_names(reference_signatures)
  n_mutations <- mutation_catalog[['count']] %>% sum
  
  signature_exposures <- run_signature_model(
    signature_model, 
    signature_names, 
    n_mutations = n_mutations,
    n_iter = n_iter
  )
  
  return(list(
    mutation_catalog = mutation_catalog,
    exposure_chain = signature_exposures %>% mutate(signature = factor(signature, levels = signature_names)),
    reference_signatures = reference_signatures,
    signature_names = signature_names,
    n_mutations = n_mutations
  ))
}


