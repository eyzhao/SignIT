#' Run the JAGS model to decipher mutation signatures by MCMC
#'
#' @param signature_model       The JAGS model obtained from get_signature_model
#'
#' @param signature_names       Signature names provided as a character vector. Must match
#'                              the number of signatures provided to the model.
#'
#' @param n_mutations           The number of somatic mutations. This is multiplied by the
#'                              fractional exposure vector to obtain the absolute number of mutations
#'                              associated with each mutation type.
#' 
#' @param n_iter                The number of iterations to run per chain.
#'
#' @return Data frame with one row per signature per iteration of the chain.
#'
#' @importFrom plyr ldply
#' @import rjags
#' @import coda
#' @import dplyr
#' @import tibble
#' @import tidyr

run_signature_model <- function(signature_model, n_iter=10000) {
  update(
    signature_model, 
    n.iter=n_iter
  )

  signature_samples <- signature_model %>% coda.samples(
    variable.names=c("exposures"),
    n.iter=n_iter
  )

  return(signature_samples)
}

extract_exposure_chain <- function(signature_model, signature_names, n_mutations) {
  signature_samples <- signature_model %>%
    lapply(function(chain) {
      chain %>%
        as_tibble %>%
        `colnames<-`(signature_names) %>%
        mutate(
          iteration = row_number()
        )
    })
  
  signature_full_chain <- signature_samples %>% 
    `names<-`(1:4) %>% 
    plyr::ldply(.id = 'chain') %>%
    as_tibble() %>% 
    mutate(chain = as.numeric(chain)) %>%
    gather(signature, exposure, -chain, -iteration) %>%
    mutate(
      signature = signature,
      exposure = exposure * n_mutations
    )
  
  return(signature_full_chain)
}

