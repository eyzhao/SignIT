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
  
  signature_model_updated <- run_signature_model(
    signature_model, 
    n_iter = n_iter
  )

  signature_exposures <- signature_model_updated %>% extract_exposure_chain(
    signature_names = signature_names, 
    n_mutations = n_mutations
  )
  
  return(list(
    sampling_tool = 'jags',
    mutation_catalog = mutation_catalog,
    exposure_chain = signature_exposures %>% mutate(signature = factor(signature, levels = signature_names)),
    reference_signatures = reference_signatures,
    signature_names = signature_names,
    n_mutations = n_mutations,
    model = signature_model_updated
  ))
}

get_stan_exposures <- function(mutation_catalog = NULL, file = NULL, reference_signatures = NULL, n_chains = 4, n_iter = 200, n_adapt = 200, n_cores = 1, stan_model = NULL) {
    if (is.null(stan_model)) {
        stan_model <- get_stan_model()
    }

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

    reference_signatures <- get_reference_signatures()

    # Reorder mutation catalog vector so that the mutation types are same order as reference_signatures
    mutation_catalog <- reference_signatures %>%
        left_join(mutation_catalog %>% dplyr::rename(count_ = count), by = 'mutation_type') %>%
        replace_na(list(count_ = 0)) %>%
        select(mutation_type, count = count_)

    stan_data = list(
        N = sum(mutation_catalog$count),
        S = reference_signatures %>% select(-mutation_type) %>% dim %>% .[2],
        R = reference_signatures %>% dim %>% .[1],
        v = mutation_catalog$count,
        ref_signatures = reference_signatures %>% select(-mutation_type) %>% as.matrix
    )

    message('Sampling')

    stan_object <- sampling(
        object = stan_model,
        data = stan_data,
        chains = n_chains,
        iter = n_iter + n_adapt,
        warmup = n_adapt,
        cores = n_cores,
        control = list(
        )
    )

    signature_names <- reference_signatures %>% select(-mutation_type) %>% colnames
    n_mutations <- sum(mutation_catalog$count)

    fit <- stan_object %>% as.array
    exposure_chain <- fit[, , grepl('exposure', dimnames(fit)$parameters)] %>% 
        plyr::adply(2, function(z) { as_tibble(z) }) %>% 
        mutate(
            iteration = row_number(), 
            chain = factor(chains) %>% as.integer
        ) %>% 
        select(-chains) %>% 
        `colnames<-`(c(signature_names, 'iteration', 'chain')) %>% 
        gather(signature, exposure, -iteration, -chain) %>% 
        as_tibble %>%
        mutate(
            signature = factor(signature, levels = signature_names),
            exposure = exposure * n_mutations
        )
 
  return(list(
    sampling_tool = 'stan',
    mutation_catalog = mutation_catalog,
    exposure_chain = exposure_chain,
    reference_signatures = reference_signatures,
    signature_names = signature_names,
    n_mutations = n_mutations,
    model = stan_object
  ))
}

get_stan_model <- function() {

    stan_model_string_vectorized <- '
        data {
          int<lower=1> S; // number of signatures
          int<lower=1> N; // number of mutations
          int<lower=1> R; // number of mutation types (vocabulary size)
          int<lower=0> v[R]; // mutation catalog vector
          matrix[R, S] ref_signatures; // reference mutation signatures (theta)
        }
        parameters {
          simplex[S] exposures; // mixing proportions
        }
        transformed parameters {
          vector[R] sim_catalog;
          simplex[R] sim_catalog_prob;
          sim_catalog = ref_signatures * exposures * N;
          sim_catalog_prob = sim_catalog / sum(sim_catalog);
        }
        model {
          v ~ multinomial(sim_catalog_prob);
          exposures ~ dirichlet(rep_vector(1, S));
        }
    '

    message('Establishing Stan model')

    stan_dso = stan_model(model_code = stan_model_string_vectorized)

}
