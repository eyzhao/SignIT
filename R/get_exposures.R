#' Compute the MCMC sample space of exposure vectors
#'
#' At the heart of SignIT is the MCMC sampling of signature exposure solutions.
#' This function provides a convenient wrapper around the MCMC sampling steps.
#' It constructs the STAN model, runs it, and then returns the result in a well formatted list.
#'
#' @param mutation_catalog      Data frame with columns mutation_type, count. The mutation types
#'                              must be equivalent to those in reference_signatures.
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
#' @import parallel
#'
#' @export

get_exposures <- function(
    mutation_catalog,
    reference_signatures = NULL, 
    n_chains = 4, 
    n_iter = 200, 
    n_adapt = 200, 
    n_cores = 1,
    stan_model = NULL,
    quiet = FALSE
) {
    if (get_os() == 'windows' && n_cores > 1) {
        stop("Multicore processing is not available on Windows. Please leave n_cores = 1")
    }

    if (is.null(stan_model)) {
        stan_model <- stanmodels$signature_model
    }

    if (! 'mutation_type' %in% names(mutation_catalog) || ! 'count' %in% names(mutation_catalog)) {
        stop('
             Mutation catalog is not properly formatted. 
             Must have two columns: mutation_type (character) and count (integer)
        ')
    }

    if (is.null(reference_signatures)) {
        reference_signatures <- get_reference_signatures()
    }

    reference_signatures <- normalize_reference_signatures(reference_signatures)

    signature_names <- reference_signatures %>% select(-mutation_type) %>% colnames
    n_mutations <- sum(mutation_catalog$count)

    if (n_mutations > 0) {
        stan_data = list(
            N = sum(mutation_catalog$count),
            S = reference_signatures %>% select(-mutation_type) %>% dim %>% .[2],
            R = reference_signatures %>% dim %>% .[1],
            v = mutation_catalog$count,
            ref_signatures = reference_signatures %>% reference_signatures_as_matrix(mutation_catalog)
        )

        message('Sampling')
        
        progress_base = interactive() && !isatty(stdout()) && !identical(Sys.getenv("RSTUDIO"), "1")

        stan_object <- sampling(
            object = stan_model,
            data = stan_data,
            chains = n_chains,
            iter = n_iter + n_adapt,
            warmup = n_adapt,
            cores = n_cores,
            control = list(
            ),
            open_progress = (! quiet) && progress_base,
            show_messages = ! quiet,
            refresh = if_else(quiet, -1, 10)
        )

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
    } else {
        warning('Your dataset contained zero mutations. Returning zero exposure vector.')
        exposure_chain <- crossing(
            chain = 1:n_chains,
            iteration = 1:n_iter,
            signature = signature_names,
            exposure = 0
        ) %>%
        mutate(
            signature = factor(signature, levels = signature_names)
        )

        stan_object <- NULL
    }
 
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


#' Watanabe-Akaike Information Criterion for Signatures
#'
#' @param mcmc_output    Output from \code{\link{get_exposures}}.
#'
#' @return WAIC value
#'
#' @import loo
#' @export

compute_signatures_waic <- function(mcmc_output) {
    log_lik <- extract_log_lik(mcmc_output$model)
    return(waic(log_lik)$waic)
}


#' Watanabe-Akaike Information Criterion for Signatures
#'
#' @param mcmc_output    Output from \code{\link{get_exposures}}.
#'
#' @return WAIC value
#'
#' @import loo
#' @export

compute_signatures_waic <- function(mcmc_output) {
    log_lik <- extract_log_lik(mcmc_output$model)
    return(waic(log_lik)$waic)
}
