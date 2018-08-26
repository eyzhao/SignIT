#' SignIT population model without mutation signatures
#'
#' Inference of mutational subpopulations without joint inference of signatures
#'
#' This method can be used to obtain full posterior estimates of mutational subpopulations.
#' It is also used under the hood by \code{\link{select_n_populations}} to perform automatic
#' estimation of the number of populations for model selection.
#'
#' @param mutation_table        Table of mutations, same as input to 
#'                              \code{\link{get_population_signatures}}.
#'
#' @param n_populations         Number of populations to fit.
#'
#' @param method                The method to use. Can be either "mcmc" for Hamiltonial 
#'                              Markov Chain Monte Carlo or "MAP" for maximum a posteriori estimate.
#'
#' @param n_chains              Number of chains to use (only for \code{method == "mcmc"})
#' @param n_iter                Number of sampling iterations per chain (only for \code{method == "mcmc"})
#' @param n_adapt               Number of adaptation iterations per chain (only for \code{method == "mcmc"})
#' @param n_cores               Number of cores for parallel sampling (by default equal to \code{n_chains}.
#'                              only for \code{method == "mcmc"})
#'
#' @return List object with Stan output, parsed data, and relevant metadata
#'
#' @import tidyr
#' @import dplyr
#' @import tibble
#' @import rstan
#'
#' @export

get_populations <- function(
    mutation_table,
    n_populations,
    method = 'mcmc',
    n_chains = 4, 
    n_iter = 300, 
    n_adapt = 200, 
    n_cores = n_chains
) {
    if (get_os() == 'windows' && n_cores > 1) {
        stop("Multicore processing is not available on Windows. Please leave n_cores = 1")
    }

    mutation_table <- mutation_table %>% distinct()
    n_mutations <- dim(mutation_table)[1]

    stopifnot(
        mutation_table$tumour_content %>% unique %>% length == 1,
        !missing(n_populations)
    )

    tumour_content <- mutation_table$tumour_content[1]

    message('Establishing Stan model')

    stan_dso = stanmodels$population_model

    stan_data = list(
        N = n_mutations,
        L = n_populations,
        x = as.numeric(mutation_table$alt_depth),
        purity = as.numeric(tumour_content),
        Ct = as.numeric(mutation_table$tumour_copy),
        Cn = as.numeric(mutation_table$normal_copy),
        d = as.numeric(mutation_table$total_depth)
    )

    message('Sampling')

    if (method == 'mcmc' || method == 'vb') {

        if (method == 'mcmc') {
            stan_fit = sampling(
                object = stan_dso,
                data = stan_data,
                chains = n_chains,
                iter = n_iter + n_adapt,
                warmup = n_adapt,
                cores = n_cores,
                init = function() {
                    list(mu = seq(0.1, 0.9, (0.9-0.1)/(stan_data$L-1)))
                }
            )
        } else {
            stan_fit = vb(
                object = stan_dso,
                data = stan_data,
                init = function() {
                    list(mu = seq(0.1, 0.9, (0.9-0.1)/(stan_data$L-1)))
                }
            )
        }

        parsed_output <- parse_stan_output(stan_fit)
        parameter_estimates <- get_population_parameter_estimates(parsed_output)

        population_mcmc_output <- list(
            mcmc_raw_output = stan_fit,
            chains = parsed_output,
            parameter_estimates = parameter_estimates
        )

    } else if (method == 'MAP') {
        estimates = optimizing(
            object = stan_dso,
            data = stan_data,
            init = function() {
                list(mu = seq(0.1, 0.9, (0.9-0.1)/(stan_data$L-1)))
            }
        )
        
        parameter_estimates <- estimates$par %>% 
            as.data.frame %>% 
            `colnames<-`('value') %>%
            rownames_to_column('parameter') %>% 
            mutate(
                parameter_name = gsub('(.*?)\\[.*', '\\1', parameter), 
                parameter_index = gsub('.*?\\[(.*?)\\]', '\\1', parameter) %>% as.numeric
            ) %>%
            arrange(parameter_name, parameter_index) %>%
            plyr::dlply('parameter_name', function(z) {z$value})

        parameter_estimates[['mu_ordered']] <- NULL
        parameter_estimates[['kappa']] <- parameter_estimates[['kappa_minus_two']] + 2
        parameter_estimates[['kappa_minus_two']] <- NULL

        population_mcmc_output <- list(
            parameter_estimates = parameter_estimates
        )
    }

    population_mcmc_output['log_likelihood'] = compute_population_model_log_likelihood(
        mutation_table, 
        population_mcmc_output
    )

    population_mcmc_output['BIC'] = compute_population_bic(mutation_table, population_mcmc_output)

    return(population_mcmc_output)
}

#' Point estimates for population Stan model output
#'
#' Parses Stan output and retrieves the mean for each parameter
#'
#' @param population_mcmc_parsed        Output from \code{\link{get_populations}} after
#'                                      processed by \code{\link{parse_stan_output}}.
#'
#' @return List of parameter estimates.
#'
#' @import dplyr
#' @import tidyr

get_population_parameter_estimates <- function(population_mcmc_parsed) {
    estimates <- population_mcmc_parsed %>%
        group_by(parameter_name, parameter_index) %>%
        summarise(
            mean = mean(value)
        ) %>%
        ungroup() %>%
        plyr::dlply(
            'parameter_name',
            function(df) {
                df %>% arrange(parameter_index) %>% .$mean
            }
        )

    estimates$kappa = estimates$kappa_minus_two + 2
    estimates[['mu_ordered']] <- NULL
    estimates[['kappa_minus_two']] <- NULL

    return(estimates)
}

#' Bayesian Information Criterion for Population Models
#'
#' @param mutation_table            Table of mutation data, same as input for 
#'                                  \code{\link{get_population_signatures}}.
#'
#' @param population_mcmc_output    Output from \code{\link{get_populations}}.
#'
#' @return BIC value

compute_population_bic <- function(mutation_table, population_mcmc_output) {
    population_parameter_estimates <- population_mcmc_output$parameter_estimates

    log_lik = compute_population_model_log_likelihood(
            mutation_table,
            population_mcmc_output
    )

    n = nrow(mutation_table)
    k = length(population_parameter_estimates$mu)

    return(log(n) * k - ( 2 * log_lik ))
}

#' Log Likelihood of Population Model
#'
#' Reports the log likelihood of data based on population model parameter estimates
#'
#' @param mutation_table            Table of mutation data, same as input for 
#'                                  \code{\link{get_population_signatures}}.
#'
#' @param population_mcmc_output    Output from \code{\link{get_populations}}.
#'
#' @return Log likelihood value associated with data, model, and parameter estimates
#'
#' @import tidyr
#' @import dplyr
#' @importFrom rmutil dbetabinom

compute_population_model_log_likelihood <- function(mutation_table, population_mcmc_output) {
    parameter_estimates <- population_mcmc_output$parameter_estimates

    successes = mutation_table$alt_depth
    trials = mutation_table$total_depth
    correction = get_vaf_correction(mutation_table)
    clone_prop = parameter_estimates$clone_prop
    mu = parameter_estimates$mu
    kappa = parameter_estimates$kappa

    tibble(
        idx = 1:length(mu), 
        mu,
        clone_prop
    ) %>% plyr::ddply(c('idx', 'clone_prop'), function(z) {
        likelihoods = dbetabinom(successes, trials, correction * z$mu, kappa)
        return(data.frame(
            mutation_id = 1:length(successes),
            likelihoods
        ))
    }) %>%
        group_by(mutation_id) %>%
        summarise(
            likelihood = log(sum(likelihoods * clone_prop))
        ) %>%
        ungroup() %>%
        .$likelihood %>%
        sum
}
