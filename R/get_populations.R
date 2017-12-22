get_vaf_correction <- function(mutation_table) {
    with(
        mutation_table,
        tumour_content / ((tumour_content * tumour_copy) + ((1 - tumour_content) * normal_copy)),
    )
}

get_populations <- function(
    mutation_table,
    n_populations,
    n_chains = 4, 
    n_iter = 300, 
    n_adapt = 200, 
    n_cores = n_chains,
    method = 'mcmc'
) {
    mutation_table <- mutation_table %>% distinct()
    n_mutations <- dim(mutation_table)[1]

    stopifnot(
        mutation_table$tumour_content %>% unique %>% length == 1,
        !missing(n_populations)
    )

    tumour_content <- mutation_table$tumour_content[1]

    message('Establishing Stan model')

    stan_dso = get_stan_model('population')

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

        parsed_output <- parse_population_model_mcmc(stan_fit)
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

parse_population_model_mcmc <- function(population_mcmc_out) {
    mcmc_chains <- population_mcmc_out %>% 
        as.array %>%
        plyr::adply(2, function(z) { as_tibble(z) }) %>%
        mutate(
            iteration = row_number(),
            chain = factor(chains) %>% as.integer
        ) %>%
        select(-chains) %>% 
        as_tibble() %>%
        gather(parameter, value, -chain, -iteration) %>%
        mutate(
            parameter_name = gsub('(.*?)\\[.*?\\].*', '\\1', parameter),
            parameter_index = gsub('.*?\\[(.*?)\\].*', '\\1', parameter) %>% as.numeric
        ) %>%
        select(-parameter)

    return(mcmc_chains)
}

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
        likelihoods = rmutil::dbetabinom(successes, trials, correction * z$mu, kappa)
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


