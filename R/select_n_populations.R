select_n_populations <- function(mutation_table) {
    message('Testing clonality models for 1-5 populations. This may take a while...')
    min_bic = Inf
    optimal_n_populations <- 0
    models <- list()

    for (n_populations in 1:5) {
        message(sprintf('Testing %s-population model', n_populations))

        population_mcmc_output <- get_populations(
            mutation_table,
            n_populations,
            method = 'MAP'
        )

        bic = population_mcmc_output$BIC

        message(sprintf('BIC for %s-population model was %s', n_populations, bic))

        if (bic < min_bic) {
            min_bic = bic
            optimal_n_populations <- n_populations
            models[n_populations] <- population_mcmc_output
        } else {
            break
        }
    }

    return(list(
        optimal_n_populations = optimal_n_populations,
        evaluated_models = models
    ))
}
