select_n_populations <- function(mutation_table, parallel = TRUE) {
    message('Testing clonality models for 1-5 populations.')

    if (parallel) {
        registerDoParallel(5)
    }

    models <- plyr::dlply(
        tibble(n_populations = 1:5),
        'n_populations',
        function(n_populations) {
            population_mcmc_output <- get_populations(
                mutation_table,
                as.integer(n_populations),
                method = 'MAP'
            )

            BIC = population_mcmc_output$BIC

            return(list(
                mcmc_output = population_mcmc_output,
                BIC = BIC
            ))
    }, .parallel = TRUE)

    BIC_vector <- sapply(models, function(z) {
        z$BIC
    })

    print(BIC_vector)

    optimal_n_populations = which(BIC_vector == min(BIC_vector))

    return(list(
        optimal_n_populations = optimal_n_populations,
        BIC_vector = BIC_vector,
        evaluated_models = models
    ))
}
