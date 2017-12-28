select_n_populations <- function(mutation_table, parallel = TRUE) {
    message('Testing clonality models for 1-5 populations.')

    if (parallel) {
        n_cores <- min(detectCores(), 5)
        registerDoParallel(n_cores)
    }

    models <- foreach(
      i = 1:5,
      .combine = c,
      .packages = c('tidyverse', 'signit')
    ) %dopar% {
            population_mcmc_output <- get_populations(
                mutation_table,
                n_populations = i,
                method = 'MAP'
            )

            BIC = population_mcmc_output$BIC

            return(list(list(
                mcmc_output = population_mcmc_output,
                BIC = BIC
            )))
    }

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
