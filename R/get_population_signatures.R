subset_reference_signatures <- function(catalog, reference_signatures) {
    n_mutations <- sum(catalog$count)

    exposures <- catalog %>% 
        get_exposures(quiet = TRUE) %>% 
        get_exposure_summary_table %>% 
        mutate(signature_present = `2.5%` > n_mutations / 2000)

    involved_signatures <- exposures %>%
        filter(signature_present) %>%
        .$signature %>%
        as.character

    reference_signatures %>% 
        gather(signature, probability, -mutation_type) %>% 
        mutate(signature = factor(signature, levels = colnames(reference_signatures))) %>%
        filter(signature %in% involved_signatures) %>% 
        spread(signature, probability)
}

get_population_signatures <- function(
    mutation_table, 
    reference_signatures = NULL, 
    subset_signatures = TRUE,
    n_clusters = NULL,
    genome = BSgenome.Hsapiens.UCSC.hg19,
    n_chains = 10,
    n_cores = n_chains,
    n_iter = 300,
    n_adapt = 200,
    method = 'mcmc'
) {
    mutation_table <- mutation_table %>%
        distinct() %>%
        mutate(
            correction = get_vaf_correction(.),
            mutation_type = get_snv_mutation_type(
                chr, pos, ref, alt, genome
            )
        )

    n_mutations <- dim(mutation_table)[1]

    stopifnot(
        mutation_table$tumour_content %>% unique %>% length == 1
    )

    tumour_content <- mutation_table$tumour_content[1]

    catalog <- mutation_table %>%
        group_by(mutation_type) %>%
        summarise(count = n()) %>%
        ungroup()

    if (is.null(n_clusters)) {
        message('No number of populations provided. Will run automatic model selection.')
        n_population_selection <- select_n_populations(mutation_table)

        n_clusters <- n_population_selection$optimal_n_populations
        mu <- n_population_selection$evaluated_models[[n_clusters]]$parameter_estimates$mu
        kappa <- n_population_selection$evaluated_models[[n_clusters]]$parameter_estimates$kappa

        message(sprintf('Engaged auto-selection of clusters. Chose a model with %s clusters.', n_clusters))
    } else {
        message(sprintf('Extracting parameters from population model with %s populations', n_clusters))

        n_population_selection <- NULL
        population_mcmc <- get_population(mutation_table, n_clusters)
        mu <- population_mcmc$parameter_estimates$mu
        kappa <- population_mcmc$parameter_estimates$kappa

        message('Done extracting population parameters.')
    }

    if (is.null(reference_signatures)) {
        reference_signatures <- get_reference_signatures() 
    }

    if (subset_signatures) {
        message('Subsetting reference signatures on request')
        reference_signatures <- subset_reference_signatures(catalog, reference_signatures)
    }

    message(sprintf(
        'Running model with provided reference signatures: %s signatures, %s mutation types, and %s clusters',
        reference_signatures %>% select(-mutation_type) %>% ncol,
        mutation_table %>% nrow,
        n_clusters
    ))

    output <- run_joint_population_signatures(
        mutation_table,
        reference_signatures,
        n_populations = n_clusters,
        n_chains = n_chains,
        n_cores = n_cores,
        n_iter = n_iter,
        n_adapt = n_adapt,
        method = method
    )

    output['model_selection_data'] = n_population_selection

    return(output)
}

run_joint_population_signatures <- function(
    mutation_table, 
    reference_signatures, 
    n_populations,
    n_chains = 10,
    n_cores = n_chains,
    n_iter = 300,
    n_adapt = 200,
    method = 'mcmc'
) {
    stan_dso = get_stan_model('joint')

    stan_data = list(
        N = nrow(mutation_table),

        # Signature-specific data
        K = reference_signatures %>% select(-mutation_type) %>% dim %>% .[2],
        R = reference_signatures %>% dim %>% .[1],
        v = factor(mutation_table$mutation_type, levels = reference_signatures$mutation_type) %>% as.numeric,
        ref_signatures = reference_signatures %>% select(-mutation_type) %>% t,

        # Clonality-specific data
        L = n_populations,
        x = as.numeric(mutation_table$alt_depth),
        d = as.numeric(mutation_table$total_depth),
        a = as.numeric(mutation_table$correction)
    )

    message('Sampling')

    parameter_dimension <- stan_data$L * stan_data$K

    if (method == 'mcmc') {
        stan_fit = sampling(
            object = stan_dso,
            data = stan_data,
            chains = n_chains,
            iter = n_iter + n_adapt,
            warmup = n_adapt,
            cores = n_cores
        )
    } else if (method == 'vb') {
        stan_fit = vb(
            object = stan_dso,
            data = stan_data
        )
    } else {
        stop('Method argument must be either "mcmc" or "vb"')
    }

    return(list(
        mutation_table = mutation_table,
        reference_signatures = reference_signatures,
        n_populations = n_populations,
        mcmc_output = stan_fit
    ))
}
