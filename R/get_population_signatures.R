#' SignIT-Pop inference of populations and signatures
#'
#' Jointly infers mutational subpopulations and their associated mutation signature exposures
#'
#' \code{get_population_signatures} is the central function which facilitates Bayesian inference
#' of mutational populations and signatures. This model infers a matrix of L x N parameters,
#' where L is the number of populations and N is the number of signatures. The posterior distribution
#' of each parameter is estimated using either automatic differentiation variational inference
#' or Hamiltonial Monte Carlo using the \code{\link[rstan]{vb}} and \code{\link[rstan]{sampling}}
#' methods respectively of the \pkg{rstan} package (an interface to the Stan probabilistic
#' programming language).
#'
#' @param mutation_table        Table of mutations, one per row. The minimum input requires the 
#'                              following columns: \itemize{
#'                                  \item \strong{total_depth}: Total number of reads covering mutated locus.
#'                                  \item \strong{alt_depth}: Total number of mutant reads covering locus.
#'                                  \item \strong{tumour_copy}: Tumour copy number at the mutated locus
#'                                  \item \strong{normal_copy}: Normal copy number at the mutated locus
#'                                  \item \strong{tumour_content}: Estimated tumour content as a fraction
#'                                                                 between 0 and 1. Must be the same value
#'                                                                 throughout the whole table.
#'                              }
#'
#' @param reference_signatures  Reference mutation signatures. This can either be from 
#'                              \code{\link{get_reference_signatures}} or a custom data frame formatted
#'                              equivalently.
#'
#' @param subset_signatures     Boolean. If TRUE (default), then \code{\link{subset_reference_signatures}}
#'                              is run to pre-select a smaller subset of signatures most likely to be active
#'                              in the cancer. This helps to reduce processing time and model complexity, but
#'                              may bias the result.
#'
#' @param n_populations         The number of populations to screen for. Must be an integer. If no value is
#'                              provided, then a model selection step is engaged to automatically estimate
#'                              the number of populations. The automatic model selection uses
#'                              \code{\link{select_n_populations}}, which performs a maximum a posteriori
#'                              estimate using the SignIT population model (without mutation signature inference).
#'
#' @param genome                A BSgenome object. This is used to determine trinucleotide contexts of mutations
#'                              to define mutation types. By default, uses BSgenome.Hsapiens.UCSC.hg19. To
#'                              define custom mutation types, simply include a column named \code{mutation_type} in
#'                              \code{mutation_table}, in which case this parameter is ignored.
#'
#' @param method                The posterior sampling method. This is a string and can either be 'vb' for
#'                              automatic variational Bayes or 'mcmc' for Hamiltonial Monte Carlo.
#'
#' @param n_chains              Number of chains to sample. Only relevant if \code{method == 'mcmc'}.
#'
#' @param n_cores               Number of cores for parallel sampling. By default this equals the number of chains.
#'                              Only relevant if \code{method == 'mcmc'}.
#'
#' @param n_iter                Number of sampling iterations per chain. These are distinct from adaptation iterations,
#'                              so the total number of iterations will be \code{n_iter + n_adapt}.
#'                              Only relevant if \code{method == 'mcmc'}.
#'
#' @param n_adapt               Number of adaptation iterations per chain. Only relevant if \code{method == 'mcmc'}.
#'
#' @return A list object with the posterior sampling of population signatures plus relevant input and metadata.
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import BSgenome
#' @import rstan
#'
#' @export

get_population_signatures <- function(
    mutation_table, 
    reference_signatures = NULL, 
    subset_signatures = TRUE,
    n_populations = NULL,
    genome = NULL,
    method = 'vb',
    n_chains = 10,
    n_cores = n_chains,
    n_iter = 300,
    n_adapt = 200
) {
    if (is.null(genome)) {
        genome = getBSgenome('BSgenome.Hsapiens.UCSC.hg19')
    }
    mutation_table <- mutation_table %>%
        distinct() %>%
        filter(
            !is.na(chr),
            !is.na(pos),
            !is.na(ref),
            !is.na(alt),
            !is.na(total_depth),
            !is.na(alt_depth),
            !is.na(normal_copy),
            !is.na(tumour_content)
        ) %>%
        mutate(
            correction = get_vaf_correction(.)
        )

    if ('mutation_type' %in% colnames(mutation_table)) {
        message("Using pre-specified mutation types.")
    } else {
        mutation_table <- mutation_table %>%
        mutate(
            mutation_type = get_snv_mutation_type(
                chr, pos, ref, alt, genome
            )
        ) %>%
        filter(
            ! grepl('N', mutation_type)
        )
    } 

    n_mutations <- dim(mutation_table)[1]

    stopifnot(
        mutation_table$tumour_content %>% unique %>% length == 1
    )

    tumour_content <- mutation_table$tumour_content[1]

    catalog <- with(
        mutation_table,
        mutations_to_catalog(chr, pos, ref, alt)
    )

    if (is.null(n_populations)) {
        message('No number of populations provided. Will run automatic model selection.')
        n_population_selection <- select_n_populations(mutation_table)

        n_populations <- n_population_selection$optimal_n_populations
        mu <- n_population_selection$evaluated_models[[n_populations]]$parameter_estimates$mu
        kappa <- n_population_selection$evaluated_models[[n_populations]]$parameter_estimates$kappa

        message(sprintf('Engaged auto-selection of clusters. Chose a model with %s clusters.', n_populations))
    } else {
        message(sprintf('Extracting parameters from population model with %s populations', n_populations))

        n_population_selection <- NULL
    }

    if (is.null(reference_signatures)) {
        reference_signatures <- get_reference_signatures() 
    }

    if (subset_signatures) {
        message('Subsetting reference signatures on request')
        reference_signatures <- subset_reference_signatures(catalog, reference_signatures)
    }

    message(sprintf(
        'Running model with provided reference signatures: %s signatures, %s mutations, and %s clusters',
        reference_signatures %>% select(-mutation_type) %>% ncol,
        mutation_table %>% nrow,
        n_populations
    ))

    stan_dso = stanmodels$signit_model_infer_populations

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

    output <- list(
        mutation_table = mutation_table,
        reference_signatures = reference_signatures,
        n_populations = n_populations,
        mcmc_output = stan_fit
    )

    output['model_selection_data'] = n_population_selection

    return(output)
}
