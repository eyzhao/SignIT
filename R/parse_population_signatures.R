#' Stan Output Parser
#'
#' Converts Stan output to a long form tidy tibble
#'
#' @param population_mcmc_output        Output from \code{\link{get_population_signatures}}.
#'
#' @return Tibble with chain information.
#'
#' @import dplyr
#' @import tidyr

parse_stan_output <- function(population_mcmc_out) {
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

#' Signature Index Mapper
#'
#' Maps parameter values to their respective mutation signatures
#'
#' @param phi_index         Parameter index in output from \code{\link{get_population_signatures}}.
#' @param signature_names   Vector of mutation signature names in order of appearance in reference signatures.
#'
#' @return Vector of signature names corresponding to parameter indices.

map_phi_to_signature <- function(phi_index, signature_names) {
  n_signatures <- length(signature_names)
  return(signature_names[((phi_index - 1) %% n_signatures) + 1])
}

#' Population Index Mapper
#'
#' Maps parameter values to their respective populations
#'
#' @param phi_index         Parameter index in output from \code{\link{get_population_signatures}}.
#' @param signature_names   Vector of mutation signature names in order of appearance in reference signatures.
#' @param n_populations     Number of populations which the model was run with
#'
#' @return Vector of population names corresponding to parameter indices.

map_phi_to_population <- function(phi_index, signature_names, n_populations) {
  return(ceiling(phi_index / length(signature_names)))
}

#' Reverse Factor Levels
#'
#' Takes a factor, reorders the levels, and reassigns values to the new levels.
#'
#' @param x     A factor
#' @return A factor with levels reversed

reverse_factor_levels <- function(x) {
  levels(x) <- rev(levels(x))
  return(x)
}

#' Reverse Factor Ranks
#'
#' Takes a factor and re-ranks its levels. None of the values change, but the level ranks are reordered.
#'
#' @param x     A factor
#' @return A factor with levels re-ranked

reverse_factor_ranks <- function(x) {
  factor(as.character(x), levels = rev(levels(x)))
}

#' Maps Population Signature Parameters to their Signatures and Population Names
#'
#' Population signature model indices each correspond with specific signature and population indices.
#'
#' @param joint_model_output        Output from \code{\link{get_population_signatures}}.
#' @return Same as input, with two new columns for signature name and population name.

map_population_signatures <- function(joint_model_output) {
  signature_names <- joint_model_output$reference_signatures %>%
    select(-mutation_type) %>%
    colnames
  n_populations <- joint_model_output$n_populations

  joint_model_output$mcmc_output %>%
    parse_stan_output %>%
    filter(parameter_name == 'phi') %>%
    mutate(
      population = map_phi_to_population(parameter_index, signature_names, n_populations),
      population = cut(
        population,
        breaks = 0:n_populations,
        labels = paste('Population', 1:n_populations)
      ) %>% reverse_factor_levels %>% reverse_factor_ranks,
      signature = factor(
        map_phi_to_signature(parameter_index, signature_names),
        levels = signature_names
      )
    )
}

#' Summarise Population Signatures
#'
#' Report summary statistics (mean, median, mode, and SD) for population signature parameters
#'
#' @param joint_model_output        Output from \code{\link{get_population_signatures}}.
#' @return Tibble of summary statistics
#'
#' @import dplyr
#' @export

summarise_population_signatures <- function(joint_model_output) {
  n_populations <- joint_model_output$n_populations

  populations <- joint_model_output$mcmc_output %>%
    parse_stan_output %>%
    filter(parameter_name == 'mu') %>%
    mutate(
      population = cut(
        parameter_index,
        breaks = 0:n_populations,
        labels = paste('Population', 1:n_populations)
      ) %>% reverse_factor_levels %>% reverse_factor_ranks
    ) %>%
    select(population, value) %>%
    group_by(population) %>%
    summarise(
        prevalence_mean = mean(value)
    ) %>%
    ungroup()

  joint_model_output %>%
    map_population_signatures %>%
    group_by(population, signature) %>%
    summarise(
      mean = mean(value),
      median = median(value),
      mode = get_density_mode(value),
      sd = sd(value)
    ) %>%
    ungroup() %>%
    group_by(
      population
    ) %>%
    mutate(
      population_proportion = sum(mean),
      mean = mean / sum(mean),
      median = median / sum(median),
      mode = mode / sum(mode),
      sd = sd / sum(mean)
    ) %>%
    left_join(populations, by = 'population')
}

#' Plot Population Signatures
#'
#' Creates a violin plot of the full posterior estimates for population signature exposures.
#'
#' @param joint_model_output    Output from \code{\link{get_population_signatures}}.
#' @return ggplot object with population mutation prevalences, proportions, and signature exposures.
#' 
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @import cowplot
#' @export

plot_population_signatures <- function(joint_model_output) {
  signature_names <- joint_model_output$reference_signatures %>%
    select(-mutation_type) %>%
    colnames

  n_populations <- joint_model_output$n_populations

  phi_df <- joint_model_output %>%
    map_population_signatures %>%
    group_by(iteration, chain, population) %>%
    mutate(phi_normalized = value / sum(value)) %>%
    ungroup()

  phi_distribution_plot <- phi_df %>%
    ggplot(aes(
      x = signature,
      y = phi_normalized,
      fill = population,
      colour = population
    )) +
    geom_violin(position = position_dodge(width = 0.6)) +
    scale_colour_brewer(palette = 'Set1') +
    scale_fill_brewer(palette = 'Set1') +
    labs(x = 'Signature', y = 'Exposure Fraction', colour='Population', fill='Population') +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    )

  clone_prop_distribution_plot <- phi_df %>%
    group_by(chain, iteration, population) %>%
    summarise(clone_prop = sum(value)) %>%
    ggplot(aes(
      x = population,
      y = clone_prop,
      fill = population,
      colour = population
    )) +
    geom_violin() +
    coord_flip() +
    scale_colour_brewer(palette = 'Set1') +
    scale_fill_brewer(palette = 'Set1') +
    labs(y = 'Proportion') +
    theme(axis.title.y = element_blank())

  mu_df <- joint_model_output$mcmc_output %>%
    parse_stan_output %>%
    filter(parameter_name == 'mu') %>%
    mutate(
      population = cut(
        parameter_index,
        breaks = 0:n_populations,
        labels = paste('Population', 1:n_populations)
      ) %>% reverse_factor_levels %>% reverse_factor_ranks
    )

  mu_distribution_plot <- mu_df %>% ggplot(aes(
    x = population,
    y = value,
    fill = population,
    colour = population
  )) +
    geom_violin() +
    coord_flip() +
    scale_colour_brewer(palette = 'Set1') +
    scale_fill_brewer(palette = 'Set1') +
    labs(y = 'Prevalence')

  clone_plots <- plot_grid(
    clone_prop_distribution_plot +
      theme(legend.position = 'none'),
    mu_distribution_plot +
      theme(
        legend.position = 'none',
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank()
      ),
    nrow = 1,
    align = 'h',
    rel_widths = c(6, 4)
  )

  plot_grid(
    phi_distribution_plot + theme(legend.position = 'none'),
    clone_plots,
    ncol = 1,
    rel_heights = c(3,1)
  ) %>%
    plot_grid(
      get_legend(phi_distribution_plot),
      nrow = 1,
      rel_widths = c(3,1)
    )
}

#' Log Likelihood of Population Signatures Model
#'
#' Reports the log likelihood of data based on population model parameter estimates
#'
#' @param mcmc_output       Output from \code{\link{get_population_signatures}}.
#'
#' @return Log likelihood value associated with data, model, and parameter estimates
#'
#' @import tidyr
#' @import dplyr

compute_population_signatures_log_likelihood <- function(mcmc_output) {
    parameter_df <- mcmc_output %>%
        summarise_population_signatures() %>%
        mutate(
            coefficient = mean * population_proportion
        ) %>%
        select(population, signature, coefficient, prevalence_mean)

    ref_signatures_long <- mcmc_output$reference_signatures %>% 
        gather(signature, mutation_type_probability, -mutation_type) %>% 
        group_by(signature) %>% 
        mutate(
            mutation_type_probability = mutation_type_probability / sum(mutation_type_probability)
        )

    kappa_estimate = 2 + (
        mcmc_output$mcmc_output %>% 
            parse_stan_output() %>% 
            filter(parameter_name == 'kappa_minus_two') %>% 
            .$value %>% 
            mean
    )

    log_likelihood <- mcmc_output$mutation_table %>% 
        mutate(mutation_id = row_number()) %>%
        crossing(parameter_df) %>%
        left_join(ref_signatures_long, by = c('mutation_type', 'signature')) %>%
        mutate(
            population_likelihood = rmutil::dbetabinom(
                alt_depth, 
                total_depth, 
                correction * prevalence_mean, 
                kappa_estimate
            ),
            signature_likelihood = mutation_type_probability,
            likelihood = population_likelihood * signature_likelihood
      ) %>%
      group_by(mutation_id) %>%
      summarise(
          log_likelihood = log(sum(coefficient * likelihood))
      ) %>%
      .$log_likelihood %>%
      sum
}

#' Bayesian Information Criterion for Population Signatures
#'
#' @param population_mcmc_output    Output from \code{\link{get_population_signatures}}.
#'
#' @return BIC value

compute_population_signatures_bic <- function(mcmc_output) {
    log_lik = compute_population_signatures_log_likelihood(
        mcmc_output
    )

    n_signatures = mcmc_output$reference_signatures %>% select(-mutation_type) %>% ncol()
    n_populations = mcmc_output$n_populations

    n = nrow(mcmc_output$mutation_table)
    k = n_signatures * n_populations

    return(log(n) * k - ( 2 * log_lik ))
}

#' Watanabe-Akaike Information Criterion for Population Signatures
#'
#' @param population_mcmc_output    Output from \code{\link{get_population_signatures}}.
#'
#' @return WAIC value
#'
#' @import loo
#' @export

compute_population_signatures_waic <- function(mcmc_output) {
    log_lik <- extract_log_lik(mcmc_output$mcmc_output)
    return(waic(log_lik)$waic)
}
