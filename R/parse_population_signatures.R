parse_population_signatures <- function(population_mcmc_out) {
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

map_phi_to_signature <- function(phi_index, signature_names) {
  n_signatures <- length(signature_names)
  return(signature_names[((phi_index - 1) %% n_signatures) + 1])
}

map_phi_to_population <- function(phi_index, signature_names, n_populations) {
  return(ceiling(phi_index / length(signature_names)))
}

plot_population_signatures <- function(joint_model_output) {
    joint_model_output <- readRDS('test_clonality.Rds')

    signature_names <- joint_model_output$reference_signatures %>% select(-mutation_type) %>% colnames
    n_populations <- joint_model_output$n_populations

    phi_df <- joint_model_output$mcmc_output %>%
      parse_population_signatures %>%
      filter(parameter_name == 'phi') %>%
      mutate(
        population = map_phi_to_population(parameter_index, signature_names, n_populations),
        population = cut(
          population,
          breaks = 0:n_populations, 
          labels = paste('Population', 1:n_populations)
        ),
        signature = factor(
          map_phi_to_signature(parameter_index, signature_names),
          levels = signature_names
        )
      ) %>%
      group_by(iteration, chain, population) %>%
      mutate(phi_normalized = value / sum(value)) %>%
      ungroup()

    phi_distribution_plot <- phi_df %>%
      ggplot(aes(
        x = signature,
        y = phi_normalized,
        fill = population
      )) +
      geom_violin(position = 'dodge') +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
      )

    clone_prop_distribution_plot <- phi_df %>%
      group_by(chain, iteration, population) %>%
      summarise(clone_prop = sum(value)) %>%
      ggplot(aes(
        x = population, 
        y = clone_prop,
        fill = population
      )) +
      geom_violin() +
      coord_flip()

    mu_df <- joint_model_output$mcmc_output %>%
      parse_population_signatures %>%
      filter(parameter_name == 'mu') %>%
      mutate(population = factor(paste('Population', parameter_index), paste('Population', 1:n_populations)))

    mu_distribution_plot <- mu_df %>% ggplot(aes(
      x = population,
      y = value,
      fill = population
    )) +
      geom_violin() +
      coord_flip()

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
