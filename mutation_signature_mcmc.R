library(tidyverse)
library(rjags)
library(purrr)
library(nnls)
library(dbscan)
library(cowplot)
library(Rtsne)

get_reference_signatures <- function(signature_set = 'cosmic_30') {
  available_reference_sets <- c(
    'cosmic_30'
  )

  if (! signature_set %in% available_reference_sets) {
    stop(sprintf('Signature set %s is not available', signature_set))
  }
  
  has_data <- function(x) { sum(!is.na(x)) > 0 }

  if (signature_set == 'cosmic_30') {  
    read_tsv('signature_matrix.txt') %>%
      select_if(has_data) %>%
      select(-`Substitution Type`, -Trinucleotide) %>%
      rename(mutation_type = `Somatic Mutation Type`) %>%
      arrange(mutation_type)
  }
}

mutation_types_match <- function(signatures_1, signatures_2, order_is_important = TRUE) {
  mutation_types_1 <- signatures_1[['mutation_type']]
  mutation_types_2 <- signatures_2[['mutation_type']]

  if (order_is_important) {
    return(identical(
      mutation_types_1,
      mutation_types_2
    ))
  } else {
    return(identical(
      mutation_types_1 %>% sort,
      mutation_types_2 %>% sort
    ))
  }
}

compute_catalogs_from_exposures <- function(reference_signatures, exposure_vector) {
  
}

reference_signatures_as_matrix <- function(reference_signatures, mutation_catalog) {
  if (! mutation_types_match(reference_signatures, mutation_catalog, order_is_important = FALSE)) {
    stop('mutation_type does not match Stratton 30-signature reference mutation types')
  }
  
  reference_signatures <- mutation_catalog %>%
    select(mutation_type) %>%
    left_join(reference_signatures, by = 'mutation_type')
  
  return(reference_signatures[, grepl('Signature ', names(reference_signatures))] %>%  as.matrix())
}

get_signature_model <- function(mutation_catalog, reference_signatures = NULL, n_chains = 4, n_adapt = 10000) {
  if (! mutation_types_match(reference_signatures, mutation_catalog, order_is_important = FALSE)) {
    stop('mutation_type does not match reference mutation types')
  }
  
  reference_signature_matrix <-  reference_signatures_as_matrix(reference_signatures, mutation_catalog)
  
  naive_nnls_solution <- nnls(reference_signature_matrix, mutation_catalog[['count']])$x
  naive_nnls_solution <- if_else(naive_nnls_solution == 0, 1, naive_nnls_solution) %>% ceiling()

  n_mutations <- sum(mutation_catalog[['count']])
  
  signature_data = list(
    mutation_catalog = mutation_catalog[['count']],
    exposure_dirichlet_prior = rep(1, length(naive_nnls_solution)),
    mutation_signatures = reference_signature_matrix,
    n_mutations = n_mutations
  )
  signature_init = list(
    exposures = naive_nnls_solution
  )
  
  signature_model = '
    model {
      mutation_catalog ~ dmulti((mutation_signatures %*% exposures) / n_mutations, n_mutations);
      exposures ~ ddirch(exposure_dirichlet_prior)
    }
  '
  
  jags.model(
    file=textConnection(signature_model), 
    data=signature_data,
    inits=signature_init,
    n.chains=n_chains, 
    n.adapt=n_adapt
  )
}

run_signature_model <- function(signature_model, signature_names, n_mutations, n_iter=10000) {
  update(
    signature_model, 
    n.iter=n_iter
  )
  
  signature_samples = coda.samples(
    signature_model, 
    variable.names=c("exposures"),
    n.iter=n_iter
  ) %>%
    lapply(function(chain) {
      chain %>%
        as_tibble %>%
        `colnames<-`(signature_names) %>%
        mutate(
          iteration = row_number()
        )
    })
  
  signature_full_chain <- do.call('rbind', signature_samples) %>%
    gather(signature, exposure, -iteration) %>%
    mutate(
      signature = signature,
      exposure = exposure * n_mutations
    )
  
  return(signature_full_chain)
}

plot_exposure_chain <- function(exposures_output) {
  exposures_output$exposure_chain %>%
    ggplot(aes(
      y = exposure,
      x = iteration
    )) +
    facet_grid(signature ~ .) +
    geom_line()
}

get_signature_names <- function(reference_signatures) {
  reference_signatures %>%
    select(-mutation_type) %>%
    colnames()
}

get_exposures <- function(mutation_catalog = NULL, file = NULL, reference_signatures = NULL, n_chains = 4, n_iter = 10000) {
  if (is.null(mutation_catalog) && is_null(file)) {
    stop('Must either provide mutation_catalog as data frame or path to TSV file')
  } else if (is.null(mutation_catalog)) {
    mutation_catalog <- read_tsv(
      file, 
      col_types = cols(
        mutation_type = col_character(),
        count = col_integer()
      )
    )
  }
  
  if (! 'mutation_type' %in% names(mutation_catalog) || ! 'count' %in% names(mutation_catalog)) {
    stop('Mutation catalog is not properly formatted. Must have two columns: mutation_type (character) and count (integer)')
  }
  
  if (is.null(reference_signatures)) {
    reference_signatures <- get_reference_signatures()
  }
  
  signature_model <- get_signature_model(
    mutation_catalog, 
    reference_signatures = reference_signatures,
    n_chains = n_chains
  )
  
  signature_names <- get_signature_names(reference_signatures)
  
  signature_exposures <- run_signature_model(
    signature_model, 
    signature_names, 
    n_mutations = mutation_catalog[['count']] %>% sum,
    n_iter = n_iter
  )
  
  return(list(
    mutation_catalog = mutation_catalog,
    exposure_chain = signature_exposures %>% mutate(signature = factor(signature, levels = signature_names)),
    reference_signatures = reference_signatures,
    signature_names = signature_names
  ))
}

plot_exposure_posteriors <- function(exposures_mcmc_output, view = 'violin') {
  plot <- exposures_mcmc_output$exposure_chain %>%
    mutate(signature = trim_signature_names(signature)) %>%
    ggplot(aes(
      x = signature %>% as.factor,
      y = exposure
    )) +
    labs(x = 'Signature', y = 'Exposure') +
    rotate_x_axis_labels()
    
  
  if (view == 'boxplot') {
    return(plot + geom_boxplot())
  } else {
    return(plot + geom_violin())
  }
}

mcmc_exposures_as_matrix <- function(exposures_mcmc_output) {
  exposures_mcmc_output$exposure_chain %>%
    select(signature, exposure) %>%
    group_by(signature) %>%
    mutate(index = row_number()) %>%
    ungroup() %>%
    spread(signature, exposure) %>%
    as.matrix()
}

density_clustering <- function(exposures_mcmc_output, minPts = 40) {
  library(dbscan)
  
  exposures_matrix <- mcmc_exposures_as_matrix(exposures_mcmc_output)
  
  # This will take a few minutes
  patient_density_clustered <- hdbscan(exposures_matrix, minPts=minPts)
  
  return(patient_density_clustered)
}

get_cluster_max_densities <- function(exposures_mcmc_output, exposure_mcmc_clusters, k_value = 400) {
  exposures_matrix <- mcmc_exposures_as_matrix(exposures_mcmc_output)
  
  max_density_per_cluster <- merge_exposures_with_clustering(exposures_mcmc_output, exposure_mcmc_clusters) %>%
    mutate(mean_knn_dist = kNNdist(exposures_matrix, k = k_value) %>% apply(1, mean)) %>%
    group_by(cluster) %>%
    filter(mean_knn_dist == max(mean_knn_dist)) %>%
    ungroup() %>%
    left_join(
      tibble(
        cluster_score = exposure_mcmc_clusters$cluster_scores
      ) %>%
        mutate(cluster = row_number()),
      by = c('cluster')
    ) %>%
    gather(signature, exposure, -index, -cluster, -mean_knn_dist, -cluster_score)
  
  return(max_density_per_cluster)
}

merge_exposures_with_clustering <- function(exposures_mcmc_output, exposure_mcmc_clusters) {
  exposures_mcmc_output %>% 
    mcmc_exposures_as_matrix %>%
    as_tibble %>%
    mutate(
      cluster = exposure_mcmc_clusters$cluster
    )
}

get_nnls_rss <- function(exposures_mcmc_output) {
  catalog <- exposures_mcmc_output$mutation_catalog %>% .$count
  exposures <- exposures_mcmc_output %>% mcmc_exposures_as_matrix() %>% as_tibble %>% select(-index)
  reference <- exposures_mcmc_output$reference_signatures %>% select(-mutation_type) %>% as.matrix
  
  apply(exposures, 1, function(e) {
    error_rss <- (catalog - (reference %*% e))^2 %>% sum
  })
}

get_cluster_best_fits <- function(exposures_mcmc_output, exposure_mcmc_clusters) {
  best_fit_per_cluster <- merge_exposures_with_clustering(exposures_mcmc_output, exposure_mcmc_clusters) %>%
    mutate(rss = get_nnls_rss(exposures_mcmc_output)) %>%
    group_by(cluster) %>%
    filter(rss == min(rss)) %>%
    ungroup() %>%
    select(-index) %>%
    distinct() %>%
    left_join(
      tibble(
        cluster_score = exposure_mcmc_clusters$cluster_scores
      ) %>%
        mutate(cluster = row_number()),
      by = c('cluster')
    ) %>%
    gather(signature, exposure, -cluster, -rss, -cluster_score)
  
  return(best_fit_per_cluster)
}

get_cluster_summary_points <- function(exposures_mcmc_output, exposure_mcmc_clusters, method = 'bestfit', k_value = 400) {
  if (method == 'mostdense') {
    max_density_per_cluster <- get_cluster_max_densities(
      exposures_mcmc_output, 
      exposure_mcmc_clusters,
      k_value
    )
  } else if (method == 'bestfit') {
    best_fit_per_cluster <- get_cluster_best_fits(exposures_mcmc_output, exposure_mcmc_clusters)
  } else {
    stop('Invalid value provided to argument "method"')
  }
}

trim_signature_names <- function(signature) {
  signature %>% `levels<-`(gsub('Signature ', '', levels(signature)))
}

plot_cluster_solutions <- function(exposures_mcmc_output, exposure_mcmc_clusters, method = 'bestfit', k_value = 400) {
  cluster_summary_points <- get_cluster_summary_points(
    exposures_mcmc_output, 
    exposure_mcmc_clusters, 
    method, 
    k_value
  ) %>%
    mutate(signature = signature %>% factor(levels = exposures_mcmc_output$signature_names))
  
  cluster_summary_points %>%
    mutate(signature = trim_signature_names(signature)) %>%
    ggplot(aes(
      x = signature,
      y = exposure
    )) +
    geom_line(aes(
      group = cluster,
      alpha = cluster_score
    )) +
    rotate_x_axis_labels()
    
}

plot_cluster_solution_errors <- function(exposures_mcmc_output, exposure_mcmc_clusters, method = 'bestfit') {
  cluster_summary_points <- get_cluster_summary_points(exposures_mcmc_output, exposure_mcmc_clusters, method)
  
  cluster_summary_points %>%
    filter(!is.na(cluster_score)) %>%
    plyr::ddply('cluster', function(cluster_table) {
      cluster_exposure <- tibble(signature = exposures_mcmc_output$signature_names) %>%
        left_join(cluster_table, by = 'signature') %>%
        .$exposure
      cluster_catalog <- exposures_mcmc_output$reference_signatures %>% 
        reference_signatures_as_matrix(patient_signature_exposures$mutation_catalog) %*% 
        cluster_exposure
      sum_of_squared_errors <- (patient_catalog$count - cluster_catalog)^2 %>% sum
      return(cluster_table[1, ] %>% select(cluster, cluster_score) %>% mutate(error = sum_of_squared_errors))
    }) %>%
    ggplot(aes(
      x = cluster_score,
      y = error
    )) +
    geom_point()
}

plot_nnls_solution <- function(exposures_mcmc_output) {
  nnls_solution <- nnls(
    exposures_mcmc_output$reference_signatures %>% select(-mutation_type) %>% as.matrix,
    exposures_mcmc_output$mutation_catalog$count
  )$x
  
  nnls_plot <- tibble(
      signature = exposures_mcmc_output$signature_names %>% factor(levels = exposures_mcmc_output$signature_names),
      nnls = nnls_solution
    ) %>%
    mutate(signature = trim_signature_names(signature)) %>%
    ggplot(aes(x = signature, y = nnls)) +
    geom_point() +
    rotate_x_axis_labels()
}

plot_two_signature_hexplot <- function(exposures_mcmc_output, signature_name_1, signature_name_2, exposure_mcmc_clusters = NULL, cluster_number = NULL, trendline = TRUE) {
  hexplot_data <- exposures_mcmc_output %>%
    mcmc_exposures_as_matrix() %>%
    as_tibble %>%
    select(-index) %>%
    `colnames<-`(exposures_mcmc_output$signature_names) 
  
  if (is.null(cluster_number) | is.null(cluster_number)) {
    title_text = 'All Clusters'
  } else {
    hexplot_data <- hexplot_data %>% 
      mutate(cluster = exposure_mcmc_clusters$cluster %>% as.factor) %>%
      filter(cluster == cluster_number)

    title_text = paste0('Cluster ', cluster_number)
  }
  
  hexplot_data <- hexplot_data[, c(
    signature_name_1 %>% as.character,
    signature_name_2 %>% as.character
  )]
  colnames(hexplot_data) <- c('signature_1', 'signature_2')
  
  spearman_rho <- with(hexplot_data, cor(signature_1, signature_2, method = 'spearman')) %>% round(3)
  
  hexplot <- hexplot_data %>%
    ggplot(aes(
      x = signature_1, 
      y = signature_2
    )) +
    labs(
      x = signature_name_1,
      y = signature_name_2,
      title = title_text
    ) +
    geom_hex(bins = 15) +
    scale_fill_distiller(palette = 'Spectral') +
    annotate(
      'text',
      label = paste0(
        'Spearman Rho = ', spearman_rho), 
      x = max(hexplot_data$signature_1), 
      y = max(hexplot_data$signature_2) + 0.2 * (max(hexplot_data$signature_2) - min(hexplot_data$signature_2)), 
      hjust=1
    )
  
  if (trendline) {
    linear_fit <- lm(signature_2 ~ signature_1, data = hexplot_data)
    lm_slope <- linear_fit$coefficients[['signature_1']]
    lm_intercept <- linear_fit$coefficients[['(Intercept)']]
    
    hexplot <- hexplot + geom_abline(slope = lm_slope, intercept = lm_intercept)
  }
  
  return(hexplot)
}

plot_density_tsne <- function(exposures_mcmc_output, patient_density_clustering) {
  tsne_data <- exposures_mcmc_output %>% mcmc_exposures_as_matrix() %>% Rtsne
  
  tsne_data$Y %>% 
    as_tibble %>%
    mutate(
      cluster = patient_density_clustering$cluster %>% as.factor %>% 
        recode(
          `0` = 'Outlier'
        )
    ) %>%
    ggplot(aes(
      x = V1,
      y = V2,
      colour = cluster
    )) +
    geom_path(colour = 'grey80') +
    geom_point()
}

get_exposure_pairwise_correlations <- function(exposures_mcmc_output) {
  exposures_mcmc_output$exposure_chain %>%
    rename(signature_1 = signature) %>%
    plyr::ddply('signature_1', function(signature_1_chain) {
      exposures_mcmc_output$exposure_chain %>%
        rename(signature_2 = signature) %>%
        plyr::ddply('signature_2', function(signature_2_chain) {
          tibble(spearman = cor(signature_1_chain$exposure, signature_2_chain$exposure, method = 'spearman'))
        })
    })
}

plot_signature_pairwise_bleed <- function(exposures_mcmc_output) {
  get_exposure_pairwise_correlations(exposures_mcmc_output) %>%
    mutate(
      signature_1 = signature_1 %>% `levels<-`(gsub('Signature ', '', levels(signature_1))),
      signature_2 = signature_2 %>% `levels<-`(gsub('Signature ', '', levels(signature_2)))
    ) %>%
    ggplot(aes(x = signature_1, y = signature_2, fill = spearman)) +
    geom_tile(width = 0.8, height = 0.8) +
    scale_fill_distiller(palette = 'Spectral') +
    rotate_x_axis_labels() +
    labs(x = 'First Signature', y = 'Second Signature', fill = 'Spearman Rho\n')
}

rotate_x_axis_labels <- function() {
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))
}