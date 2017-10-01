#' Plots a hexplot between two signatures
#'
#' The MCMC sample space allows one to examine for multicollinearity between
#' mutation signatures. For more details on signature bleed, see plot_signature_pairwise_bleed.
#'
#' This function provides a closer view on correlated signatures in the
#' MCMC sample space. It plots the relationship between two signatures as a
#' hexplot, and reports the Spearman correlation coefficient.
#'
#' @param exposures_mcmc_output     Output from get_exposures.
#'
#' @param signature_name_1          String denoting signature name (must match one of those in exposures_mcmc_output)
#'
#' @param signature_name_2          String denoting signature name (must match one of those in exposures_mcmc_output)
#'
#' @param exposure_mcmc_clusters    Optionally, you can provide the output from density_clustering, which allows
#'                                  you to visualize the comparison plot for a single cluster.
#'
#' @param cluster_number            If you wish to visualize a single cluster, provide the cluster number here.
#'
#' @param trendline                 Boolean value specifying whether to plot a linear trendline.
#'
#' @return A ggplot object showing the correlation of MCMC samples between two signatures.
#'
#' @import tibble
#' @import dplyr
#' @import ggplot
#'
#' @export

plot_two_signature_hexplot <- function(exposures_mcmc_output, signature_name_1, signature_name_2, exposure_mcmc_clusters = NULL, cluster_number = NULL, trendline = TRUE) {
  hexplot_data <- exposures_mcmc_output %>%
    mcmc_exposures_as_matrix() %>%
    as_tibble %>%
    select(-iteration, -chain) %>%
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


