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
#' @param trendline                 Boolean value specifying whether to plot a linear trendline.
#'
#' @return A ggplot object showing the correlation of MCMC samples between two signatures.
#'
#' @import tibble
#' @import dplyr
#' @import ggplot2
#'
#' @export

plot_two_signature_hexplot <- function(exposures_mcmc_output, signature_name_1, signature_name_2, trendline = TRUE) {
  hexplot_data <- exposures_mcmc_output %>%
    mcmc_exposures_as_matrix() %>%
    as_tibble %>%
    select(-iteration, -chain) %>%
    `colnames<-`(exposures_mcmc_output$signature_names) 
  
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
      y = signature_name_2
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


