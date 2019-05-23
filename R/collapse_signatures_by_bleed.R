#' Collapse signatures using bleed information
#' 
#' This is a convenience function provided. It will return all levels of bleed
#' collapsing attempted as a single large list. For a more convenient
#' implementation, use collapse_signatures_by_bleed, which auto-selects the
#' optimal solution using the WAIC statistic.
#' 
#' @param exposure_mcmc_output Output from get_exposures
#' @param bleed_threshold The higher this is set, the fewer levels of collapsing will be attempted
#' @importFrom purrr accumulate
#' @export
#' @return List with all attempted collapsed exposure objects

collapse_signatures_by_bleed_all_solutions <- function(exposure_mcmc_output, bleed_threshold = 0.4) {
  accumulate(
    c(list(exposure_mcmc_output), (length(exposure_mcmc_output$signature_names)-1) : 3),
    function(previous_exposures, n_signatures) {
    
    if (is.null(previous_exposures)) {
      return(NULL)
    }
      
    ref <- previous_exposures$reference_signatures
  
    bleed_table <- get_exposure_pairwise_correlations(previous_exposures) %>%
      filter(as.numeric(signature_2) > as.numeric(signature_1))
    
    if (min(bleed_table$spearman) < -bleed_threshold) {
      top_bleed <- bleed_table %>%
        arrange(spearman) %>%
        filter(row_number() == 1) %>%
        gather(id, name, signature_1, signature_2) %>%
        .$name
      
      message(sprintf('Max bleed was %s. Merging %s.', -min(round(bleed_table$spearman, 3)), paste(top_bleed, collapse=' and ')))
      
      if (!is.null(previous_exposures$merged_signatures)) {
        signatures_to_merge %>%
          mutate(
            signature %in% previous_exposures$merged_signatures$merged_signature
          )
      } else {
        signatures_to_merge <- get_exposure_summary_table(previous_exposures) %>%
          filter(signature %in% top_bleed) %>%
          mutate(coefficient = mean_exposure / sum(mean_exposure)) %>%
          select(signature, coefficient) %>%
          mutate(
            signature = as.character(signature),
            merged_signature = paste(signature, collapse='/')
          )
      }
      
      merged_signature <- signatures_to_merge %>%
        left_join(
          ref %>% gather(signature, probability, -mutation_type) %>% mutate(signature = as.character(signature)),
          by = 'signature'
        ) %>%
        group_by(mutation_type) %>%
        summarise(
          signature = unique(merged_signature),
          probability = sum(probability * coefficient)
        ) %>%
        ungroup()
      
      new_ref_signatures <- ref %>%
        gather(signature, probability, -mutation_type) %>%
        filter(! signature %in% top_bleed) %>%
        bind_rows(merged_signature) %>%
        spread(signature, probability)
      
      capture.output(new_exposures <- get_exposures(previous_exposures$mutation_catalog, new_ref_signatures, quiet=T))
      
      message(sprintf('%s signature solution obtained.', new_ref_signatures %>% select(-mutation_type) %>% names() %>% length))

      new_exposures$waic = compute_signatures_waic(new_exposures)
      new_exposures$signature_names <- c(
        previous_exposures$signature_names[previous_exposures$signature_names %in% new_exposures$signature_names],
        new_exposures$signature_names[! new_exposures$signature_names %in% previous_exposures$signature_names]
      ) %>% unique

      new_exposures$exposure_chain$signature = factor(new_exposures$exposure_chain$signature, levels = new_exposures$signature_names)
      new_exposures$reference_signatures <- new_exposures$reference_signatures %>%
         gather(signature, exposure, -mutation_type) %>%
         mutate(signature = factor(signature, levels = new_exposures$signature_names)) %>%
         arrange(signature, mutation_type) %>%
         spread(signature, exposure)
      
      return(new_exposures)
    } else {
      return(NULL)
    }
  })
}

#' Collapse signatures using bleed information
#' 
#' @param exposure_mcmc_output Output from get_exposures
#' @param bleed_threshold The higher this is set, the fewer levels of collapsing will be attempted
#' @export
#' @return New exposure object corresponding to the collapsed solution with
#'         greatest WAIC

collapse_signatures_by_bleed <- function(exposure_mcmc_output, bleed_threshold=0.4) {
  new_exposures = collapse_signatures_by_bleed_all_solutions(exposure_mcmc_output, bleed_threshold=bleed_threshold)
  
  message('Choosing signature set with maximum WAIC')
  
  waic_values <- sapply(new_exposures, function(z) {
    if (!is.null(z)) {
      z$waic
    }
  }) %>% unlist
  
  exposure_index <- which(waic_values == min(waic_values))
  
  output <- new_exposures[[exposure_index]]
  output[['waic_values']] = waic_values
  output[['bleed_threshold']] = 0.4
  output[['exposure_index']] = exposure_index
  output[['n_signatures']] = length(new_exposures[[exposure_index]]$signature_names)
  message(sprintf('%s signature solution chosen', output[['n_signatures']]))
  return(output)
}
