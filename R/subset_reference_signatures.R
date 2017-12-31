#' Automatic subsetting of reference mutation signatures
#'
#' Subsets mutation signatures based on which ones are present in a somatic mutation catalog
#'
#' @param catalog               Mutation catalog. A data-frame or tibble with columns
#'                              named mutation_type and count.
#'
#' @param reference_signatures  Reference mutation signatures, formatted as in those
#'                              obtained by \code{\link{get_reference_signatures}}.
#'
#' @param threshold             A numeric value specifying the minimum lower credible interval
#'                              (at the 95% confidence level). Any signatures with a 2.5%
#'                              interval limit above this threshold will be included in the subset.
#'
#' @return A new reference signature dataframe with a subset of the initial signatures
#'
#' @import tidyr
#' @import dplyr

subset_reference_signatures <- function(catalog, reference_signatures, threshold = NULL) {
    n_mutations <- sum(catalog$count)
    
    if (is.null(threshold)) {
        threshold = n_mutations / 2000
    }

    exposures <- catalog %>% 
        get_exposures(quiet = TRUE) %>% 
        get_exposure_summary_table %>% 
        mutate(signature_present = `2.5%` > threshold)

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
