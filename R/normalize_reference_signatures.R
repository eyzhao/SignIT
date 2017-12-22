#' Normalization of reference signatures.
#'
#' Ensures that all columns of a reference signature matrix add up to 1,
#' and that all non-zero values are replaced with a very small number
#' (.Machine$double.eps).
#'
#' @param reference_signature_df        Reference signature data frame, such as the
#'                                      output of get_reference_signatures()
#'
#' @return Reference signatures in the same format, but normalized and non-zero
#'
#' @import tidyr
#' @import dplyr

normalize_reference_signatures <- function(reference_signature_df) {
    # Ensure reference signature columns add to 1
    reference_signature_df %>% 
        gather(signature, probability, -mutation_type) %>% 
        group_by(signature) %>% 
        mutate(
            probability = if_else(probability == 0, .Machine$double.eps, probability),
            probability = probability / sum(probability)
        ) %>% 
        ungroup() %>% 
        spread(signature, probability)
}