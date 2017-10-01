#' Retrieves the names of mutation signatures from reference signatures
#'
#' @param reference_signatures    Reference mutation signatures (i.e. from get_reference_signatures)
#'
#' @return Character vector of mutation signature names
#'
#' @import dplyr

get_signature_names <- function(reference_signatures) {
  reference_signatures %>%
    select(-mutation_type) %>%
    colnames()
}
