#' Convert reference signatures to matrix
#'
#' Converts a data frame of reference signature data to the matrix form.
#'
#' When generating the matrix, a mutation catalog is required. Because the matrix output
#' does not have a mutation_type column, the matrix rows are arranged so that mutation types
#' match the order of the provide mutation catalog. This way, the matrix output can immediately
#' in correspondence with the provided mutation catalog.
#'
#' @param reference_signatures Reference mutation signature as output by get_reference_signatures()
#'
#' @param mutation_catalog A mutation catalog with mutation types matching those of the reference
#'                         signatures. This catalog is used only to order the rows of the output
#'                         matrix such that they correspond with the mutation type order in mutation_catalog.
#'
#' @return A matrix with named columns, one for each signature.
#'
#' @import dplyr
#'
#' @export

reference_signatures_as_matrix <- function(reference_signatures, mutation_catalog) {
  if (! mutation_types_match(reference_signatures, mutation_catalog, order_is_important = FALSE)) {
    stop('mutation_type does not match Stratton 30-signature reference mutation types')
  }

  catalog_columns <- colnames(mutation_catalog)
  
  reference_signatures <- mutation_catalog %>%
    select(mutation_type) %>%
    left_join(reference_signatures, by = 'mutation_type') %>%
    select(-mutation_type)

  reference_signatures <- reference_signatures[! reference_signatures %in% catalog_columns]
  
  return(reference_signatures %>% as.matrix())
}
