#' Check Whether Mutation Types Match
#'
#' Returns boolean value indicating whether two catalogs/reference signatures have matching mutation types
#'
#' Given two mutation catalogs, check that their mutation types match.
#' This is an important safeguard against potential errors which can arise
#' when catalogs have mismatched mutation types.
#'
#' @param catalogs_1 Mutation catalog object (data frame with cols: mutation_type, count)
#'
#' @param catalogs_2 Mutation catalog object (data frame with cols: mutation_type, count)
#'
#' @param order_is_important Boolean. If FALSE, then the catalogs are sorted by mutation
#'                           type before they are compared.
#'
#' @return Boolean value. TRUE if mutation types match. FALSE otherwise.
#'
#' @import magrittr
#'
#' @export

mutation_types_match <- function(catalogs_1, catalogs_2, order_is_important = TRUE) {
  mutation_types_1 <- catalogs_1[['mutation_type']]
  mutation_types_2 <- catalogs_2[['mutation_type']]

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
