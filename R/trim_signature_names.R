#' Trims the names of mutation signatures
#'
#' Any signature names starting with "Signature" get trimmed.
#'
#' @param signature     A factor containing signature names.
#' @return A factor where all names starting with "Signature" are trimmed.
#'
#' @import magrittr

trim_signature_names <- function(signature, replace_string='Signature') {
  signature %>% 
    `levels<-`(
        trimws(
            gsub(replace_string, '', levels(signature)
            )
        )
    )
}
