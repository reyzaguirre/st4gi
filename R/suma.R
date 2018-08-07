#' Compute sum of two vectors
#'
#' Compute the sum of two vectors elementwise. Missing values do not propagate.
#' @param v1 Name of vector 1.
#' @param v2 Name of vector 2.
#' @details Missing values do not propagate. If NA is present for both vectors then NA
#' is applied to the sum.
#' @return It returns the sum of the two vectors elementwise.
#' @author Raul Eyzaguirre.
#' @examples
#' suma(pjpz09$trw, pjpz09$vw)
#' @export

suma <- function(v1, v2) {
  
  v3 <- apply(cbind(v1, v2), 1, sum, na.rm = TRUE)
  v3[is.na(v1) & is.na(v2)] <- NA
  
  # Return
  
  v3
  
}