#' Compute sum of two traits
#'
#' Compute the sum of two traits. Missing values do not propagate.
#' @param a Name of trait 1 to sum
#' @param b Name of trait 2 to sum.
#' @details Missing values do not propagate. If NA is present for both traits then NA
#' is applied to the sum.
#' @return It returns the sum of the two traits.
#' @author Raul Eyzaguirre.
#' @examples
#' ## Compute total biomass as the sum of trw and vw
#' suma(pjpz09$trw, pjpz09$vw)
#' @export

suma <- function(a, b) {
  
  s <- apply(cbind(a, b), 1, sum, na.rm = TRUE)
  s[is.na(a) & is.na(b)] <- NA
  s
  
}