#' Factorization
#'
#' Decompose a number in its prime factors.
#' 
#' @param n Number for factorization.
#' 
#' @details It decomposes a number in its prime factors.
#' 
#' @return It returns the prime factors of a number.
#' 
#' @author Raul Eyzaguirre.
#' @examples
#' factors(18)
#' 
#' @export

factors <- function(n) {
  
  print(paste("The factors of",n,"are:"))
  for(i in 2:n) {
    if((n %% i) == 0) {
      print(i)
    }
  }
  
}
