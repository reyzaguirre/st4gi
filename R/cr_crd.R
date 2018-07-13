#' Completely Randomized Design
#'
#' This function creates the fieldbook and fieldplan for a CRD.
#' @param geno The list of genotypes.
#' @param nrep Number of replications.
#' @param nc Number of columns.
#' @author Raul Eyzaguirre.
#' @details The genotypes are randomly allocated on a field following a CRD.
#' @return It returns the fieldbook and fieldplan.
#' @examples
#' cr.crd(1:20, 3, 12)
#' cr.crd(1:20, 2, 7)
#' @export

cr.crd <- function(geno, nrep, nc) {
  
  # Error messages
  
  if (nrep < 2)
    stop("Include at least 2 replications.")

  ng <- length(geno)

  if (ng < 2)
    stop("Include at least 2 genotypes.")

  # Number of rows
  
  nr <- ceiling(ng * nrep / nc)

  # Fieldplan array
  
  plan <- array(dim = c(nr, nc))

  rownames(plan) <- paste("row", 1:nr)
  colnames(plan) <- paste("col", 1:nc)
  
  # Include genotypes at random
  
  geno <- sample(rep(geno, nrep))
  
  k <- 1
  
  for (i in 1:nr)
    for (j in 1:nc) {
      plan[i, j] <- geno[k]
      k <- k + 1
    }
  
  # Create fielbook
  
  row <- as.integer(gl(nr, nc))
  col <- rep(1:nc, nr)

  book <- data.frame(plot = 1:(nr * nc), row, col, geno = c(t(plan)),
                     stringsAsFactors = F)
  book <- book[!is.na(book$geno), ]
  
  # Return
  
  list(plan = plan, book = book)
  
}
