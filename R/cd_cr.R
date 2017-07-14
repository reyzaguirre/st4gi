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
#' dcr(1:20, 3, 12)
#' dcr(1:20, 2, 7)
#' @export

dcr <- function(geno, nrep, nc) {
  
  # Dimensions
  
  ng <- length(geno)          # Number of genotypes
  nr <- ceiling(ng * nrep / nc) # Number of rows
  
  # Fieldplan array
  
  plan <- array(dim = c(nr, nc))
  
  # Include genotypes at random
  
  geno <- sample(rep(geno, nrep))
  
  k <- 1
  
  for (i in 1:nr)
    for (j in 1:nc) {
      plan[i, j] <- geno[k]
      k <- k + 1
    }
  
  # Row and column names
  
  rownames(plan) <- paste("row", 1:nr)
  colnames(plan) <- paste("col", 1:nc)
  
  # Create fielbook
  
  row <- as.integer(gl(nr, nc))
  col <- rep(1:nc, nr)
  book <- data.frame(plot = 1:(nr * nc), row, col,
                     geno = c(t(plan)), stringsAsFactors = F)
  book <- book[!is.na(book$geno), ]
  
  # Return
  
  list(plan = plan, book = book)
}
