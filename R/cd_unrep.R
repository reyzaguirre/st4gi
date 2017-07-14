#' Unreplicated design
#'
#' This function creates the fieldbook and fieldplan for an unreplicated design
#' with genotypes randomly allocated on a field.
#' @param geno The list of genotypes.
#' @param nc Number of columns.
#' @author Raul Eyzaguirre.
#' @details The genotypes are randomly allocated on a field.
#' @return It returns the fieldbook and fieldplan.
#' @examples
#' dunrep(1:100, 5)
#' dunrep(1:100, 7)
#' @export

dunrep <- function(geno, nc) {
  
  # Dimensions
  
  ng <- length(geno)      # Number of genotypes
  nr <- ceiling(ng / nc)  # Number of rows
  
  # Fieldplan array
  
  plan <- array(dim = c(nr, nc))

  # Include genotypes at random
  
  geno <- sample(geno)

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
  
  book <- data.frame(plot = 1:(nr * nc), row, col, geno = c(t(plan)),
                     stringsAsFactors = F)
  book <- book[!is.na(book$geno), ]

  # Return
  
  list(plan = plan, book = book)
}
