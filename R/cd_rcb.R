#' Randomized Complete Block Design
#'
#' This function creates the fieldbook and fieldplan for a RCBD.
#' @param geno The list of genotypes.
#' @param nb Number of blocks.
#' @param nc Number of columns.
#' @author Raul Eyzaguirre.
#' @details The genotypes are randomly allocated on a field following a RCBD.
#' The blocks are disposed alongside the rows.
#' @return It returns the fieldbook and fieldplan.
#' @examples
#' cd.rcb(1:20, 3, 12)
#' cd.rcb(1:20, 2, 7)
#' @export

cd.rcb <- function(geno, nb, nc) {
  
  # Dimensions
  
  ng <- length(geno)          # Number of genotypes
  nr <- ceiling(ng * nb / nc) # Number of rows
  
  # Fieldplan array
  
  plan <- array(dim = c(nr, nc))
  
  # Include genotypes at random
  
  sg <- sample(geno) # sorted genotypes
  for (i in 2:nb)
    sg <- c(sg, sample(geno))

  k <- 1
  
  for (i in 1:nr)
    for (j in 1:nc) {
      plan[i, j] <- sg[k]
      k <- k + 1
    }
  
  # Row and column names
  
  rownames(plan) <- paste("row", 1:nr)
  colnames(plan) <- paste("col", 1:nc)

  # Create fielbook

  row <- as.integer(gl(nr, nc))
  col <- rep(1:nc, nr)
  block <- as.integer(gl(nb, ng))
  length(block) <- nr * nc
  book <- data.frame(plot = 1:(nr * nc), row, col, block,
                     geno = c(t(plan)), stringsAsFactors = F)
  book <- book[!is.na(book$geno), ]

  # Return
  
  list(plan = plan, book = book)
}
