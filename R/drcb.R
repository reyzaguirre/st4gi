#' Randomized Complete Block Design
#'
#' This function creates the fieldbook and fieldplan for a RCBD.
#' @param geno The list of genotypes.
#' @param nb Number of blocks.
#' @param nr Number of rows.
#' @param byrow Logical. If \code{"TRUE"} (the default) the blocks are
#' allocated alongside the rows.
#' @author Raul Eyzaguirre.
#' @details The genotypes are randomly allocated on a field following a RCBD.
#' @return It returns the fieldbook and fieldplan.
#' @examples
#' drcb(1:20, 3, 12, TRUE)
#' drcb(1:20, 2, 7, FALSE)
#' @export

drcb <- function(geno, nb, nr, byrow = TRUE) {
  
  # Dimensions
  
  ng <- length(geno)      # Number of genotypes
  nc <- ceiling(ng * nb / nr)  # Number of columns
  
  # Fieldplan array
  
  plan <- array(dim = c(nr, nc))
  
  # Include genotypes at random
  
  sg <- sample(geno) # sorted genotypes
  for (i in 2:nb)
    sg <- c(sg, sample(geno))

  k <- 1
  
  if (byrow == TRUE) {
    for (i in 1:nr)
      for (j in 1:nc) {
        plan[i, j] <- sg[k]
        k <- k + 1
      }
  } else {
    for (i in 1:nc)
      for (j in 1:nr) {
        plan[j, i] <- sg[k]
        k <- k + 1
      }
  }
  
  # Row and column names
  
  rownames(plan) <- paste("row", 1:nr)
  colnames(plan) <- paste("col", 1:nc)

  # Create fielbook
  
  if (byrow == TRUE) {
    row <- as.integer(gl(nr, nc))
    col <- rep(1:nc, nr)
    genoplan <- t(plan)
  } else {
    row <- rep(1:nr, nc)
    col <- as.integer(gl(nc, nr))
    genoplan <- plan
  }

  block <- as.integer(gl(nb, ng))
  length(block) <- nr * nc

  book <- data.frame(plot = 1:(nr * nc),
                     row, col, block, geno = c(genoplan), stringsAsFactors = F)
  book <- book[!is.na(book$geno), ]

  # Return
  
  list(plan = plan, book = book)
}
