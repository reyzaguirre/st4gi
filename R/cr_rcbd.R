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
#' cr.rcbd(1:20, 3, 10)
#' cr.rcbd(1:20, 2, 7)
#' @export

cr.rcbd <- function(geno, nb, nc) {
  
  # As character
  
  geno <- as.character(geno)

  # Error messages
  
  if (nb < 2)
    stop("Include at least 2 blocks.")
  
  ng <- length(geno)
  
  if (ng < 2)
    stop("Include at least 2 genotypes.")
  
  # Number of rows for each block
  
  nr <- ceiling(ng / nc)
  
  # Fieldplan array
  
  plan <- array(dim = c(nr, nc, nb))
  
  rownames(plan) <- paste("row", 1:nr)
  colnames(plan) <- paste("col", 1:nc)
  dimnames(plan)[[3]] <- paste("block", 1:nb)

  # Include genotypes at random

  for (i in 1:nb) {
    sg <- sample(geno)
    k <- 1
    for (j in 1:nr)
      for (l in 1:nc) {
        plan[j, l, i] <- sg[k]
        k <- k + 1
      }
  }
  
  # Create fielbook

  block <- as.integer(gl(nb, nr * nc))
  row <- rep(as.integer(gl(nr, nc)), nb)
  col <- rep(rep(1:nc, nr), nb)
  
  geno <- NULL

  for (i in 1:nb)
    geno <- c(geno, c(t(plan[, , i])))
  
  book <- data.frame(block, row, col, geno, stringsAsFactors = F)
  book <- book[!is.na(book$geno), ]
  book$plot <- 1:dim(book)[1]
  book <- book[, c(5, 1, 2, 3, 4)]
  rownames(book) <- 1:dim(book)[1]

  # Return
  
  list(plan = plan, book = book)
}
