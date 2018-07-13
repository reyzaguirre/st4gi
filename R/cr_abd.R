#' Augmented Block Design
#'
#' This function creates the fieldbook and fieldplan for an ABD.
#' @param geno The list of genotypes.
#' @param checks The list of checks.
#' @param nb Number of blocks.
#' @param nc Number of columns.
#' @author Raul Eyzaguirre.
#' @details The genotypes are randomly allocated on a field following an ABD.
#' The blocks are disposed alongside the rows.
#' @return It returns the fieldbook and fieldplan.
#' @examples
#' checks <- paste("ch", 1:4, sep = "_")
#' cr.abd(1:40, checks, 4, 10)
#' cr.abd(1:50, checks, 3, 7)
#' @export

cr.abd <- function(geno, checks, nb, nc) {
  
  # Error messages

  if (nb < 2)
    stop("Include at least 2 blocks.")
    
  nch <- length(checks)
  
  if (nch < 2)
    stop("Include at least 2 checks.")
  
  ng <- length(geno)
  
  if (ng < nb)
    stop(paste("Include at least", nb, "genotypes."))

  # Create blocks with checks
  
  blocks <- list()
  
  for (i in 1:nb)
    blocks[[i]] <- checks
  
  # Add genotypes at random to the blocks
  
  sg <- sample(geno)

  for (i in 1:ng) {
    j <- ((i - 1) %% nb) + 1
    blocks[[j]] <- c(blocks[[j]], sg[i])
  }
  
  # Sort full blocks
  
  blocks <- lapply(blocks, sample)
  
  # Number of rows for each block
  
  nr <- ceiling(max(unlist(lapply(blocks, length))) / nc)
  
  # Fieldplan array
  
  plan <- array(dim = c(nr, nc, nb))
  
  rownames(plan) <- paste("row", 1:nr)
  colnames(plan) <- paste("col", 1:nc)
  dimnames(plan)[[3]] <- paste("block", 1:nb)
  
  # Add genotypes and checks to the fieldplan

  for (i in 1:nb) {
    k <- 1
    for (j in 1:nr)
      for (l in 1:nc) {
        plan[j, l, i] <- blocks[[i]][k]
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
