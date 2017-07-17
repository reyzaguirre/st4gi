#' Augmented Block Design
#'
#' This function creates the fieldbook and fieldplan for an ABD
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
#' cd.ab(1:40, checks, 4, 10)
#' cd.ab(1:40, checks, 4, 7)
#' @export

cd.ab <- function(geno, checks, nb, nc) {
  
  # Error messages
    
  if (nb < 2)
    stop("Include at least 2 blocks.")
    
  nch <- length(checks)
  if (nch < 2)
    stop("Include at least 2 checks.")
  
  ng <- length(geno)
  if (ng < nb)
    stop(paste("Include at least", nb, "genotypes."))
    
  # Number of rows
  
  nr <- ceiling((ng + nch * nb) / nc)
  
  # Fieldplan array
  
  plan <- array(dim = c(nr, nc))
  
  # Create blocks and add genotypes at random
  
  blocks <- list()
  
  for (i in 1:nb)
    blocks[[i]] <- checks
  
  sg <- sample(geno)
  
  for (i in 1:ng) {
    j <- i%%4
    if (j == 0) j <- 4
    blocks[[j]] <- c(blocks[[j]], sg[i])
  }
  
  blocks <- lapply(blocks, sample)
  
  # Add genotypes and checks to the fieldplan

  sgc <-unlist(blocks)
  
  k <- 1
  
  for (i in 1:nr)
    for (j in 1:nc) {
      plan[i, j] <- sgc[k]
      k <- k + 1
    }
  
  # Row and column names
  
  rownames(plan) <- paste("row", 1:nr)
  colnames(plan) <- paste("col", 1:nc)

  # Create fielbook

  row <- as.integer(gl(nr, nc))
  col <- rep(1:nc, nr)
  block <- rep(1:nb, sapply(blocks, length))
  length(block) <- nr * nc
  book <- data.frame(plot = 1:(nr * nc), row, col, block,
                     geno = c(t(plan)), stringsAsFactors = F)
  book <- book[!is.na(book$geno), ]

  # Return
  
  list(plan = plan, book = book)
}
