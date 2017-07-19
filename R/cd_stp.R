#' Strip-Split-Plot Design
#'
#' This function creates the fieldbook and fieldplan for a Strip-Split-Plot design.
#' @param A The levels of factor A.
#' @param B The levels of factor B.
#' @param nrep Number of replications (blocks).
#' @author Raul Eyzaguirre.
#' @details The levels of the factors are randomly allocated on a field
#' following a Strip-Split-Plot design. Row and column numbers are specific
#' to each replication.
#' @return It returns the fieldbook and fieldplan.
#' @examples
#' A <- paste("a", 1:4, sep = "")
#' B <- paste("b", 1:3, sep = "")
#' cd.str(A, B, 3)
#' @export

cd.str <- function(A, B, nrep) {
  
  # Error messages
  
  if (nrep < 2)
    stop("Include at least 2 replications.")

  na <- length(A)
  if (na < 2)
    stop("Include at least 2 levels for factor A.")

  nb <- length(B)
  if (nb < 2)
    stop("Include at least 2 levels for factor B.")

  # Fieldplan array
  
  plan <- array(dim = c(nb, na, nrep))
  
  # Include treatments at random
  
  ta <- as.integer(gl(length(A), length(B)))
  ta <- A[ta]
  tb <- rep(B, length(A))
  tab <- paste(ta, tb, sep = "_")
  
  ran <- sample(1:nt)
  
  for (i in 2:nb)
    ran <- c(ran, sample(1:nt))
  
  ta <- ta[ran]
  tb <- tb[ran]
  tab <- tab[ran]
  
  k <- 1
  
  for (i in 1:nr)
    for (j in 1:nc) {
      plan[i, j] <- tab[k]
      k <- k + 1
    }
  
  # Row and column names
  
  rownames(plan) <- paste("row", 1:nr)
  colnames(plan) <- paste("col", 1:nc)
  
  # Create fielbook
  
  row <- as.integer(gl(nr, nc))
  col <- rep(1:nc, nr)
  block <- as.integer(gl(nb, nt))
  
  length(ta) <- nr * nc
  length(tb) <- nr * nc
  length(block) <- nr * nc
  
  book <- data.frame(plot = 1:(nr * nc), row, col, block, A = ta, B = tb,
                     treat = c(t(plan)), stringsAsFactors = F)
  book <- book[!is.na(book$treat), ]
  
  # Return
  
  list(plan = plan, book = book)
}
