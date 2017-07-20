#' Completely Randomized Design for a two-factor factorial
#'
#' This function creates the fieldbook and fieldplan for a two-factor
#' factorial with a CRD.
#' @param A The levels of factor A.
#' @param B The levels of factor B.
#' @param nrep Number of replications.
#' @param nc Number of columns.
#' @author Raul Eyzaguirre.
#' @details The treatments are randomly allocated on a field following a CRD.
#' @return It returns the fieldbook and fieldplan.
#' @examples
#' A <- paste("a", 1:5, sep = "")
#' B <- paste("b", 1:3, sep = "")
#' cd.2fcr(A, B, 3, 12)
#' cd.2fcr(A, B, 2, 7)
#' @export

cd.2fcr <- function(A, B, nrep, nc) {
  
  # Error messages
  
  if (nrep < 2)
    stop("Include at least 2 replications.")

  nla <- length(A)
  if (nla < 2)
    stop("Include at least 2 levels for factor A.")

  nlb <- length(B)
  if (nlb < 2)
    stop("Include at least 2 levels for factor B.")

  # Number of rows
  
  nt <- nla * nlb
  nr <- ceiling(nt * nrep / nc)

  # Fieldplan array
  
  plan <- array(dim = c(nr, nc))

  rownames(plan) <- paste("row", 1:nr)
  colnames(plan) <- paste("col", 1:nc)
  
  # Include treatments at random
  
  ta <- as.integer(gl(nla, nlb))
  ta <- rep(A[ta], nrep)
  tb <- rep(rep(B, nla), nrep)
  tab <- paste(ta, tb, sep = "_")
  
  ran <- sample(1:(nt * nrep))
  
  ta <- ta[ran]
  tb <- tb[ran]
  tab <- tab[ran]
  
  k <- 1
  
  for (i in 1:nr)
    for (j in 1:nc) {
      plan[i, j] <- tab[k]
      k <- k + 1
    }
  
  # Create fielbook
  
  row <- as.integer(gl(nr, nc))
  col <- rep(1:nc, nr)
  
  length(ta) <- nr * nc
  length(tb) <- nr * nc
  
  book <- data.frame(plot = 1:(nr * nc), row, col, A = ta, B = tb,
                     treat = c(t(plan)), stringsAsFactors = F)
  book <- book[!is.na(book$treat), ]
  
  # Return
  
  list(plan = plan, book = book)
}
