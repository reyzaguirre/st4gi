#' Randomized Complete Block Design for a two-factor factorial
#'
#' This function creates the fieldbook and fieldplan for a two-factor
#' factorial with a RCBD.
#' @param A The levels of factor A.
#' @param B The levels of factor B.
#' @param nb Number of blocks.
#' @param nc Number of columns.
#' @author Raul Eyzaguirre.
#' @details The treatments are randomly allocated on a field following a RCBD.
#' The blocks are disposed alongside the rows.
#' @return It returns the fieldbook and fieldplan.
#' @examples
#' A <- paste("a", 1:5, sep = "")
#' B <- paste("b", 1:3, sep = "")
#' cd.2frcb(A, B, 3, 9)
#' cd.2frcb(A, B, 2, 7)
#' @export

cd.2frcb <- function(A, B, nb, nc) {
  
  # Error messages
  
  if (nb < 2)
    stop("Include at least 2 blocks.")
  
  na <- length(A)
  if (na < 2)
    stop("Include at least 2 levels for factor A.")
  
  nb <- length(B)
  if (nb < 2)
    stop("Include at least 2 levels for factor B.")
  
  # Number of rows
  
  nt <- length(A) * length(B)
  nr <- ceiling(nt * nb / nc)

  # Fieldplan array
  
  plan <- array(dim = c(nr, nc))

  rownames(plan) <- paste("row", 1:nr)
  colnames(plan) <- paste("col", 1:nc)
  
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
