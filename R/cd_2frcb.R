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
  
  # As character
  
  A <- as.character(A)
  B <- as.character(B)

  # Error messages
  
  nla <- length(A)
  nlb <- length(B)

  if (nb < 2)
    stop("Include at least 2 blocks.")
  
  if (nla < 2)
    stop("Include at least 2 levels for factor A.")
  
  if (nlb < 2)
    stop("Include at least 2 levels for factor B.")
  
  # Number of rows for each block
  
  nt <- nla * nlb
  nr <- ceiling(nt / nc)
  
  # Fieldplan array
  
  plan <- array(dim = c(nr, nc, nb))
  
  rownames(plan) <- paste("row", 1:nr)
  colnames(plan) <- paste("col", 1:nc)
  dimnames(plan)[[3]] <- paste("block", 1:nb)
  
  # Include treatments at random
  
  ta <- as.integer(gl(nla, nlb))
  ta <- A[ta]
  tb <- rep(B, nla)
  tab <- paste(ta, tb, sep = "_")
  
  ran <- array(dim = c(nr * nc, nb))

  for (i in 1:nb) {
    st <- sample(1:nt)
    length(st) <- nr * nc
    k <- 1
    for (j in 1:nr)
      for (l in 1:nc) {
        plan[j, l, i] <- tab[st[k]]
        k <- k + 1
      }
    ran[, i] <- st
  }
  
  # Create fielbook

  block <- as.integer(gl(nb, nr * nc))
  row <- rep(as.integer(gl(nr, nc)), nb)
  col <- rep(rep(1:nc, nr), nb)
  
  sta <- NULL
  stb <- NULL
  stab <- NULL
  
  for (i in 1:nb) {
    sta <- c(sta, ta[ran[, i]])
    stb <- c(stb, tb[ran[, i]])
    stab <- c(stab, tab[ran[, i]])
  }

  book <- data.frame(block, row, col, A = sta, B = stb, treat = stab,
                     stringsAsFactors = F)
  book <- book[!is.na(book$treat), ]
  book$plot <- 1:dim(book)[1]
  book <- book[, c(7, 1, 2, 3, 4, 5, 6)]
  rownames(book) <- 1:dim(book)[1]

  # Return
  
  list(plan = plan, book = book)
}
