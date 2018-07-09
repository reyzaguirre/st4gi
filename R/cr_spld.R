#' Split-Plot Design
#'
#' This function creates the fieldbook and fieldplan for a Split-Plot design.
#' @param A The levels of factor A (row factor to whole plots).
#' @param B The levels of factor B (column factor to subplots).
#' @param nrep Number of replications (or blocks).
#' @param nc Number of columns in each replication. Default is the
#' number of levels of factor B.
#' @author Raul Eyzaguirre.
#' @details The levels of the factors are randomly allocated on a field
#' following a Split-Plot design. Row and column numbers are specific
#' to each replication and identify subplots. Each replication is a
#' complete block for factor A.
#' @return It returns the fieldbook and fieldplan.
#' @examples
#' A <- paste0("a", 1:4)
#' B <- paste0("b", 1:3)
#' cr.spld(A, B, 3)
#' @export

cr.spld <- function(A, B, nrep, nc = NULL) {
  
  # As character
  
  A <- as.character(A)
  B <- as.character(B)

  # Error messages
  
  nla <- length(A)
  nlb <- length(B)
  
  if (nrep < 2)
    stop("Include at least 2 replications.")

  if (nla < 2)
    stop("Include at least 2 levels for factor A.")

  if (nlb < 2)
    stop("Include at least 2 levels for factor B.")
  
  # Number of rows for each plot
  
  if (is.null(nc))
    nc <- nlb

  nr <- ceiling(nlb / nc)

  # Fieldplan array
  
  plan <- array(dim = c(nla * nr, nc, nrep))

  rownames(plan) <- paste("row", 1:(nla * nr))
  colnames(plan) <- paste("col", 1:nc)
  dimnames(plan)[[3]] <- paste("rep", 1:nrep)
  
  # Include treatments at random

  rana <- array(dim = c(nla, nrep))
  ranb <- array(dim = c(nlb, nla, nrep))
  
  for (i in 1:nrep) {
    rana[, i] <- sample(1:nla)
    for (j in 1:nla) {
      ranb[, j, i] <- sample(1:nlb)
      stab <- paste(A[rana[j, i]], B[ranb[, j, i]], sep = "_")
      k <- 1
      for (l in ((j - 1) * nr + 1):(j * nr))
        for (m in 1:nc) {
          plan[l, m, i] <- stab[k]
          k <- k + 1
        }
    }
  }
  
  # Create fielbook
  
  plot <- as.integer(gl(nla * nrep, nr * nc))
  block <- as.integer(gl(nrep, nla * nr * nc))
  row <- rep(as.integer(gl(nla * nr, nc)), nrep)
  col <- rep(rep(1:nc, nla * nr), nrep)
  
  sta <- NULL
  stb <- NULL
  stab <- NULL
  
  for (i in 1:nrep) {
    sta <- c(sta, c(sapply(A[rana[, i]], rep, nr * nc)))
    for (j in 1:nla) {
      temp <- B[c(ranb[, j, i])]
      length(temp) <- nr * nc
      stb <- c(stb, temp)
    }
    stab <- c(stab, c(t(plan[, , i])))
  }
  
  book <- data.frame(plot, block, row, col, A = sta, B = stb,
                     treat = stab, stringsAsFactors = F)
  book <- book[!is.na(book$treat), ]
  book$subplot <- 1:dim(book)[1]
  book <- book[, c(1, 8, 2, 3, 4, 5, 6, 7)]
  rownames(book) <- 1:dim(book)[1]
  
  # Return
  
  list(plan = plan, book = book)
}
