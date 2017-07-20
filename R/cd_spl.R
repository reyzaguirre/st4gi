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
#' A <- paste("a", 1:4, sep = "")
#' B <- paste("b", 1:3, sep = "")
#' cd.spl(A, B, 3)
#' @export

cd.spl <- function(A, B, nrep, nc = NULL) {
  
  # As character
  
  A <- as.character(A)
  B <- as.character(B)

  # Error messages
  
  if (nrep < 2)
    stop("Include at least 2 replications.")

  na <- length(A)
  if (na < 2)
    stop("Include at least 2 levels for factor A.")

  nb <- length(B)
  if (nb < 2)
    stop("Include at least 2 levels for factor B.")
  
  # Number of rows for each plot
  
  if (is.null(nc))
    nc <- nb

  nr <- ceiling(nb / nc)

  # Fieldplan array
  
  plan <- array(dim = c(na * nr, nc, nrep))

  rownames(plan) <- paste("row", 1:(na * nr))
  colnames(plan) <- paste("col", 1:nc)
  dimnames(plan)[[3]] <- paste("rep", 1:nrep)
  
  # Include treatments at random

  for (i in 1:nrep) {
    ta <- sample(A)
    for (j in 1:na) {
      tb <- sample(B)
      tab <- paste(ta[j], tb, sep = "_")
      k <- 1
      for (l in ((j - 1) * nr + 1):(j * nr))
        for (m in 1:nc) {
          plan[l, m, i] <- tab[k]
          k <- k + 1
        }
    }
  }
  
  # Create fielbook
  
  plot <- as.integer(gl(na * nrep, nr * nc))
  block <- as.integer(gl(nrep, na * nr * nc))
  row <- rep(as.integer(gl(na * nr, nc)), nrep)
  col <- rep(rep(1:nc, na * nr), nrep)
  
  tab <- c(t(plan[, , 1]))
  
  for (i in 2:nrep)
    tab <- c(tab, c(t(plan[, , i])))
  
  for (i in 1:length(tab)) {
    ta[i] <- unlist(strsplit(tab[i], "_"))[1]
    tb[i] <- unlist(strsplit(tab[i], "_"))[2]
  }
  
  book <- data.frame(plot, block, row, col, A = ta, B = tb,
                     treat = tab, stringsAsFactors = F)
  book <- book[!is.na(book$treat), ]
  book$subplot <- 1:dim(book)[1]
  book <- book[, c(1, 8, 2, 3, 4, 5, 6, 7)]
  
  # Return
  
  list(plan = plan, book = book)
}
