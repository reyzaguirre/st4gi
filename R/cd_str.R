#' Strip-Split-Plot Design
#'
#' This function creates the fieldbook and fieldplan for a Strip-Split-Plot design.
#' @param A The levels of factor A (row factor).
#' @param B The levels of factor B (column factor).
#' @param nrep Number of replications (or blocks).
#' @author Raul Eyzaguirre.
#' @details The levels of the factors are randomly allocated on a field
#' following a Strip-Split-Plot design. Row and column numbers are specific
#' to each replication. Each replication is a complete block for factor A
#' and for factor B.
#' @return It returns the fieldbook and fieldplan.
#' @examples
#' A <- paste("a", 1:4, sep = "")
#' B <- paste("b", 1:3, sep = "")
#' cd.str(A, B, 3)
#' @export

cd.str <- function(A, B, nrep) {
  
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

  # Fieldplan array
  
  plan <- array(dim = c(na, nb, nrep))

  rownames(plan) <- paste("row", 1:na)
  colnames(plan) <- paste("col", 1:nb)
  dimnames(plan)[[3]] <- paste("rep", 1:nrep)
  
  # Random order for A and B levels
  
  for (i in 1:nrep) {
    ta <- sample(A)
    tb <- sample(B)
    plan[, , i] <- outer(ta, tb, paste, sep = "_")
  }
   
  # Create fielbook
  
  row <- rep(as.integer(gl(na, nb)), nrep)
  col <- rep(rep(1:nb, na), nrep)
  block <- as.integer(gl(nrep, na * nb))
  
  tab <- c(t(plan[, , 1]))

  for (i in 2:nrep)
    tab <- c(tab, c(t(plan[, , i])))
  
  for (i in 1:length(tab)) {
    ta[i] <- unlist(strsplit(tab[i], "_"))[1]
    tb[i] <- unlist(strsplit(tab[i], "_"))[2]
  }
  
  book <- data.frame(plot = 1:(na * nb * nrep), block, row, col,
                     A = ta, B = tb, treat = tab, stringsAsFactors = F)

  # Return
  
  list(plan = plan, book = book)
}
