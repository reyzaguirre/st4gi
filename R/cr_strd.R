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
#' A <- paste0("a", 1:4)
#' B <- paste0("b", 1:3)
#' cr.strd(A, B, 3)
#' @export

cr.strd <- function(A, B, nrep) {
  
  # Error messages

  nla <- length(A)
  nlb <- length(B)
  
  if (nrep < 2)
    stop("Include at least 2 replications.")

  if (nla < 2)
    stop("Include at least 2 levels for factor A.")

  if (nlb < 2)
    stop("Include at least 2 levels for factor B.")

  # Fieldplan array
  
  plan <- array(dim = c(nla, nlb, nrep))

  rownames(plan) <- paste("row", 1:nla)
  colnames(plan) <- paste("col", 1:nlb)
  dimnames(plan)[[3]] <- paste("rep", 1:nrep)
  
  # Random order for A and B levels
  
  rana <- array(dim = c(nla, nrep))
  ranb <- array(dim = c(nlb, nrep))
  
  for (i in 1:nrep) {
    rana[, i] <- sample(1:nla)
    ranb[, i] <- sample(1:nlb)
    plan[, , i] <- outer(A[rana[, i]], B[ranb[, i]], paste, sep = "_")
  }
   
  # Create fielbook
  
  row <- rep(as.integer(gl(nla, nlb)), nrep)
  col <- rep(rep(1:nlb, nla), nrep)
  block <- as.integer(gl(nrep, nla * nlb))
  
  sta <- NULL
  stb <- NULL
  stab <- NULL

  for (i in 1:nrep) {
    sta <- c(sta, c(sapply(A[rana[, i]], rep, nlb)))
    stb <- c(stb, rep(B[ranb[, i]], nla))
    stab <- c(stab, c(t(plan[, , i])))
  }
  
  book <- data.frame(plot = 1:(nla * nlb * nrep), block, row, col,
                     A = sta, B = stb, treat = stab, stringsAsFactors = F)

  # Return
  
  list(plan = plan, book = book)
  
}
