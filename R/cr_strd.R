#' Strip-Split-Plot Design
#'
#' This function creates the fieldbook and fieldplan for a Strip-Split-Plot design.
#' @param A The levels of factor A (row factor).
#' @param B The levels of factor B (column factor).
#' @param nrep Number of replications (or blocks).
#' @param serpentine \code{"yes"} or \code{"no"}, default \code{"yes"}.
#' @return It returns the fieldbook and fieldplan.
#' @author Raul Eyzaguirre.
#' @examples
#' A <- paste0("a", 1:4)
#' B <- paste0("b", 1:3)
#' cr.strd(A, B, 3)
#' @export

cr.strd <- function(A, B, nrep, serpentine = c("yes", "no")) {
  
  # Match arguments
  
  serpentine <- match.arg(serpentine)

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
  
  plan.id <- t(array(1:(nlb*nla), dim = c(nlb, nla)))
  
  if (serpentine == 'yes')
    for (i in seq(2, nla, 2))
      plan.id[i, ] <- sort(plan.id[i, ], decreasing = TRUE)
  
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
  plot.num <- NULL

  for (i in 1:nrep) {
    sta <- c(sta, c(sapply(A[rana[, i]], rep, nlb)))
    stb <- c(stb, rep(B[ranb[, i]], nla))
    stab <- c(stab, c(t(plan[, , i])))
    plot.num <- c(plot.num, c(t(plan.id)) + nla * nlb * (i - 1))
  }
  
  book <- data.frame(plot.num, block, row, col,
                     A = sta, B = stb, treat = stab, stringsAsFactors = F)
  
  # Sort by plot number
  
  if (serpentine == 'yes')
    book <- book[sort(book$plot.num, index.return = TRUE)$ix, ]
  
  rownames(book) <- 1:dim(book)[1]
  
  # Return
  
  list(plan = plan, book = book)
  
}
