#' Split-Plot Design
#'
#' This function creates the fieldbook and fieldplan for a Split-Plot design.
#' @param A The levels of factor A (row factor to whole plots).
#' @param B The levels of factor B (column factor to subplots).
#' @param nrep Number of replications (or blocks).
#' @param nc Number of columns in each replication. Default is the
#' number of levels of factor B.
#' @param serpentine \code{"yes"} or \code{"no"}, default \code{"yes"}.
#' @return It returns the fieldbook and fieldplan.
#' @author Raul Eyzaguirre.
#' @examples
#' A <- paste0("a", 1:4)
#' B <- paste0("b", 1:3)
#' cr.spld(A, B, 3)
#' @export

cr.spld <- function(A, B, nrep, nc = NULL, serpentine = c("yes", "no")) {
  
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
  
  # Number of rows for each plot
  
  if (is.null(nc))
    nc <- nlb

  nr <- ceiling(nlb / nc)

  # Fieldplan array
  
  plan.id <- fp(nr, nc, serpentine)
  
  plan <- array(dim = c(nla * nr, nc, nrep))

  rownames(plan) <- paste("row", 1:(nla * nr))
  colnames(plan) <- paste("col", 1:nc)
  dimnames(plan)[[3]] <- paste("rep", 1:nrep)
  
  # Include treatments at random

  for (k in 1:nrep) {
    rana <- sample(A)
    for (l in 1:nla) {
      ranb <- sample(B)
      stab <- paste(rana[l], ranb, sep = "_")
      for (i in 1:nr)
        for (j in 1:nc)
          plan[i + (l - 1) * nr, j, k] <- stab[plan.id[i, j]]
    }
  }
  
  # Create fielbook
  
  row <- rep(as.integer(gl(nla * nr, nc)), nrep)
  col <- rep(rep(1:nc, nla * nr), nrep)
  block <- as.integer(gl(nrep, nla * nr * nc))
  
  treat <- NULL
  subplot.num <- NULL
  
  for (k in 1:nrep) {
    treat <- c(treat, c(t(plan[, , k])))
    for (l in 1:nla)
      subplot.num <- c(subplot.num, c(t(plan.id)) + (k - 1) * nlb * nla + nlb * (l - 1))
  }
  
  book <- data.frame(subplot.num, block, row, col, treat, stringsAsFactors = FALSE)
  book <- book[!is.na(book$treat), ]
  
  # Falta generar niveles de A y B
  
  book$A <- c(data.frame(sapply(book$treat, strsplit, "_"),
                         stringsAsFactors = FALSE)[1, ])
  book$B <- c(data.frame(sapply(book$treat, strsplit, "_"),
                         stringsAsFactors = FALSE)[2, ])

  # Reorder columns
  
  book <- book[, c(1, 2, 3, 4, 6, 7, 5)]

  # Sort by plot number
  
  if (serpentine == 'yes' & nr > 1)
    book <- book[sort(book$subplot.num, index.return = TRUE)$ix, ]
  
  rownames(book) <- 1:dim(book)[1]
  
  # Return
  
  list(plan = plan, book = book)
  
}
