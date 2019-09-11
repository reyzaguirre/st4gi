#' Strip-Split-Plot Design
#'
#' This function creates the fieldbook and fieldplan for a Strip-Split-Plot design.
#' @param A The levels of factor A (row factor).
#' @param B The levels of factor B (column factor).
#' @param nb Number of blocks.
#' @param serpentine \code{"yes"} or \code{"no"}, default \code{"yes"}.
#' @param alongside \code{"no"} for independent blocks, or \code{"rows"}
#' or \code{"columns"} if blocks are together alongside rows or columns.
#' @return It returns the fieldbook and fieldplan.
#' @author Raul Eyzaguirre.
#' @examples
#' A <- paste0("a", 1:4)
#' B <- paste0("b", 1:3)
#' cr.strd(A, B, 3)
#' @export

cr.strd <- function(A, B, nb, serpentine = c("yes", "no"),
                    alongside = c("no", "rows", "columns")) {
  
  # Match arguments
  
  serpentine <- match.arg(serpentine)
  alongside <- match.arg(alongside)
  
  # Error messages

  nla <- length(A)
  nlb <- length(B)
  
  if (nb < 2)
    stop("Include at least 2 blocks.")

  if (nla < 2)
    stop("Include at least 2 levels for factor A.")

  if (nlb < 2)
    stop("Include at least 2 levels for factor B.")

  # Fieldplan array
  
  plan.id <- fp(nla, nlb, serpentine)

  plan <- array(dim = c(nla, nlb, nb))

  rownames(plan) <- paste("row", 1:nla)
  colnames(plan) <- paste("col", 1:nlb)
  dimnames(plan)[[3]] <- paste("block", 1:nb)
  
  # Random order for A and B levels
  
  rana <- array(dim = c(nla, nb))
  ranb <- array(dim = c(nlb, nb))
  
  for (i in 1:nb) {
    rana[, i] <- sample(1:nla)
    ranb[, i] <- sample(1:nlb)
    plan[, , i] <- outer(A[rana[, i]], B[ranb[, i]], paste, sep = ":-p")
  }
   
  # Create fielbook
  
  row <- rep(as.integer(gl(nla, nlb)), nb)
  col <- rep(rep(1:nlb, nla), nb)
  block <- as.integer(gl(nb, nla * nlb))
  
  sta <- NULL
  stb <- NULL
  stab <- NULL
  plot.num <- NULL

  for (i in 1:nb) {
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
  
  # Replace characters for treatment names
  
  plan <- gsub(":-p", "_", plan)
  book$treat <- gsub(":-p", "_", book$treat)

  # Change row and column numbers if required
  
  if (alongside == "rows") {
    plan <- t(apply(plan, 1, rbind))
    colnames(plan) <- paste("col", 1:dim(plan)[2])
    book$col <- book$col + (book$block - 1) * nlb
  }
  
  if (alongside == "columns") {
    plan <- apply(plan, 2, rbind)
    rownames(plan) <- paste("row", 1:dim(plan)[1])
    book$row <- book$row + (book$block - 1) * nla
  }
  
  # Return
  
  list(plan = plan, book = book)
  
}
