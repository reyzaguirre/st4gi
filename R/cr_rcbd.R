#' Randomized Complete Block Design
#'
#' This function creates the fieldbook and fieldplan for a RCBD.
#' @param geno The list of genotypes.
#' @param nb Number of blocks.
#' @param nc Number of available columns on the field for each block.
#' @param serpentine \code{"yes"} or \code{"no"}, default \code{"yes"}.
#' @param alongside \code{"no"} for independent blocks, or \code{"rows"}
#' or \code{"columns"} if blocks are together alongside rows or columns.
#' @return It returns the fieldbook and fieldplan. 
#' @author Raul Eyzaguirre.
#' @examples
#' cr.rcbd(1:20, 3, 10)
#' cr.rcbd(1:20, 2, 7)
#' @export

cr.rcbd <- function(geno, nb, nc = NULL, serpentine = c("yes", "no"),
                    alongside = c("no", "rows", "columns")) {
  
  # Match arguments
  
  serpentine <- match.arg(serpentine)
  alongside <- match.arg(alongside)

  # Error messages
  
  if (nb < 2)
    stop("Include at least 2 blocks.")
  
  ng <- length(geno)
  
  if (ng < 2)
    stop("Include at least 2 genotypes.")
  
  # Number of rows and columns for each block
  
  if (is.null(nc))
    nc <- gnc(ng)

  nr <- ceiling(ng / nc)
  
  # Fieldplan array for each block
  
  plan.id <- fp(nr, nc, serpentine)
  
  # Create fieldplan

  plan <- array(dim = c(nr, nc, nb))
  
  rownames(plan) <- paste("row", 1:nr)
  colnames(plan) <- paste("col", 1:nc)
  dimnames(plan)[[3]] <- paste("block", 1:nb)

  # Include genotypes at random

  for (k in 1:nb) {
    sg <- sample(geno)
    plan[, , k] <- array(sg[plan.id], c(nr, nc))
  }
  
  # Create fielbook

  block <- as.integer(gl(nb, nr * nc))
  row <- rep(as.integer(gl(nr, nc)), nb)
  col <- rep(rep(1:nc, nr), nb)
  
  geno <- NULL
  plot.num <- NULL

  for (i in 1:nb) {
    geno <- c(geno, c(t(plan[, , i])))
    plot.num <- c(plot.num, c(t(plan.id)) + ng * (i - 1))
  }
  
  book <- data.frame(plot.num, block, row, col, geno, stringsAsFactors = FALSE)
  book <- book[!is.na(book$geno), ]
  
  # Sort by plot number
  
  if (serpentine == 'yes' & nr > 1)
    book <- book[sort(book$plot.num, index.return = TRUE)$ix, ]

  rownames(book) <- 1:dim(book)[1]

  # Change row and column numbers if required
  
  if (alongside == "rows") {
    plan <- t(apply(plan, 1, rbind))
    colnames(plan) <- paste("col", 1:dim(plan)[2])
    book$col <- book$col + (book$block - 1) * nc
  }
  
  if (alongside == "columns") {
    plan <- apply(plan, 2, rbind)
    rownames(plan) <- paste("row", 1:dim(plan)[1])
    book$row <- book$row + (book$block - 1) * nr
  }  

  # Return
  
  list(plan = plan, book = book)
  
}
