#' Augmented Block Design
#'
#' This function creates the fieldbook and fieldplan for an ABD.
#' @param geno The list of genotypes.
#' @param checks The list of checks.
#' @param nb Number of blocks.
#' @param nc Number of available columns on the field.
#' @param serpentine \code{"yes"} or \code{"no"}, default \code{"yes"}.
#' @param alongside \code{"no"} for independent blocks, or \code{"rows"}
#' or \code{"columns"} if blocks are together alongside rows or columns.
#' @return It returns the fieldbook and fieldplan.
#' @author Raul Eyzaguirre.
#' @examples
#' checks <- paste("ch", 1:4, sep = "_")
#' cr.abd(1:40, checks, 4, 10)
#' cr.abd(1:50, checks, 3, 7)
#' @export

cr.abd <- function(geno, checks, nb, nc = NULL, serpentine = c("yes", "no"),
                   alongside = c("no", "rows", "columns")) {
  
  # Match arguments
  
  serpentine <- match.arg(serpentine)
  alongside <- match.arg(alongside)
  
  # Error messages

  if (nb < 2)
    stop("Include at least 2 blocks.")
    
  nch <- length(checks)
  
  if (nch < 2)
    stop("Include at least 2 checks.")
  
  ng <- length(geno)
  
  if (ng < nb)
    stop(paste("Include at least", nb, "genotypes."))

  # Create blocks with checks
  
  blocks <- list()
  
  for (i in 1:nb)
    blocks[[i]] <- checks
  
  # Add genotypes at random to the blocks
  
  sg <- sample(geno)

  for (i in 1:ng) {
    j <- ((i - 1) %% nb) + 1
    blocks[[j]] <- c(blocks[[j]], sg[i])
  }
  
  # Sort full blocks
  
  blocks <- lapply(blocks, sample)
  
  # Get maximum block size

  bs <- max(unlist(lapply(blocks, length)))
  
  # Number of rows and columns for each block

  if (is.null(nc))
    nc <-  gnc(bs)
  
  nr <- ceiling(bs / nc)
  
  # Fieldplan array
  
  plan.id <- fp(nr, nc, serpentine)

  plan <- array(dim = c(nr, nc, nb))
  
  rownames(plan) <- paste("row", 1:nr)
  colnames(plan) <- paste("col", 1:nc)
  dimnames(plan)[[3]] <- paste("block", 1:nb)
  
  # Add genotypes and checks to the fieldplan

  for (k in 1:nb)
    plan[, , k] <- array(blocks[[k]][plan.id], c(nr, nc))
  
  # Create fielbook
  
  block <- as.integer(gl(nb, nr * nc))
  row <- rep(as.integer(gl(nr, nc)), nb)
  col <- rep(rep(1:nc, nr), nb)
  
  geno <- NULL
  plot.num <- NULL
  to.add <- 0

  for (i in 1:nb) {
    geno <- c(geno, c(t(plan[, , i])))
    plot.num <- c(plot.num, c(t(plan.id)) + to.add)
    to.add <- to.add + sum(!is.na(plan[, , i]))
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
