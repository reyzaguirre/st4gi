#' Unreplicated experiment
#'
#' This function creates the fieldbook and fieldplan for an unreplicated
#' experiment with genotypes randomly allocated on a field.
#' @param geno The list of genotypes.
#' @param nc Number of available columns on the field.
#' @param serpentine \code{"yes"} or \code{"no"}, default \code{"yes"}.
#' @return It returns the fieldbook and fieldplan.
#' @author Raul Eyzaguirre.
#' @examples
#' cr.ur(1:100, 5)
#' cr.ur(1:100, 7)
#' @export

cr.ur <- function(geno, nc, serpentine = c("yes", "no")) {
  
  # Match arguments
  
  serpentine <- match.arg(serpentine)

  # Error messages
  
  ng <- length(geno)
  
  if (ng < 2)
    stop("Include at least 2 genotypes.")
  
  # Number of rows
  
  nr <- ceiling(ng / nc)
  
  # Fieldplan array
  
  plan.id <- fp(nr, nc, serpentine)
  
  plan <- array(dim = c(nr, nc))

  rownames(plan) <- paste("row", 1:nr)
  colnames(plan) <- paste("col", 1:nc)

  # Include genotypes at random
  
  geno <- sample(geno)

  for (i in 1:nr)
    for (j in 1:nc)
      plan[i, j] <- geno[plan.id[i, j]]
  
  # Rows and columns numbers
  
  row <- as.integer(gl(nr, nc))
  col <- rep(1:nc, nr)

  # Create fielbook
  
  book <- data.frame(plot.num = c(t(plan.id)), row, col, geno = c(t(plan)),
                     stringsAsFactors = F)
  book <- book[!is.na(book$geno), ]

  # Sort by plot number
  
  if (serpentine == 'yes' & nr > 1) {
    book <- book[sort(book$plot.num, index.return = TRUE)$ix, ]
    rownames(book) <- 1:dim(book)[1]
  }
  
  # Return
  
  list(plan = plan, book = book)
  
}
