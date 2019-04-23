#' Completely Randomized Design
#'
#' This function creates the fieldbook and fieldplan for a CRD.
#' @param geno The list of genotypes.
#' @param nrep Number of replications.
#' @param nc Number of available columns on the field.
#' @param serpentine \code{"yes"} or \code{"no"}, default \code{"yes"}.
#' @return It returns the fieldbook and fieldplan.
#' @author Raul Eyzaguirre.
#' @examples
#' cr.crd(1:20, 3, 12)
#' cr.crd(1:20, 2, 7)
#' @export

cr.crd <- function(geno, nrep, nc = NULL, serpentine = c("yes", "no")) {
  
  # Match arguments
  
  serpentine <- match.arg(serpentine)

  # Error messages
  
  ng <- length(geno)

  if (ng < 2)
    stop("Include at least 2 genotypes.")

  # Number of rows and columns
  
  if (is.null(nc))
    nc <- gnc(ng * nrep)
  
  nr <- ceiling(ng * nrep / nc)

  # Fieldplan array
  
  plan.id <- fp(nr, nc, serpentine)

  # Sort genotypes
  
  geno <- sample(rep(geno, nrep))
  
  # Create fieldplan
  
  plan <- array(geno[plan.id], c(nr, nc))
      
  rownames(plan) <- paste("row", 1:nr)
  colnames(plan) <- paste("col", 1:nc)
      
  # Rows and columns numbers
  
  row <- as.integer(gl(nr, nc))
  col <- rep(1:nc, nr)
  
  # Create fielbook
    
  book <- data.frame(plot.num = c(t(plan.id)), row, col, geno = c(t(plan)),
                     stringsAsFactors = FALSE)
  book <- book[!is.na(book$geno), ]
  
  # Sort by plot number
  
  if (serpentine == 'yes' & nr > 1) {
    book <- book[sort(book$plot.num, index.return = TRUE)$ix, ]
    rownames(book) <- 1:dim(book)[1]
  }
    
  # Return
  
  list(plan = plan, book = book)
  
}
