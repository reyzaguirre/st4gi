#' Completely Randomized Design
#'
#' This function creates the fieldbook and fieldplan for a CRD.
#' @param geno The list of genotypes.
#' @param nrep Number of replications.
#' @param nc Number of columns in the field.
#' @param fillby Allocate the plots by \code{"rows"} or \code{"columns"},
#' default \code{"rows"}.
#' @param serpentine \code{"yes"} or \code{"no"}, default \code{"yes"}.
#' @details If \code{"nc"} is \code{"NULL"}, the function will look for the
#' most squared rectangular field.
#' @return It returns the fieldbook and fieldplan.
#' @author Raul Eyzaguirre.
#' @examples
#' ds <- cr.crd(1:20, 3, 12)
#' ds$plan
#' ds$book
#' plot(ds)
#' ds <- cr.crd(1:20, 2, 7)
#' ds$plan
#' ds$book
#' plot(ds)
#' @export

cr.crd <- function(geno, nrep, nc = NULL,
                   fillby = c('rows', 'columns'),
                   serpentine = c("yes", "no")) {
  
  # Match arguments
  
  fillby <- match.arg(fillby)
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
  
  plan.id <- fp(nr, nc, fillby, serpentine)

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
    
  book <- data.frame(plot = c(t(plan.id)), row, col, geno = c(t(plan)),
                     stringsAsFactors = FALSE)
  book <- book[!is.na(book$geno), ]
  
  # Sort by plot number
  
  book <- book[order(book$plot), ]
  rownames(book) <- 1:dim(book)[1]
  
  # Return
  
  output <- list(plan = plan, book = book, design = 'CRD')
  class(output) <- "st4gi_design"
  invisible(output)
  
}
