#' Collapse design
#'
#' Collapse a blocks design to remove missing plots.
#' @param book The fieldbook design with columns \code{plot}, \code{block}
#' and \code{geno}.
#' @param nc Number of columns in the field.
#' @param fillby Allocate the plots by \code{"rows"} or \code{"columns"},
#' default \code{"rows"}.
#' @param serpentine \code{"yes"} or \code{"no"}, default \code{"yes"}.
#' @details If \code{"nc"} is \code{"NULL"}, the function will look for the
#' most squared rectangular field.
#' @return It returns the fieldbook and fieldplan.
#' @author Raul Eyzaguirre.
#' @examples
#' checks <- paste("ch", 1:4, sep = "_")
#' genos <- paste("g", 1:20, sep = "_")
#' book <- cr.abd(genos, checks, 4, 5)$book
#' cds(book, 5, 'rows', 'yes')
#' @export

cds <- function(book, nc = NULL,
                fillby = c('rows', 'columns'),
                serpentine = c("yes", "no")) {
  
  # Match arguments
  
  fillby <- match.arg(fillby)
  serpentine <- match.arg(serpentine)
  
  # Sort by plots
  
  book <- book[order(book$plot), ]
  
  # Number of plots
  
  np <- dim(book)[1]

  # Number of rows and columns
  
  if (is.null(nc))
    nc <- gnc(np)
  
  nr <- ceiling(np / nc)

  # Fieldplan array
  
  plan.id <- fp(nr, nc, fillby, serpentine)

  # Create fieldplan
  
  geno <- book$geno
  block <- book$block
  plan <- array(geno[plan.id], c(nr, nc))
  blockplan <- array(block[plan.id], c(nr, nc))
      
  rownames(plan) <- paste("row", 1:nr)
  colnames(plan) <- paste("col", 1:nc)
      
  # Rows and columns numbers
  
  row <- as.integer(gl(nr, nc))
  col <- rep(1:nc, nr)
  
  # Create fielbook with new rows and columns
    
  book <- data.frame(plot = c(t(plan.id)), block = c(t(blockplan)),
                     row, col, geno = c(t(plan)),
                     stringsAsFactors = FALSE)
  book <- book[!is.na(book$geno), ]
  
  # Sort by plot number
  
  book <- book[order(book$plot), ]
  rownames(book) <- 1:dim(book)[1]
  
  # Return
  
  list(plan = plan, book = book)
  
}
