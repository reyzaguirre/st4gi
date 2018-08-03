#' Check row and column positions
#'
#' This function checks that there is only one genotype in each row and column position.
#' This is a wrapper for \code{ck.pos} function.
#' @param row Label for rows.
#' @param col Label for columns.
#' @param rep Label for replications.
#' @param dfr The name of the data frame.
#' @return For each replication a list of plots (unique row and column position)
#' with more than one genotype.
#' @author Raul Eyzaguirre.
#' @examples
#' # Create a design
#' dfr <- cr.rcbd(1:20, 3, 10)
#' dfr <- dfr$book
#' # Check positions
#' check.pos('row', 'col', 'block', dfr)
#' @export

check.pos <- function(row, col, rep, dfr) {
  
  out <- ck.pos(row, col, rep, dfr)
  
  for (i in 1:out$nrep) {
    
    # Print list of plots with problems
    
    cat('------------------------------\n')
    cat('Replication', out$lrep[i], '\n')
    cat('------------------------------\n')
    if (out$nplot[i] > 0) {
      cat('More than one genotype in the same position: \n')
      print(out$lplot[[i]])
      cat('\n')
    } else {
      cat('OK \n')
      cat('\n')
    }
  }
  
}
