#' Check row and column positions
#'
#' This function checks that there is only one genotype in each row and column position.
#' This is a wrapper for \code{ck.pos} function.
#' @param dfr The name of the data frame.
#' @param row The name of the column that identifies the rows.
#' @param col The name of the column that identifies the columns.
#' @param rep The name of the column that identifies the replications.
#' @return For each replication a list of plots (unique row and column position)
#' with more than one genotype.
#' @author Raul Eyzaguirre.
#' @examples
#' # Create a design
#' dfr <- cr.rcbd(1:20, 3, 10)
#' dfr <- dfr$book
#' # Check positions
#' check.pos(dfr, 'row', 'col', 'block')
#' @export

check.pos <- function(dfr, row = 'row', col = 'col', rep = 'rep') {
  
  out <- ck.pos(dfr, row, col, rep)
  
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
