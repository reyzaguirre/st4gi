#' Check row and column positions
#'
#' This function checks that there is only one genotype in each row and column position.
#' @param dfr The name of the data frame.
#' @param row The name of the column that identifies the rows.
#' @param col The name of the column that identifies the columns.
#' @param rep The name of the column that identifies the replications.
#' @details If \code{rep} is not \code{NULL}, then it checks
#' positions for each replication.
#' @return A list of plots (unique row and column position)
#' with more than one genotype.
#' @author Raul Eyzaguirre.
#' @examples
#' # Create a design
#' dfr <- cr.rcbd(1:20, 3, 10)
#' dfr <- dfr$book
#' # Check positions
#' check.pos(dfr, 'row', 'col', 'block')
#' @export

check.pos <- function(dfr, row, col, rep = NULL) {
  
  out <- ck.pos(dfr, row, col, rep)
  
  if (out$nrep == 1) {
    
    cat('------------------------------\n')
    cat('Full experiment ',           '\n')
    cat('------------------------------\n')
    if (out$nplot > 0) {
      cat('More than one genotype in the same position: \n')
      print(out$lplot[[1]])
      cat('\n')
    } else {
      cat('OK \n')
      cat('\n')
    }
    
  } else {
    
    for (i in 1:out$nrep) {
      
      cat('------------------------------\n')
      cat('Replication', out$lrep[i],   '\n')
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
  
}

ck.pos <- function(dfr, row, col, rep = NULL) {
  
  # Define units to check
  
  if (is.null(rep)) {
    dfr[, "rep"] <- 1
    rep <- "rep"
  }
  
  # Check and remove rows with missing values for factors
  
  out <- ck.fs(dfr, c(row, col), rep)
  dfr <- out$dfr
  nrep <- out$nrep
  lrep <- out$lrep
  
  # Number and list of plots with more than one genotype
  
  nplot <- NULL
  lplot <- list()
  
  # Check row and column
  
  for (i in 1:nrep) {
    
    # Compute frequencies
    
    tmp <- dfr[dfr[, rep] == lrep[i], ]
    ttt <- as.data.frame(table(tmp[, row], tmp[, col]))
    colnames(ttt) <- c('Row', 'Column', 'Freq')
    
    # Number of plots with more than one genotype
    
    nplot[i] <- dim(ttt[ttt$Freq > 1, ])[1]
    
    # List of plots with more than one genotype if any
    
    if (nplot[i] > 0)
      lplot[[i]] <- ttt[ttt$Freq > 1, ] 
  }
  
  # Return
  
  list(nplot = nplot, lplot = lplot, nrep = nrep, lrep = lrep)
  
}