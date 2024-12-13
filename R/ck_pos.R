#' Check row and column positions
#'
#' This function checks that there is only one genotype in each row and column position
#' for each replication.
#' @param dfr The name of the data frame.
#' @param row The name of the column that identifies the rows.
#' @param col The name of the column that identifies the columns.
#' @param rep The name of the column that identifies the replications.
#' @return For each replication, the number (\code{nplot}) and list (\code{lplot})
#' of plots (unique row and column position) with more than one genotype, and the 
#' number (\code{nrep}) and list (\code{lrep}) of replications.
#' @author Raul Eyzaguirre.
#' @examples
#' # Create a design
#' dfr <- cr.rcbd(1:20, 3, 10)
#' dfr <- dfr$book
#' # Check positions
#' ck.pos(dfr, 'row', 'col', 'block')
#' @export

ck.pos <- function(dfr, row, col, rep = NULL) {
  
  # Define replications
  
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
    
    temp <- dfr[dfr[, rep] == lrep[i], ]
    ttt <- as.data.frame(table(temp[, row], temp[, col]))
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
