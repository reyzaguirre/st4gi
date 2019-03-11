#' Check row and column positions
#'
#' This function checks that there is only one genotype in each row and column position
#' for each replication.
#' @param row Label for rows.
#' @param col Label for columns.
#' @param rep Label for replications.
#' @param dfr The name of the data frame.
#' @return For each replication, the number (\code{nplot}) and list (\code{lplot})
#' of plots (unique row and column position) with more than one genotype, and the 
#' number (\code{nrep}) and list (\code{lrep}) of replications.
#' @author Raul Eyzaguirre.
#' @examples
#' # Create a design
#' dfr <- cr.rcbd(1:20, 3, 10)
#' dfr <- dfr$book
#' # Check positions
#' ck.pos('row', 'col', 'block', dfr)
#' @export

ck.pos <- function(row, col, rep = NULL, dfr) {
  
  # Define replications
  
  if (is.null(rep)) {
    dfr[, "rep"] <- 1
    rep <- "rep"
  }
  
  # Check and remove rows with missing values for factors
  
  out <- ck.fs(c(row, col), rep, dfr)
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
