#' Remove factors' levels without data
#' 
#' This function removes rows for factors' levels without data. 
#' @param trait The trait.
#' @param factors The factors.
#' @param dfr The name of the data frame.
#' @return The number of rows in the data frame with factor levels without data
#' (\code{nmis.lev}), and the data frame after removal of all these rows.
#' @author Raul Eyzaguirre.
#' @examples
#' # Create a design
#' dfr <- cr.crd(1:10, 3, 10)
#' dfr <- dfr$book
#' # Create some random data
#' dfr$y <- rnorm(30)
#' # Delete values for some factor levels
#' dfr[dfr$geno == 5, 'y'] <- NA
#' # Check missing values
#' rm.lna('y', 'geno', dfr)
#' @export

rm.lna <- function(trait, factors, dfr) {
  
  # Remove rows with missing values for factors and get a copy
  
  dfr <- rm.fna(factors, dfr)$dfr
  temp <- dfr
  
  # To factor to preserve levels
  
  for (i in 1:length(factors))
    temp[, factors[i]] <- factor(temp[, factors[i]])
  
  # Number of rows before
    
  nrb <- dim(temp)[1]
      
  # Get a subset with data
    
  temp <- temp[!is.na(temp[, trait]), ]
    
  # Check factor levels and remove levels without data
    
  for (i in 1:length(factors)) {
    tfreq <- data.frame(table(temp[, factors[i]]))
    nodata <- as.character(tfreq[tfreq$Freq == 0, 'Var1'])
    dfr <- dfr[!(dfr[, factors[i]] %in% nodata), ]
  }
    
  # Number of rows after
    
  nra <- dim(dfr)[1]
    
  # Number of rows deleted
    
  nmis.lev <- nrb - nra
    
  # Return
  
  list(nmis.lev = nmis.lev, dfr = dfr)
  
}
