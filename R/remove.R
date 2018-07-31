#' Remove factors with NA values
#' 
#' This function removes rows for factors with missing values.
#' 
#' @param factors The factors
#' @param dfr The name of the data frame
#' @return The number of rows in the data frame with missing values for factors
#' (\code{nmis.fac}), and the data frame after removal of all these rows.
#' @author Raul Eyzaguirre.
#' @export

rm.fna <- function(factors, dfr) {
  
  # Check missing values for factors
  
  cond <- apply(data.frame(dfr[, c(factors)]), 1, function(x) sum(is.na(x)) > 0)
  
  # Number of missing values for factors
  
  nmis.fac <- sum(cond)
  
  # Remove rows with missing values for factors
  
  if (nmis.fac > 0)
    dfr <- dfr[!cond, ]
  
  # Return
  
  list(nmis.fac = nmis.fac, dfr = dfr)
  
}

#' Remove factors' levels without data
#' 
#' This function removes rows for factors' levels without data.
#' 
#' @param trait The trait
#' @param factors The factors
#' @param dfr The name of the data frame
#' @return The number of rows in the data frame with factor levels without data
#' (\code{nmis.lev}), and the data frame after removal of all these rows.
#' @author Raul Eyzaguirre.
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
    
  # Check factor levels
    
  for (i in 1:length(factors)) {
      
    tfreq <- data.frame(table(temp[, factors[i]]))
      
    # Identify and remove levels without data
      
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
