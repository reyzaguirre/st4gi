#' Remove factors with NA values
#' 
#' This function removes rows for factors with missing values.
#' 
#' @param factors The factors.
#' @param dfr The name of the data frame.
#' @return The number of rows in the data frame with missing values for factors
#' (\code{nmis.fac}), and the data frame after removal of all these rows.
#' @author Raul Eyzaguirre.
#' @examples
#' # Create a design
#' dfr <- cr.crd(1:10, 3, 10)
#' dfr <- dfr$book
#' # Delete some values for classification factors
#' dfr[5, 2] <- NA
#' dfr[7, 3] <- NA
#' dfr[11, 4] <- NA
#' # Check missing values
#' rm.fna(c('row', 'col', 'geno'), dfr)
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
