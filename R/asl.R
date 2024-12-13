#' Add non significant lines
#'
#' This function adds columns to a datra frame with averages to identify
#' non significant pairwise comparisons.
#' @param dfr The name of the data frame.
#' @param geno The name of the column that identifies the genotypes.
#' @param aver The name of the column with the averages.
#' @param lsd Least significant value.
#' @param symbol If \code{NULLL} it uses etters by default. See \code{details}
#' @details The data frame is sorted descending and columns are added to
#' identify non singnificant pairwise comparisons. By default it works with
#' letters unless a different symbol is specified in \code{symbol}.
#' @return It returns a data frame with averages and non significant lines.
#' @author Raul Eyzaguirre.
#' @examples 
#' # Adjust model
#' model <- aov.rcbd(pjpz09, "trw", "geno", "rep")
#' # LSD values
#' lsd <- sqrt(model[3, 3]) * qt(.975, model[3, 1])
#' # Compute means
#' ests <- docomp(pjpz09, 'mean', 'trw', 'geno')
#' # Add non significant lines
#' ests <- asl(ests, 'geno', 'trw', lsd)
#' @export

asl <- function(dfr, geno, aver, lsd, symbol = NULL) {
  
  # Define symbols
  
  if (is.null(symbol))
    symbol = letters
  
  nsym <- length(symbol)

  # Sort data frame
  
  dfr <- dfr[order(dfr[, aver], decreasing = TRUE), ]
  
  # First line
  
  i <- 1 # count for genotypes
  j <- 1 # count for column names
  k <- 1 # count for symbols
  
  lsup <- dfr[i, aver]
  linf <- dfr[i, aver] - lsd
  
  # Column name
  
  col.name <- paste0('line.', j)
  
  # Mark first line
  
  dfr[, col.name] <- NA
  dfr[dfr[, aver] <= lsup & dfr[, aver] >= linf , col.name] <- symbol[k]
  
  # Detect last value
  
  last.value <- min(dfr[dfr[, aver] > linf, aver])

  # Additional lines
  
  nrows <- dim(dfr)[1]

  while (linf > dfr[nrows, aver]) {
    
    i <- i + 1
    
    # Define limits
    
    lsup <- dfr[i, aver]
    linf <- dfr[i, aver] - lsd
    
    # Check if new line is longer than previous line
    
    new.last.value <- min(dfr[dfr[, aver] > linf, aver])
    
    if (new.last.value < last.value) {
      
      j <- j + 1
      k <- k + 1
        
      # Column name
          
      col.name <- paste0('line.', j)
      
      # Recycle simbols
      
      if (k > nsym)
        k <- k - nsym
          
      # Mark line
          
      dfr[, col.name] <- NA
      dfr[dfr[, aver] <= lsup & dfr[, aver] >= linf , col.name] <- symbol[k]
      
      # Last value
      
      last.value <- new.last.value
      
    }
    
  }
  
  # Return
  
  dfr
  
}
