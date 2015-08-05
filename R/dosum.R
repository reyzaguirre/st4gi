#' Compute sum over some factors
#'
#' Compute the sum for several traits for some specific factors.
#' @param traits List of traits to compute the sum. 
#' @param factors List of factors.
#' @param addcol Additional columns to keep.
#' @param data The name of the data frame containing the data.
#' @author Raul Eyzaguirre
#' @details This function computes the sum for all the \code{traits}
#' for each level's combination of the \code{factors}. Additional columns can be
#' kept if specified in \code{addcol}.
#' @return It returns a data frame with the sum.
#' @examples
#' # The data
#' head(spg)
#' str(spg)
#'
#' # Compute the sum for all the traits across the two replications
#' # for each genotype and location.
#' traits <- c('rytha', 'bc', 'dm', 'star', 'nocr')
#' factors <- c('geno', 'loc')
#' dosum(traits, factors, data=spg)
#' 
#' # Save the output in a data.frame with name 'output1'
#' # and compute the sum for each genotype across the two locations.
#' output1 <- dosum(traits, factors, data=spg)
#' dosum(traits, 'geno', data=output1)
#' @export

dosum <- function(traits, factors, addcol = NULL, data){
  
  # Create data.frame
  
  data$temp <- seq(1:dim(data)[1])
    
  dataout <- data[, c("temp", factors, addcol)]
  dataout$dup <- duplicated(dataout[, factors])
  dataout <- subset(dataout, dup==F)
  dataout <- dataout[, -dim(dataout)[2]]
  
  # Number of factors and traits
  
  nt <- length(traits)
  nf <- length(factors)
  
  # Create column to index sum
    
  data$x <- data[, factors[1]]
  dataout$x <- dataout[, factors[1]]
  
  if (nf > 1)
    for (i in 2:nf){
      data$x <- paste(data$x, data[, factors[i]])
      dataout$x <- paste(dataout$x, dataout[, factors[i]])
    }
  
  # Compute sum
  
  for(i in 1:nt){
    for (j in 1:dim(dataout)[1])
      dataout[j, traits[i]] <- sum(subset(data, x == dataout$x[j])[, traits[i]], na.rm=T)
  }
  
  # Remove x and temp
  
  dataout <- dataout[, !(names(dataout) %in% c("temp", "x"))]
  
  # return data.frame with maxima
    
  return(dataout)
}
