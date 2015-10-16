#' Count number of valid values
#'
#' Count the number of valid values for several traits for some specific factors.
#' @param traits List of traits to count for. 
#' @param factors List of factors.
#' @param addcol Additional columns to keep.
#' @param data The name of the data frame containing the data.
#' @author Raul Eyzaguirre
#' @details This function counts the number of valid values for all the \code{traits}
#' for each level's combination of the \code{factors}. Additional columns can be
#' kept if specified in \code{addcol}.
#' @return It returns a data frame with the counts.
#' @examples
#' # The data
#' head(spg)
#' str(spg)
#'
#' # Count the number of values for all the traits across the two replications
#' # and the two locations for each genotype.
#' traits <- c("rytha", "bc", "dm", "star", "nocr")
#' dosum(traits, "geno", data = spg)
#' @export

docount <- function(traits, factors, addcol = NULL, data) {
  
  # Create data.frame
  
  data$temp <- seq(1:dim(data)[1])
    
  dataout <- data[, c("temp", factors, addcol)]
  dataout$dup <- duplicated(dataout[, factors])
  dataout <- subset(dataout, dataout$dup == FALSE)
  dataout <- dataout[, -dim(dataout)[2]]
  
  # Number of factors and traits
  
  nt <- length(traits)
  nf <- length(factors)
  
  # Create column to index count
    
  data$x <- data[, factors[1]]
  dataout$x <- dataout[, factors[1]]
  
  if (nf > 1)
    for (i in 2:nf) {
      data$x <- paste(data$x, data[, factors[i]])
      dataout$x <- paste(dataout$x, dataout[, factors[i]])
    }
  
  # Do count
  
  for (i in 1:nt) {
    for (j in 1:dim(dataout)[1])
      dataout[j, traits[i]] <- sum(!is.na(subset(data, data$x == dataout$x[j])[, traits[i]]))
  }
  
  # Remove x and temp
  
  dataout <- dataout[, !(names(dataout) %in% c("temp", "x"))]
  
  # return data.frame with counts
    
  return(dataout)
}
