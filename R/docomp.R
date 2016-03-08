#' Do computations over some factors
#'
#' Do computations for several traits for some specific factors.
#' @param do The computation to perform. Implemented options are \code{count},
#' \code{max}, \code{mean}, \code{min}, and \code{sum}.
#' @param traits List of traits. 
#' @param factors List of factors.
#' @param addcol Additional columns to keep.
#' @param data The name of the data frame containing the data.
#' @author Raul Eyzaguirre
#' @details This function do a specific computation for all the \code{traits}
#' for each level's combination of the \code{factors}. Additional columns can be
#' kept if specified in \code{addcol}.
#' @return It returns a data frame with the computations.
#' @examples
#' # The data
#' head(spg)
#' str(spg)
#'
#' # Compute means for all the traits across the two replications
#' # for each genotype and location and then for each genotype 
#' # across the two locations.
#' traits <- c("rytha", "bc", "dm", "star", "nocr")
#' factors <- c("geno", "loc")
#' output1 <- docomp("mean", traits, factors, data = spg)
#' docomp("mean", traits, "geno", data = output1)
#' 
#' # Compute maxima for all the traits across the two replications
#' # for each genotype and location.
#' docomp("max", traits, factors, data = spg)
#' @export

docomp <- function(do, traits, factors, addcol = NULL, data) {

  # Create data.frame
  
  dataout <- data.frame(data[, c(factors, addcol)])
  colnames(dataout) <- c(factors, addcol)
  dataout <- subset(dataout, duplicated(dataout[, factors]) == F)

  # Number of factors and traits
  
  nt <- length(traits)
  nf <- length(factors)
  
  # Create matching vars
    
  idin <- data[, factors[1]]
  idout <- dataout[, factors[1]]
  
  if (nf > 1)
    for (i in 2:nf) {
      idin <- paste(idin, data[, factors[i]])
      idout <- paste(idout, dataout[, factors[i]])
    }
  
  # Do computations
  
  for (i in 1:nt) {
    for (j in 1:dim(dataout)[1]){
      if (do == "count")
        dataout[j, traits[i]] <- sum(!is.na(data[idin == idout[j], traits[i]]))
      if (do == "max")
        dataout[j, traits[i]] <- max(data[idin == idout[j], traits[i]], na.rm = TRUE)
      if (do == "mean")
        dataout[j, traits[i]] <- mean(data[idin == idout[j], traits[i]], na.rm = TRUE)
      if (do == "min")
        dataout[j, traits[i]] <- min(data[idin == idout[j], traits[i]], na.rm = TRUE)
      if (do == "sum")
        dataout[j, traits[i]] <- sum(data[idin == idout[j], traits[i]], na.rm = TRUE)
    }
  }
  
  # return data.frame with maxima
    
  dataout
}


