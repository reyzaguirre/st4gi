#' Do computations over some factors
#'
#' Do computations for several traits for some specific factors.
#' @param do The computation to perform. Implemented options are \code{count},
#' and standard functions like \code{mean}, \code{median}, \code{min}, \code{max},
#' \code{sd}, \code{var}, \code{sum}, etc.
#' @param traits List of traits. 
#' @param factors List of factors.
#' @param addcol Additional columns to keep.
#' @param data The name of the data frame containing the data.
#' @author Raul Eyzaguirre
#' @details This function do a specific computation for all the \code{traits}
#' for each level's combination of the \code{factors}. Additional columns can be
#' kept if specified in \code{addcol}. All \code{factors} and \code{addcol} values
#' are converted to character. \code{do = "count"} counts the number
#' of valid cases (excluding missing values).
#' @return It returns a data frame with the computations.
#' @examples
#' ## Compute means across replications and then across locations for each genotype
#' traits <- c("rytha", "bc", "dm", "star", "nocr")
#' factors <- c("geno", "loc")
#' output1 <- docomp("mean", traits, factors, data = spg)
#' docomp("mean", traits, "geno", data = output1)
#' 
#' ## Compute maxima across replications for each genotype and location.
#' docomp("max", traits, factors, data = spg)
#' @export

docomp <- function(do, traits, factors, addcol = NULL, data) {

  # Convert to character
  
  n <- dim(data[, c(factors, addcol)])[2]
  for (i in 1:n)
    data[, c(factors, addcol)][, i] <- as.character(data[, c(factors, addcol)][, i])
  
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
      else
        dataout[j, traits[i]] <- eval(parse(text = do))(data[idin == idout[j], traits[i]], na.rm = TRUE)
    }
  }
  
  # return data.frame with maxima
    
  dataout
}
