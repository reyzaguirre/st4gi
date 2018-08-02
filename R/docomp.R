#' Do computations over some factors
#'
#' Do computations for several traits for some specific factors.
#' @param do The computation to perform. Implemented options are \code{count},
#' and standard functions like \code{mean}, \code{median}, \code{min}, \code{max},
#' \code{sd}, \code{var}, \code{sum}, etc.
#' @param traits List of traits. 
#' @param factors List of factors.
#' @param adc Additional columns to keep.
#' @param dfr The name of the data frame containing the data.
#' @author Raul Eyzaguirre
#' @details This function do a specific computation for all the \code{traits}
#' for each level's combination of the \code{factors}. Additional columns can be
#' kept if specified in \code{adc}. All \code{factors} and \code{adc} values
#' are converted to character. \code{do = "count"} counts the number
#' of valid cases (excluding missing values).
#' @return It returns a data frame with the computations.
#' @examples
#' # Compute means across replications and then across locations for each genotype
#' traits <- c("rytha", "bc", "dm", "star", "nocr")
#' factors <- c("geno", "loc")
#' output1 <- docomp("mean", traits, factors, dfr = spg)
#' docomp("mean", traits, "geno", dfr = output1) 
#' # Compute maxima across replications for each genotype and location.
#' docomp("max", traits, factors, dfr = spg)
#' @export

docomp <- function(do, traits, factors, adc = NULL, dfr) {

  # Create data.frame
  
  dfrout <- data.frame(dfr[, c(factors, adc)])
  colnames(dfrout) <- c(factors, adc)
  dfrout <- subset(dfrout, !duplicated(dfrout[, factors]))
  
  # Number of factors and traits
  
  nt <- length(traits)
  nf <- length(factors)
  
  # Create matching vars
    
  idin <- dfr[, factors[1]]
  idout <- dfrout[, factors[1]]
  
  if (nf > 1)
    for (i in 2:nf) {
      idin <- paste(idin, dfr[, factors[i]])
      idout <- paste(idout, dfrout[, factors[i]])
    }
  
  # Do computations
  
  for (i in 1:nt) {
    for (j in 1:dim(dfrout)[1]){
      if (do == "count")
        dfrout[j, traits[i]] <- sum(!is.na(dfr[idin == idout[j], traits[i]]))
      else
        dfrout[j, traits[i]] <- eval(parse(text = do))(dfr[idin == idout[j], traits[i]], na.rm = TRUE)
    }
  }
  
  # return data.frame
    
  dfrout
  
}
