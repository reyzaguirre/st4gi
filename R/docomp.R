#' Do computations over some factors
#'
#' Do computations for several traits for some specific factors.
#' @param do The computation to perform. Implemented options are \code{count},
#' \code{mode}, and standard functions like \code{mean}, \code{median},
#' \code{min}, \code{max}, \code{sd}, \code{var}, \code{sum}, etc.
#' @param traits List of traits. 
#' @param factors List of factors.
#' @param keep Additional columns to keep.
#' @param dfr The name of the data frame.
#' @param method Use \code{fast} or \code{slow} method. 
#' @details This function do a specific computation for all the \code{traits}
#' for each level's combination of the \code{factors}. Additional columns can be
#' kept if specified in \code{keep}. All \code{factors} and \code{keep} values
#' are converted to character. \code{do = "count"} counts the number
#' of valid cases (excluding missing values).
#' @return It returns a data frame with the computations.
#' @author Raul Eyzaguirre.
#' @examples
#' # Compute means across replications and then across locations for each genotype
#' traits <- c("rytha", "bc", "dm", "star", "nocr")
#' factors <- c("geno", "loc")
#' output1 <- docomp("mean", traits, factors, dfr = spg)
#' docomp("mean", traits, "geno", dfr = output1) 
#' # Compute maxima across replications for each genotype and location.
#' docomp("max", traits, factors, dfr = spg)
#' @export

docomp <- function(do, traits, factors, keep = NULL, dfr, method = c("fast", "slow")) {

  # Match arguments
  
  method <- match.arg(method)
  
  # Internal mode function
  
  getmode <- function(v) {
    x <- names(which.max(table(unlist(v))))
    if (is.null(x))
      x <- NA
    x
  }
  
  # Create data.frame
  
  dfr.out <- data.frame(dfr[, c(factors, keep)])
  colnames(dfr.out) <- c(factors, keep)
  dfr.out <- subset(dfr.out, !duplicated(dfr.out[, factors]))
  
  # Number of factors and traits
  
  nt <- length(traits)
  nf <- length(factors)
  
  # Create matching vars
    
  idin <- dfr[, factors[1]]
  idout <- dfr.out[, factors[1]]
  
  if (nf > 1)
    for (i in 2:nf) {
      idin <- paste(idin, dfr[, factors[i]])
      idout <- paste(idout, dfr.out[, factors[i]])
    }
  
  # Do computations
  
  if (method == "slow") {

    for (i in 1:nt) {
      for (j in 1:dim(dfr.out)[1]){
        if (do == "count")
          dfr.out[j, traits[i]] <- sum(!is.na(dfr[idin == idout[j], traits[i]]))
        if (do == 'mode')
          dfr.out[j, traits[i]] <- getmode(dfr[idin == idout[j], traits[i]])
        if (!do %in% c('count', 'mode')) {
          if (sum(!is.na(dfr[idin == idout[j], traits[i]])) == 0)
            dfr.out[j, traits[i]] <- NA
          else
            dfr.out[j, traits[i]] <- eval(parse(text = do))(dfr[idin == idout[j], traits[i]], na.rm = TRUE)
        }
      }
    }
    
    row.names(dfr.out) <- 1:dim(dfr.out)[1]
    
  }

  if (method == "fast") {
    
    if (do == "count")
      foo <- function(x) sum(!is.na(x))
    if (do == 'mode')
      foo <- getmode
    if (!do %in% c('count', 'mode'))
      foo <- function(x) {
        if (sum(!is.na(x)) == 0)
          NA
        else
          eval(parse(text = do))(x, na.rm = TRUE)
      }
    
    dfr.out[, "id.to.merge"] <- idout
        
    for (i in 1:nt) {
        
      dfr.trait <- data.frame(tapply(dfr[, traits[i]], idin, foo))
      dfr.trait[, 'id.to.merge'] <- rownames(dfr.trait)
      colnames(dfr.trait)[1] <- traits[i]
      dfr.out <- merge(dfr.out, dfr.trait, sort = FALSE)
    
    }

    dfr.out <- dfr.out[, colnames(dfr.out) != 'id.to.merge']
        
  }
  
  # return data.frame
    
  dfr.out
  
}
