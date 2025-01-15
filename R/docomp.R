#' Do computations over some factors
#'
#' Do computations for several variables for some specific factors.
#' @param dfr The name of the data frame.
#' @param do The computation to perform. Implemented options are \code{count},
#' \code{mode}, and standard functions like \code{mean}, \code{median},
#' \code{min}, \code{max}, \code{sd}, \code{var}, \code{sum}, etc.
#' @param vars The names of the columns for the variables. 
#' @param factors The names of the columns for the factors.
#' @param keep The names of additional columns to keep.
#' @param method Use \code{fast} or \code{slow} method. 
#' @details This function do a specific computation for all the \code{vars}
#' for each level's combination of the \code{factors}. Additional columns can be
#' kept if specified in \code{keep}. All \code{factors} and \code{keep} values
#' are converted to character. \code{do = "count"} counts the number
#' of valid cases (excluding missing values).
#' @return It returns a data frame with the computations.
#' @author Raul Eyzaguirre.
#' @examples
#' # Compute means across replications and then across locations for each genotype
#' vars <- c("rytha", "bc", "dm", "star", "nocr")
#' factors <- c("geno", "loc")
#' output1 <- docomp(spg, "mean", vars, factors)
#' docomp("mean", vars, "geno", dfr = output1) 
#' # Compute maxima across replications for each genotype and location.
#' docomp(spg, "max", vars, factors)
#' @export

docomp <- function(dfr, do, vars, factors, keep = NULL, method = c("fast", "slow")) {

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
  
  # Number of factors and variables
  
  nt <- length(vars)
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
          dfr.out[j, vars[i]] <- sum(!is.na(dfr[idin == idout[j], vars[i]]))
        if (do == 'mode')
          dfr.out[j, vars[i]] <- getmode(dfr[idin == idout[j], vars[i]])
        if (!do %in% c('count', 'mode')) {
          if (sum(!is.na(dfr[idin == idout[j], vars[i]])) == 0)
            dfr.out[j, vars[i]] <- NA
          else
            dfr.out[j, vars[i]] <- eval(parse(text = do))(dfr[idin == idout[j], vars[i]], na.rm = TRUE)
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
        
      dfr.y <- data.frame(tapply(dfr[, vars[i]], idin, foo))
      dfr.y[, 'id.to.merge'] <- rownames(dfr.y)
      colnames(dfr.y)[1] <- vars[i]
      dfr.out <- merge(dfr.out, dfr.y, sort = FALSE)
    
    }

    dfr.out <- dfr.out[, colnames(dfr.out) != 'id.to.merge']
        
  }
  
  # return data.frame
    
  dfr.out
  
}
