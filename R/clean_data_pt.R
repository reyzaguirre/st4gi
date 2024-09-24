#' Clean potato data
#'
#' This is a wrapper for functions \code{setna.pt} and \code{setzero.pt} and
#' applies both in that order.
#' 
#' @param dfr The name of the data frame.
#' @param f Factor for extreme values detection.
#' 
#' @details The data frame must use the labels (lower or upper case) listed
#' in function \code{check.names.pt}. Then functions \code{setna.pt} and
#' \code{setzero.pt} are applied to the data.
#' 
#' @return It returns the data frame with all impossible values set to \code{NA},
#' some values set to \code{0} and a list of warnings with all the rows that have
#' been modified.
#' @author Raul Eyzaguirre.
#' @importFrom stats IQR quantile
#' @export

clean.data.pt <- function(dfr, f = 10) {
  
  .Deprecated("clean.data")
  
  dfr <- setna.pt(dfr, f)
  
  dfr <- setzero.pt(dfr)
  
  # Return
  
  dfr
  
}
