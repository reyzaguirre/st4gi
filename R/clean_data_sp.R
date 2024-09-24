#' Clean sweetpotato data
#'
#' This is a wrapper for functions \code{setna.sp} and \code{setzero.sp} and
#' applies both in that order.
#' 
#' @param dfr The name of the data frame.
#' @param f Factor for extreme values detection.
#' 
#' @details The data frame must use the labels (lower or upper case) listed
#' in function \code{check.names.sp}. Then functions \code{setna.sp} and
#' \code{setzero.sp} are applied to the data.
#' 
#' @return It returns the data frame with all impossible values set to \code{NA},
#' some values set to \code{0} and a list of warnings with all the rows that have
#' been modified.
#' @author Raul Eyzaguirre.
#' @importFrom stats IQR quantile
#' @export

clean.data.sp <- function(dfr, f = 10) {
  
  .Deprecated("clean.data")
  
  dfr <- setna.sp(dfr, f)
  
  dfr <- setzero.sp(dfr)
  
  # Return
  
  dfr
  
}
