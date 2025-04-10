#' Clean potato and sweetpotato data
#'
#' This is a wrapper for functions \code{setna} and \code{setzero} and
#' applies both in that order.
#' @param dfr The name of the data frame.
#' @param f Factor for extreme values detection.
#' @param crop \code{"auto"} for autodetection or \code{"pt"} for potato and \code{"sp"} for sweetpotato.
#' @param checknames Logical indicating if column names should be checked, default \code{TRUE}.
#' @details The data frame must use the labels (lower or upper case) listed
#' in functions \code{ptont()} and \code{spont()}. Then functions \code{setna} and
#' \code{setzero} are applied to the data.
#' @return It returns the data frame with all impossible values set to \code{NA},
#' some values set to \code{0} and a list of warnings with all the rows that have
#' been modified.
#' @author Raul Eyzaguirre.
#' @examples
#' dfr <- data.frame(mtwp = c(2.2, 5.0, 3.6, 12, 1600, -4, 0),
#'                   dm = c(21, 23, 105, 24, -3, 30, NA),
#'                   nmtp = c(1.3, 10, 11, NA, 2, 5, NA))
#' clean.data(dfr)
#' dfr <- data.frame(crw = c(2.2, 5.0, 3.6, 12, 1600, -4, 0),
#'                   dm = c(21, 23, 105, 24, -3, 30, NA),
#'                   nocr = c(1.3, 10, 11, NA, 2, 5, NA),
#'                   scol = c(1, 0, 15, 5, 4, 7, NA),
#'                   fcol.cc = c(1, 15, 12, 24, 55, 20, NA))
#' clean.data(dfr)
#' @importFrom stats IQR quantile
#' @export

clean.data <- function(dfr, f = 10, crop = c('auto', 'pt', 'sp'), checknames = TRUE) {
  
  # Match arguments
  
  crop = match.arg(crop)
  
  if (crop == 'auto') {
    crop <- detect.crop(dfr)
    warning(crop, " crop detected", call. = FALSE)
  }
  
  if (checknames)
    dfr <- check.names(dfr, crop = crop)
  
  dfr <- setna(dfr, f, crop = crop, checknames = FALSE)
  
  dfr <- setzero(dfr, crop = crop, checknames = FALSE)
  
  # Return
  
  dfr
  
}
