#' Area Under the Disease Progress Curve
#'
#' Computes the absolute (AUDPC) or relative (rAUDPC) area under the disease progress curve.
#' @param evals The names of the columns for the evaluations.
#' @param dates A vector with the dates of evaluation.
#' @param dfr The name of the data frame.
#' @details Evaluations are subjective ranging from 0 (no disease) to 100 (full disease).
#' Dates should follow the format YYYY-MM-DD.
#' \code{audpc} is expressed in days x % and \code{raudpc} as a proportion.
#' @return It returns the original data frame with added columns for AUDPC and rAUDPC.
#' @author Raul Eyzaguirre.
#' @examples
#' head(lbb3c3)
#' evals <- c("lb1", "lb2", "lb3", "lb4", "lb5", "lb6")
#' dates <- c("2017-11-20", "2017-11-27", "2017-12-04",
#'            "2017-12-11", "2017-12-18", "2017-12-25")
#' tmp <- audpc(evals, dates, lbb3c3)
#' head(tmp)              
#' @export

audpc <- function(evals, dates, dfr) {

  # Check number of evals and dates
  
  ne <- length(evals)
  nd <- length(dates)
  
  if (ne != nd)
    stop("Number of evaluations and number of dates are different.")
  
  # Check dates

  dates <- as.Date(dates, format = "%Y-%m-%d")
  
  for (i in 2:nd) {
    if (dates[i] <= dates[i - 1])
      stop("Dates must be ordered.")
  }
    
  # Compute AUDPC
  
  dfr$audpc <- 0
  
  for (i in 2:nd) {
    
    days <- as.numeric(dates[i] - dates[i - 1])
    dfr$audpc <- dfr$audpc + (dfr[, evals[i - 1]] + dfr[, evals[i]]) / 2 * days
    
  }
  
  total.days <- as.numeric(dates[nd] - dates[1])
  
  dfr$raudpc <- dfr$audpc / total.days / 100
  
  # result

  dfr
  
}
