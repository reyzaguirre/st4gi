#' Area Under the Disease Progress Curve
#'
#' Computes the absolute (AUDPC) or relative (rAUDPC) area under the disease progress curve.
#' @param dfr The name of the data frame.
#' @param evals The names of the columns for the evaluations.
#' @param dates A vector with the dates of evaluation.
#' @param na Logical, if \code{TRUE}, \code{NA} values are replace with the 
#' contiguous values. Default is \code{FALSE}.
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
#' tmp <- audpc(lbb3c3, evals, dates)
#' head(tmp)              
#' @export

audpc <- function(dfr, evals, dates, na = FALSE) {

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
  
  # Check if there are missing values and report
  
  if (sum(is.na(dfr[, evals])) > 0) {
    
    report.table <- data.frame(evaluation = evals,
                               date = dates,
                               n.missing = NA,
                               p.missing = NA)
    
    n.total <- dim(dfr)[1]
    
    for (i in 1:nd) {
      
      report.table$n.missing[i] <- sum(is.na(dfr[, evals[i]]))
      report.table$p.missing[i] <- round(report.table$n.missing[i] / n.total, 4)
      
    }
    
    cat('----------------------------------------\n')
    cat('There are missing values in evaluations','\n')
    cat('----------------------------------------\n')
    cat('\n')
    print(report.table)
    cat('\n')

  }
  
  # Compute AUDPC
  
  dfr$audpc <- 0
  
  for (i in 2:nd) {
    
    days <- as.numeric(dates[i] - dates[i - 1])
    
    # Compute audpc without NAs
    
    if (!na)
      dfr$audpc <- dfr$audpc + (dfr[, evals[i - 1]] + dfr[, evals[i]]) / 2 * days
    
    # Compute audpc with NAs
    
    if (na)
      dfr$audpc <- dfr$audpc + apply(dfr[, evals[(i - 1):i]], 1, mean, na.rm = TRUE) * days
      
  }
  
  if (sum(is.na(dfr$audpc)) > 0) {
    
    cat(sum(is.na(dfr$audpc)), 'entries cannot be computed\n')
    
  }
    
  total.days <- as.numeric(dates[nd] - dates[1])
  
  dfr$raudpc <- dfr$audpc / total.days / 100
  
  # result

  dfr
  
}
