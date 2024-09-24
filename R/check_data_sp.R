#' Check consistency for sweetpotato experimental data
#'
#' Set of rules to check for consistency of sweetpotato experimental data.
#' Data labels must be defined as specified in the Procedures for the evaluation
#' and analysis of sweetpotato trials, ISBN 978-92-9060-522-5.
#' @param dfr The name of the data frame.
#' @param f Factor for extreme values detection. See details.
#' @param out.mod Statistical model for outliers' detection. See details.
#' @param out.max Threshold for outliers' detection.
#' @param add Additional quantitative traits.
#' @param print.text Logical, if \code{TRUE} the output is printed on screen.
#' @details The data frame must use the labels (lower or upper case) listed in
#' function \code{check.names.sp}.
#' 
#' Extreme low and high values are detected using the interquartile range.
#' The rule is to detect any value out of the interval 
#' \eqn{[Q_1 - f \times IQR; Q_3 + f \times IQR]}. By default \code{f = 5}.
#' If \code{f = 0}, the detection of extreme values is not executed.
#' 
#' Outliers are detected based on standardized residuals for some statistical
#' models. Options are \code{"rcbd"} and \code{"met"} for a randomized complete
#' block design and a multienvironment trial with RCBD in each environment.
#' By default the threshold value is \code{out.max = 4}.
#' 
#' @return It returns:
#' \itemize{
#' \item \code{$Inconsist.List}, a \code{data.frame} with a list of all the
#' rows with some kind of inconsistency.
#' \item \code{$Inconsist.Matrix}, a \code{data.frame} with the positions
#' in the fieldbook data frame where inconsistencies occur. These are coded
#' with: (1) for inconsistencies among traits, (2) for out of range values,
#' (3) for extreme values or outliers.
#' }
#' @author Raul Eyzaguirre.
#' @examples
#' check.data.sp(pjpz09)
#' @importFrom stats IQR quantile rstandard
#' @export

check.data.sp <- function(dfr, f = 5, out.mod = c("none", "rcbd", "met"),
                          out.max = 4, add = NULL, print.text = TRUE) {
  
  # .Deprecated("check.data")

  # Match arguments
  
  out.mod = match.arg(out.mod)

  # Check names
  
  dfr <- check.names.sp(dfr, add)
  if (!is.null(add))
    add <- tolower(add)

  # Transform dmvf, dmvd, dmf, and dmd to kilograms
  
  if (exists("dmvf", dfr))
    dfr$dmvf <- dfr$dmvf / 1000
  if (exists("dmvd", dfr))
    dfr$dmvd <- dfr$dmvd / 1000
  if (exists("dmf", dfr))
    dfr$dmf <- dfr$dmf / 1000
  if (exists("dmd", dfr))
    dfr$dmd <- dfr$dmd / 1000
  
  # Create inconsistencies matrix
  
  im <- structure(rep(0, dim(dfr)[1] * dim(dfr)[2]),
                  .Dim = c(dim(dfr)[1], dim(dfr)[2]),
                  .Dimnames = list(NULL, names(dfr)))
  
  # Run check
  
  output <- rules.sp(dfr, im, f, out.mod, out.max, add, print.text)
  
  if (dim(output$Inconsist.List)[1] > 0) {
    rownames(output$Inconsist.List) <- 1:dim(output$Inconsist.List)[1]
  }
  
  # Output
  
  class(output) <- "st4gi_dc"
  invisible(output)
  
}
