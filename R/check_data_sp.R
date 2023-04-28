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
#' @param format Output format as \code{"plain.text"} or \code{"data.frame"}.
#' @details The data frame must use the labels (lower or upper case) listed in
#' function \code{check.names.sp}.
#' 
#' Extreme low and high values are detected using the interquartile range.
#' The rule is to detect any value out of the interval 
#' \eqn{[Q_1 - f \times IQR; Q_3 + f \times IQR]}. By default \code{f = 5}.
#' 
#' Outliers are detected based on standardized residuals for some statistical
#' models. Options are \code{"rcbd"} and \code{"met"} for a randomized complete
#' block design and a multi environment trial with RCBD in each environment.
#' By default the threshold value is \code{out.max = 4}.
#' 
#' @return It returns all rows with some kind of inconsistency or outliers.
#' @author Raul Eyzaguirre.
#' @examples
#' check.data.sp(pjpz09)
#' @importFrom stats IQR quantile rstandard
#' @export

check.data.sp <- function(dfr, f = 5, out.mod = c("none", "rcbd", "met"),
                          out.max = 4, add = NULL,
                          format = c("plain.text", "data.frame")) {
  
  # Match arguments
  
  out.mod = match.arg(out.mod)
  format = match.arg(format)
  
  # Run check
  
  dfr.out <- rules.sp(dfr, f, out.mod, out.max, add, format)
  
  if (format == 'data.frame') {
    rownames(dfr.out) <- 1:dim(dfr.out)[1]
    dfr.out
  }
  
}
