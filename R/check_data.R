#' Check consistency for potato and sweetpotato experimental data
#'
#' Set of rules to check for consistency of potato and sweetpotato experimental data.
#' Data labels must be defined as specified in
#' \url{http://www.cropontology.org/ontology/CO_330/Potato} 
#' and in the Procedures for the evaluation and analysis of sweetpotato trials,
#' ISBN 978-92-9060-522-5, for potato and sweetpotato respectively.
#' @param dfr The name of the data frame.
#' @param crop \code{'pt'} for potato or \code{'sp'} for sweetpotato.
#' @param f Factor for extreme values detection. See details.
#' @param out.mod Statistical model for outliers' detection. See details.
#' @param out.max Threshold for outliers' detection.
#' @param add Additional quantitative traits.
#' @details The data frame must use the labels (lower or upper case) listed in function
#' \code{check.names.pt} for potato and in function \code{check.names.sp} for sweetpotato.
#' See \code{?check.names.pt} and \code{?check.names.sp} for details.
#' 
#' Extreme low and high values are detected using the interquartile range.
#' The rule is to detect any value out of the interval 
#' \eqn{[Q_1 - f \times IQR; Q_3 + f \times IQR]}. By default \code{f = 3}.
#' 
#' Outliers are detected based on standardized residuals for some statistical models.
#' Options are \code{'rcbd'} and \code{'met'} for a randomized complete block design
#' and a multi environment trial with RCBD in each environment. By default the threshold
#' value is \code{out.max = 4}.
#' @return It returns all rows with some kind of inconsistency or outliers.
#' @author Raul Eyzaguirre, Johan Ninanya.
#' @examples
#' check.data(potatoyield, "pt")
#' check.data(pjpz09, "sp")
#' @importFrom stats IQR quantile rstandard
#' @export

check.data <- function(dfr, crop, f = 3, out.mod = c("none", "rcbd", "met"),
                       out.max = 4, add = NULL) {
  
  # Match arguments

  out.mod = match.arg(out.mod)
  
  # Run check
  
  if (crop == 'pt')
    check.data.pt(dfr, f, out.mod, out.max, add)
  
  if (crop == 'sp')
    check.data.sp(dfr, f, out.mod, out.max, add)

}
