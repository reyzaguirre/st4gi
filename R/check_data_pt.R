#' Check consistency for potato experimental data
#'
#' Set of rules to check for consistency of potato experimental data.
#' Data labels must be defined as specified in
#' \url{http://www.cropontology.org/ontology/CO_330/Potato}.
#' @param dfr The name of the data frame.
#' @param f Factor for extreme values detection. See details.
#' @param out.mod Statistical model for outliers' detection. See details.
#' @param out.max Threshold for outliers' detection.
#' @param add Additional quantitative traits.
#' @details The data frame must use the labels (lower or upper case) listed in
#' function \code{check.names.pt}.
#' 
#' Extreme low and high values are detected using the interquartile range.
#' The rule is to detect any value out of the interval 
#' \eqn{[Q_1 - f \times IQR; Q_3 + f \times IQR]}. By default \code{f = 5}.
#' 
#' Outliers are detected based on standardized residuals for some statistical
#' models. Options are \code{"rcbd"} and \code{"met"} for a randomized complete
#' block design and a multi environment trial with RCBD in each environment.
#' By default the threshold value is \code{out.max = 4}.
#' @return It returns all rows with some kind of inconsistency or outliers.
#' @author Johan Ninanya, Raul Eyzaguirre.
#' @examples
#' check.data.pt(potatoyield)
#' @importFrom stats IQR quantile rstandard
#' @export

check.data.pt <- function(dfr, f = 5, out.mod = c("none", "rcbd", "met"),
                          out.max = 4, add = NULL) {
  
  # Match arguments
  
  out.mod = match.arg(out.mod)
  
  # Run check
  
  rules.pt(dfr, f, out.mod, out.max, add)

}
