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
#' @param print.text Logical, if \code{TRUE} the output is printed on screen.
#' @details The data frame must use the labels (lower or upper case) listed in
#' function \code{check.names.pt}.
#' 
#' Extreme low and high values are detected using the interquartile range.
#' The rule is to detect any value out of the interval 
#' \eqn{[Q_1 - f \times IQR; Q_3 + f \times IQR]}. By default \code{f = 5}.
#' If \code{f = 0}, the detection of extreme values is not executed.
#' 
#' Outliers are detected based on standardized residuals for some statistical
#' models. Options are \code{"rcbd"} and \code{"met"} for a randomized complete
#' block design and a multi environment trial with RCBD in each environment.
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
#' @importFrom stats IQR quantile rstandard
#' @export

check.data.pt <- function(dfr, f = 5, out.mod = c("none", "rcbd", "met"),
                          out.max = 4, add = NULL, print.text = TRUE) {

  .Deprecated("check.data")
  
  # Match arguments
  
  out.mod = match.arg(out.mod)

  # Check names
  
  dfr <- check.names.pt(dfr, add)
  if (!is.null(add))
    add <- tolower(add)
  
  # Create inconsistencies matrix
  
  im <- structure(rep(0, dim(dfr)[1] * dim(dfr)[2]),
                  .Dim = c(dim(dfr)[1], dim(dfr)[2]),
                  .Dimnames = list(NULL, names(dfr)))
  
  # Run check
  
  output <- rules.pt(dfr, im, f, out.mod, out.max, add, print.text)
  
  if (dim(output$Inconsist.List)[1] > 0) {
    rownames(output$Inconsist.List) <- 1:dim(output$Inconsist.List)[1]
  }
  
  # Output

  class(output) <- "st4gi_dc"
  invisible(output)
  
}
