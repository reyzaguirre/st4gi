#' Set impossible values to NA
#'
#' Detect impossible values and set them to missing value (NA).
#' @param dfr The name of the data frame.
#' @param crop \code{"pt"} for potato or \code{"sp"} for sweetpotato.
#' @param f Factor for extreme values detection. See details.
#' @param add.con Additional continuous positive traits.
#' @param add.cat Additional 1 to k categorical traits.
#' @param k The \code{k} values.
#' @details The data frame must use the labels (lower or upper case) listed in function
#' \code{check.names.pt} for potato and in function \code{check.names.sp} for sweetpotato.
#' See \code{?check.names.pt} and \code{?check.names.sp} for details.
#' Rules are:
#' \itemize{
#'  \item Continuous positive traits cannot be negative.
#'  \item 0 to 100 percentage values cannot be less than 0 or greater than 100.
#'  \item Discrete positive traits cannot be negative and should take integer values.
#'  \item 1 to 9 categorical data cannot take out of scale values.
#'  \item 1 to k categorical data cannot take out of scale values for k defined by user.
#'  \item Beta carotene values determined by RHS color charts cannot take values
#'  different from the possible values in the RHS color chart.
#'  \item Extreme low and high values are detected using the interquartile range.
#'  The rule is to detect any value out of the interval
#'  \eqn{[Q_1 - f \times (IQR + m); Q_3 + f \times (IQR + m)]}
#'  where \code{m} is the mean. By default \code{f = 10} and cannot be less than 10.
#' }
#' @return It returns the data frame with all impossible values set to \code{NA}.
#' @author Raul Eyzaguirre.
#' @examples
#' mad.na(pjpz09, "sp")
#' @importFrom stats IQR quantile rstandard
#' @export

mad.na <- function(dfr, crop, f = 10, add.con = NULL, add.cat = NULL, k = NULL) {
  
  # Check f
  
  if (f < 10)
    f <- 10
  
  # Run fix
  
  if (crop == 'pt')
    mad.na.pt(dfr, f, add.con, add.cat, k)
  
  if (crop == 'sp')
    mad.na.sp(dfr, f, add.con, add.cat, k)

}
