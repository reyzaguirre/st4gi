#' Data consistency visualization
#'
#' This function produces a visualization of data consistency.
#' @param x An object of class \code{checkdata}.
#' @param ... Additional plot arguments.
#' @details It produces a visual matrix representation of data consistency problems.
#' @return It returns a plot.
#' @author Raul Eyzaguirre.
#' @examples
#' checks <- check.data.sp(pjpz09)
#' plot(checks)
#' @importFrom graphics image
#' @export

plot.st4gi_dc <- function(x, ...) {
  
  x <- x$Inconsist.Matrix
  
  image(t(x), col = c("grey", "red"), xaxt = "n", yaxt = "n")
  
}
