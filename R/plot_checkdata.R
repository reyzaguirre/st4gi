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
#' checks <- check.data.pt(potatoyield)
#' plot(checks)
#' @importFrom graphics image
#' @export

plot.st4gi_dc <- function(x, ...) {
  
  x <- x$Inconsist.Matrix
  
  nr <- nrow(x)
  nc <- ncol(x)
  nr1 <- nr - 1
  nc1 <- nc - 1
  x <- x[nr:1, ]
  i <- 0:nr1
  j <- 0:nc1
  
  image(t(x), col = c("grey", "red"), xlab = 'Variable', ylab = 'Plot',
        axes = FALSE, xaxt = "n", yaxt = "n")
  
  axis(2, at = 1 - i/nr1, labels = as.character(1:nr), col = 'gray', las = 1)
  axis(3, at = j/nc1, labels = colnames(x), col = 'gray', las = 2)

}
