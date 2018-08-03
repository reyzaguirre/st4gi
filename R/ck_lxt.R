#' Check data for a line by tester genetic design
#'
#' This function checks the data structure for a line by tester genetic design.
#' @param line The lines.
#' @param tester The testers.
#' @param dfr The name of the data frame.
#' @return Three control values, the number of lines (\code{nlin}),
#' and the number of testers (\code{ntes}). The control values are:
#' \itemize{
#'  \item \code{c1}: TRUE if all lines appear as parents.
#'  \item \code{c2}: TRUE if all testers appear as parents.
#'  \item \code{c3}: TRUE if all lines x testers are in the crosses.
#' }
#' @author Raul Eyzaguirre.
#' @examples
#' ck.lxt("line", "tester", lxt)
#' @export

ck.lxt <- function(line, tester, dfr) {
  
  # Get names and number of lines, testers, and crosses
  
  out <- ck.fs(c("line", "tester"), NULL, dfr)
  
  nlin <- out$nl[1]
  ntes <- out$nl[2]
  nlxt <- out$nt
  llin <- out$lf[[1]]
  ltes <- out$lf[[2]]

  # Check that all lines appear as parents
  
  temp <- dfr[is.na(dfr[, tester]), ]
  c1 <- identical(llin, unique(temp[!is.na(temp[, line]), line]))
  
  # Check that all testers appear as parents
  
  temp <- dfr[is.na(dfr[, line]), ]
  c2 <- identical(ltes, unique(temp[!is.na(temp[, tester]), tester]))
  
  # Check that there are nl x nt crosses
  
  c3 <- as.numeric(nlin * ntes) == nlxt
  
  # Return
  
  list(c1 = c1, c2 = c2, c3 = c3, nlin = nlin, ntes = ntes)
  
}
