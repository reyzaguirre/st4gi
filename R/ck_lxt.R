#' Check data for a line by tester genetic design
#'
#' This function checks the data structure for a line by tester genetic design.
#' @param line The lines.
#' @param tester The testers.
#' @param data The name of the data frame.
#' @return Three control values (\code{c1}, \code{c2}, and \code{c3},
#' the number of lines (\code{nl}), and the number of testers (\code{nt}).
#' @author Raul Eyzaguirre.
#' @details The control values are:
#' \itemize{
#'  \item \code{c1}: TRUE if all lines appear as parents.
#'  \item \code{c2}: TRUE if all testers appear as parents.
#'  \item \code{c3}: TRUE if all lines x testers are in the crosses.
#' }
#' @export

ck.lxt <- function(line, tester, data) {
  
  # Get a list of lines and testers
  
  ll <- unique(data[!is.na(data[, line]), line])
  lt <- unique(data[!is.na(data[, tester]), tester])
  
  # Get number of lines and testers
  
  nl <- length(ll)
  nt <- length(lt)
  
  # Get number of crosses
  
  temp <- data[!is.na(data[, line]) & !is.na(data[, tester]), ]
  temp$lxt <- paste(temp[, line], temp[, tester])
  
  nlxt <- length(unique(temp$lxt))
  
  # Check that all lines appear as parents
  
  temp <- data[is.na(data[, tester]), ]
  c1 <- identical(ll, unique(temp[!is.na(temp[, line]), line]))
  
  # Check that all testers appear as parents
  
  temp <- data[is.na(data[, line]), ]
  c2 <- identical(lt, unique(temp[!is.na(temp[, tester]), tester]))
  
  # Check that there are nl x nt crosses
  
  c3 <- nl * nt == nlxt
  
  # Return
  
  list(c1 = c1, c2 = c2, c3 = c3, nl = nl, nt = nt)
  
}
