#' Set values to \code{NA}.
#'
#' Set values to \code{NA} for selected traits.
#' @param traits List of traits (at least two).
#' @param add List of additional traits.
#' @param dfr The name of the data frame.
#' @details If all traits in \code{traits} are zero or NA, then all traits
#' in \code{traits} and \code{add} are set to \code{NA}.
#' @return It returns a data frame and a list of warnings with all the rows
#' that have been modified.
#' @author Raul Eyzaguirre
#' @examples
#' dfr <- data.frame(noph = c(2, 0, 0, 2, 1),
#'                   trt1 = c(1, 0, 0, 0, NA),
#'                   trt2 = c(3, 0, 1, NA, 2),
#'                   trt3 = c(4, 0, 5, 0, 0),
#'                   trt4 = c(1, NA, 2, 4, 5))
#' setna(c('trt1', 'trt2', 'trt3'), 'trt4', dfr)
#' @export

setna <- function(traits, add = NULL, dfr) {
  
  # Number of traits
  
  ntr <- length(traits)
  
  if (ntr == 1)
    stop("Enter at least two traits")
  
  # Condition
  
  cond <- apply(dfr[, traits] == 0 | is.na(dfr[, traits]), 1, sum) == ntr

  # Replace
  
  dfr[cond, c(traits, add)] <- NA

  if (sum(cond) > 0)
    warning("Rows set to NA: ", paste0(rownames(dfr)[cond], " "), call. = FALSE)
  
  # Return data.frame
  
  dfr
  
}
