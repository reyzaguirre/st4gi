#' Set values to NA or zero.
#'
#' Set values to NA or zero for selected traits.
#' @param traits List of traits.
#' @param add List of additional traits.
#' @param dfr The name of the data frame.
#' @details This function sets values to NA or zero with the following rules:
#' \itemize{
#'  \item If all traits in \code{traits} are zero or NA, then all traits in \code{traits}
#'  and \code{add} are set to NA, and if \code{noph} exists and it is zero or NA,
#'  then it is set to zero.
#'  \item If there is at least one trait in \code{traits} with data and \code{noph} is zero,
#'  then \code{noph} is set to NA.
#'  }
#' @return It returns a data frame.
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

  # 1. All traits = 0 or NA
  
  cond <- apply(dfr[, traits] == 0 | is.na(dfr[, traits]), 1, sum) == ntr

  dfr[cond, c(traits, add)] <- NA
  
  if (exists('noph', dfr))
    dfr[cond & is.na(dfr[, 'noph']), 'noph'] <- 0
  
  # 2. At least one trait with data
  
  if (exists('noph', dfr))
    dfr[!cond & dfr[, 'noph'] == 0 & !is.na(dfr[, 'noph']), 'noph'] <- NA
  
  # return data.frame
    
  dfr

}
