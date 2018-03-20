#' Set values to NA or zero.
#'
#' Set values to NA or zero for selected traits.
#' @param traits List of traits.
#' @param fb The name of the fieldbook data frame.
#' @author Raul Eyzaguirre
#' @details This function sets values to NA or zero with the following rules:
#' \itemize{
#'  \item If all traits are zero or NA, then all traits are set to NA, and if \code{noph}
#'  exists and it is zero or NA, then it is set to zero.
#'  \item If there is at least one trait with data and \code{noph} is zero, then
#'  it is set to NA.
#'  }
#' @return It returns a data frame.
#' @examples
#' ddd <- data.frame(noph = c(2, 0, 0, 2, 1),
#'                   trt1 = c(1, 0, 0, 0, NA),
#'                   trt2 = c(3, 0, 1, NA, 2),
#'                   trt3 = c(4, 0, 5, 0, 0))
#' setna(c('trt1', 'trt2', 'trt3'), ddd)
#' @export

setna <- function(traits, fb) {
  
  # Number of traits
  
  nt <- length(traits)

  ## 1. All traits = 0 or NA
  
  cond <- apply(fb[, traits] == 0 | is.na(fb[, traits]), 1, sum) == nt

  fb[cond, traits] <- NA
  
  if (exists('noph', fb))
    fb[cond & is.na(fb[, 'noph']), 'noph'] <- 0
  
  ## 2. At least one trait with data
  
  if (exists('noph', fb))
    fb[!cond & fb[, 'noph'] == 0 & !is.na(fb[, 'noph']), 'noph'] <- NA
  
  # return data.frame
    
  fb
}
