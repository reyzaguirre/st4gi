#' Set values to NA or zero.
#'
#' Set values to NA or zero for selected traits.
#' @param traits List of traits.
#' @param vines Trait for vines.
#' @param roots Trait for roots. 
#' @param fb The name of the fieldbook data frame.
#' @author Raul Eyzaguirre
#' @details This function sets values to NA or zero with the following rules:
#' \itemize{
#'  \item If \code{noph = 0} then all traits are set to \code{NA}.
#'  \item If \code{noph > 0} then all traits with \code{NA} are set to zero.
#'  \item If \code{is.na(noph)} and \code{vines = 0} and \code{roots = 0} then set
#'  \code{noph} to zero and all traits to \code{NA}.
#'  \item If \code{is.na(noph)} and \code{vines > 0} or \code{roots > 0} then set
#'  all traits to zero.
#'  }
#' @return It returns a data frame.
#' @examples
#' traits <- c('nocr', 'nonc', 'crw', 'ncrw', 'trw', 'vw')
#' setna(traits, 'vw', 'trw', pjpz09)
#' @export

setna <- function(traits, vines = NULL, roots = NULL, fb) {
  
  # Number of traits
  
  nt <- length(traits)

  # Check there is noph on data frame
  
  if (!exists("noph", fb))
    stop("Number of plants harvested, noph, is missing.")
  
  # 1. noph = 0
  
  cond <- fb[, 'noph'] == 0 & !is.na(fb[, 'noph'])
  fb[cond, traits] <- NA
  
  # 2. noph > 0
  
  for (i in 1:nt) {
    cond <- fb[, 'noph'] > 0 & !is.na(fb[, 'noph']) & is.na(fb[, traits[i]])
    fb[cond, traits[i]] <- 0
  }
  
  # 3. is.na(noph)
  
  if (!is.null(vines) & !is.null(roots)) {
    
    c1 <- is.na(fb[, 'noph'])

    # 3.1. no vines no roots
    
    c2 <- fb[, vines] == 0 | is.na(fb[, vines])
    c3 <- fb[, roots] == 0 | is.na(fb[, roots])
    fb[c1 & c2 & c3, 'noph'] <- 0
    fb[c1 & c2 & c3, traits] <- NA
    
    # 3.2. vines or roots

    c2 <- fb[, vines] > 0 & !is.na(fb[, vines])
    c3 <- fb[, roots] > 0 & !is.na(fb[, roots])
    for (i in 1:nt)
      fb[c1 & (c2 | c3) & is.na(fb[, traits[i]]), traits[i]] <- 0
    
  }
    
  # return data.frame
    
  fb
}
