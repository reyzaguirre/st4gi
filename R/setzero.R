#' Set values to \code{0}.
#'
#' Set values to \code{0} for harvested traits.
#' @param dfr The name of the data frame.
#' @param on Logic values to turn on or off the rules for setting values
#' to \code{0} (see details). Default is \code{TRUE}.
#' @details This function sets values to \code{0} for all traits at harvest
#' (\code{vw}, \code{nocr}, \code{nonc}, \code{crw}, and \code{ncrw}) which
#' are \code{NA} according to the following rules:
#' \itemize{
#'  \item If \code{noph == 0}, then all traits are set to \code{0}.
#'  \item If \code{nopr == 0}, then all traits with exception of \code{vw}
#'  are set to \code{0}.
#'  \item If \code{nocr == 0}, then \code{crw} is set to \code{0}.
#'  \item If \code{nonc == 0}, then \code{ncrw} is set to \code{0}.
#'  \item If there is no information for \code{noph} and all traits are
#'  \code{NA}, then all traits are set to \code{0}. By default, this only
#'  rule is set to \code{FALSE}.
#' }
#' @return It returns a data frame and a list of warnings with all the rows
#' that have been modified.
#' @author Raul Eyzaguirre
#' @examples
#' dfr <- data.frame(noph = c(NA,  0,  3,  3,  3, 3, 3),
#'                   nopr = c(NA,  0,  0,  2,  3, 3, 3),
#'                   vw   = c(NA, NA,  3,  2,  3, 3, 6),
#'                   crw  = c(NA, NA, NA, NA,  8, 2, 0),
#'                   ncrw = c(NA, NA, NA,  2, NA, 2, 4),
#'                   nocr = c(NA, NA, NA,  0,  6, 2, 0),
#'                   nonc = c(NA, NA, NA,  4,  0, 3, 5))
#' setzero(dfr)
#' @export

setzero <- function(dfr, on = c(TRUE, TRUE, TRUE, TRUE, FALSE)) {
  
  # Check noph
  
  if (!exists("noph", dfr))
    dfr$noph.tmp <- NA
  else
    dfr$noph.tmp <- dfr$noph
    
  # Harvest traits
  
  har <- c("vw", "nocr", "nonc", "crw", "ncrw")
  
  # Subset in fielddook and number of traits
  
  har <- har[har %in% colnames(dfr)]
  ntr <- length(har)

  # Rule 1 and Rule 5
  
  if (ntr > 0) {
  
    # Rule 5: No information for noph and all traits are NA, then all traits are set to 0
    
    if (on[1]) {
      
      cond <- apply(is.na(dfr[, c("noph.tmp", har)]), 1, sum) == ntr + 1
      
      dfr[cond, har] <- 0
      
      if (sum(cond) > 0)
        warning("Rule 1: Rows with NA replaced with 0 for all traits: ",
                paste0(rownames(dfr)[cond], " "), call. = FALSE)
      
    }
    
    # Rule 1: If noph == 0, then all traits are set to 0
    
    if (on[2] & exists("noph", dfr)) {
      
      for (i in 1:length(har)) {
        
        cond <- dfr[, "noph"] == 0 & !is.na(dfr[, "noph"]) & is.na(dfr[, har[i]])
        
        dfr[cond, har[i]] <- 0
        
        if (sum(cond) > 0)
          warning("Rule 2: Rows with NA replaced with 0 for trait ",
                  har[i], ": ", paste0(rownames(dfr)[cond], " "), call. = FALSE)
        
      }
    }
    
  }
  
  # Subset in fieldbook and number of traits without vw
  
  har <- har[har != "vw"]
  ntr <- length(har)
  
  # Rule 2: If nopr == 0, then all traits with exception of vw are set to 0
  
  if (ntr > 0) {
    
    if (on[3] & exists("nopr", dfr)) {
      
      for (i in 1:length(har)) {
        
        cond <- dfr[, "nopr"] == 0 & !is.na(dfr[, "nopr"]) & is.na(dfr[, har[i]])
        
        dfr[cond, har[i]] <- 0
        
        if (sum(cond) > 0)
          warning("Rule 3: Rows with NA replaced with 0 for trait ",
                  har[i], ": ", paste0(rownames(dfr)[cond], " "), call. = FALSE)
        
      }
    }
    
  }
    
  # Rule 3: If nocr == 0 then crw is set to 0
  
  if (on[4] & exists("nocr", dfr) & exists("crw", dfr)) {
    
    cond <- dfr[, "nocr"] == 0 & !is.na(dfr[, "nocr"]) & is.na(dfr[, "crw"])
    
    dfr[cond, "crw"] <- 0
    
    if (sum(cond) > 0)
      warning("Rule 4: Rows with NA replaced with 0 for crw: ",
              paste0(rownames(dfr)[cond], " "), call. = FALSE)
    
  }
  
  # Rule 4: If nonc == 0 then ncrw is set to 0
  
  if (on[5] & exists("nonc", dfr) & exists("ncrw", dfr)) {
    
    cond <- dfr[, "nonc"] == 0 & !is.na(dfr[, "nonc"]) & is.na(dfr[, "ncrw"])
    
    dfr[cond, "ncrw"] <- 0
    
    if (sum(cond) > 0)
      warning("Rule 5: Rows with NA replaced with 0 for ncrw: ",
              paste0(rownames(dfr)[cond], " "), call. = FALSE)
    
  }
  
  # Return data frame
  
  if (exists("noph.tmp", dfr))
    dfr <- dfr[, colnames(dfr) != "noph.tmp"]
  
  dfr
  
}
