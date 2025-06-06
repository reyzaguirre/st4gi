#' Set values to \code{0} for potato and sweetpotato data.
#'
#' Set values to \code{0} for variables at harvest.
#' @param dfr The name of the data frame.
#' @param crop \code{"auto"} for autodetection or \code{"pt"} for potato,
#' \code{"sp"} for sweetpotato and \code{"uk"} for unknown.
#' @param checknames Logical indicating if column names should be checked, default \code{FALSE}.
#' @details This function sets values to \code{0} for all variables at harvest
#' according to some rules.
#' 
#' For potato, variables (\code{nph}, \code{nmtp}, \code{nnomtp}, \code{mtwp},
#' and \code{nomtwp}) which are \code{NA} according to the following rules:
#' \itemize{
#'  \item If \code{npe > 0} and \code{nph == 0}, then all variables are set to \code{0}.
#'  \item If \code{nmtp == 0}, then \code{mtwp} is set to \code{0}.
#'  \item If \code{mtwp == 0}, then \code{nmtp} is set to \code{0}.
#'  \item If \code{nnomtp == 0}, then \code{nomtwp} is set to \code{0}.
#'  \item If \code{nomtwp == 0}, then \code{nnomtp} is set to \code{0}.
#' }
#' For sweetpotato, variables (\code{noph}, \code{nopr}, \code{vw}, \code{nocr},
#' \code{nonc}, \code{crw}, and \code{ncrw}) which are \code{NA} according to
#' the following rules:
#' \itemize{
#'  \item If \code{nope > 0} and \code{noph == 0}, then all variables are set to \code{0}.
#'  \item If all variables are \code{0}, then \code{noph} is set to \code{0}.
#'  \item If \code{nopr == 0}, then all variables with exception of \code{vw}
#'  are set to \code{0}.
#'  \item If all variables with exception of \code{vw} are \code{0}, then
#'  \code{nopr} is set to \code{0}.
#'  \item If \code{nocr == 0}, then \code{crw} is set to \code{0}.
#'  \item If \code{crw == 0}, then \code{nocr} is set to \code{0}.
#'  \item If \code{nonc == 0}, then \code{ncrw} is set to \code{0}.
#'  \item If \code{ncrw == 0}, then \code{nonc} is set to \code{0}.
#' }
#' @return It returns a data frame and a list of warnings with all the rows
#' that have been modified.
#' @author Raul Eyzaguirre
#' @examples
#' dfr <- data.frame(nph = c(NA,  0,  3,  3,  3, 3, 3),
#'                   mtwp  = c(NA, NA, NA, NA,  8, 2, 0),
#'                   nomtwp = c(NA, NA, NA,  2, NA, 2, 4),
#'                   nmtp = c(NA, NA, NA,  0,  6, 2, 0),
#'                   nnomtp = c(NA, NA, NA,  4,  0, 3, 5))
#' setzero(dfr)
#' dfr <- data.frame(noph = c(NA,  0,  3,  3,  3, 3, 3),
#'                   nopr = c(NA,  0,  0,  2,  3, 3, 3),
#'                   vw   = c(NA, NA,  3,  2,  3, 3, 6),
#'                   crw  = c(NA, NA, NA, NA,  8, 2, 0),
#'                   ncrw = c(NA, NA, NA,  2, NA, 2, 4),
#'                   nocr = c(NA, NA, NA,  0,  6, 2, 0),
#'                   nonc = c(NA, NA, NA,  4,  0, 3, 5))
#' setzero(dfr)
#' @export

setzero <- function(dfr, crop = c('auto', 'pt', 'sp', 'uk'), checknames = FALSE) {

  # Match arguments
  
  crop = match.arg(crop)
  
  if (crop == 'auto')
    crop <- detect.crop(dfr)

  # Check names
  if (checknames)
    dfr <- check.names(dfr, crop = crop)
  
  #--------------------
  # Run potato
  #--------------------
  
  if (crop == 'pt') {
    
    # Harvest variables
    
    har <- c("nmtp", "nnomtp", "mtwp", "nomtwp")
    
    # Subset in fielddook and number of variables
    
    har <- har[har %in% colnames(dfr)]
    ntr <- length(har)
    
    # If nph == 0, then all variables are set to 0
    
    if (ntr > 0) {
      if (exists("nph", dfr)) {
        for (i in 1:length(har)) {
          if (exists("npe", dfr)) {
            cond <- dfr[, "npe"] > 0 & !is.na(dfr[, "npe"]) &
              dfr[, "nph"] == 0 & !is.na(dfr[, "nph"]) & is.na(dfr[, har[i]])
          } else {
            cond <- dfr[, "nph"] == 0 & !is.na(dfr[, "nph"]) & is.na(dfr[, har[i]])
          }
          dfr[cond, har[i]] <- 0
          if (sum(cond) > 0)
            warning("Rows with NA replaced with 0 for variable ",
                    har[i], ": ", paste0(rownames(dfr)[cond], " "), call. = FALSE)
        }
      }
    }
    
    # If nmtp == 0 then mtwp is set to 0, and if mtwp == 0 then nmtp is set to 0
    
    if (exists("nmtp", dfr) & exists("mtwp", dfr)) {
      cond <- dfr[, "nmtp"] == 0 & !is.na(dfr[, "nmtp"]) & is.na(dfr[, "mtwp"])
      dfr[cond, "mtwp"] <- 0
      if (sum(cond) > 0)
        warning("Rows with NA replaced with 0 for mtwp: ",
                paste0(rownames(dfr)[cond], " "), call. = FALSE)
      cond <- dfr[, "mtwp"] == 0 & !is.na(dfr[, "mtwp"]) & is.na(dfr[, "nmtp"])
      dfr[cond, "nmtp"] <- 0
      if (sum(cond) > 0)
        warning("Rows with NA replaced with 0 for nmtp: ",
                paste0(rownames(dfr)[cond], " "), call. = FALSE)
    }
    
    # If nnomtp == 0 then nomtwp is set to 0, and if nomtwp == 0 then nnomtp is set to 0
    
    if (exists("nnomtp", dfr) & exists("nomtwp", dfr)) {
      cond <- dfr[, "nnomtp"] == 0 & !is.na(dfr[, "nnomtp"]) & is.na(dfr[, "nomtwp"])
      dfr[cond, "nomtwp"] <- 0
      if (sum(cond) > 0)
        warning("Rows with NA replaced with 0 for nomtwp: ",
                paste0(rownames(dfr)[cond], " "), call. = FALSE)
      cond <- dfr[, "nomtwp"] == 0 & !is.na(dfr[, "nomtwp"]) & is.na(dfr[, "nnomtp"])
      dfr[cond, "nnomtp"] <- 0
      if (sum(cond) > 0)
        warning("Rows with NA replaced with 0 for nnomtp: ",
                paste0(rownames(dfr)[cond], " "), call. = FALSE)
    }
   
  }

  #--------------------
  # Run sweetpotato
  #--------------------
  
  if (crop == 'sp') {

    # Harvest variables
    
    har <- c("vw", "nocr", "nonc", "crw", "ncrw")
    
    # Subset in fielddook and number of variables
    
    har <- har[har %in% colnames(dfr)]
    ntr <- length(har)
    
    # noph = NA and nopr = NA conditions for nocr and crw
    
    if (exists("nocr", dfr) & !exists("crw", dfr))
      cr.cond <- dfr[, "nocr"] == 0 & !is.na(dfr[, "nocr"])
    
    if (!exists("nocr", dfr) & exists("crw", dfr))
      cr.cond <- dfr[, "crw"] == 0 & !is.na(dfr[, "crw"])
    
    if (exists("nocr", dfr) & exists("crw", dfr))
      cr.cond <- dfr[, "nocr"] == 0 & !is.na(dfr[, "nocr"]) & dfr[, "crw"] == 0 & !is.na(dfr[, "crw"])
    
    # noph = NA and nopr = NA conditions for nonc and ncrw
    
    if (exists("nonc", dfr) & !exists("ncrw", dfr))
      ncr.cond <- dfr[, "nonc"] == 0 & !is.na(dfr[, "nonc"])
    
    if (!exists("nonc", dfr) & exists("ncrw", dfr))
      ncr.cond <- dfr[, "ncrw"] == 0 & !is.na(dfr[, "ncrw"])
    
    if (exists("nonc", dfr) & exists("ncrw", dfr))
      ncr.cond <- dfr[, "nonc"] == 0 & !is.na(dfr[, "nonc"]) & dfr[, "ncrw"] == 0 & !is.na(dfr[, "ncrw"])
    
    # If noph == 0, then all variables are set to 0
    
    if (ntr > 0) {
      if (exists("noph", dfr)) {
        for (i in 1:length(har)) {
          if (exists('nope', dfr)) {
            cond <- dfr[, "nope"] > 0 & !is.na(dfr[, "nope"]) &
              dfr[, "noph"] == 0 & !is.na(dfr[, "noph"]) & is.na(dfr[, har[i]])
          } else {
            cond <- dfr[, "noph"] == 0 & !is.na(dfr[, "noph"]) & is.na(dfr[, har[i]])
          }
          dfr[cond, har[i]] <- 0
          if (sum(cond) > 0)
            warning("Rows with NA replaced with 0 for variable ",
                    har[i], ": ", paste0(rownames(dfr)[cond], " "), call. = FALSE)
        }
      }
    }
    
    # If all == 0, then noph is set to 0
    
    if (exists("noph", dfr) & (exists("nocr", dfr) | exists("crw", dfr)) &
        (exists("nonc", dfr) | exists("ncrw", dfr)) & exists("vw", dfr)) {
      cond <- is.na(dfr[, "noph"]) & dfr[, "vw"] == 0 & !is.na(dfr[, "vw"]) & cr.cond & ncr.cond
      dfr[cond, "noph"] <- 0
      if (sum(cond) > 0)
        warning("Rows with NA replaced with 0 for noph: ",
                paste0(rownames(dfr)[cond], " "), call. = FALSE)
    }
    
    # Subset in fieldbook and number of variables without vw
    
    har <- har[har != "vw"]
    ntr <- length(har)
    
    # If nopr == 0, then all variables with exception of vw are set to 0
    
    if (ntr > 0) {
      if (exists("nopr", dfr)) {
        for (i in 1:length(har)) {
          cond <- dfr[, "nopr"] == 0 & !is.na(dfr[, "nopr"]) & is.na(dfr[, har[i]])
          dfr[cond, har[i]] <- 0
          if (sum(cond) > 0)
            warning("Rows with NA replaced with 0 for variable ",
                    har[i], ": ", paste0(rownames(dfr)[cond], " "), call. = FALSE)
        }
      }
    }
    
    # If all == 0, then nopr is set to 0
    
    if (exists("nopr", dfr) & (exists("nocr", dfr) | exists("crw", dfr)) &
        (exists("nonc", dfr) | exists("ncrw", dfr))) {
      cond <- is.na(dfr[, "nopr"]) & cr.cond & ncr.cond
      dfr[cond, "nopr"] <- 0
      if (sum(cond) > 0)
        warning("Rows with NA replaced with 0 for nopr: ",
                paste0(rownames(dfr)[cond], " "), call. = FALSE)
    }
    
    # If nocr == 0 then crw is set to 0, and if crw == 0 then nocr is set to 0
    
    if (exists("nocr", dfr) & exists("crw", dfr)) {
      cond <- dfr[, "nocr"] == 0 & !is.na(dfr[, "nocr"]) & is.na(dfr[, "crw"])
      dfr[cond, "crw"] <- 0
      if (sum(cond) > 0)
        warning("Rows with NA replaced with 0 for crw: ",
                paste0(rownames(dfr)[cond], " "), call. = FALSE)
      cond <- dfr[, "crw"] == 0 & !is.na(dfr[, "crw"]) & is.na(dfr[, "nocr"])
      dfr[cond, "nocr"] <- 0
      if (sum(cond) > 0)
        warning("Rows with NA replaced with 0 for nocr: ",
                paste0(rownames(dfr)[cond], " "), call. = FALSE)
    }
    
    # If nonc == 0 then ncrw is set to 0, and if ncrw == 0 then nonc is set to 0
    
    if (exists("nonc", dfr) & exists("ncrw", dfr)) {
      cond <- dfr[, "nonc"] == 0 & !is.na(dfr[, "nonc"]) & is.na(dfr[, "ncrw"])
      dfr[cond, "ncrw"] <- 0
      if (sum(cond) > 0)
        warning("Rows with NA replaced with 0 for ncrw: ",
                paste0(rownames(dfr)[cond], " "), call. = FALSE)
      cond <- dfr[, "ncrw"] == 0 & !is.na(dfr[, "ncrw"]) & is.na(dfr[, "nonc"])
      dfr[cond, "nonc"] <- 0
      if (sum(cond) > 0)
        warning("Rows with NA replaced with 0 for nonc: ",
                paste0(rownames(dfr)[cond], " "), call. = FALSE)
    }
    
  }
  
  # Return data frame
  
  dfr
  
}
