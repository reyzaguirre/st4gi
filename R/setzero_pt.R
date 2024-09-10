#' Set values to \code{0} for potato data.
#'
#' Set values to \code{0} for harvested traits.
#' @param dfr The name of the data frame.
#' @details This function sets values to \code{0} for all traits at harvest
#' (\code{nph}, \code{nmtp}, \code{nnomtp}, \code{mtwp}, and \code{nomtwp})
#' which are \code{NA} according to the following rules:
#' \itemize{
#'  \item If \code{npe > 0} and \code{nph == 0}, then all traits are set to \code{0}.
#'  \item If \code{nmtp == 0}, then \code{mtwp} is set to \code{0}.
#'  \item If \code{mtwp == 0}, then \code{nmtp} is set to \code{0}.
#'  \item If \code{nnomtp == 0}, then \code{nomtwp} is set to \code{0}.
#'  \item If \code{nomtwp == 0}, then \code{nnomtp} is set to \code{0}.
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
#' setzero.pt(dfr)
#' @export

setzero.pt <- function(dfr) {
  
  # Harvest traits
  
  har <- c("nmtp", "nnomtp", "mtwp", "nomtwp")
  
  # Subset in fielddook and number of traits
  
  har <- har[har %in% colnames(dfr)]
  ntr <- length(har)

  # If nph == 0, then all traits are set to 0
  
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
          warning("Rows with NA replaced with 0 for trait ",
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
  
  # Return data frame
  
  dfr
  
}
