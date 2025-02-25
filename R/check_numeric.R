#' Check fieldbook variables are numeric
#'
#' Check that fieldbook variables are stored as numeric.
#' @param dfr The name of the data frame.
#' @param crop \code{"auto"} for autodetection or \code{"pt"} for potato and
#' \code{"sp"} for sweetpotato.
#' @details It checks that all variables recognized by \code{check.names}
#' are stored as numeric. Non-numeric columns are transformed to numeric
#' and the NAs values introduced by coercion are listed.
#' @return It returns the data frame with all non-numeric variables transformed
#' to numeric, and a list of all NAs values introduced by coercion.
#' @author Raul Eyzaguirre.
#' @examples
#' tmp <- check.names(potatoyield)
#' check.numeric(tmp)
#' tmp <- check.names(pjpz09)
#' check.numeric(tmp)
#' @export

check.numeric <- function(dfr, crop = c('auto', 'pt', 'sp')) {
  
  # Match arguments
  
  crop = match.arg(crop)
  
  if (crop == 'auto') {
    crop <- detect.crop(dfr)
    warning(crop, " crop detected", call. = FALSE)
  }
  
  # Check names
  
  dfr <- check.names(dfr, crop = crop)
  
  # Valid names for variables
  
  if (crop == 'pt')
    vars <- ptont$Label
  
  if (crop == 'sp')
    vars <- spont$Label
  
  # Check variables are numeric
  
  nonumeric.list <- NULL
  
  nonumeric.nas <- list()
  
  column.class <- unlist(lapply(dfr, class))
  
  for(i in colnames(dfr)) {
    if(i %in% vars & column.class[i] != "numeric") {
      tmp <- dfr[, i] 
      dfr[, i] <- suppressWarnings(as.numeric(as.character(dfr[, i])))
      nonumeric.list <- c(nonumeric.list, i)
      tmp <- data.frame(tmp, dfr[, i])
      tmp <- tmp[!is.na(tmp[, 1]) & is.na(tmp[, 2]), ]
      colnames(tmp)[1] <- i
      nonumeric.nas[[i]] <- tmp[, 1]
    }
  }
  
  if (!is.null(nonumeric.list)) {
    print('Non-numeric values detected:')
    print(nonumeric.nas)
  }
  
  # Return
  
  dfr
  
}
