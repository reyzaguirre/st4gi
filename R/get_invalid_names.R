#' Get column names not defined in crop ontology
#' 
#' Run \code{get_invalid_names()} after running \code{check.names()}
#'
#' Check that fieldbook factors and variables' names correspond with the names defined
#' in crop ontology \url{https://cropontology.org} and in the potato and sweetpotato
#' CIP protocols.
#' @param dfr The name of the data frame.
#' @param add Additional variables. See details.
#' @param crop \code{"auto"} for autodetection or \code{"pt"} for potato and \code{"sp"} for sweetpotato.
#' @details Type \code{pt.ont()} or \code{sp.ont()} to see the list of variables and
#' corresponding short labels and CO numbers.
#' Additional variables are checked for extreme values only.
#' @return A character vector of invalid column names
#' @author Raul Eyzaguirre.
#' @examples
#' \dontrun{
#' check.names(potatoyield) |> get_invalid_names()
#' check.names(pjpz09) |> get_invalid_names()
#' }
#' @export

get_invalid_names <- function(dfr, add = NULL, crop = c('auto', 'pt', 'sp')) {
  
  crop <- match.arg(crop)
  
  if (crop == 'auto') {
    crop <- detect.crop(dfr)
    warning(crop, " crop detected", call. = FALSE)
  }
  
  factors <- c("plot", "row", "col", "rep", "block",
               "loc", "year", "season", "env",
               "geno", 'type', 'is_a_control',
               "treat", "harvest")
  
  if(crop == "pt")
    colnames.valid <- c(factors, ptont$Label)
  
  if(crop == "sp")
    colnames.valid <- c(factors, spont$Label)
  
  names.not.valid <- !(colnames(dfr) %in% colnames.valid)
  
  if (max(names.not.valid) == 1) {
    
    # Return
    
    colnames(dfr[,names.not.valid])
    
  } else {
    
    message("No invalid names")
    
  }

}
