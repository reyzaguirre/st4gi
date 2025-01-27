#' get column names not defined in crop ontology
#' 
#'run \code{get_invalid_names()} after running \code{check.names()}
#'
#' Check that fieldbook factors and variables' names correspond with the names defined
#' in crop ontology \url{https://cropontology.org} and in the potato and sweetpotato
#' CIP protocols. It also checks that all variables are stored as numeric.
#' @param dfr The name of the data frame.
#' @param add Additional variables. See details.
#' @param crop \code{"auto"} for autodetection or \code{"pt"} for potato and \code{"sp"} for sweetpotato.
#' @details Type \code{pt.ont()} or \code{sp.ont()} to see the list of variables and
#' corresponding short labels and CO numbers.
#' Additional variables are checked for extreme values only.
#' @return a character vector of invalid column names
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
  crop <- detect.names(dfr)
  warning(crop, " crop detected", call. = FALSE)
}

factors <-
  c(
    "loc",
    "year",
    "season",
    "env",
    "geno",
    'type',
    "rep",
    "block",
    "treat",
    "harvest",
    'is_a_control',
    'plot',
    "row",
    "col"
  )

if(crop == "pt"){
  colnames.valid <- c(
    factors,
    pt.ont()$Label
  )
}else{
  colnames.valid <- c(
    factors,
    sp.ont()$Label
  )
}
  
  names.not.valid <- !(colnames(dfr) %in% colnames.valid)
  
  if (max(names.not.valid) == 1){
    
    
    # Return
    
    colnames(dfr[,names.not.valid])
    
  }else{
    message("No invalid names")
  }

}



