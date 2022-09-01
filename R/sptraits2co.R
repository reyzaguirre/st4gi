#' Sweetpotato CO numbers
#'
#' Converts sweetpotato short labels to Crop Ontology variable numbers.
#' @param dfr The name of the data frame.
#' @details All labels (lower or upper case) listed in function \code{check.names.sp}
#' are converted to CO variable numbers.
#' @return It returns a data frame with all short labels for traits converted to
#' Crop Ontology variable numbers.
#' @author Raul Eyzaguirre.
#' @examples
#' spt2co(pjpz09)
#' @export

spt2co <- function(dfr) {
  
  # Check names
  
  dfr <- check.names.sp(dfr)
    
  # Convert
  
  id <- match(colnames(dfr), spont$Label)
  colnames(dfr)[!is.na(id)] <- spont$Variable.ID[id[!is.na(id)]]
  
  # Return
  
  dfr
  
}
