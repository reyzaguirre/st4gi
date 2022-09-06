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
#' convert.co.sp(pjpz09)
#' @export

convert.co.sp <- function(dfr) {
  
  # Check names
  
  dfr <- check.names.sp(dfr)
    
  # Convert
  
  id <- match(colnames(dfr), spont$Label)
  colnames(dfr)[!is.na(id)] <- spont$Variable.ID[id[!is.na(id)]]
  
  # Return
  
  dfr
  
}

#' Potato CO numbers
#'
#' Converts potato short labels to Crop Ontology variable numbers.
#' @param dfr The name of the data frame.
#' @details All labels (lower or upper case) listed in function \code{check.names.pt}
#' are converted to CO variable numbers.
#' @return It returns a data frame with all short labels for traits converted to
#' Crop Ontology variable numbers.
#' @author Raul Eyzaguirre.
#' @examples
#' convert.co.pt(potatoyield)
#' @export

convert.co.pt <- function(dfr) {
  
  # Check names
  
  dfr <- check.names.pt(dfr)
  
  # Convert
  
  id <- match(colnames(dfr), ptont$Label)
  colnames(dfr)[!is.na(id)] <- ptont$Variable.ID[id[!is.na(id)]]
  
  # Return
  
  dfr
  
}