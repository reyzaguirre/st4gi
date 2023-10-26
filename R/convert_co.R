#' Sweetpotato CO numbers
#'
#' Interchanges sweetpotato short labels with Crop Ontology variable numbers.
#' @param dfr The name of the data frame.
#' @param direction \code{labels.to.co} or \code{co.to.labels}
#' @details All labels (lower or upper case) listed in function \code{check.names.sp}
#' are interchanged with CO_331 variable numbers.
#' @return It returns a data frame with all short labels for traits interchanged
#' with Crop Ontology variable numbers.
#' @author Raul Eyzaguirre.
#' @examples
#' convert.co.sp(pjpz09)
#' @export

convert.co.sp <- function(dfr, direction = c('labels.to.co', 'co.to.labels')) {
  
  direction <- match.arg(direction)
  
  if (direction == 'labels.to.co') {
    
    # Check names
    
    dfr <- check.names.sp(dfr)
    
    # Convert
    
    id <- match(colnames(dfr), spont$Label)
    colnames(dfr)[!is.na(id)] <- spont$Variable.ID[id[!is.na(id)]]
    
  }
  
  if (direction == 'co.to.labels') {
    
    # Remove extra text
    
    colnames(dfr) <- gsub('.*CO_331.', 'CO_331:', colnames(dfr))
    colnames(dfr) <- gsub('.*COMP.', 'COMP:', colnames(dfr))

    # Obsolete values
    
    colnames(dfr)[colnames(dfr) == 'CO_331:2000036'] <- 'CO_331:0006024'
    
    # Convert
    
    id <- match(colnames(dfr), spont$Variable.ID)
    colnames(dfr)[!is.na(id)] <- spont$Label[id[!is.na(id)]]
    
    # Factors
    
    colnames(dfr)[colnames(dfr) == 'studyName'] <- 'trial'
    colnames(dfr)[colnames(dfr) == 'locationName'] <- 'loc'
    colnames(dfr)[colnames(dfr) == 'germplasmName'] <- 'geno'
    colnames(dfr)[colnames(dfr) == 'replicate'] <- 'rep'
    colnames(dfr)[colnames(dfr) == 'blockNumber'] <- 'block'
    colnames(dfr)[colnames(dfr) == 'plotNumber'] <- 'plot'
    colnames(dfr)[colnames(dfr) == 'rowNumber'] <- 'row'
    colnames(dfr)[colnames(dfr) == 'colNumber'] <- 'col'
    colnames(dfr)[colnames(dfr) == 'entryType'] <- 'type'

  }
  
  # Return
  
  dfr
  
}

#' Potato CO numbers
#'
#' Interchanges potato short labels with Crop Ontology variable numbers.
#' @param dfr The name of the data frame.
#' @param direction \code{labels.to.co} or \code{co.to.labels}
#' @details All labels (lower or upper case) listed in function \code{check.names.pt}
#' are interchanged with CO_330 variable numbers.
#' @return It returns a data frame with all short labels for traits interchanged
#' with Crop Ontology variable numbers.
#' @author Raul Eyzaguirre.
#' @examples
#' convert.co.pt(potatoyield)
#' @export

convert.co.pt <- function(dfr, direction = c('labels.to.co', 'co.to.labels')) {
  
  direction <- match.arg(direction)

  if (direction == 'labels.to.co') {
    
    # Check names
    
    dfr <- check.names.pt(dfr)
    
    # Convert
    
    id <- match(colnames(dfr), ptont$Label)
    colnames(dfr)[!is.na(id)] <- ptont$Variable.ID[id[!is.na(id)]]
    
  }
  
  if (direction == 'co.to.labels') {
    
    # Remove extra text
    
    colnames(dfr) <- gsub('.*CO_330.', 'CO_330:', colnames(dfr))
    colnames(dfr) <- gsub('.*COMP.', 'COMP:', colnames(dfr))
    
    # Convert
    
    id <- match(colnames(dfr), ptont$Variable.ID)
    colnames(dfr)[!is.na(id)] <- ptont$Label[id[!is.na(id)]]
    
    # Factors
    
    colnames(dfr)[colnames(dfr) == 'studyName'] <- 'trial'
    colnames(dfr)[colnames(dfr) == 'locationName'] <- 'loc'
    colnames(dfr)[colnames(dfr) == 'germplasmName'] <- 'geno'
    colnames(dfr)[colnames(dfr) == 'replicate'] <- 'rep'
    colnames(dfr)[colnames(dfr) == 'blockNumber'] <- 'block'
    colnames(dfr)[colnames(dfr) == 'plotNumber'] <- 'plot'
    colnames(dfr)[colnames(dfr) == 'rowNumber'] <- 'row'
    colnames(dfr)[colnames(dfr) == 'colNumber'] <- 'col'
    colnames(dfr)[colnames(dfr) == 'entryType'] <- 'type'
    
  }
  
  # Return
  
  dfr
  
}