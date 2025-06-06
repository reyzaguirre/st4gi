#' CO numbers
#'
#' Interchanges potato and sweetpotato short labels with Crop Ontology variable numbers.
#' @param dfr The name of the data frame.
#' @param direction \code{labels.to.co} or \code{co.to.labels}
#' @param crop \code{"auto"} for autodetection or \code{"pt"} for potato,
#' \code{"sp"} for sweetpotato and \code{"uk"} for unknown.
#' @param checknames Logical indicating if column names should be checked, default \code{FALSE}.
#' @details All labels (lower or upper case) listed in functions \code{ptont()}
#' and \code{spont()} are interchanged with CO variable numbers.
#' @return It returns a data frame with all short labels for variables interchanged
#' with Crop Ontology variable numbers.
#' @author Raul Eyzaguirre.
#' @examples
#' convert.co(potatoyield)
#' convert.co(pjpz09)
#' @export

convert.co <- function(dfr, direction = c('labels.to.co', 'co.to.labels'),
                       crop = c('auto', 'pt', 'sp', 'uk'), checknames = FALSE) {
  
  # Match arguments
  
  direction <- match.arg(direction)

  crop = match.arg(crop)
  
  if (crop == 'auto')
    crop <- detect.crop(dfr)

  if (direction == 'labels.to.co') {
    
    # Check names
    
    if (checknames)
      dfr <- check.names(dfr, crop = crop)
    
    # Convert
    
    if (crop == 'pt') {
      id <- match(colnames(dfr), pt_ont$Label)
      colnames(dfr)[!is.na(id)] <- pt_ont$ID[id[!is.na(id)]]
    }
      
    if (crop == 'sp') {
      id <- match(colnames(dfr), sp_ont$Label)
      colnames(dfr)[!is.na(id)] <- sp_ont$ID[id[!is.na(id)]]
    }
    
    # Factors
    
    colnames(dfr)[colnames(dfr) == 'plot'] <- 'plot_number'
    colnames(dfr)[colnames(dfr) == 'geno'] <- 'accession_name'
    colnames(dfr)[colnames(dfr) == 'rep'] <- 'rep_number'
    colnames(dfr)[colnames(dfr) == 'block'] <- 'block_number'
    colnames(dfr)[colnames(dfr) == 'row'] <- 'row_number'
    colnames(dfr)[colnames(dfr) == 'col'] <- 'col_number'

  }
  
  if (direction == 'co.to.labels') {
    
    # Remove extra text

    colnames(dfr) <- sub('.*CO_330.([0-9]+).*', 'CO_330:\\1', colnames(dfr))
    colnames(dfr) <- sub('.*CO_331.([0-9]+).*', 'CO_331:\\1', colnames(dfr))
    colnames(dfr) <- sub('.*COMP.([0-9]+).*', 'COMP:\\1', colnames(dfr))

    # Obsolete values
    
    colnames(dfr)[colnames(dfr) == 'CO_331:2000036'] <- 'CO_331:0006024'
    
    # Convert
    
    if (crop == 'pt') {
      id <- match(colnames(dfr), pt_ont$ID)
      colnames(dfr)[!is.na(id)] <- pt_ont$Label[id[!is.na(id)]]
    }
    
    if (crop == 'sp') {
      id <- match(colnames(dfr), sp_ont$ID)
      colnames(dfr)[!is.na(id)] <- sp_ont$Label[id[!is.na(id)]]
    }

    
    # Factors
    
    colnames(dfr)[colnames(dfr) == 'studyYear'] <- 'year'
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
