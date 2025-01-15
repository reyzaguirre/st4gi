#' Create potato and sweetpotato fieldbook
#'
#' Creates a potato and sweetpotato fieldbook with short labels or CO numbers.
#' @param design The name of the design data frame.
#' @param crop \code{"pt"} for potato and \code{"sp"} for sweetpotato.
#' @param label Use \code{standard} for standard sweetpotato short labels or
#' \code{CO} for CO numbers. Default is \code{standard}.
#' @param minimal Logical, if \code{TRUE}, a minimal list of variables is included.
#' Default is \code{TRUE}.
#' @param add Additional variables to include. Only if \code{minimal = TRUE}.
#' @param computation Logical, if \code{TRUE}, computed variables are included.
#' Only if \code{minimal = FALSE}. Default is \code{FALSE}.
#' @details Only labels listed in function \code{check.names} are valid.
#' Uppercase labels are converted to lowercase.
#' @return It returns a data frame with fieldbook design and variables.
#' @author Raul Eyzaguirre.
#' @examples
#' book <- cr.rcbd(1:20, 3, 10)$book
#' # Get fieldbook with minimal set of variables for potato
#' create.fb(book, 'pt')
#' # Add additional variables
#' create.fb(book, 'pt', add = c('fedw', 'zndw'))
#' # Get the fieldbook with CO numbers
#' create.fb(book, 'pt', 'CO', add = c('dm', 'fedw', 'zndw'))
#' # Get fieldbook with minimal set of variables for sweetpotato
#' create.fb(book, 'sp')
#' # Add additional variables
#' create.fb(book, 'sp', add = c('bc', 'fe', 'zn'))
#' # Get the fieldbook with CO numbers
#' create.fb(book, 'sp', 'CO', add = c('bc', 'fe', 'zn'))
#' @export

create.fb <- function(design, crop = c('pt', 'sp'), label = c("standard", "CO"),
                      minimal = TRUE, add = NULL, computation = FALSE) {
  
  # Match arguments
  
  label <- match.arg(label)

  crop = match.arg(crop)
  
  if (crop == 'pt') {
  
    # Minimal list of variables
    
    minimal.list <- c("ntp", "npe", "nph", "nmtp", "nnomtp", "mtwp", "nomtwp",
                      "plant_unif", "plant_vigor", "tuber_apper", "tub_unif")
    
    # Ontology
    
    ont <- ptont
    Method <- 'PotatoMethod'
    
  }
  
  if (crop == 'sp') {

    # Minimal list of variables
    
    minimal.list <- c("nops", "nope", "noph", "vir", "alt", "vv", "vw", "nopr",
                      "nocr", "nonc", "crw", "ncrw", "scol", "fcol", "rs", "rf")
    
    # Ontology
    
    ont <- spont
    Method <- 'SweetpotatoMethod'
    
  }
  
  # Additional variables
  
  if (!is.null(add)) {
    
    add <- tolower(add)
    
    minimal.list <- c(minimal.list, add[add %in% ont$Label])
  
    if (prod(add %in% ont$Label) == 0)
      warning("Some invalid names for labels: ", list(add[!add %in% ont$Label]), call. = FALSE)
    
  }
  
  # Add variables
  
  if (minimal == TRUE)
    add <- ont[ont$Label %in% minimal.list, "Label"]
  
  if (minimal == FALSE) 
    if (computation == FALSE)
      add <- ont[ont[, Method] != "Computation", "Label"]
  
  if (minimal == FALSE) 
    if (computation == TRUE)
      add <- ont[, "Label"]
  
  design[, add] <- NA

  # Choose labels

  if (label == 'CO')
    suppressWarnings(design <- convert.co(design, 'labels.to.co'))
  
  # Return

  design
  
}
