#' Create potato and sweetpotato fieldbook
#'
#' Creates a potato and sweetpotato fieldbook with short labels or CO numbers.
#' @param design The name of the design data frame.
#' @param label Use \code{standard} for standard sweetpotato short labels or
#' \code{CO} for CO numbers. Default is \code{standard}.
#' @param minimal Logical, if \code{TRUE}, a minimal list of traits is included.
#' Default is \code{TRUE}.
#' @param add Additional traits to include. Only if \code{minimal = TRUE}.
#' @param computation Logical, if \code{TRUE}, computed traits are included.
#' Only if \code{minimal = FALSE}. Default is \code{FALSE}.
#' @param crop \code{"pt"} for potato and \code{"sp"} for sweetpotato.
#' @details Only labels listed in function \code{check.names} are valid.
#' Uppercase labels are converted to lowercase.
#' @return It returns a data frame with fieldbook design and traits.
#' @author Raul Eyzaguirre.
#' @examples
#' book <- cr.rcbd(1:20, 3, 10)$book
#' # Get fieldbook with minimal set of traits for potato
#' create.fb(book, crop = 'pt')
#' # Add additional traits
#' create.fb(book, add = c('fedw', 'zndw'), crop = 'pt')
#' # Get the fieldbook with CO numbers
#' create.fb(book, label = 'CO', add = c('dm', 'fedw', 'zndw'), crop = 'pt')
#' # Get fieldbook with minimal set of traits for sweetpotato
#' create.fb(book, crop = 'sp')
#' # Add additional traits
#' create.fb(book, add = c('bc', 'fe', 'zn'), crop = 'sp')
#' # Get the fieldbook with CO numbers
#' create.fb(book, label = 'CO', add = c('bc', 'fe', 'zn'), crop = 'sp')
#' @export

create.fb <- function(design, label = c("standard", "CO"), minimal = TRUE,
                      add = NULL, computation = FALSE, crop = c('pt', 'sp')) {
  
  # Match arguments
  
  label <- match.arg(label)

  crop = match.arg(crop)
  
  if (crop == 'pt') {
  
    # Minimal list of traits
    
    minimal.list <- c("ntp", "npe", "nph", "nmtp", "nnomtp", "mtwp", "nomtwp",
                      "plant_unif", "plant_vigor", "tuber_apper", "tub_unif")
    
    # Ontology
    
    ont <- ptont
    Method <- 'PotatoMethod'
    
  }
  
  if (crop == 'sp') {

    # Minimal list of traits
    
    minimal.list <- c("nops", "nope", "noph", "vir", "alt", "vv", "vw", "nopr",
                      "nocr", "nonc", "crw", "ncrw", "scol", "fcol", "rs", "rf")
    
    # Ontology
    
    ont <- spont
    Method <- 'SweetpotatoMethod'
    
  }
  
  # Additional traits
  
  if (!is.null(add)) {
    
    add <- tolower(add)
    
    minimal.list <- c(minimal.list, add[add %in% ont$Label])
  
    if (prod(add %in% ont$Label) == 0)
      warning("Some invalid names for labels: ", list(add[!add %in% ont$Label]), call. = FALSE)
    
  }
  
  # Add traits
  
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
