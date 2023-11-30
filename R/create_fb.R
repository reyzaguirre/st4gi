#' Create sweetpotato fieldbook
#'
#' Creates a sweetpotato fieldbook with short labels.
#' @param design The name of the design data frame.
#' @param label Use \code{standard} for standard sweetpotato short labels or
#' \code{CO} for CO numbers. Default is \code{standard}.
#' @param minimal Logical, if \code{TRUE}, a minimal list of traits is included.
#' Default is \code{TRUE}.
#' @param add Additional traits to include. Only if \code{minimal = TRUE}.
#' @param computation Logical, if \code{TRUE}, computed traits are included.
#' Only if \code{minimal = FALSE}. Default is \code{FALSE}. 
#' @details Only labels listed in function \code{check.names.sp} are valid.
#' Uppercase labels are converted to lowercase.
#' @return It returns a data frame with fieldbook design and traits.
#' @author Raul Eyzaguirre.
#' @examples
#' book <- cr.rcbd(1:20, 3, 10)$book
#' # Get fieldbook with minimal set of traits
#' create.fb.sp(book)
#' # Add additional traits
#' create.fb.sp(book, add = c('bc', 'fe', 'zn'))
#' # Get the fieldbook with CO numbers
#' create.fb.sp(book, label = 'CO', add = c('bc', 'fe', 'zn'))
#' @export

create.fb.sp <- function(design, label = c("standard", "CO"),
                         minimal = TRUE, add = NULL, computation = FALSE) {
  
  # Match arguments
  
  label <- match.arg(label)
  
  # Minimal list of traits
  
  minimal.list <- c("nops", "nope", "noph", "vir", "alt", "vv", "vw", "nopr",
                    "nocr", "nonc", "crw", "ncrw", "scol", "fcol", "rs", "rf")
  
  # Additional traits
  
  if (!is.null(add)) {
    
    add <- tolower(add)
    
    minimal.list <- c(minimal.list, add[add %in% spont$Label])
  
    if (prod(add %in% spont$Label) == 0)
      warning("Some invalid names for labels: ", list(add[!add %in% spont$Label]), call. = FALSE)
  }
  
  # Add traits
  
  if (minimal == TRUE)
    add <- spont[spont$Label %in% minimal.list, "Label"]
  
  if (minimal == FALSE) 
    if (computation == FALSE)
      add <- spont[spont$SweetpotatoMethod != "Computation", "Label"]
  
  if (minimal == FALSE) 
    if (computation == TRUE)
      add <- spont[, "Label"]
  
  design[, add] <- NA

  # Choose labels

  if (label == 'CO')
    suppressWarnings(design <- convert.co.sp(design, 'labels.to.co'))
  
  # Return

  design
  
}

#' Create potato fieldbook
#'
#' Creates a potato fieldbook with short labels.
#' @param design The name of the design data frame.
#' @param label Use \code{standard} for standard potato short labels or
#' \code{CO} for CO numbers. Default is \code{standard}. 
#' @param minimal Logical, if \code{TRUE}, a minimal list of traits is included.
#' Default is \code{TRUE}.
#' @param add Additional traits to include. Only if \code{minimal = TRUE}.
#' @param computation Logical, if \code{TRUE}, computed traits are included.
#' Only if \code{minimal = FALSE}. Default is \code{FALSE}. 
#' @details Only labels listed in function \code{check.names.pt} are valid.
#' Uppercase labels are converted to lowercase.
#' @return It returns a data frame with fieldbook design and traits.
#' @author Raul Eyzaguirre.
#' @examples
#' book <- cr.rcbd(1:20, 3, 10)$book
#' # Get fieldbook with minimal set of traits
#' create.fb.pt(book)
#' # Add additional traits
#' create.fb.pt(book, add = c('bc', 'fe', 'zn'))
#' # Get the fieldbook with CO numbers
#' create.fb.pt(book, label = 'CO', add = c('dm', 'fedw', 'zndw'))
#' @export

create.fb.pt <- function(design, label = c("standard", "CO"),
                         minimal = TRUE, add = NULL, computation = FALSE) {
  
  # Match arguments
  
  label <- match.arg(label)
  
  # Minimal list of traits
  
  minimal.list <- c("ntp", "npe", "nph", "nmtp", "nnomtp", "mtwp", "nomtwp",
                    "plant_unif", "plant_vigor", "tuber_apper", "tub_unif")
  
  # Additional traits
  
  if (!is.null(add)) {
    
    add <- tolower(add)
  
    minimal.list <- c(minimal.list, add[add %in% ptont$Label])
    
    if (prod(add %in% ptont$Label) == 0)
      warning("Some invalid names for labels: ", list(add[!add %in% ptont$Label]), call. = FALSE)
  }

  # Add traits
  
  if (minimal == TRUE)
    add <- ptont[ptont$Label %in% minimal.list, "Label"]
  
  if (minimal == FALSE) 
    if (computation == FALSE)
      add <- ptont[ptont$PotatoMethod != "Computation", "Label"]
  
  if (minimal == FALSE) 
    if (computation == TRUE)
      add <- ptont[, "Label"]
  
  design[, add] <- NA
  
  # Choose labels
  
  if (label == 'CO')
    suppressWarnings(design <- convert.co.pt(design, 'labels.to.co'))

  # Return
  
  design
  
}
