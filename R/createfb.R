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
#' @details Labels listed in function \code{check.names.sp} are included.
#' @return It returns a data frame with fieldbook design and traits.
#' @author Raul Eyzaguirre.
#' @examples
#' book <- cr.rcbd(1:20, 3, 10)$book
#' createfb.sp(book)
#' @export

createfb.sp <- function(design, label = c("standard", "CO"),
                        minimal = TRUE, add = NULL, computation = FALSE) {
  
  # Match arguments
  
  label <- match.arg(label)
  
  # Minimal list of traits
  
  minimal.list <- c("nops", "nope", "noph", "vir", "alt", "vv", "vw", "nopr",
                    "nocr", "nonc", "crw", "ncrw", "scol", "fcol", "rs", "rf")
  
  # Additional traits
  
  if (!is.null(add)) {
    
    minimal.list <- c(minimal.list, add[add %in% spont$Label])
  
    if (prod(add %in% spont$Label) == 0)
      warning("Some invalid names for labels: ", list(add[!add %in% spont$Label]), call. = FALSE)
  }
  
  # Choose labels
  
  if(label == "standard")
    labeltouse <- "Label" else
      labeltouse <- "Variable.ID"

  # Add traits
  
  if (minimal == TRUE)
    add <- spont[spont$Label %in% minimal.list, labeltouse]
  
  if (minimal == FALSE) 
    if (computation == FALSE)
      add <- spont[spont$SweetpotatoMethod != "Computation", labeltouse]
  
  if (minimal == FALSE) 
    if (computation == TRUE)
      add <- spont[, labeltouse]
  
  design[, add] <- NA

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
#' @details Labels listed in function \code{check.names.pt} are included.
#' @return It returns a data frame with fieldbook design and traits.
#' @author Raul Eyzaguirre.
#' @examples
#' book <- cr.rcbd(1:20, 3, 10)$book
#' createfb.pt(book)
#' @export

createfb.pt <- function(design, label = c("standard", "CO"),
                        minimal = TRUE, add = NULL, computation = FALSE) {
  
  # Match arguments
  
  label <- match.arg(label)
  
  # Minimal list of traits
  
  minimal.list <- c("ntp", "npe", "nph", "nmtp", "nnomtp", "mtwp", "nomtwp",
                    "plant_unif", "plant_vigor", "tuber_apper", "tub_unif")
  
  # Additional traits
  
  if (!is.null(add)) {
    
    minimal.list <- c(minimal.list, add[add %in% ptont$Label])
    
    if (prod(add %in% ptont$Label) == 0)
      warning("Some invalid names for labels: ", list(add[!add %in% ptont$Label]), call. = FALSE)
  }

  # Choose labels
  
  if(label == "standard")
    labeltouse <- "Label" else
      labeltouse <- "Variable.ID"
  
  # Add traits
  
  if (minimal == TRUE)
    add <- ptont[ptont$Label %in% minimal.list, labeltouse]
  
  if (minimal == FALSE) 
    if (computation == FALSE)
      add <- ptont[ptont$PotatoMethod != "Computation", labeltouse]
  
  if (minimal == FALSE) 
    if (computation == TRUE)
      add <- ptont[, labeltouse]
  
  design[, add] <- NA
  
  # Return
  
  design
  
}
