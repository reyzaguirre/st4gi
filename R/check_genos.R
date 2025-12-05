#' Check genotype names for potato and sweetpotato
#'
#' Check that genotype's names correspond to CIP numbers instead of variety names
#' or breeder codes.
#' @param dfr The name of the data frame.
#' @param geno The name of the column in \code{dfr} that identifies the genotypes.
#' @param crop \code{"auto"} for autodetection or \code{"pt"} for potato,
#' \code{"sp"} for sweetpotato and \code{"uk"} for unknown.
#' @details It checks for genotypes with cultivar names or breeder codes that
#' have a CIP number.
#' @return It returns a list with all the genotypes with cultivar name or
#' breeder code.
#' @author Raul Eyzaguirre.
#' @examples
#' check.genos(pjpz09, "geno")
#' @export

check.genos <- function(dfr, geno, crop = c('auto', 'pt', 'sp', 'uk')) {
  
  # Match arguments
  
  crop = match.arg(crop)

  # Count number of coincidences
  
  tmp <- germcat[germcat$crop == 'pt', ]
  pt.ind <- tolower(dfr[, geno]) %in% tolower(c(tmp$cultivar_name, tmp$breeder_code))
  names.pt <- sum(pt.ind)
  tmp <- germcat[germcat$crop == 'sp', ]
  sp.ind <- tolower(dfr[, geno]) %in% tolower(c(tmp$cultivar_name, tmp$breeder_code))
  names.sp <- sum(sp.ind)
  
  if (crop == 'auto' & names.pt + names.sp > 0) {
    
    if (names.pt == names.sp) {
      names.pt <- names.pt / dim(germcat[germcat$crop == 'potato', ])[1]
      names.sp <- names.sp / dim(germcat[germcat$crop == 'sweetpotato', ])[1]
    }
      
    if (names.pt > names.sp) {
      crop <- 'pt'
      warning("pt crop detected", call. = FALSE)
    }
    if (names.pt < names.sp) {
      crop <- 'sp'
      warning("sp crop detected", call. = FALSE)
    }
    if (names.pt == names.sp) {
      crop <- 'uk'
      warning("Unknown crop", call. = FALSE)
    }
      
  }
  
  if (names.pt + names.sp > 0) {
    if (crop == 'pt')
      index <- pt.ind
    if (crop == 'sp')
      index <- sp.ind
    if (crop == 'uk')
      index <- pt.ind | sp.ind
    message("Some genotypes with cultivar name or breeder code: ", list(unique(dfr[index, geno])))
  } else {
    message("No genotypes detected with cultivar name or breeder code.")
  }

}
