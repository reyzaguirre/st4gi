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
#' @return It returns a data frame with the CIP number, cultivar name and
#' breeder code of all the genotypes identified.
#' @author Raul Eyzaguirre.
#' @examples
#' check.genos(pjpz09, "geno")
#' @export

check.genos <- function(dfr, geno, crop = c('auto', 'pt', 'sp', 'uk')) {
  
  # Match arguments
  
  crop = match.arg(crop)

  # Count number of coincidences
  
  tmp <- gcat[gcat$crop == 'pt', ]
  pt.ind <- sapply(dfr[, geno], function(x) sum(grepl(x, c(tmp$cultivar_name, tmp$breeder_code), ignore.case = TRUE)) > 0)
  names.pt <- sum(pt.ind)
  tmp <- gcat[gcat$crop == 'sp', ]
  sp.ind <- sapply(dfr[, geno], function(x) sum(grepl(x, c(tmp$cultivar_name, tmp$breeder_code), ignore.case = TRUE)) > 0)
  names.sp <- sum(sp.ind)
  
  if (crop == 'auto' & names.pt + names.sp > 0) {
    
    if (names.pt == names.sp) {
      names.pt <- names.pt / dim(gcat[gcat$crop == 'potato', ])[1]
      names.sp <- names.sp / dim(gcat[gcat$crop == 'sweetpotato', ])[1]
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
    
    if (crop == 'pt') {
      genos <- unique(dfr[pt.ind, geno])
      output <- gcat[gcat$crop == 'pt', ]
    }
    
    if (crop == 'sp') {
      genos <- unique(dfr[sp.ind, geno])
      output <- gcat[gcat$crop == 'sp', ]
    }
    
    if (crop == 'uk') {
      genos <- unique(dfr[pt.ind | sp.ind, geno])
      output <- gcat
    }
    
    tmp1 <- sapply(genos, function(x) grepl(x, output$cultivar_name, ignore.case = TRUE))
    tmp2 <- sapply(genos, function(x) grepl(x, output$breeder_code, ignore.case = TRUE))
    tmp <- cbind(tmp1, tmp2)
    tmp <- apply(tmp, 1, function(x) sum(x) > 0)
    
    output[tmp, ]
    
  } else {
    
    message("No genotypes detected with cultivar name or breeder code.")
    
  }

}
