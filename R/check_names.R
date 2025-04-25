#' Check fieldbook factors and variables names for potato and sweetpotato
#'
#' Check that fieldbook factors and variables' names correspond with the names defined
#' in crop ontology \url{https://cropontology.org} and in the potato and sweetpotato
#' CIP protocols.
#' @param dfr The name of the data frame.
#' @param crop \code{"auto"} for autodetection or \code{"pt"} for potato and
#' \code{"sp"} for sweetpotato.
#' @details The following list of factors are recognized:
#' \itemize{
#' \item \code{plot}: The plot number.
#' \item \code{row}: The row number in the field.
#' \item \code{col}: The columnt number in the field.
#' \item \code{rep}: The replication number.
#' \item \code{block}: The block number.
#' \item \code{loc}: The location name.
#' \item \code{year}: The year.
#' \item \code{season}: The season.
#' \item \code{env}: The environment name.
#' \item \code{geno}: The genotype name.
#' \item \code{is_a_control}: Identification variable for control checks (1 for 
#' checks and NA for test genotypes).
#' \item \code{type}: The genotype type (clone, entry, parent, control, etc).
#' \item \code{treat}: The treatment name.
#' }
#' Type \code{ptont()} or \code{spont()} to see the list of variables and
#' corresponding short labels and CO numbers.
#' @return It returns:
#' \itemize{
#' \item The fieldbook data frame with all column names in lowercase and
#' with some possible modifications in the names.
#' \item A list of warnings for all the column names that have been changed.
#' \item A list of warnings for all the column names not recognized.
#' }
#' @author Raul Eyzaguirre.
#' @examples
#' check.names(potatoyield)
#' check.names(pjpz09)
#' @export

check.names <- function(dfr, crop = c('auto', 'pt', 'sp')) {
  
  # Match arguments
  
  crop = match.arg(crop)
  
  if (crop == 'auto') {
    crop <- detect.crop(dfr)
    warning(crop, " crop detected", call. = FALSE)
  }
    
  # Valid names for factors
  
  factors <- c("plot", "row", "col", "rep", "block",
               "loc", "year", "season", "env",
               "geno", 'is_a_control', 'type', "treat")
  
  # Valid names for variables
  
  if (crop == 'pt')
    vars <- pt_ont$Label
  
  if (crop == 'sp')
    vars <- sp_ont$Label
  
  # Valid names for factors and variables
  
  colnames.valid <- c(factors, vars)
  
  # Factors and variables in fieldbook (original names)
    
  colnames.fb <- colnames(dfr)

  # Convert all fieldbook names to lower case (except CO numbers)
  
  cond <- !substring(colnames(dfr), 1, 7) %in% c('CO_330:', 'CO_331:') & substring(colnames(dfr), 1, 5) != 'COMP:'
  colnames(dfr)[cond] <- tolower(colnames(dfr))[cond]
  
  if (sum(colnames.fb != colnames(dfr)) > 0)
    warning("Some labels converted to lower case", call. = FALSE)

  # Solve synonyms for factors

  old.geno <- c("accession_name", "cipno", "cipn", "cip.number", "genotype", "instn", "clon", "clone")
  new.geno <- rep('geno', length(old.geno))
  
  old.rep <- c("replication", "rep_number")
  new.rep <- rep('rep', length(old.rep))
  
  old.row <- c("row_number", "fila")
  new.row <- rep('row', length(old.row))
  
  old.col <- c("col_number", "column", "columna", "range")
  new.col <- rep('col', length(old.col))

  old.names.factors <- c('plot_number', 'location', "block_number", old.geno, old.rep, old.row, old.col)
  new.names.factors <- c('plot',        'loc',      "block",        new.geno, new.rep, new.row, new.col)

  changed.names.factors <- NULL
  
  for (i in 1:length(old.names.factors)) {
    if (exists(old.names.factors[i], dfr) & !exists(new.names.factors[i], dfr)) {
      changed.names.factors <- c(changed.names.factors, old.names.factors[i])
      colnames(dfr)[colnames(dfr) == old.names.factors[i]] <- new.names.factors[i]
    }
  }  
  
  if (!is.null(changed.names.factors)) {
    cond <- old.names.factors %in% changed.names.factors
    warning("Factors' names ", list(old.names.factors[cond]), " changed to ", list(new.names.factors[cond]), call. = FALSE)
  }
  
  # Solve synonyms for variables
  
  if (crop == 'pt') {
    old.names.vars <- c("mwt", "mwmt", "stfw", "stdw", "pdm", 'avdm', "protein", 'chipping')
    new.names.vars <- c("atw", "atmw", "sfw",  "sdw",  "dm",  'dm',   "pro",     'chip_color')
  }
  
  if (crop == 'sp') {
    old.names.vars <- c("trwd",  "biomd",  "cythaaj",  "rythaaj",  "dmryaj",  'vwd',  "fythaaj",  "dmvyaj",  "bythaaj",  "dmbyaj")
    new.names.vars <- c("trw.d", "biom.d", "cytha.aj", "rytha.aj", "dmry.aj", 'vw.d', "fytha.aj", "dmvy.aj", "bytha.aj", "dmby.aj")
  }

  changed.names.vars <- NULL 
  
  for (i in 1:length(old.names.vars)) {
    if (exists(old.names.vars[i], dfr) & !exists(new.names.vars[i], dfr)) {
      changed.names.vars <- c(changed.names.vars, old.names.vars[i])
      colnames(dfr)[colnames(dfr) == old.names.vars[i]] <- new.names.vars[i]
    }
  }  
  
  if (!is.null(changed.names.vars)) {
    cond <- old.names.vars %in% changed.names.vars
    warning("Variables' names ", list(old.names.vars[cond]), " changed to ", list(new.names.vars[cond]), call. = FALSE)
  }
  
  # Names not valid
  
  cond <- !(colnames(dfr) %in% colnames.valid)

  if (max(cond) == 1)
    warning("Some columns with invalid names: ", list(colnames(dfr)[cond]), call. = FALSE)

  # Return
  
  dfr
  
}

#' Potato ontology
#'
#' Lists all variables in potato ontology used by \code{st4gi} package.
#' 
#' @return It returns a data frame.
#' @author Raul Eyzaguirre.
#' @examples
#' ptont()
#' @export
 
ptont <- function() {
  pt_ont[, c('ID', 'Label', 'Name')]
}

#' Sweetpotato ontology
#'
#' Lists all variables in sweetpotato ontology used by \code{st4gi} package.
#' 
#' @return It returns a data frame.
#' @author Raul Eyzaguirre.
#' @examples
#' spont()
#' @export

spont <- function() {
  sp_ont[, c('ID', 'Label', 'Name')]
}

#' Get column names not defined in crop ontology
#' 
#' Run \code{get_invalid_names()} after running \code{check.names()}
#'
#' Check that fieldbook factors and variables' names correspond with the names defined
#' in crop ontology \url{https://cropontology.org} and in the potato and sweetpotato
#' CIP protocols.
#' @param dfr The name of the data frame.
#' @param add Additional variables.
#' @param crop \code{"auto"} for autodetection or \code{"pt"} for potato and \code{"sp"} for sweetpotato.
#' @details Type \code{ptont()} or \code{spont()} to see the list of variables and
#' corresponding short labels and CO numbers.
#' @return A character vector of invalid column names
#' @author Raul Eyzaguirre.
#' @examples
#' \dontrun{
#' tmp <- check.names(potatoyield)
#' get.invalid.names(tmp)
#' tmp <- check.names(pjpz09)
#' get.invalid.names(tmp)
#' }
#' @export

get.invalid.names <- function(dfr, add = NULL, crop = c('auto', 'pt', 'sp')) {
  
  crop <- match.arg(crop)
  
  if (crop == 'auto') {
    crop <- detect.crop(dfr)
    warning(crop, " crop detected", call. = FALSE)
  }
  
  factors <- c("plot", "row", "col", "rep", "block",
               "loc", "year", "season", "env",
               "geno", 'is_a_control', 'type', "treat")
  
  if(crop == "pt")
    colnames.valid <- c(factors, pt_ont$Label)
  
  if(crop == "sp")
    colnames.valid <- c(factors, sp_ont$Label)
  
  names.not.valid <- !(colnames(dfr) %in% colnames.valid)
  
  if (max(names.not.valid) == 1) {
    
    # Return
    
    colnames(dfr[,names.not.valid])
    
  } else {
    
    message("No invalid names")
    
  }
  
}

# Detect crop automatically
  
detect.crop <- function(dfr) {
    
  # Remove extra text
    
  colnames(dfr) <- gsub('.*CO_330.', 'CO_330:', colnames(dfr))
  colnames(dfr) <- gsub('.*CO_331.', 'CO_331:', colnames(dfr))
  colnames(dfr) <- gsub('.*COMP.', 'COMP:', colnames(dfr))
    
  # Count number of coincidences
    
  names.pt <- sum(tolower(colnames(dfr)) %in% tolower(c(pt_ont$Label, pt_ont$ID)))
  names.sp <- sum(tolower(colnames(dfr)) %in% tolower(c(sp_ont$Label, sp_ont$ID)))
    
  if (names.pt == names.sp) {
    names.pt <- names.pt / length(c(pt_ont$Label, pt_ont$ID))
    names.sp <- names.sp / length(c(sp_ont$Label, sp_ont$ID))
  }
    
  if (names.pt > names.sp)
    crop <- 'pt'
  if (names.pt < names.sp)
    crop <- 'sp'
    
  return(crop)
    
}
