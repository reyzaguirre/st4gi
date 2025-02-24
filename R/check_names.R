#' Check fieldbook factors and variables names for potato and sweetpotato
#'
#' Check that fieldbook factors and variables' names correspond with the names defined
#' in crop ontology \url{https://cropontology.org} and in the potato and sweetpotato
#' CIP protocols. It also checks that all variables are stored as numeric.
#' @param dfr The name of the data frame.
#' @param crop \code{"auto"} for autodetection or \code{"pt"} for potato and \code{"sp"} for sweetpotato.
#' @details Type \code{pt.ont()} or \code{sp.ont()} to see the list of variables and
#' corresponding short labels and CO numbers.
#' Additional variables are checked for extreme values only.
#' @return It returns:
#' \itemize{
#' \item The fieldbook data frame with all column names in lowercase and
#' with some possible modifications in the names. Variables that are stored
#' with a non-numeric class are transformed to numeric.
#' \item A list of warnings for all the column names that have been changed.
#' \item A list of warnings for all the column names not recognized.
#' \item A list of warnings for all the column variables that have been changed to numeric.
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
               "geno", 'type', 'is_a_control',
               "treat", "harvest")
  
  # Valid names for variables
  
  if (crop == 'pt')
    vars <- c(ptont$Label, "nmtci", "nmtcii", "mtwci", "mtwcii", "fwts", "dwts",
              "fwts1", "fwts2", "dwts1", "dwts2", "dm1", "dm2")
  
  if (crop == 'sp')
    vars <- c(spont$Label, 'fcol.cc')
  
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

  old.geno <- c("accession_name", "cipno", "cip.number", 'genotype', "instn")
  new.geno <- rep('geno', length(old.geno))
  
  old.names.factors <- c('plot_number', 'location', 'replication', "rep_number", "block_number", "row_number", "col_number", old.geno)
  new.names.factors <- c('plot',        'loc',      'rep',         "rep",        "block",        "row",        "col",        new.geno)

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
    old.names.vars <- c("mwt", "mwmt", "stfw", "stdw", "pdm", 'avdm', 'dm_oven', 'dm_liof', 'dm_hyd', "protein", 'chipping')
    new.names.vars <- c("atw", "atmw", "sfw",  "sdw",  "dm",  'dm',   'dm',      'dm',      'dm',     "pro",     'chip_color')
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

  # Check variables are numeric
  
  nonumeric.list <- NULL
  
  column.class <- unlist(lapply(dfr, class))
  
  for(i in colnames(dfr)) {
    if(i %in% vars & column.class[i] != "numeric") {
      dfr[, i] <- suppressWarnings(as.numeric(as.character(dfr[, i])))
      nonumeric.list <- c(nonumeric.list, i)
    }
  }
  
  if (!is.null(nonumeric.list))
    warning("Some variables converted to numeric: ", list(nonumeric.list), call. = FALSE)
  
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
#' pt.ont()
#' @export
 
pt.ont <- function() {
  ptont[, c('Label', 'Name', 'ID')]
}

#' Sweetpotato ontology
#'
#' Lists all variables in sweetpotato ontology used by \code{st4gi} package.
#' 
#' @return It returns a data frame.
#' @author Raul Eyzaguirre.
#' @examples
#' sp.ont()
#' @export

sp.ont <- function() {
  spont[, c('Label', 'Name', 'ID')]
}

# Detect crop automatically

detect.crop <- function(dfr) {
  
  # Remove extra text
  
  colnames(dfr) <- gsub('.*CO_330.', 'CO_330:', colnames(dfr))
  colnames(dfr) <- gsub('.*CO_331.', 'CO_331:', colnames(dfr))
  colnames(dfr) <- gsub('.*COMP.', 'COMP:', colnames(dfr))

  # Count number of coincidences
  
  names.pt <- sum(tolower(colnames(dfr)) %in% tolower(c(ptont$Label, ptont$ID)))
  names.sp <- sum(tolower(colnames(dfr)) %in% tolower(c(spont$Label, spont$ID)))
  
  if (names.pt == names.sp) {
    names.pt <- names.pt / length(c(ptont$Label, ptont$ID))
    names.sp <- names.sp / length(c(spont$Label, spont$ID))
  }
  
  if (names.pt > names.sp)
    crop <- 'pt'
  if (names.pt < names.sp)
    crop <- 'sp'
  
  return(crop)
  
}
