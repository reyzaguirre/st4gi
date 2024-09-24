#' Check fieldbook factors and variables names for potato and sweetpotato
#'
#' Check that fieldbook factors and variables' names correspond with the names defined
#' in crop ontology \url{https://cropontology.org} and in the potato and sweetpotato
#' CIP protocols. It also checks that all variables are stored as numeric.
#' @param dfr The name of the data frame.
#' @param add Additional variables. See details.
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

check.names <- function(dfr, add = NULL, crop = c('auto', 'pt', 'sp')) {
  
  # Match arguments
  
  crop = match.arg(crop)
  
  if (crop == 'auto') {
    crop <- detect.names(dfr)
    warning(crop, " crop detected", call. = FALSE)
  }
    
  # Valid names for factors
  
  plot.id <- c("plot", "row", "col")
  
  factors <- c("loc", "year", "season", "env", "geno", 'type', "rep", "block",
               "treat", "harvest", 'is_a_control')
  
  # Valid names for variables
  
  if (crop == 'pt')
    vars <- c(ptont$Label, "nmtci", "nmtcii", "mtwci", "mtwcii",
                "fwts", "dwts", "fwts1", "fwts2", "dwts1", "dwts2",
                "dm1", "dm2", tolower(add))
  
  if (crop == 'sp')
    vars <- c(spont$Label, 'fcol.cc', tolower(add))
  
  # Valid names for factors and variables
  
  colnames.valid <- c(plot.id, factors, vars)
  
  # Factors and variables in field book (original names)
    
  colnames.fb <- colnames(dfr)

  # Convert all fieldbook names to lower case (except CO numbers)
  
  cond <- !substring(colnames(dfr), 1, 7) %in% c('CO_330:', 'CO_331:') & substring(colnames(dfr), 1, 5) != 'COMP:'
  colnames(dfr)[cond] <- tolower(colnames(dfr))[cond]
  
  if (sum(colnames.fb != colnames(dfr)) > 0)
    warning("Some labels converted to lower case", call. = FALSE)

  # Solve synonyms for factors

  old.geno <- c("accession_name", "cipno", "cip.number", 'genotype', "instn")
  new.geno <- rep('geno', length(old.geno))
  
  old.names.f <- c('plot_number', 'location', 'replication', "rep_number", "block_number", "row_number", "col_number", old.geno)
  new.names.f <- c('plot',        'loc',      'rep',         "rep",        "block",        "row",        "col",        new.geno)

  change.names.f <- NULL
  
  for (i in 1:length(old.names.f)) {
    if (exists(old.names.f[i], dfr) & !exists(new.names.f[i], dfr)) {
      change.names.f <- c(change.names.f, old.names.f[i])
      colnames(dfr)[colnames(dfr) == old.names.f[i]] <- new.names.f[i]
    }
  }  
  
  if (!is.null(change.names.f)) {
    change.names.list <- old.names.f %in% change.names.f
    warning("Factors' names ", list(old.names.f[change.names.list]), " changed to ", list(new.names.f[change.names.list]), call. = FALSE)
  }
  
  # Solve synonyms for variables
  
  if (crop == 'pt') {
    old.names.t <- c("mwt", "mwmt", "stfw", "stdw", "pdm", 'avdm', "protein", 'chipping')
    new.names.t <- c("atw", "atmw", "sfw",  "sdw",  "dm",  'dm',   "pro",     'chip_color')
  }
  
  if (crop == 'sp') {
    old.names.t <- c("trwd",  "biomd",  "cythaaj",  "rythaaj",  "dmryaj",  'vwd',  "fythaaj",  "dmvyaj",  "bythaaj",  "dmbyaj")
    new.names.t <- c("trw.d", "biom.d", "cytha.aj", "rytha.aj", "dmry.aj", 'vw.d', "fytha.aj", "dmvy.aj", "bytha.aj", "dmby.aj")
  }

  change.names.t <- NULL 
  
  for (i in 1:length(old.names.t)) {
    if (exists(old.names.t[i], dfr) & !exists(new.names.t[i], dfr)) {
      change.names.t <- c(change.names.t, old.names.t[i])
      colnames(dfr)[colnames(dfr) == old.names.t[i]] <- new.names.t[i]
    }
  }  
  
  if (!is.null(change.names.t)) {
    change.names.list <- old.names.t %in% change.names.t
    warning("Variables' names ", list(old.names.t[change.names.list]), " changed to ", list(new.names.t[change.names.list]), call. = FALSE)
  }
  
  # Names not valid
  
  names.not.valid <- !(colnames(dfr) %in% colnames.valid)

  if (max(names.not.valid) == 1)
    warning("Some columns with invalid names: ", list(colnames(dfr)[names.not.valid]), call. = FALSE)

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
  ptont[, c('Label', 'Full.Name', 'Variable.ID')]
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
  spont[, c('Label', 'Full.Name', 'Variable.ID')]
}

# Detect names automatically

detect.names <- function(dfr) {
  
  names.pt <- sum(tolower(colnames(dfr)) %in% tolower(c(ptont$Label, ptont$Variable.ID)))
  names.sp <- sum(tolower(colnames(dfr)) %in% tolower(c(spont$Label, spont$Variable.ID)))

  if (names.pt == names.sp) {
    names.pt <- names.pt / length(c(ptont$Label, ptont$Variable.ID))
    names.sp <- names.sp / length(c(spont$Label, spont$Variable.ID))
  }
  
  if (names.pt > names.sp)
    crop <- 'pt'
  if (names.pt < names.sp)
    crop <- 'sp'
  
  return(crop)
  
}
