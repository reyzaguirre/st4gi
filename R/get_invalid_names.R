#' get column names not defined in crop ontology
#' 
#'run \code{get_invalid_names()} after running \code{check.names()}
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
#' @return a character vector of invalid column names
#' @author Raul Eyzaguirre.
#' @examples
#' \dontrun{
#' check.names(potatoyield) |> get_invalid_names()
#' check.names(pjpz09) |> get_invalid_names()
#' }
#' @export

get_invalid_names <- function(dfr, add = NULL, crop = c('auto', 'pt', 'sp')) {
  
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
  
  # colnames.fb <- colnames(dfr)
  
  # Convert all fieldbook names to lower case (except CO numbers)
  
  # cond <- !substring(colnames(dfr), 1, 7) %in% c('CO_330:', 'CO_331:') & substring(colnames(dfr), 1, 5) != 'COMP:'
  # colnames(dfr)[cond] <- tolower(colnames(dfr))[cond]
  # 
  # if (sum(colnames.fb != colnames(dfr)) > 0)
    # warning("Some labels converted to lower case", call. = FALSE)
  
  # Solve synonyms for factors
  
  # old.geno <- c("accession_name", "cipno", "cip.number", 'genotype', "instn")
  # new.geno <- rep('geno', length(old.geno))
  # 
  # old.names.f <- c('plot_number', 'location', 'replication', "rep_number", "block_number", "row_number", "col_number", old.geno)
  # new.names.f <- c('plot',        'loc',      'rep',         "rep",        "block",        "row",        "col",        new.geno)
  # 
  # change.names.f <- NULL
  # 
  # for (i in 1:length(old.names.f)) {
  #   if (exists(old.names.f[i], dfr) & !exists(new.names.f[i], dfr)) {
  #     change.names.f <- c(change.names.f, old.names.f[i])
  #     colnames(dfr)[colnames(dfr) == old.names.f[i]] <- new.names.f[i]
  #   }
  # }  
  # 
  # if (!is.null(change.names.f)) {
  #   change.names.list <- old.names.f %in% change.names.f
  #   # warning("Factors' names ", list(old.names.f[change.names.list]), " changed to ", list(new.names.f[change.names.list]), call. = FALSE)
  # }
  
  # Solve synonyms for variables
  
  if (crop == 'pt') {
    old.names.t <- c("mwt", "mwmt", "stfw", "stdw", "pdm", 'avdm', "protein", 'chipping')
    new.names.t <- c("atw", "atmw", "sfw",  "sdw",  "dm",  'dm',   "pro",     'chip_color')
  }
  
  if (crop == 'sp') {
    old.names.t <- c("trwd",  "biomd",  "cythaaj",  "rythaaj",  "dmryaj",  'vwd',  "fythaaj",  "dmvyaj",  "bythaaj",  "dmbyaj")
    new.names.t <- c("trw.d", "biom.d", "cytha.aj", "rytha.aj", "dmry.aj", 'vw.d', "fytha.aj", "dmvy.aj", "bytha.aj", "dmby.aj")
  }
  
  # change.names.t <- NULL 
  # 
  # for (i in 1:length(old.names.t)) {
  #   if (exists(old.names.t[i], dfr) & !exists(new.names.t[i], dfr)) {
  #     change.names.t <- c(change.names.t, old.names.t[i])
  #     colnames(dfr)[colnames(dfr) == old.names.t[i]] <- new.names.t[i]
  #   }
  # }  
  # 
  # if (!is.null(change.names.t)) {
  #   change.names.list <- old.names.t %in% change.names.t
  #   # warning("Variables' names ", list(old.names.t[change.names.list]), " changed to ", list(new.names.t[change.names.list]), call. = FALSE)
  # }
  
  
  # Names not valid
  colnames.valid <- c(colnames.valid, new.names.t)
  
  names.not.valid <- !(colnames(dfr) %in% colnames.valid)
  
  if (max(names.not.valid) == 1){
    
    
    # Return
    
    colnames(dfr[,names.not.valid])
    
  }else{
    message("No invalid names")
  }

}



