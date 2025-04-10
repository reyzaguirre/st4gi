#' Set values to \code{NA} for potato and sweetpotato data.
#'
#' Detect impossible values for potato and sweetpotato data and set them to
#' missing value (\code{NA}) according to some rules.
#' @param dfr The name of the data frame.
#' @param f Factor for extreme values detection. See details.
#' @param crop \code{"auto"} for autodetection or \code{"pt"} for potato and \code{"sp"} for sweetpotato.
#' @param checknames Logical indicating if column names should be checked, default \code{TRUE}.
#' @details The data frame must use the labels (lower or upper case) listed in
#' functions \code{ptont()} and \code{spont()}.
#' 
#' Values are set to \code{NA} with the following rules:
#' \itemize{
#'  \item All values out of the variable scale are set to \code{NA}.
#'  \item Extreme low and high values are detected using the interquartile range.
#'  The rule is to detect any value out of the interval
#'  \eqn{[Q_1 - f \times (m/3 + IQR); Q_3 + f \times (m/3 + IQR)]} where \code{m}
#'  is the mean. By default \code{f = 10} and if less than 10 a warning is shown.
#'  Values out of this range are set to \code{NA}.
#'  \item If the number of plants established is \code{0} and there is some data for
#'  any variable, then the number of plants established is set to \code{NA}.
#'  \item If the number of plants established is \code{0} and there is no data but
#'  some variables that are \code{0}, then all those variables are set to \code{NA}.
#'  \item If the number of plants harvested is \code{0} and there is some data for
#'  any non-pre-harvest variable, then the number of plants harvested is set to \code{NA}.
#'  \item If the number of plants with roots or tubers is \code{0} and there is
#'  some data for any variable evaluated with roots or tubers, then the number of
#'  plants with roots or tubers is set to \code{NA}.
#'  \item If the number of plants with roots or tubers is \code{> 0} and all
#'  variables at harvest with roots are 0, then the non-commercial yields are
#'  set to \code{NA}.  
#'  \item If the number of commercial roots or tubers is \code{0} and the
#'  commercial root or tuber weight is \code{> 0}, then the number of commercial
#'  roots or tubers is set to \code{NA}.
#'  \item If the number of commercial roots or tubers is \code{> 0} and the
#'  commercial root or tuber weight is \code{0}, then the commercial root or 
#'  tuber weight is set to \code{NA}.
#'  \item If the number of non-commercial roots or tubers is \code{0} and the
#'  non-commercial root or tuber weight is \code{> 0}, then the number of commercial
#'  roots or tubers is set to \code{NA}.
#'  \item If the number of non-commercial roots or tubers is \code{> 0} and the
#'  non-commercial root or tuber weight is \code{0}, then the non-commercial root
#'  or tuber weight is set to \code{NA}.
#' }
#' @return It returns the data frame with all impossible values set to \code{NA}
#' and a list of warnings with all the rows that have been modified.
#' @author Raul Eyzaguirre.
#' @examples
#' dfr <- data.frame(mtwp = c(2.2, 5.0, 3.6, 12, 1600, -4, 0),
#'                   dm = c(21, 23, 105, 24, -3, 30, NA),
#'                   nmtp = c(1.3, 10, 11, NA, 2, 5, NA))
#' setna(dfr)
#' dfr <- data.frame(trw = c(2.2, 5.0, 3.6, 12, 1600, -4),
#'                   dm = c(21, 23, 105, 24, -3, 30),
#'                   tnr = c(1.3, 10, 11, NA, 2, 5),
#'                   scol = c(1, 0, 15, 5, 4, 7),
#'                   fcol.cc = c(1, 15, 12, 24, 55, 20))
#' setna(dfr)
#' @importFrom stats IQR quantile
#' @export

setna <- function(dfr, f = 10, crop = c('auto', 'pt', 'sp'), checknames = TRUE) {
  
  # Match arguments
  
  crop = match.arg(crop)
  
  if (crop == 'auto') {
    crop <- detect.crop(dfr)
    warning(crop, " crop detected", call. = FALSE)
  }
  
  # Check f
  
  if (f < 10)
    warning("f < 10 can lead to delete true values", call. = FALSE)

  # Check names
  
  if (checknames)
    dfr <- check.names(dfr, crop = crop)
  
  # March arguments for both crops

  if (crop == 'pt') {
    ont <- pt_ont
    nope <- 'npe'
    noph <- 'nph'
    nopr <- 'npt'
    nocr <- 'nmtp'
    nonc <- 'nnomtp'
    crw <- 'mtwp'
    ncrw <- 'nomtwp'
  }

  if (crop == 'sp') {
    ont <- sp_ont
    nope <- 'nope'
    noph <- 'noph'
    nopr <- 'nopr'
    nocr <- 'nocr'
    nonc <- 'nonc'
    crw <- 'crw'
    ncrw <- 'ncrw'
  }

  #----------------------------------------------------------------------------
  # Impossible values (out of scale)
  #----------------------------------------------------------------------------
  
  # Impossible values for numerical variables
    
  all.vars <- ont[ont$Scale %in% c('continuous', 'integer'), ]
    
  for (i in 1:dim(all.vars)[1])
      
    if (exists(all.vars$Label[i], dfr)) {
        
      # Check minimum
      cond <- dfr[, all.vars$Label[i]] < all.vars$Minimum[i] & !is.na(dfr[, all.vars$Label[i]])
        
      # Check maximum
      if (!is.na(all.vars$Maximum[i]))
        cond <- cond | dfr[, all.vars$Label[i]] > all.vars$Maximum[i] & !is.na(dfr[, all.vars$Label[i]])
        
      # Check non integer
      if (all.vars$Scale[i] == 'integer')
        cond <- cond | dfr[, all.vars$Label[i]] %% 1 > 0 & !is.na(dfr[, all.vars$Label[i]])
        
      # Replace and print warning
      dfr[cond, all.vars$Label[i]] <- NA
      if (sum(cond) > 0)
        warning("Rows with out of scale values replaced with NA for variable ",
                all.vars$Label[i], ": ", paste0(rownames(dfr)[cond], " "), call. = FALSE)
        
    }

  # Impossible values for categorical variables
    
  all.vars <- ont[ont$Scale == 'categorical', ]
    
  for (i in 1:dim(all.vars)[1])
      
    if (exists(all.vars$Label[i], dfr)) {
        
      values <- all.vars$Values[i]
      values <- as.numeric(strsplit(values, '/')[[1]])
        
      # Check values
      cond <- !(dfr[, all.vars$Label[i]] %in% values) & !is.na(dfr[, all.vars$Label[i]])
        
      # Replace and print warning
      dfr[cond, all.vars$Label[i]] <- NA
      if (sum(cond) > 0)
        warning("Rows with out of scale values replaced with NA for variable ",
                all.vars$Label[i], ": ", paste0(rownames(dfr)[cond], " "), call. = FALSE)
        
    }
    
  #----------------------------------------------------------------------------
  # Extreme values (almost impossible)
  #----------------------------------------------------------------------------
  
  all.vars <- ont[ont$Scale %in% c('continuous', 'integer') & ont$Evaluation != 'plant', ]

  for (i in 1:dim(all.vars)[1])
      
    if (exists(all.vars$Label[i], dfr)) {
        
      m <- mean(dfr[dfr[, all.vars$Label[i]] != 0, all.vars$Label[i]], na.rm = TRUE)
      if (!is.nan(m)) {
        q1 <- quantile(dfr[, all.vars$Label[i]], 0.25, na.rm = TRUE)
        q3 <- quantile(dfr[, all.vars$Label[i]], 0.75, na.rm = TRUE)
        tol <- (m / 3 + IQR(dfr[, all.vars$Label[i]], na.rm = TRUE))
        cond1 <- dfr[, all.vars$Label[i]] < q1 - f * tol & !is.na(dfr[, all.vars$Label[i]])
        cond2 <- dfr[, all.vars$Label[i]] > q3 + f * tol & !is.na(dfr[, all.vars$Label[i]])
        cond <- cond1 | cond2
        dfr[cond, all.vars$Label[i]] <- NA
        if (sum(cond) > 0)
          warning("Rows with extreme values replaced with NA for variable ",
                  all.vars$Label[i], ": ", paste0(rownames(dfr)[cond], " "), call. = FALSE)
      }
        
    }
    
  #----------------------------------------------------------------------------
  # 1. If nope = 0 and there is some data for any variable,
  # then nope is set to NA
  # 2. If nope = 0 and there is no data but some variables that are 0,
  # then all those variables are set to NA
  #----------------------------------------------------------------------------
  
  all.vars <- ont[ont$Evaluation %in% c('plant', 'root', 'vine'), ]
  all.vars <- all.vars[all.vars$Label %in% colnames(dfr), ]
    
  if (dim(all.vars)[1] > 0 & exists(nope, dfr)) {
      
    # nope = 0 and some data
      
    if (dim(all.vars)[1] == 1)
      cond <- dfr[, all.vars$Label] > 0 & !is.na(dfr[, all.vars$Label]) & dfr[, nope] == 0 & !is.na(dfr[, nope])
    if (dim(all.vars)[1] > 1)
      cond <- apply(dfr[, all.vars$Label] > 0 & !is.na(dfr[, all.vars$Label]), 1, sum) > 0 & dfr[, nope] == 0 & !is.na(dfr[, nope])
    
    dfr[cond, nope] <- NA
    if (sum(cond) > 0)
      warning("Rows with 0 replaced with NA for variable ", nope, ": ",
              paste0(rownames(dfr)[cond], " "), call. = FALSE)
    
    # nope = 0 and no data
      
    if (dim(all.vars)[1] == 1)
      cond <- dfr[, all.vars$Label] == 0 & !is.na(dfr[, all.vars$Label]) & dfr[, nope] == 0 & !is.na(dfr[, nope])
    if (dim(all.vars)[1] > 1)
      cond <- apply(dfr[, all.vars$Label] == 0 & !is.na(dfr[, all.vars$Label]), 1, sum) > 0 & dfr[, nope] == 0 & !is.na(dfr[, nope])
      
    if (sum(cond) > 0)
      for (i in 1:dim(all.vars)[1]) {
        cond.tmp <- dfr[, all.vars$Label[i]] == 0 & !is.na(dfr[, all.vars$Label[i]])
        cond2 <- cond & cond.tmp
        dfr[cond2, all.vars$Label[i]] <- NA
        if (sum(cond2) > 0)
          warning("Rows with 0 replaced with NA for variable ",
                  all.vars$Label[i], ": ", paste0(rownames(dfr)[cond2], " "), call. = FALSE)
      }
      
  }
  
  #----------------------------------------------------------------------------
  # If noph = 0 and there is some data for any variable,
  # then noph is set to NA
  #----------------------------------------------------------------------------
  
  all.vars <- ont[ont$Evaluation %in% c(nopr, 'root', 'vine'), ]
  all.vars <- all.vars[all.vars$Label %in% colnames(dfr), ]
  
  if (dim(all.vars)[1] > 0 & exists(noph, dfr)) {
    if (dim(all.vars)[1] == 1)
      cond <- dfr[, all.vars$Label] > 0 & !is.na(dfr[, all.vars$Label]) & dfr[, noph] == 0 & !is.na(dfr[, noph])
    if (dim(all.vars)[1] > 1)
      cond <- apply(dfr[, all.vars$Label] > 0 & !is.na(dfr[, all.vars$Label]), 1, sum) > 0 & dfr[, noph] == 0 & !is.na(dfr[, noph])
    dfr[cond, noph] <- NA
    if (sum(cond) > 0)
      warning("Rows with 0 replaced with NA for variable ", noph, ": ",
              paste0(rownames(dfr)[cond], " "), call. = FALSE)
  }
  
  #----------------------------------------------------------------------------
  # If nopr = 0 and there is some data for any variable,
  # then nopr is set to NA
  #----------------------------------------------------------------------------
  
  all.vars <- ont[ont$Evaluation %in% 'root', ]
  all.vars <- all.vars[all.vars$Label %in% colnames(dfr), ]
  
  if (dim(all.vars)[1] > 0 & exists(nopr, dfr)) {
    if (dim(all.vars)[1] == 1)
      cond <- dfr[, all.vars$Label] > 0 & !is.na(dfr[, all.vars$Label]) & dfr[, nopr] == 0 & !is.na(dfr[, nopr])
    if (dim(all.vars)[1] > 1)
      cond <- apply(dfr[, all.vars$Label] > 0 & !is.na(dfr[, all.vars$Label]), 1, sum) > 0 & dfr[, nopr] == 0 & !is.na(dfr[, nopr])
    dfr[cond, nopr] <- NA
    if (sum(cond) > 0)
      warning("Rows with 0 replaced with NA for variable ", nopr, ": ",
              paste0(rownames(dfr)[cond], " "), call. = FALSE)
  }
    
  #----------------------------------------------------------------------------
  # Conditions for nocr and crw and for nonc and ncrw
  #----------------------------------------------------------------------------
  
  # Commercial roots
  
  if (exists(nocr, dfr) & !exists(crw, dfr))
    cr.cond <- dfr[, nocr] == 0 & !is.na(dfr[, nocr])
    
  if (!exists(nocr, dfr) & exists(crw, dfr))
    cr.cond <- dfr[, crw] == 0 & !is.na(dfr[, crw])
    
  if (exists(nocr, dfr) & exists(crw, dfr))
    cr.cond <- dfr[, nocr] == 0 & !is.na(dfr[, nocr]) & dfr[, crw] == 0 & !is.na(dfr[, crw])
    
  # Non-commercial roots
    
  if (exists(nonc, dfr) & !exists(ncrw, dfr)) {
    ncr.cond <- dfr[, nonc] == 0 & !is.na(dfr[, nonc])
    ncr.variables <- nonc
  }
    
  if (!exists(nonc, dfr) & exists(ncrw, dfr)) {
    ncr.cond <- dfr[, ncrw] == 0 & !is.na(dfr[, ncrw])
    ncr.variables <- ncrw
  }
    
  if (exists(nonc, dfr) & exists(ncrw, dfr)) {
    ncr.cond <- dfr[, nonc] == 0 & !is.na(dfr[, nonc]) & dfr[, ncrw] == 0 & !is.na(dfr[, ncrw])
    ncr.variables <- c(nonc, ncrw)
  }
    
  #----------------------------------------------------------------------------
  # If nopr > 0 and all variables (commercial and noncommercial) are 0, 
  # then non-commercial yields are set to NA
  # If some root > 0 and all variables (commercial and noncommercial) are 0, 
  # then non-commercial yields are set to NA
  #----------------------------------------------------------------------------
  
  if (exists(nopr, dfr) & (exists(nocr, dfr) | exists(crw, dfr)) & (exists(nonc, dfr) | exists(ncrw, dfr))) {
    cond <- dfr[, nopr] > 0 & !is.na(dfr[, nopr]) & cr.cond & ncr.cond
    dfr[cond, ncr.variables] <- NA
    if (sum(cond) > 0)
      warning("Rows with 0 replaced with NA for variables ", nonc, " and ", ncrw, ": ",
              paste0(rownames(dfr)[cond], " "), call. = FALSE)
  }

  all.vars <- ont[ont$Evaluation %in% 'root' & !ont$Evaluation %in% c(nocr, crw, nonc, ncrw), ]
  all.vars <- all.vars[all.vars$Label %in% colnames(dfr), ]
  
  if (dim(all.vars)[1] > 0 & (exists(nonc, dfr) | exists(ncrw, dfr))) {
    if (dim(all.vars)[1] == 1)
      cond <- dfr[, all.vars$Label] > 0 & !is.na(dfr[, all.vars$Label]) & ncr.cond
    if (dim(all.vars)[1] > 1)
      cond <- apply(dfr[, all.vars$Label] > 0 & !is.na(dfr[, all.vars$Label]), 1, sum) > 0 & ncr.cond
    dfr[cond, ncr.variables] <- NA
    if (sum(cond) > 0)
      warning("Rows with 0 replaced with NA for variables ", nonc, " and ", ncrw, ": ",
              paste0(rownames(dfr)[cond], " "), call. = FALSE)
  }

  #----------------------------------------------------------------------------
  # 1. If nocr = 0 and crw > 0, then nocr is set to NA
  # 2. If nocr > 0 and crw = 0, then crw is set to NA
  #----------------------------------------------------------------------------

  if (exists(nocr, dfr) & exists(crw, dfr)) {
    
    cond <- dfr[, nocr] == 0 & !is.na(dfr[, nocr]) & dfr[, crw] > 0 & !is.na(dfr[, crw])
    dfr[cond, nocr] <- NA
    if (sum(cond) > 0)
      warning("Rows with 0 replaced with NA for variable ", nocr, ": ",
              paste0(rownames(dfr)[cond], " "), call. = FALSE)
    
    cond <- dfr[, nocr] > 0 & !is.na(dfr[, nocr]) & dfr[, crw] == 0 & !is.na(dfr[, crw])
    dfr[cond, crw] <- NA
    if (sum(cond) > 0)
      warning("Rows with 0 replaced with NA for variable ", crw, ": ",
              paste0(rownames(dfr)[cond], " "), call. = FALSE)
    
  }
    
  #----------------------------------------------------------------------------
  # 1. If nonc = 0 and ncrw > 0, then nonc is set to NA
  # 2. If nonc > 0 and ncrw = 0, then ncrw is set to NA
  #----------------------------------------------------------------------------
    
  if (exists(nonc, dfr) & exists(ncrw, dfr)) {
    
    cond <- dfr[, nonc] == 0 & !is.na(dfr[, nonc]) & dfr[, ncrw] > 0 & !is.na(dfr[, ncrw])
    dfr[cond, nonc] <- NA
    if (sum(cond) > 0)
      warning("Rows with 0 replaced with NA for variable ", nonc, ": ",
              paste0(rownames(dfr)[cond], " "), call. = FALSE)
    
    cond <- dfr[, nonc] > 0 & !is.na(dfr[, nonc]) & dfr[, ncrw] == 0 & !is.na(dfr[, ncrw])
    dfr[cond, ncrw] <- NA
    if (sum(cond) > 0)
      warning("Rows with 0 replaced with NA for variable ", ncrw, ": ",
              paste0(rownames(dfr)[cond], " "), call. = FALSE)
  }
    
  # Return data frame
  
  dfr
  
}
