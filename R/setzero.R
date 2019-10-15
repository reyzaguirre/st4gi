#' Set values to \code{0}.
#'
#' Set values to \code{0} for harvested traits.
#' @param dfr The name of the data frame.
#' @details This function sets values to \code{0} for all traits at harvest
#' (\code{vw}, \code{nocr}, \code{nonc}, \code{crw}, and \code{ncrw}) which are
#' \code{NA} if at least one of these traits has a value greater than \code{0}.
#' @return It returns a data frame and a list of warnings with all the rows
#' that have been modified.
#' @author Raul Eyzaguirre
#' @examples
#' dfr <- data.frame(noph = c(2, 3, 3, 2, 1),
#'                   crw = c(1, 0, 0, 0, NA),
#'                   ncrw = c(3, 0, 1, NA, 2),
#'                   nocr = c(4, 0, 5, 0, 0),
#'                   nonc = c(1, NA, 2, 4, 5))
#' setzero(dfr)
#' @export

setzero <- function(dfr) {
  
  # Harvest traits
  
  har <- c("vw", "nocr", "nonc", "crw", "ncrw")
  
  # Subset in fieldook
  
  har <- har[har %in% colnames(dfr)]

  # Condition
  
  if (length(har) > 0) {
    
    if (length(har) == 1)
      cond <- dfr[, har] > 0 & !is.na(dfr[, har])
    
    if (length(har) > 1)
      cond <- apply(dfr[, har] > 0 & !is.na(dfr[, har]), 1, sum) > 0
    
    # Replace
    
    for (i in 1:length(har)) {
      to.rep <- is.na(dfr[, har[i]]) & cond
      dfr[to.rep, har[i]] <- 0
      if (sum(to.rep) > 0)
        warning("Rows with NA replaced with 0 for trait ",
                  har[i], ": ", paste0(rownames(dfr)[to.rep], " "), call. = FALSE)
      }
  }
  
  # Return data.frame
  
  dfr
  
}
