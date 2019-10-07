#' Set values to \code{NA} or \code{0}
#'
#' Set values to \code{NA} or \code{0} according to some rules.
#' @param dfr The name of the data frame.
#' @details Consider the following groups of traits:
#' \itemize{
#'  \item \code{rhar} (traits evaluated with roots on harvest): \code{nopr},
#'  \code{nocr}, \code{nonc}, \code{crw}, \code{ncrw}, \code{tnr}, and \code{trw}.
#'
#'  \item \code{rpos} (traits evaluated with roots post-harvest): \code{scol},
#'  \code{fcol}, \code{fcol.cc}, \code{rs}, \code{rf}, \code{damr}, \code{rspr},
#'  \code{wed}, \code{dmf}, \code{dmd}, \code{dm}, \code{fraw}, \code{fraw1},
#'  \code{suraw}, \code{suraw1}, \code{straw}, \code{straw1}, \code{coof},
#'  \code{coof1}, \code{coosu}, \code{coosu1}, \code{coost}, \code{coost1},
#'  \code{coot}, \code{coot1}, \code{cooap}, \code{cooap1}, \code{fraw2},
#'  \code{suraw2}, \code{straw2}, \code{coof2}, \code{coosu2}, \code{coost2},
#'  \code{coot2}, \code{cooap2}, \code{prot}, \code{fe}, \code{zn}, \code{ca},
#'  \code{mg}, \code{bc}, \code{bc.cc}, \code{tc}, \code{star}, \code{fruc},
#'  \code{gluc}, \code{sucr}, and \code{malt}.
#'  
#'  \item \code{vpre} (traits evaluated with vines pre-harvest): \code{vir},
#'  \code{vir1}, \code{vir2}, \code{alt}, \code{alt1}, \code{alt2}, and \code{vv}.
#'  
#'  \item \code{vpos} (traits evaluated with vines post-harvest): \code{dmvf},
#'  \code{dmvd}, and \code{dmv}.
#' }
#' Values are set to \code{NA} or \code{0} with the following rules:
#' \itemize{
#' 
#'  \item If \code{nope == 0} and there is some data for \code{vpre}, then
#'  \code{nope} is set to \code{NA}.
#'  
#'  \item If \code{noph == 0} and there is some data for \code{rhar}, then
#'  \code{noph} is set to \code{NA}.
#' }
#' @return It returns a data frame.
#' @author Raul Eyzaguirre
#' @examples
#' dfr <- data.frame(nope = c(2, 0, 2, 2, 1),
#'                   noph = c(1, 0, 2, 0, NA),
#'                   vv = c(3, 2, 1, NA, 2),
#'                   crw = c(4, 0, 5, 3, 0),
#'                   trt4 = c(1, NA, 2, 4, 5))
#' setna(dfr)
#' @export

setna <- function(dfr) {
  
  # Groups of traits
  
  rhar <- c("nopr", "nocr", "nonc", "crw", "ncrw", "tnr", "trw")
  
  rpos <- c("scol", "fcol", "fcol.cc", "rs", "rf", "damr", "rspr", "wed", "dmf",
            "dmd", "dm", "fraw", "fraw1", "suraw", "suraw1", "straw", "straw1",
            "coof", "coof1", "coosu", "coosu1", "coost", "coost1", "coot", "coot1",
            "cooap", "cooap1", "fraw2", "suraw2", "straw2", "coof2", "coosu2",
            "coost2", "coot2", "cooap2", "prot", "fe", "zn", "ca", "mg", "bc",
            "bc.cc", "tc", "star", "fruc", "gluc", "sucr", "malt")
  
  vpre <- c("vir", "vir1", "vir2", "alt", "alt1", "alt2", "vv")
  
  vpos <- c("dmvf", "dmvd", "dmv")
  
  # Subset in fieldook
  
  rhar <- rhar[rhar %in% colnames(dfr)]
  rpos <- rpos[rpos %in% colnames(dfr)]
  vpre <- vpre[vpre %in% colnames(dfr)]
  vpos <- vpos[vpos %in% colnames(dfr)]
  
  # nope and vpre
  
  if (length(vpre) > 0 & exists("nope", dfr)) {
    if (length(vpre) == 1)
      cond <- dfr[, vpre] %in% 1:9 & !is.na(dfr[, vpre])
    if (length(vpre) > 1)
      cond <- apply(dfr[, vpre] %in% 1:9 & !is.na(dfr[, vpre]), 1, sum) > 0
    dfr[cond & dfr[, 'nope'] == 0 & !is.na(dfr[, 'nope']), 'nope'] <- NA
  }
    
  # noph and rhar
  
  if (length(rhar) > 0 & exists("noph", dfr)) {
    if (length(rhar) == 1)
      cond <- dfr[, rhar] > 0 & !is.na(dfr[, rhar])
    if (length(rhar) > 1)
      cond <- apply(dfr[, rhar] > 0 & !is.na(dfr[, rhar]), 1, sum) > 0
    dfr[cond & dfr[, 'noph'] == 0 & !is.na(dfr[, 'noph']), 'noph'] <- NA
  }

  # return data.frame
    
  dfr

}
