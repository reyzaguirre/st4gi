#' Clean data
#'
#' Detect impossible values for sweetpotato data and set them to missing value
#' (\code{NA}) according to some rules.
#' 
#' @param dfr The name of the data frame.
#' @param f Factor for extreme values detection. See details.
#' 
#' @details The data frame must use the labels (lower or upper case) listed in
#' function \code{check.names.sp}; see \code{?check.names.sp} for details.
#' Consider the following groups of traits:
#' \itemize{
#'  \item \code{pre} (traits evaluated pre-harvest): \code{vir},
#'  \code{vir1}, \code{vir2}, \code{alt}, \code{alt1}, \code{alt2}, and \code{vv}.
#'  
#'  \item \code{nrt} (traits evaluated without roots): \code{vir}, \code{vir1},
#'  \code{vir2}, \code{alt}, \code{alt1}, \code{alt2}, \code{vv}, \code{vw},
#'  \code{dmvf}, \code{dmvd}, \code{dmv}, \code{vw.d}, \code{fytha}, \code{fytha.aj},
#'  \code{dmvy}, \code{dmvy.aj}, \code{vpp}, \code{vpsp}, and \code{shi}.
#'  
#'  \item \code{cnn} (continuos non-negative traits): \code{vw}, \code{crw},
#'  \code{ncrw}, \code{trw}, \code{trw.d}, \code{biom}, \code{biom.d}, \code{cytha},
#'  \code{cytha.aj}, \code{rytha}, \code{rytha.aj}, \code{dmry}, \code{dmry.aj},
#'  \code{vw.d}, \code{fytha}, \code{fytha.aj}, \code{dmvy}, \code{dmvy.aj},
#'  \code{bytha}, \code{bytha.aj}, \code{dmby}, \code{dmby.aj}, \code{nrpp},
#'  \code{nrpsp}, \code{ncrpp}, \code{ncrpsp}, \code{ypp}, \code{ypsp}, \code{vpp},
#'  \code{vpsp}, and \code{rfr}.
#'  
#'  \item \code{cpo} (continuous positive traits): \code{dmf}, \code{dmd},
#'  \code{dmvf}, \code{dmvd}, \code{fe}, \code{zn}, \code{ca}, \code{mg},
#'  \code{bc}, \code{tc}, \code{acrw}, \code{ancrw}, and \code{atrw}.
#'  
#'  \item \code{pnn} (percentage non-negative traits): \code{ci}, \code{hi},
#'  and \code{shi}.
#'  
#'  \item \code{ppo} (percentage positive traits): \code{dm}, \code{dmv},
#'  \code{prot}, \code{star}, \code{fruc}, \code{gluc}, \code{sucr}, and \code{malt}.
#'
#'  \item \code{dnn} (discrete non-negative traits): \code{nops}, \code{nope},
#'  \code{noph}, \code{nopr}, \code{nocr}, \code{nonc}, and \code{tnr}.
#'  
#'  \item \code{ctg} {categorical 1 to 9 traits}: \code{vir}, \code{vir1},
#'  \code{vir2}, \code{alt}, \code{alt1}, \code{alt2}, \code{vv}, \code{scol},
#'  \code{fcol}, \code{rs}, \code{rf}, \code{damr}, \code{rspr}, \code{wed},
#'  \code{fraw}, \code{fraw1}, \code{suraw}, \code{suraw1}, \code{straw},
#'  \code{straw1}, \code{coof}, \code{coof1}, \code{coosu}, \code{coosu1},
#'  \code{coost}, \code{coost1}, \code{coot}, \code{coot1}, \code{cooap},
#'  \code{cooap1}, \code{fraw2}, \code{suraw2}, \code{straw2}, \code{coof2},
#'  \code{coosu2}, \code{coost2}, \code{coot2}, and \code{cooap2}.
#' }
#' Values are set to \code{NA} or \code{0} with the following rules:
#' \itemize{
#'  \item \code{cnn} traits with negative values are set to \code{NA}.
#'  \item \code{cpo} traits with non-positive values are set to \code{NA}.
#'  \item \code{pnn} traits with values out of the [0, 100] interval are set to \code{NA}.  
#'  \item \code{ppo} with values out of the (0, 100] interval are set to \code{NA}.
#'  \item \code{dnn} traits with negative and non-integer values are set to \code{NA}.
#'  \item \code{ctg} traits with out of scale values are set to \code{NA}.
#'  \item Beta carotene values determined by RHS color charts with values different from
#'  the possible values in the RHS color chart are set to \code{NA}.
#'  \item Extreme low and high values are detected using the interquartile range.
#'  The rule is to detect any value out of the interval
#'  \eqn{[Q_1 - m - f \times IQR; Q_3 + m + f \times IQR]} where \code{m} is the
#'  median. By default \code{f = 10} and cannot be less than 10. Values out of this
#'  range are set to \code{NA}.
#'  \item If \code{nope == 0} and there is some data for any trait,
#'  then \code{nope} is set to \code{NA}.  
#'  \item If \code{noph == 0} and there is some data for any non-pre-harvest trait,
#'  then \code{noph} is set to \code{NA}.
#'  \item If \code{nopr == 0} and there is some data for any trait evaluated with roots,
#'  then \code{nopr} is set to \code{NA}.
#'  \item If \code{nocr == 0} and \code{crw > 0}, then \code{nocr} is set to \code{NA}.
#'  \item If \code{nocr > 0} and \code{crw == 0}, then \code{crw} is set to \code{NA}.
#'  \item If \code{nonc == 0} and \code{ncrw > 0}, then \code{nonc} is set to \code{NA}.
#'  \item If \code{nonc > 0} and \code{ncrw == 0}, then \code{ncrw} is set to \code{NA}.
#' }
#' @return It returns the data frame with all impossible values set to \code{NA}
#' or \code{0} and a list of warnings with all the rows that have been modified.
#' @author Raul Eyzaguirre.
#' @examples
#' dfr <- data.frame(trw = c(2.2, 5.0, 3.6, 12, 1600, -4),
#'                   dm = c(21, 23, 105, 24, -3, 30),
#'                   tnr = c(1.3, 10, 11, NA, 2, 5),
#'                   scol = c(1, 0, 15, 5, 4, 7),
#'                   fcol.cc = c(1, 15, 12, 24, 55, 20))
#' clean.data(dfr)
#' @importFrom stats IQR quantile median
#' @export

clean.data <- function(dfr, f = 10) {
  
  #############################################################################
  # Preliminary settings
  #############################################################################

  # Check f
  
  if (f < 10)
    f <- 10
  
  # Check names
  
  dfr <- check.names.sp(dfr)

  # Pre-harvest traits
  
  pre <- c("vir", "vir1", "vir2", "alt", "alt1", "alt2", "vv")
  
  # Traits evaluated without roots
  
  nrt <- c("vir", "vir1", "vir2", "alt", "alt1", "alt2", "vv", "vw", "dmvf",
           "dmvd", "dmv", "vw.d", "fytha", "fytha.aj", "dmvy", "dmvy.aj",
           "vpp", "vpsp", "shi")
  
  # Continuous non-negative traits
  
  cnn <- c("vw", "crw", "ncrw", "trw", "trw.d", "biom", "biom.d", "cytha",
           "cytha.aj", "rytha", "rytha.aj", "dmry", "dmry.aj", "vw.d", "fytha",
           "fytha.aj", "dmvy", "dmvy.aj", "bytha", "bytha.aj", "dmby", "dmby.aj",
           "nrpp", "nrpsp", "ncrpp", "ncrpsp", "ypp", "ypsp", "vpp", "vpsp", "rfr")
  
  # Continuous positive traits
  
  cpo <- c("dmf", "dmd", "dmvf", "dmvd", "fe", "zn", "ca", "mg", "bc", "tc",
           "acrw", "ancrw", "atrw")
  
  # Percentage non-negative traits
  
  pnn <- c("ci", "hi", "shi")
  
  # Percentage positive traits
  
  ppo <- c("dm", "dmv", "prot", "star", "fruc", "gluc", "sucr", "malt")
  
  # Discrete non-negative traits
  
  dnn <- c("nops", "nope", "noph", "nopr", "nocr", "nonc", "tnr")
  
  # Categorical 1 to 9 traits
  
  ctg <- c("vir", "vir1", "vir2", "alt", "alt1", "alt2", "vv", "scol", "fcol",
           "rs", "rf", "damr", "rspr", "wed", "fraw", "fraw1", "suraw", "suraw1",
           "straw", "straw1", "coof", "coof1", "coosu", "coosu1", "coost", "coost1",
           "coot", "coot1", "cooap", "cooap1", "fraw2", "suraw2", "straw2", "coof2",
           "coosu2", "coost2", "coot2", "cooap2")

  # Special traits
  
  bc.cc <- "bc.cc"
  fcol.cc <- "fcol.cc"
  
  #############################################################################
  # Impossible values
  #############################################################################

  # Impossible values for continuous non-negative traits
  
  for (i in 1:length(cnn))
    if (exists(cnn[i], dfr)) {
      cond <- dfr[, cnn[i]] < 0 & !is.na(dfr[, cnn[i]])
      dfr[cond, cnn[i]] <- NA
      if (sum(cond) > 0)
        warning("- Rows with negative values replaced with NA for trait ",
                cnn[i], ": ", paste0(rownames(dfr)[cond], " "), call. = FALSE)
    }
  
  # Impossible values for continuous positive traits
  
  for (i in 1:length(cpo))
    if (exists(cpo[i], dfr)) {
      cond <- dfr[, cpo[i]] < 0 & !is.na(dfr[, cpo[i]])
      dfr[cond, cpo[i]] <- NA
      if (sum(cond) > 0)
        warning("- Rows with non-positive values replaced with NA for trait ",
                cpo[i], ": ", paste0(rownames(dfr)[cond], " "), call. = FALSE)
    }

  # Impossible values for percentage non-negative traits
  
  for (i in 1:length(pnn))
    if (exists(pnn[i], dfr)) {
      cond1 <- dfr[, pnn[i]] < 0 & !is.na(dfr[, pnn[i]])
      cond2 <- dfr[, pnn[i]] > 100 & !is.na(dfr[, pnn[i]])
      cond <- cond1 | cond2
      dfr[cond, pnn[i]] <- NA
      if (sum(cond) > 0)
        warning("- Rows with values out of [0-100] replaced with NA for trait ",
                pnn[i], ": ", paste0(rownames(dfr)[cond], " "), call. = FALSE)
    }
  
  # Impossible values for percentage positive traits
  
  for (i in 1:length(ppo))
    if (exists(ppo[i], dfr)) {
      cond1 <- dfr[, ppo[i]] <= 0 & !is.na(dfr[, ppo[i]])
      cond2 <- dfr[, ppo[i]] > 100 & !is.na(dfr[, ppo[i]])
      cond <- cond1 | cond2
      dfr[cond, ppo[i]] <- NA
      if (sum(cond) > 0)
        warning("- Rows with values out of (0-100] replaced with NA for trait ",
                ppo[i], ": ", paste0(rownames(dfr)[cond], " "), call. = FALSE)
    }

  # Impossible values for discrete non-negative traits
  
  for (i in 1:length(dnn))
    if (exists(dnn[i], dfr)) {
      cond1 <- dfr[, dnn[i]] < 0 & !is.na(dfr[, dnn[i]])
      cond2 <- dfr[, dnn[i]] %% 1 > 0 & !is.na(dfr[, dnn[i]])
      cond <- cond1 | cond2
      dfr[cond, dnn[i]] <- NA
      if (sum(cond) > 0)
        warning("- Rows with negative or non integer values replaced with NA for trait ",
                dnn[i], ": ", paste0(rownames(dfr)[cond], " "), call. = FALSE)
    }
  
  # Impossible values for 1 to 9 categorical traits
  
  for (i in 1:length(ctg))
    if (exists(ctg[i], dfr)) {
      cond <- !(dfr[, ctg[i]] %in% 1:9) & !is.na(dfr[, ctg[i]])
      dfr[cond, ctg[i]] <- NA
      if (sum(cond) > 0)
        warning("- Rows with values out of 1-9 integer scale replaced with NA for trait ",
                ctg[i], ": ", paste0(rownames(dfr)[cond], " "), call. = FALSE)
    }
  
  # Impossible values for bc.cc and fcol.cc
  
  if (exists(bc.cc, dfr)) {
    bc.cc.values <- c(0.03, 0, 0.12, 0.02, 0.15, 1.38, 1.65, 1.5, 1.74, 1.76,
                      0.69, 1.17, 1.32, 1.04, 4.41, 4.92, 6.12, 5.46, 3.96, 5.49,
                      3.03, 3.76, 4.61, 7.23, 7.76, 10.5, 11.03, 12.39, 14.37)
    cond <- !(dfr[, bc.cc] %in% bc.cc.values) & !is.na(dfr[, bc.cc])
    dfr[cond, bc.cc] <- NA
    if (sum(cond) > 0)
      warning("- Rows with values out of scale replaced with NA for trait ",
              bc.cc, ": ", paste0(rownames(dfr)[cond], " "), call. = FALSE)
  }

  if (exists(fcol.cc, dfr)) {
    cond <- !(dfr[, fcol.cc] %in% 1:30) & !is.na(dfr[, fcol.cc])
    dfr[cond, fcol.cc] <- NA
    if (sum(cond) > 0)
      warning("- Rows with values out of integer scale replaced with NA for trait ",
              fcol.cc, ": ", paste0(rownames(dfr)[cond], " "), call. = FALSE)
  }
  
  # Extreme values (almost impossible)
  
  t.all <- c(cnn, cpo, pnn, ppo, dnn)
  
  for (i in 1:length(t.all))
    if (exists(t.all[i], dfr)) {
      tol <- IQR(dfr[, t.all[i]], na.rm = TRUE)
      m <- median(dfr[, t.all[i]], na.rm = TRUE)
      cond1 <- dfr[, t.all[i]] < quantile(dfr[, t.all[i]], 0.25, na.rm = TRUE) -
        m - f * tol & !is.na(dfr[, t.all[i]])
      cond2 <- dfr[, t.all[i]] > quantile(dfr[, t.all[i]], 0.75, na.rm = TRUE) + 
        m + f * tol & !is.na(dfr[, t.all[i]])
      cond <- cond1 | cond2
      dfr[cond, t.all[i]] <- NA
      if (sum(cond) > 0)
        warning("- Rows with extreme values replaced with NA for trait ",
                t.all[i], ": ", paste0(rownames(dfr)[cond], " "), call. = FALSE)
    }
  
  #############################################################################
  # nope, noph, nopr consistency
  #############################################################################
  
  # Subset in fieldook all traits
  
  t.all <- c(cnn, cpo, pnn, ppo, dnn, ctg, bc.cc, fcol.cc)
  t.all <- t.all[t.all %in% colnames(dfr)]
  t.all <- t.all[!(t.all %in% c("nops", "nope"))]
  
  # Subset in fieldook all non-pre-harvest traits
  
  t.pos <- t.all[!(t.all %in% pre)]
  t.pos <- t.pos[t.pos != "noph"]

  # Subset in fieldook all traits evaluated with roots
  
  t.rot <- t.pos[!(t.pos %in% nrt)]
  t.rot <- t.rot[t.rot != "nopr"]
  
  # Rule 1 for nope
  
  if (length(t.all) > 0 & exists("nope", dfr)) {
    if (length(t.all) == 1)
      cond <- dfr[, t.all] > 0 & !is.na(dfr[, t.all])
    if (length(t.all) > 1)
      cond <- apply(dfr[, t.all] > 0 & !is.na(dfr[, t.all]), 1, sum) > 0
    dfr[cond & dfr[, 'nope'] == 0 & !is.na(dfr[, 'nope']), 'nope'] <- NA
    if (sum(cond) > 0)
      warning("- Rows replaced with NA for trait nope:",
              paste0(rownames(dfr)[cond], " "), call. = FALSE)
  }
  
  # Rule 2 for noph
  
  if (length(t.pos) > 0 & exists("noph", dfr)) {
    if (length(t.pos) == 1)
      cond <- dfr[, t.pos] > 0 & !is.na(dfr[, t.pos])
    if (length(t.pos) > 1)
      cond <- apply(dfr[, t.pos] > 0 & !is.na(dfr[, t.pos]), 1, sum) > 0
    dfr[cond & dfr[, 'noph'] == 0 & !is.na(dfr[, 'noph']), 'noph'] <- NA
    if (sum(cond) > 0)
      warning("- Rows replaced with NA for trait noph:",
              paste0(rownames(dfr)[cond], " "), call. = FALSE)
  }

  # Rule 3 for nopr
  
  if (length(t.rot) > 0 & exists("nopr", dfr)) {
    if (length(t.rot) == 1)
      cond <- dfr[, t.rot] > 0 & !is.na(dfr[, t.rot])
    if (length(t.rot) > 1)
      cond <- apply(dfr[, t.rot] > 0 & !is.na(dfr[, t.rot]), 1, sum) > 0
    dfr[cond & dfr[, 'nopr'] == 0 & !is.na(dfr[, 'nopr']), 'nopr'] <- NA
    if (sum(cond) > 0)
      warning("- Rows replaced with NA for trait nopr:",
              paste0(rownames(dfr)[cond], " "), call. = FALSE)
  }
  
  #############################################################################
  # Other traits consistency
  #############################################################################
  
  # nocr and crw
  
  if (exists("nocr", dfr) & exists("crw", dfr)) {
    cond <- dfr[, "nocr"] == 0 & !is.na(dfr[, "nocr"]) & dfr[, "crw"] > 0 & !is.na(dfr[, "crw"])
    dfr[cond, 'nocr'] <- NA
    if (sum(cond) > 0)
      warning("- Rows replaced with NA for trait nocr:",
              paste0(rownames(dfr)[cond], " "), call. = FALSE)
    cond <- dfr[, "nocr"] > 0 & !is.na(dfr[, "nocr"]) & dfr[, "crw"] == 0 & !is.na(dfr[, "crw"])
    dfr[cond, 'crw'] <- NA
    if (sum(cond) > 0)
      warning("- Rows replaced with NA for trait crw:",
              paste0(rownames(dfr)[cond], " "), call. = FALSE)
  }
  
  # nonc and ncrw
  
  if (exists("nonc", dfr) & exists("ncrw", dfr)) {
    cond <- dfr[, "nonc"] == 0 & !is.na(dfr[, "nonc"]) & dfr[, "ncrw"] > 0 & !is.na(dfr[, "ncrw"])
    dfr[cond, 'nonc'] <- NA
    if (sum(cond) > 0)
      warning("- Rows replaced with NA for trait nonc:",
              paste0(rownames(dfr)[cond], " "), call. = FALSE)
    cond <- dfr[, "nonc"] > 0 & !is.na(dfr[, "nonc"]) & dfr[, "ncrw"] == 0 & !is.na(dfr[, "ncrw"])
    dfr[cond, 'ncrw'] <- NA
    if (sum(cond) > 0)
      warning("- Rows replaced with NA for trait ncrw:",
              paste0(rownames(dfr)[cond], " "), call. = FALSE)
  }
  
  # Return data frame
  
  dfr
  
}