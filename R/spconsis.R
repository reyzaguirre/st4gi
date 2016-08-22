#' Check consistency for sweetpotato experimental data
#'
#' Set of rules to check for consistency of sweetpotato experimental data.
#' Data labels must be defined as specified in the PROCEDURES FOR THE EVALUATION
#' AND ANALYSIS OF SWEETPOTATO TRIALS document.
#' @param fb The name of the fieldbook data frame.
#' @param f Factor for extreme values detection. See details.
#' @param out.mod Statistical model for outliers' detection. See details.
#' @param out.max Threshold for outliers' detection.
#' @param aqt Additional quantitative traits.
#' @param width Number of columns for the output.
#' @param file Logigal, if TRUE the output goes to a file with name output.txt.
#' @details The data frame must use the labels (lower or upper case) listed in function \code{checknames}.
#' See \code{?checknames} for details.
#' 
#' Extreme low and high values are detected using the interquartile range.
#' The rule is to detect any value out of the interval 
#' \eqn{[Q_1 - f \times IQR; Q_3 + f \times IQR]}. By default \code{f = 3}.
#' 
#' Outliers are detected based on standardized residuals for some statistical models.
#' Options are \code{'rcbd'} and \code{'met'} for a randomized complete block design
#' and a multi environment trial with RCBD in each environment. By default the threshold
#' value is \code{out.max = 4}.
#' @return If \code{file = TRUE} it returns a file with name checks.txt with a list of
#' all rows with some kind of inconsistency and all rows with outliers. If \code{file = FALSE}
#' the output is shown in the R console.
#' @author Raul Eyzaguirre.
#' @examples
#' spconsis(pjpz09)
#' @importFrom stats IQR quantile rstandard
#' @export

spconsis <- function(fb, f = 3, out.mod = c("none", "rcbd", "met"),
                     out.max = 4, aqt = NULL, width = 240, file = FALSE) {

  options(width = width)
  
  out.mod = match.arg(out.mod)
  
  fb <- checknames(fb, aqt)
  if (!is.null(aqt))
    aqt <- toupper(aqt)

  if (file == TRUE) sink(paste(getwd(), "/checks.txt", sep = ""))

  # Inconsistencies for NOPS > NOPE > NOPH > NOPR.

  spc01(fb, 1, "NOPE", "NOPS", "- Number of plants established (NOPE) is greater than number of plants sowed (NOPS):")
  spc01(fb, 1, "NOPH", "NOPS", "- Number of plants harvested (NOPH)  is greater than number of plants sowed (NOPS):")
  spc01(fb, 1, "NOPR", "NOPS", "- Number of plants with roots (NOPR) is greater than number of plants sowed (NOPS):")
  spc01(fb, 1, "NOPH", "NOPE", "- Number of plants harvested (NOPH) is greater than number of plants established (NOPE):")
  spc01(fb, 1, "NOPR", "NOPE", "- Number of plants with roots (NOPR) is greater than number of plants established (NOPE):")
  spc01(fb, 1, "NOPR", "NOPH", "- Number of plants with roots (NOPR) is greater than number of plants harvested (NOPH):")

  # Inconsistencies for NOPE and dependencies.

  spc01(fb, 2, "NOPE", "VIR", "- Number of plants established (NOPE) is zero or NA but there is data for virus symptoms (VIR):")
  spc01(fb, 2, "NOPE", "VIR1", "- Number of plants established (NOPE) is zero or NA but there is data for virus symptoms first evaluation (VIR1):")
  spc01(fb, 2, "NOPE", "VIR2", "- Number of plants established (NOPE) is zero or NA but there is data for virus symptoms second evaluation (VIR2):")
  spc01(fb, 2, "NOPE", "ALT", "- Number of plants established (NOPE) is zero or NA but there is data for alternaria symptoms (ALT):")
  spc01(fb, 2, "NOPE", "ALT1", "- Number of plants established (NOPE) is zero or NA but there is data for alternaria symptoms first evaluation (ALT1):")
  spc01(fb, 2, "NOPE", "ALT2", "- Number of plants established (NOPE) is zero or NA but there is data for alternaria symptoms second evaluation (ALT2):")
  spc01(fb, 2, "NOPE", "VV", "- Number of plants established (NOPE) is zero or NA but there is data for vine vigor (VV):")

  # NOPH and VW
  
  spc01(fb, 3, "NOPH", "VW", "- Number of plants harvested (NOPH) is zero or NA but vine weight (VW) is greater than zero:") 
  spc01(fb, 3, "VW", "NOPH", "- Vine weight (VW) is zero or NA but the number of plants harvested (NOPH) is greater than zero:") 
  spc01(fb, 3, "NOPH", "FYTHA", "- Number of plants harvested (NOPH) is zero or NA but foliage yield in tons per hectare (FYTHA) is greater than zero:") 
  spc01(fb, 3, "FYTHA", "NOPH", "- Foliage yield in tons per hectare (FYTHA) is zero or NA but the number of plants harvested (NOPH) is greater than zero:") 
  spc01(fb, 3, "NOPH", "FYTHA.AJ", "- Number of plants harvested (NOPH) is zero or NA but foliage yield in tons per hectare (FYTHA.AJ) is greater than zero:") 
  spc01(fb, 3, "FYTHA.AJ", "NOPH", "- Foliage yield in tons per hectare (FYTHA.AJ) is zero or NA but the number of plants harvested (NOPH) is greater than zero:") 
  
  # VW and dependencies
  
  spc01(fb, 3, "VW", "DMVF", "- Vine weight (VW) is zero or NA but there is fresh weight vines for dry matter assessment (DMVF):") 
  spc01(fb, 3, "VW", "DMVD", "- Vine weight (VW) is zero or NA but there is dry weight vines for dry matter assessment (DMVD):") 
  spc01(fb, 1, "DMVD", "DMVF", "- Dry weight vines for dry matter assessment (DMVD) is greater than fresh weight vines for dry matter assessment (DBVF):")

  # NOPR and roots
  
  spc02(fb, 1, "NOPR", "NOCR", "NONC", "- Number of plants with roots (NOPR) is zero or NA but number of roots (NOCR + NONC) is greater than zero:")
  spc02(fb, 2, "NOPR", "NOCR", "NONC", "- Number of roots (NOCR + NONC) is zero or NA but number of plants with roots (NOPR) is greater than zero:")
  spc02(fb, 1, "NOPR", "CRW", "NCRW", "- Number of plants with roots (NOPR) is zero or NA but root weight (CRW + NCRW) is greater than zero:")
  spc02(fb, 2, "NOPR", "CRW", "NCRW", "- Root weight (CRW + NCRW) is zero or NA but number of plants with roots (NOPR) is greater than zero:")
  spc01(fb, 3, "NOPR", "TRW", "- Number of plants with roots (NOPR) is zero or NA but total root weight (TRW) is greater than zero::")
  spc01(fb, 3, "TRW", "NOPR", "- Total root weight (TRW) is zero or NA but number of plants with roots (NOPR) is greater than zero:")
  spc01(fb, 3, "NOPR", "RYTHA", "- Number of plants with roots (NOPR) is zero or NA but root yield in tons per hectare (RYTHA) is greater than zero::")
  spc01(fb, 3, "RYTHA", "NOPR", "- Root yield in tons per hectare (RYTHA) is zero or NA but number of plants with roots (NOPR) is greater than zero:")
  spc01(fb, 3, "NOPR", "RYTHA.AJ", "- Number of plants with roots (NOPR) is zero or NA but root yield in tons per hectare (RYTHA.AJ) is greater than zero::")
  spc01(fb, 3, "RYTHA.AJ", "NOPR", "- Root yield in tons per hectare (RYTHA.AJ) is zero or NA but number of plants with roots (NOPR) is greater than zero:")
  spc01(fb, 3, "NOPR", "CYTHA", "- Number of plants with roots (NOPR) is zero or NA but commercial root yield in tons per hectare (CYTHA) is greater than zero::")
  spc01(fb, 3, "NOPR", "CYTHA.AJ", "- Number of plants with roots (NOPR) is zero or NA but commercial root yield in tons per hectare (CYTHA.AJ) is greater than zero::")
  spc01(fb, 2, "NOPR", "WED", "- Number of plants with roots (NOPR) is zero or NA but there is data for weevil damage (WED):")
  
  # Number of roots and root weight
  
  spc01(fb, 3, "NOCR", "CRW", "- Number of commercial roots (NOCR) is zero or NA but the commercial root weight (CRW) is greater than zero:") 
  spc01(fb, 3, "CRW", "NOCR", "- Commercial root weight (CRW) is zero or NA but the number of commercial roots (NOCR) is greater than zero:") 
  spc01(fb, 3, "NOCR", "CYTHA", "- Number of commercial roots (NOCR) is zero or NA but the commercial root yield in tons per hectare (CYTHA) is greater than zero:") 
  spc01(fb, 3, "CYTHA", "NOCR", "- Commercial root yield in tons per hectare (CYTHA) is zero or NA but the number of commercial roots (NOCR) is greater than zero:") 
  spc01(fb, 3, "NOCR", "CYTHA.AJ", "- Number of commercial roots (NOCR) is zero or NA but the commercial root yield in tons per hectare (CYTHA.AJ) is greater than zero:") 
  spc01(fb, 3, "CYTHA.AJ", "NOCR", "- Commercial root yield in tons per hectare (CYTHA.AJ) is zero or NA but the number of commercial roots (NOCR) is greater than zero:") 
  spc01(fb, 3, "NONC", "NCRW", "- Number of non commercial roots (NONC) is zero or NA but the non commercial root weight (NCRW) is greater than zero:")
  spc01(fb, 3, "NCRW", "NONC", "- Non commercial root weight (NCRW) is zero or NA but the number of non commercial roots (NONC) is greater than zero:")
  spc02(fb, 3, "TRW", "NOCR", "NONC", "- Total root weight (TRW) is zero or NA but number of roots (NOCR + NONC) is greater than zero:")
  spc02(fb, 2, "TRW", "NOCR", "NONC", "- Number of roots (NOCR + NONC) is zero or NA but total root weight (TRW) is greater than zero:")
  spc02(fb, 3, "RYTHA", "NOCR", "NONC", "- Root yield in tons per hectare (RYTHA) is zero or NA but number of roots (NOCR + NONC) is greater than zero:")
  spc02(fb, 2, "RYTHA", "NOCR", "NONC", "- Number of roots (NOCR + NONC) is zero or NA but root yield in tons per hectare (RYTHA) is greater than zero:")
  spc02(fb, 3, "RYTHA.AJ", "NOCR", "NONC", "- Root yield in tons per hectare (RYTHA.AJ) is zero or NA but number of roots (NOCR + NONC) is greater than zero:")
  spc02(fb, 2, "RYTHA.AJ", "NOCR", "NONC", "- Number of roots (NOCR + NONC) is zero or NA but root yield in tons per hectare (RYTHA.AJ) is greater than zero:")
  
  # Roots and dependencies

  temp <- array(FALSE, dim(fb)[1])
  do <- FALSE
  
  if (exists("NOPR", fb)) {
    temp <- temp | (!is.na(fb$NOPR) & fb$NOPR > 0)
    do <- TRUE
  }
  if (exists("NOCR", fb)) {
    temp <- temp | (!is.na(fb$NOCR) & fb$NOCR > 0)
    do <- TRUE
  }
  if (exists("NONC", fb)) {
    temp <- temp | (!is.na(fb$NONC) & fb$NONC > 0)
    do <- TRUE
  }
  if (exists("CRW", fb)) {
    temp <- temp | (!is.na(fb$CRW) & fb$CRW > 0)
    do <- TRUE
  }
  if (exists("NCRW", fb)) {
    temp <- temp | (!is.na(fb$NCRW) & fb$NCRW > 0)
    do <- TRUE
  }
  if (exists("TRW", fb)) {
    temp <- temp | (!is.na(fb$TRW) & fb$TRW > 0)
    do <- TRUE
  }
  if (exists("RYTHA", fb)) {
    temp <- temp | (!is.na(fb$RYTHA) & fb$RYTHA > 0)
    do <- TRUE
  }
  if (exists("RYTHA.AJ", fb)) {
    temp <- temp | (!is.na(fb$RYTHA.AJ) & fb$RYTHA.AJ > 0)
    do <- TRUE
  }

  spc03(fb, temp, do, "RFCP", "- There are no roots but there is data for root primary flesh color (RFCP):")
  spc03(fb, temp, do, "RFCS", "- There are no roots but there is data for root secondary flesh color (RFCS):")
  spc03(fb, temp, do, "SCOL", "- There are no roots but there is data for storage root skin color (SCOL):")
  spc03(fb, temp, do, "FCOL", "- There are no roots but there is data for storage root flesh color (FCOL):")
  spc03(fb, temp, do, "RS", "- There are no roots but there is data for root size (RS):")
  spc03(fb, temp, do, "RF", "- There are no roots but there is data for root form (RF):")
  spc03(fb, temp, do, "DAMR", "- There are no roots but there is data for root defects (DAMR):")
  spc03(fb, temp, do, "RSPR", "- There are no roots but there is data for root sprouting (RSPR):")
  spc03(fb, temp, do, "WED", "- There are no roots but there is data for weevil damage (WED):")
  spc03(fb, temp, do, "DMF", "- There are no roots but there is data for fresh weight of roots for dry matter assessment (DMF):")
  spc03(fb, temp, do, "DMD", "- There are no roots but there is data for dry weight of roots for dry matter assessment (DMD):")
  spc01(fb, 1, "DMD", "DMF", "- Dry weight of roots for dry matter assessment (DMD) is greater than fresh weight of roots for dry matter assessment (DMF):")
  spc03(fb, temp, do, "FRAW", "- There are no roots but there is data for root fiber (FRAW):")
  spc03(fb, temp, do, "SURAW", "- There are no roots but there is data for root sugar (SURAW):")
  spc03(fb, temp, do, "STRAW", "- There are no roots but there is data for root starch (STRAW):")
  spc03(fb, temp, do, "COOF", "- There are no roots but there is data for cooked fiber (COOF):")
  spc03(fb, temp, do, "COOSU", "- There are no roots but there is data for cooked sugars (COOSU):")
  spc03(fb, temp, do, "COOST", "- There are no roots but there is data for cooked starch (COOST):")
  spc03(fb, temp, do, "COOT", "- There are no roots but there is data for cooked taste (COOT):")
  spc03(fb, temp, do, "COOAP", "- There are no roots but there is data for cooked appearance (COOAP):")
  spc03(fb, temp, do, "FRAW1", "- There are no roots but there is data for root fiber first determination (FRAW1):")
  spc03(fb, temp, do, "SURAW1", "- There are no roots but there is data for root sugar first determination (SURAW1):")
  spc03(fb, temp, do, "STRAW1", "- There are no roots but there is data for root starch first determination (STRAW1):")
  spc03(fb, temp, do, "COOF1", "- There are no roots but there is data for cooked fiber first evaluation (COOF1):")
  spc03(fb, temp, do, "COOSU1", "- There are no roots but there is data for cooked sugars first evaluation (COOSU1):")
  spc03(fb, temp, do, "COOST1", "- There are no roots but there is data for cooked starch first evaluation (COOST1):")
  spc03(fb, temp, do, "COOT1", "- There are no roots but there is data for cooked taste first evaluation (COOT1):")
  spc03(fb, temp, do, "COOAP1", "- There are no roots but there is data for cooked appearance first evaluation (COOAP1):")
  spc03(fb, temp, do, "FRAW2", "- There are no roots but there is data for root fiber second determination (FRAW2):")
  spc03(fb, temp, do, "SURAW2", "- There are no roots but there is data for root sugar second determination (SURAW2):")
  spc03(fb, temp, do, "STRAW2", "- There are no roots but there is data for root starch second determination (STRAW2):")
  spc03(fb, temp, do, "COOF2", "- There are no roots but there is data for cooked fiber second evaluation (COOF2):")
  spc03(fb, temp, do, "COOSU2", "- There are no roots but there is data for cooked sugars second evaluation (COOSU2):")
  spc03(fb, temp, do, "COOST2", "- There are no roots but there is data for cooked starch second evaluation (COOST2):")
  spc03(fb, temp, do, "COOT2", "- There are no roots but there is data for cooked taste second evaluation (COOT2):")
  spc03(fb, temp, do, "COOAP2", "- There are no roots but there is data for cooked appearance second evaluation (COOAP2):")
  spc03(fb, temp, do, "PROT", "- There are no roots but there is data for protein (PROT):")
  spc03(fb, temp, do, "FE", "- There are no roots but there is data for iron (FE):")
  spc03(fb, temp, do, "ZN", "- There are no roots but there is data for zinc (ZN):")
  spc03(fb, temp, do, "CA", "- There are no roots but there is data for calcium (CA):")
  spc03(fb, temp, do, "MG", "- There are no roots but there is data for magnesium (MG):")
  spc03(fb, temp, do, "BC", "- There are no roots but there is data for beta-carotene (BC):")
  spc03(fb, temp, do, "BC.CC", "- There are no roots but there is data for beta-carotene (BC.CC):")
  spc03(fb, temp, do, "TC", "- There are no roots but there is data for total carotenoids (TC):")
  spc03(fb, temp, do, "STAR", "- There are no roots but there is data for starch (STAR):")
  spc03(fb, temp, do, "FRUC", "- There are no roots but there is data for fructose (FRUC):")
  spc03(fb, temp, do, "GLUC", "- There are no roots but there is data for glucose (GLUC):")
  spc03(fb, temp, do, "SUCR", "- There are no roots but there is data for sucrose (SUCR):")
  spc03(fb, temp, do, "MALT", "- There are no roots but there is data for maltose (MALT):")
  
  # Extreme values detection and values out of range for field data

  spc04(fb, c(0:100, NA), "NOPE", "- Out of range values for number of plants established (NOPE):")
  spc04(fb, c(1:9, NA), "VIR", "- Out of range values for virus symptoms (VIR):")
  spc04(fb, c(1:9, NA), "VIR1", "- Out of range values for virus symptoms first evaluation (VIR1):")
  spc04(fb, c(1:9, NA), "VIR2", "- Out of range values for virus symptoms second evaluation (VIR2):")
  spc04(fb, c(1:9, NA), "ALT", "- Out of range values for alternaria symptoms (ALT):")
  spc04(fb, c(1:9, NA), "ALT1", "- Out of range values for alternaria symptoms first evaluation (ALT1):")
  spc04(fb, c(1:9, NA), "ALT2", "- Out of range values for alternaria symptoms second evaluation (ALT2):")
  spc04(fb, c(1:9, NA), "VV", "- Out of range values for vine vigor (VV):")
  spc05(fb, "lower", "VW", "- Out of range values for vine weight (VW):")
  spc06(fb, f, "low", "VW", "- Extreme low values for vine weight (VW):")
  spc06(fb, f, "high", "VW", "- Extreme high values for vine weight (VW):")
  spc05(fb, "lower", "NOPH", "- Out of range values for number of plants harvested (NOPH):")
  spc05(fb, "lower", "NOPR", "- Out of range values for number of plants with roots (NOPR):")
  spc05(fb, "lower", "NOCR", "- Out of range values for number of commercial roots (NOCR):")
  spc06(fb, f, "low", "NOCR", "- Extreme low values for number of commercial roots (NOCR):")
  spc06(fb, f, "high", "NOCR", "- Extreme high values for number of commercial roots (NOCR):")
  spc05(fb, "lower", "NONC", "- Out of range values for number of non commercial roots (NONC):")
  spc06(fb, f, "low", "NONC", "- Extreme low values for number of non commercial roots (NONC):")
  spc06(fb, f, "high", "NONC", "- Extreme high values for number of non commercial roots (NONC):")
  spc05(fb, "lower", "CRW", "- Out of range values for commercial root weight (CRW):")
  spc06(fb, f, "low", "CRW", "- Extreme low values for commercial root weight (CRW):")
  spc06(fb, f, "high", "CRW", "- Extreme high values for commercial root weight (CRW):")
  spc05(fb, "lower", "NCRW", "- Out of range values for non commercial root weight (NCRW):")
  spc06(fb, f, "low", "NCRW", "- Extreme low values for non commercial root weight (NCRW):")
  spc06(fb, f, "high", "NCRW", "- Extreme high values for non commercial root weight (NCRW):")
  spc04(fb, c(1:9, NA), "SCOL", "- Out of range values for storage root skin color (SCOL):")
  spc04(fb, c(1:9, NA), "FCOL", "- Out of range values for storage root flesh color (FCOL):")
  spc04(fb, c(1:9, NA), "RS", "- Out of range values for root size (RS):")
  spc04(fb, c(1:9, NA), "RF", "- Out of range values for root form (RF):")
  spc04(fb, c(1:9, NA), "DAMR", "- Out of range values for root defects (DAMR):")
  spc04(fb, c(1:9, NA), "RSPR", "- Out of range values for root sprouting (RSPR):")
  spc04(fb, c(1:9, NA), "WED", "- Out of range values for weevil damage (WED):")

  # Extreme values detection and values out of range for DM data
  
  spc05(fb, "lower", "DMF", "- Out of range values for fresh weight of roots for dry matter assessment (DMF):")
  spc06(fb, f, "low", "DMF", "- Extreme low values for fresh weight of roots for dry matter assessment (DMF):")
  spc06(fb, f, "high", "DMF", "- Extreme high values for fresh weight of roots for dry matter assessment (DMF):")
  spc05(fb, "lower", "DMD", "- Out of range values for dry weight of roots for dry matter assessment (DMD):")
  spc06(fb, f, "low", "DMD", "- Extreme low values for dry weight of roots for dry matter assessment (DMD):")
  spc06(fb, f, "high", "DMD", "- Extreme high values for dry weight of roots for dry matter assessment (DMD):")
  spc05(fb, "lower", "DMVF", "- Out of range values for fresh weight vines for dry matter assessment (DMVF):")
  spc06(fb, f, "low", "DMVF", "- Extreme low values for fresh weight of vines for dry matter assessment (DMVF):")
  spc06(fb, f, "high", "DMVF", "- Extreme high values for fresh weight of vines for dry matter assessment (DMVF):")
  spc05(fb, "lower", "DMVD", "- Out of range values for dry weight of vines for dry matter assessment (DMVD):")
  spc06(fb, f, "low", "DMVD", "- Extreme low values for dry weight of vines for dry matter assessment (DMVD):")
  spc06(fb, f, "high", "DMVD", "- Extreme high values for dry weight of vines for dry matter assessment (DMVD):")
  spc05(fb, "lower", "DM", "- Out of range values for storage root dry matter content (DM):")
  spc06(fb, f, "low", "DM", "- Extreme low values for storage root dry matter content (DM):")
  spc06(fb, f, "high", "DM", "- Extreme high values for storage root dry matter content (DM):")
  spc05(fb, "lower", "DMV", "- Out of range values for vine dry matter content (DMV):")
  spc06(fb, f, "low", "DMV", "- Extreme low values for vine dry matter content (DMV):")
  spc06(fb, f, "high", "DMV", "- Extreme high values for vine dry matter content (DMV):")
  spc05(fb, "lower", "DMFY", "- Out of range values for dry matter foliage yield (DMFY):")
  spc06(fb, f, "low", "DMFY", "- Extreme low values for dry matter foliage yield (DMFY):")
  spc06(fb, f, "high", "DMFY", "- Extreme high values for dry matter foliage yield (DMFY):")
  spc05(fb, "lower", "DMRY", "- Out of range values for dry matter root yield (DMRY):")
  spc06(fb, f, "low", "DMRY", "- Extreme low values for dry matter root yield (DMRY):")
  spc06(fb, f, "high", "DMRY", "- Extreme high values for dry matter root yield (DMRY):")
  
  # Extreme values detection and values out of range for cooked traits
  
  spc04(fb, c(1:9, NA), "FRAW", "- Out of range values for root fiber (FRAW):")
  spc04(fb, c(1:9, NA), "SURAW", "- Out of range values for root sugar (SURAW):")
  spc04(fb, c(1:9, NA), "STRAW", "- Out of range values for root starch (STRAW):")
  spc04(fb, c(1:9, NA), "COOF", "- Out of range values for cooked fiber (COOF):")
  spc04(fb, c(1:9, NA), "COOSU", "- Out of range values for cooked sugars (COOSU):")
  spc04(fb, c(1:9, NA), "COOST", "- Out of range values for cooked starch (COOST):")
  spc04(fb, c(1:9, NA), "COOT", "- Out of range values for cooked taste (COOT):")
  spc04(fb, c(1:9, NA), "COOAP", "- Out of range values for cooked appearance (COOAP):")
  spc04(fb, c(1:9, NA), "FRAW1", "- Out of range values for root fiber first determination (FRAW1):")
  spc04(fb, c(1:9, NA), "SURAW1", "- Out of range values for root sugar first determination (SURAW1):")
  spc04(fb, c(1:9, NA), "STRAW1", "- Out of range values for root starch first determination (STRAW1):")
  spc04(fb, c(1:9, NA), "COOF1", "- Out of range values for cooked fiber first evaluation (COOF1):")
  spc04(fb, c(1:9, NA), "COOSU1", "- Out of range values for cooked sugars first evaluation (COOSU1):")
  spc04(fb, c(1:9, NA), "COOST1", "- Out of range values for cooked starch first evaluation (COOST1):")
  spc04(fb, c(1:9, NA), "COOT1", "- Out of range values for cooked taste first evaluation (COOT1):")
  spc04(fb, c(1:9, NA), "COOAP1", "- Out of range values for cooked appearance first evaluation (COOAP1):")
  spc04(fb, c(1:9, NA), "FRAW2", "- Out of range values for root fiber second determination (FRAW2):")
  spc04(fb, c(1:9, NA), "SURAW2", "- Out of range values for root sugar second determination (SURAW2):")
  spc04(fb, c(1:9, NA), "STRAW2", "- Out of range values for root starch second determination (STRAW2):")
  spc04(fb, c(1:9, NA), "COOF2", "- Out of range values for cooked fiber second evaluation (COOF2):")
  spc04(fb, c(1:9, NA), "COOSU2", "- Out of range values for cooked sugars second evaluation (COOSU2):")
  spc04(fb, c(1:9, NA), "COOST2", "- Out of range values for cooked starch second evaluation (COOST2):")
  spc04(fb, c(1:9, NA), "COOT2", "- Out of range values for cooked taste second evaluation (COOT2):")
  spc04(fb, c(1:9, NA), "COOAP2", "- Out of range values for cooked appearance second evaluation (COOAP2):")
  
  # Extreme values detection and values out of range for lab data
  
  spc05(fb, "lower", "PROT", "- Out of range values for protein (PROT):")
  spc06(fb, f, "low", "PROT", "- Extreme low values for protein (PROT):")
  spc06(fb, f, "high", "PROT", "- Extreme high values for protein (PROT):")
  spc05(fb, "lower", "FE", "- Out of range values for iron (FE):")
  spc06(fb, f, "low", "FE", "- Extreme low values for iron (FE):")
  spc06(fb, f, "high", "FE", "- Extreme high values for iron (FE):")
  spc05(fb, "lower", "ZN", "- Out of range values for zinc (ZN):")
  spc06(fb, f, "low", "ZN", "- Extreme low values for zinc (ZN):")
  spc06(fb, f, "high", "ZN", "- Extreme high values for zinc (ZN):")
  spc05(fb, "lower", "CA", "- Out of range values for calcium (CA):")
  spc06(fb, f, "low", "CA", "- Extreme low values for calcium (CA):")
  spc06(fb, f, "high", "CA", "- Extreme high values for calcium (CA):")
  spc05(fb, "lower", "MG", "- Out of range values for magnesium (MG):")
  spc06(fb, f, "low", "MG", "- Extreme low values for magnesium (MG):")
  spc06(fb, f, "high", "MG", "- Extreme high values for magnesium (MG):")
  spc05(fb, "lower", "BC", "- Out of range values for beta-carotene (BC):")
  spc06(fb, f, "low", "BC", "- Extreme low values for beta-carotene (BC):")
  spc06(fb, f, "high", "BC", "- Extreme high values for beta-carotene (BC):")
  bc.cc.values <- c(0.03, 0, 0.12, 0.02, 0.15, 1.38, 1.65, 1.5, 1.74, 1.76, 0.69, 1.17, 1.32,
                    1.04, 4.41, 4.92, 6.12, 5.46, 3.96, 5.49, 3.03, 3.76, 4.61, 7.23, 7.76,
                    10.5, 11.03, 12.39, 14.37, NA)
  spc04(fb, bc.cc.values, "BC.CC", "- Out of range values for beta-carotene (BC.CC):")
  spc05(fb, "lower", "TC", "- Out of range values for total carotenoids (TC):")
  spc06(fb, f, "low", "TC", "- Extreme low values for total carotenoids (TC):")
  spc06(fb, f, "high", "TC", "- Extreme high values for total carotenoids (TC):")
  spc05(fb, "lower", "STAR", "- Out of range values for starch (STAR):")
  spc06(fb, f, "low", "STAR", "- Extreme low values for starch (STAR):")
  spc06(fb, f, "high", "STAR", "- Extreme high values for starch (STAR):")
  spc05(fb, "lower", "FRUC", "- Out of range values for fructose (FRUC):")
  spc06(fb, f, "low", "FRUC", "- Extreme low values for fructose (FRUC):")
  spc06(fb, f, "high", "FRUC", "- Extreme high values for fructose (FRUC):")
  spc05(fb, "lower", "GLUC", "- Out of range values for glucose (GLUC):")
  spc06(fb, f, "low", "GLUC", "- Extreme low values for glucose (GLUC):")
  spc06(fb, f, "high", "GLUC", "- Extreme high values for glucose (GLUC):")
  spc05(fb, "lower", "SUCR", "- Out of range values for sucrose (SUCR):")
  spc06(fb, f, "low", "SUCR", "- Extreme low values for sucrose (SUCR):")
  spc06(fb, f, "high", "SUCR", "- Extreme high values for sucrose (SUCR):")
  spc05(fb, "lower", "MALT", "- Out of range values for maltose (MALT):")
  spc06(fb, f, "low", "MALT", "- Extreme low values for maltose (MALT):")
  spc06(fb, f, "high", "MALT", "- Extreme high values for maltose (MALT):")
  
  # Extreme values detection and values out of range for derived variables

  spc05(fb, "lower", "TRW", "- Out of range values for total root weight (TRW):")
  spc06(fb, f, "low", "TRW", "- Extreme low values for total root weight (TRW):")
  spc06(fb, f, "high", "TRW", "- Extreme high values for total root weight (TRW):")
  spc05(fb, "lower", "CYTHA", "- Out of range values for commercial root yield in tons per hectare (CYTHA):")
  spc06(fb, f, "low", "CYTHA", "- Extreme low values for commercial root yield in tons per hectare (CYTHA):")
  spc06(fb, f, "high", "CYTHA", "- Extreme high values for commercial root yield in tons per hectare (CYTHA):")
  spc05(fb, "lower", "CYTHA.AJ", "- Out of range values for commercial root yield in tons per hectare (CYTHA.AJ):")
  spc06(fb, f, "low", "CYTHA.AJ", "- Extreme low values for commercial root yield in tons per hectare (CYTHA.AJ):")
  spc06(fb, f, "high", "CYTHA", "- Extreme high values for commercial root yield in tons per hectare (CYTHA.AJ):")
  spc05(fb, "lower", "RYTHA", "- Out of range values for total root yield in tons per hectare (RYTHA):")
  spc06(fb, f, "low", "RYTHA", "- Extreme low values for total root yield in tons per hectare (RYTHA):")
  spc06(fb, f, "high", "RYTHA", "- Extreme high values for total root yield in tons per hectare (RYTHA):")
  spc05(fb, "lower", "RYTHA.AJ", "- Out of range values for total root yield in tons per hectare (RYTHA.AJ):")
  spc06(fb, f, "low", "RYTHA.AJ", "- Extreme low values for total root yield in tons per hectare (RYTHA.AJ):")
  spc06(fb, f, "high", "RYTHA.AJ", "- Extreme high values for total root yield in tons per hectare (RYTHA.AJ):")
  spc05(fb, "lower", "ACRW", "- Out of range values for average commercial root weight (ACRW):")
  spc06(fb, f, "low", "ACRW", "- Extreme low values for average commercial root weight (ACRW):")
  spc06(fb, f, "high", "ACRW", "- Extreme high values for average commercial root weight (ACRW):")
  spc05(fb, "lower", "NRPP", "- Out of range values for number of roots per plant (NRPP):")
  spc06(fb, f, "low", "NRPP", "- Extreme low values for number of roots per plant (NRPP):")
  spc06(fb, f, "high", "NRPP", "- Extreme high values for number of roots per plant (NRPP):")
  spc05(fb, "lower", "YPP", "- Out of range values for yield per plant (YPP):")
  spc06(fb, f, "low", "YPP", "- Extreme low values for yield per plant (YPP):")
  spc06(fb, f, "high", "YPP", "- Extreme high values for yield per plant (YPP):")
  spc05(fb, "both", "CI", "- Out of range values for commercial index (CI):")
  spc06(fb, f, "low", "CI", "- Extreme low values for commercial index (CI):")
  spc06(fb, f, "high", "CI", "- Extreme high values for commercial index (CI):")
  spc05(fb, "both", "HI", "- Out of range values for harvest index (HI):")
  spc06(fb, f, "low", "HI", "- Extreme low values for harvest index (HI):")
  spc06(fb, f, "high", "HI", "- Extreme high values for harvest index (HI):")
  spc05(fb, "both", "SHI", "- Out of range values for harvest sowing index (SHI):")
  spc06(fb, f, "low", "SHI", "- Extreme low values for harvest sowing index (SHI):")
  spc06(fb, f, "high", "SHI", "- Extreme high values for harvest sowing index (SHI):")
  spc05(fb, "lower", "BIOM", "- Out of range values for biomass yield (BIOM):")
  spc06(fb, f, "low", "BIOM", "- Extreme low values for biomass yield (BIOM):")
  spc06(fb, f, "high", "BIOM", "- Extreme high values for biomass yield (BIOM):")
  spc05(fb, "lower", "BIOM.AJ", "- Out of range values for biomass yield (BIOM.AJ):")
  spc06(fb, f, "low", "BIOM.AJ", "- Extreme low values for biomass yield (BIOM.AJ):")
  spc06(fb, f, "high", "BIOM.AJ", "- Extreme high values for biomass yield (BIOM.AJ):")
  spc05(fb, "lower", "FYTHA", "- Out of range values for foliage total yield in tons per hectare (FYTHA):")
  spc06(fb, f, "low", "FYTHA", "- Extreme low values for foliage total yield in tons per hectare (FYTHA):")
  spc06(fb, f, "high", "FYTHA", "- Extreme high values for foliage total yield in tons per hectare (FYTHA):")
  spc05(fb, "lower", "FYTHA.AJ", "- Out of range values for foliage total yield in tons per hectare (FYTHA.AJ):")
  spc06(fb, f, "low", "FYTHA.AJ", "- Extreme low values for foliage total yield in tons per hectare (FYTHA.AJ):")
  spc06(fb, f, "high", "FYTHA.AJ", "- Extreme high values for foliage total yield in tons per hectare (FYTHA.AJ):")
  spc05(fb, "both", "RFR", "- Out of range values for root foliage ratio (RFR):")
  spc06(fb, f, "low", "RFR", "- Extreme low values for root foliage ratio (RFR):")
  spc06(fb, f, "high", "RFR", "- Extreme high values for root foliage ratio (RFR):")

  # Extreme values detection and values out of range for additional traits
  
  if (!is.null(aqt)) {
    for (i in 1:length(aqt)) {
      spc06(fb, f, "low", aqt[i], paste("- Extreme low values for (", aqt[i], "):", sep = ""))
      spc06(fb, f, "high", aqt[i], paste("- Extreme high values for (", aqt[i], "):", sep = ""))
    }
  }
  
  # Outliers' detection
  
  oc <- 0
  
  if (out.mod == "rcbd") {
    oc <- 1
    if (exists('G', fb)) {
      geno <- fb[, 'G']
    } else {
      if (exists('GENO', fb)) {
        geno <- fb[, 'GENO']
      } else {
        if (exists('NAME', fb)) {
          geno <- fb[, 'NAME']
        } else {
          oc <- 0
          warning("Genotypes are not defined. Use G or GENO as labels.", call. = FALSE)  
        }
      }
    }
    if (exists('R', fb)) {
      rep <- fb[, 'R']
    } else {
      if (exists('REP', fb)) {
        rep <- fb[, 'REP']
      } else {
        oc <- 0
        warning('Blocks are not defined. Use R or REP as labels.', call. = FALSE)
      }
    }
    env <- NULL
  }
  
  if (out.mod == "met") {
    oc <- 1
    if (exists('G', fb)) {
      geno <- fb[, 'G']
    } else {
      if (exists('GENO', fb)) {
        geno <- fb[, 'GENO']
      } else {
        if (exists('NAME', fb)) {
          geno <- fb[, 'NAME']
        } else {
          oc <- 0
          warning("Genotypes are not defined. Use G or GENO as labels.", call. = FALSE)  
        }
      }
    }
    if (exists('E', fb)) {
      env <- fb[, 'E']
    } else {
      if (exists('ENV', fb)) {
        env <- fb[, 'ENV']
      } else {
        oc <- 0
        warning('Environments are not defined. Use E or ENV as labels.', call. = FALSE)
      }
    }
    if (exists('R', fb)) {
      rep <- fb[, 'R']
    } else {
      if (exists('REP', fb)) {
        rep <- fb[, 'REP']
      } else {
        oc <- 0
        warning('Blocks are not defined. Use R or REP as labels.', call. = FALSE)
      }
    }
  }
  
  if (oc == 1) {
    spc07(fb, geno, env, rep, "VW", out.mod, out.max, "- Outliers for vine weight (VW):")
    spc07(fb, geno, env, rep, "NOCR", out.mod, out.max, "- Outliers for number of commercial roots (NOCR):")
    spc07(fb, geno, env, rep, "NONC", out.mod, out.max, "- Outliers for number of non commercial roots (NONC):")
    spc07(fb, geno, env, rep, "CRW", out.mod, out.max, "- Outliers for commercial root weight (CRW):")
    spc07(fb, geno, env, rep, "NCRW", out.mod, out.max, "- Outliers for non commercial root weight (NCRW):")
    spc07(fb, geno, env, rep, "DM", out.mod, out.max, "- Outliers for storage root dry matter content (DM):")
    spc07(fb, geno, env, rep, "DMRY", out.mod, out.max, "- Outliers for dry matter root yield (DMRY):")
    spc07(fb, geno, env, rep, "DMV", out.mod, out.max, "- Outliers for vines dry matter content (DMV):")
    spc07(fb, geno, env, rep, "DMFY", out.mod, out.max, "- Outliers for dry matter foliage yield (DMFY):")
    spc07(fb, geno, env, rep, "PROT", out.mod, out.max, "- Outliers for protein (PROT):")
    spc07(fb, geno, env, rep, "FE", out.mod, out.max, "- Outliers for iron (FE):")
    spc07(fb, geno, env, rep, "ZN", out.mod, out.max, "- Outliers for zinc (ZN):")
    spc07(fb, geno, env, rep, "CA", out.mod, out.max, "- Outliers for calcium (CA):")
    spc07(fb, geno, env, rep, "MG", out.mod, out.max, "- Outliers for magnesium (MG):")
    spc07(fb, geno, env, rep, "BC", out.mod, out.max, "- Outliers for beta-carotene (BC):")
    spc07(fb, geno, env, rep, "BC.CC", out.mod, out.max, "- Outliers for beta-carotene (BC.CC):")
    spc07(fb, geno, env, rep, "TC", out.mod, out.max, "- Outliers for total carotenoids (TC):")
    spc07(fb, geno, env, rep, "STAR", out.mod, out.max, "- Outliers for starch (STAR):")
    spc07(fb, geno, env, rep, "FRUC", out.mod, out.max, "- Outliers for fructose (FRUC):")
    spc07(fb, geno, env, rep, "GLUC", out.mod, out.max, "- Outliers for glucose (GLUC):")
    spc07(fb, geno, env, rep, "SUCR", out.mod, out.max, "- Outliers for sucrose (SUCR):")
    spc07(fb, geno, env, rep, "MALT", out.mod, out.max, "- Outliers for maltose (MALT):")
    spc07(fb, geno, env, rep, "TRW", out.mod, out.max, "- Outliers for total root weight (TRW):")
    spc07(fb, geno, env, rep, "CYTHA", out.mod, out.max, "- Outliers for commercial root yield in tons per hectare (CYTHA):")
    spc07(fb, geno, env, rep, "CYTHA.AJ", out.mod, out.max, "- Outliers for commercial root yield in tons per hectare (CYTHA.AJ):")
    spc07(fb, geno, env, rep, "RYTHA", out.mod, out.max, "- Outliers for total root yield in tons per hectare (RYTHA):")
    spc07(fb, geno, env, rep, "RYTHA.AJ", out.mod, out.max, "- Outliers for total root yield in tons per hectare (RYTHA.AJ):")
    spc07(fb, geno, env, rep, "ACRW", out.mod, out.max, "- Outliers for average commercial root weight (ACRW):")
    spc07(fb, geno, env, rep, "NRPP", out.mod, out.max, "- Outliers for number of roots per plant (NRPP):")
    spc07(fb, geno, env, rep, "YPP", out.mod, out.max, "- Outliers for yield per plant (YPP):")
    spc07(fb, geno, env, rep, "CI", out.mod, out.max, "- Outliers for commercial index (CI):")
    spc07(fb, geno, env, rep, "HI", out.mod, out.max, "- Outliers for harvest index (HI):")
    spc07(fb, geno, env, rep, "SHI", out.mod, out.max, "- Outliers for harvest sowing index (SHI):")
    spc07(fb, geno, env, rep, "BIOM", out.mod, out.max, "- Outliers for biomass yield (BIOM):")
    spc07(fb, geno, env, rep, "BIOM.AJ", out.mod, out.max, "- Outliers for biomass yield (BIOM.AJ):")
    spc07(fb, geno, env, rep, "FYTHA", out.mod, out.max, "- Outliers for foliage total yield in tons per hectare (FYTHA):")
    spc07(fb, geno, env, rep, "FYTHA.AJ", out.mod, out.max, "- Outliers for foliage total yield in tons per hectare (FYTHA.AJ):")
    spc07(fb, geno, env, rep, "RFR", out.mod, out.max, "- Outliers for root foliage ratio (RFR):")
    
    # Outliers' detection for additional traits
    
    if (!is.null(aqt))
      for (i in 1:length(aqt))
        spc07(fb, geno, env, rep, aqt[i], out.mod, out.max, paste("- Outlilers for (", aqt[i], "):", sep = ""))
    
  }
  
  if (file == TRUE) sink()
}

# Print results
output <- function(fb, cond, tx) {
  if (sum(cond, na.rm = TRUE) > 0) {
    cat("\n", tx, "\n", sep = "")
    print(subset(fb, cond))
  }
}

# Two traits conditions
spc01 <- function(fb, type, t1, t2, tx) {
  if (exists(t1, fb) & exists(t2, fb)) {
    if (type == 1)
      cond <- fb[, t1] > fb[, t2]
    if (type == 2)
      cond <- (fb[, t1] == 0 | is.na(fb[, t1])) & !is.na(fb[, t2])
    if (type == 3)
      cond <- (fb[, t1] == 0 | is.na(fb[, t1])) & fb[, t2] > 0
    output(fb, cond, tx)
  }
}

# Three traits conditions
spc02 <- function(fb, type, t1, t2, t3, tx) {
  if (exists(t1, fb) & exists(t2, fb) & exists(t3, fb)) {
    if (type == 1)
      cond <- (fb[, t1] == 0 | is.na(fb[, t1])) & (fb[, t2] > 0 | fb[, t3] > 0)
    if (type == 2)
      cond <- fb[, t1] > 0 & (suma(fb[, t2], fb[, t3]) == 0 | is.na(suma(fb[, t2], fb[, t3])))
    if (type == 3)
      cond <- (fb[, t1] == 0 | is.na(fb[, t1])) & suma(fb[, t2], fb[, t3]) > 0
    output(fb, cond, tx)
  }
}

# Roots and dependencies
spc03 <- function(fb, temp, do, t1, tx) {
  if (exists(t1, fb) & do == TRUE) {
    cond <- (temp == FALSE) & (fb[, t1] > 0)
    output(fb, cond, tx)
  }
}

# Detect out of discrete range
spc04 <- function(fb, vv, t1, tx) {
  if (exists(t1, fb)) {
    cond <- !(fb[, t1] %in% vv)
    output(fb, cond, tx)
  }
}

# Detect out of continous range
spc05 <- function(fb, ex, t1, tx) { 
  if (exists(t1, fb)) {
    if (ex == "lower")
      cond <- fb[, t1] < 0
    if (ex == "both")
      cond <- fb[, t1] < 0 | fb[, t1] > 100
    output(fb, cond, tx)
  }
}

# Extreme values
spc06 <- function(fb, f, ex, t1, tx) { 
  if (exists(t1, fb)) {
    if (ex == "low")
      cond <- fb[, t1] < quantile(fb[, t1], 0.25, na.rm = TRUE) - f * IQR(fb[, t1], na.rm = TRUE)
    if (ex == "high")
      cond <- fb[, t1] > quantile(fb[, t1], 0.75, na.rm = TRUE) + f * IQR(fb[, t1], na.rm = TRUE)
    output(fb, cond, tx)
  } 
}

# Outliers' detection
spc07 <- function(fb, geno, env, rep, t1, out.mod, out.max, tx) {
  if (exists(t1, fb)) {
    fb$id <- as.numeric(row.names(fb))
    if (out.mod == "rcbd")
      model <- aov(fb[, t1] ~ geno + rep)
    if (out.mod == "met")
      model <- aov(fb[, t1] ~ geno + env + rep %in% env + geno:env)
    res <- data.frame(residual = rstandard(model))
    res$id <- as.numeric(row.names(res))
    fb <- merge(fb, res, all = T)[, -1]
    cond <- abs(fb[, 'residual']) > out.max
    output(fb, cond, tx)
  }
}
