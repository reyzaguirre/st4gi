#' Check consistency for sweetpotato experimental data
#'
#' Set of rules to check for consistency of sweetpotato experimental data.
#' Data labels must be defined as specified in the PROCEDURES FOR THE EVALUATION
#' AND ANALYSIS OF SWEETPOTATO TRIALS document.
#' @param fb The name of the fieldbook data frame.
#' @param plot.size Plot size in square meters.
#' @param f Factor for extreme values detection. See details.
#' @param width Number of columns for the output.
#' @param file Logigal, if TRUE the output goes to a file with name output.txt.
#' @details The data frame must use the labels (lower or upper case) listed in function \code{checknames}.
#' See \code{?checknames} for details.
#' 
#' \code{Plot.size} must be specified to check if traits expressed in tons per
#' hectare have been correctly computed from kilograms per plot. If \code{NULL} this
#' is not verified.
#' 
#' Extreme values are detected using the interquartile range.
#' The rule is to detect any value out of the interval 
#' \eqn{[Q_1 - f \times IQR; Q_3 + f \times IQR]}. By default \code{f = 3}.
#' @return If \code{file = TRUE} it returns a file with name checks.txt with a list of
#' all rows with some kind of inconsistency and all rows with outliers. If \code{file = FALSE}
#' the output is shown in the R console.
#' @author Raul Eyzaguirre.
#' @examples
#' spconsis(pjpz09)
#' @importFrom stats IQR quantile
#' @export

spconsis <- function(fb, plot.size = NULL, f = 3, width = 240, file = FALSE) {

  options(width = width)
  
  fb <- checknames(fb)

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
  spc01(fb, 2, "NOPE", "VV1", "- Number of plants established (NOPE) is zero or NA but there is data for vine vigor first evaluation (VV1):")
  
  # NOPH and VW
  
  spc01(fb, 3, "NOPH", "VW", "- Number of plants harvested (NOPH) is zero or NA but vine weight (VW) is greater than zero:") 
  spc01(fb, 3, "VW", "NOPH", "- Vine weight (VW) is zero or NA but the number of plants harvested (NOPH) is greater than zero:") 
  spc01(fb, 3, "NOPH", "FYTHA", "- Number of plants harvested (NOPH) is zero or NA but foliage yield in tons per hectare (FYTHA) is greater than zero:") 
  spc01(fb, 3, "FYTHA", "NOPH", "- Foliage yield in tons per hectare (FYTHA) is zero or NA but the number of plants harvested (NOPH) is greater than zero:") 
  
  # VW and dependencies
  
  spc01(fb, 3, "VW", "DMVF", "- Vine weight (VW) is zero or NA but there is fresh weight vines for dry matter assessment (DMVF):") 
  spc01(fb, 3, "VW", "DMVD", "- Vine weight (VW) is zero or NA but there is dry weight vines for dry matter assessment (DMVD):") 
  spc01(fb, 1, "DMVD", "DMVF", "- Dry weight vines for dry matter assessment (DMVD) is greater than fresh weight vines for dry matter assessment (DBVF):")
  spc01(fb, 2, "VW", "VV2", "- Vine weight (VW) is zero or NA but there is data for vine vigor second evaluation (VV2):")
  spc01(fb, 2, "VW", "VIR3", "- Vine weight (VW) is zero or NA but there is data for virus symptoms third evaluation (VIR3):")
  
  # NOPR and roots
  
  spc02(fb, 1, "NOPR", "NOCR", "NONC", "- Number of plants with roots (NOPR) is zero or NA but number of roots (NOCR + NONC) is greater than zero:")
  spc02(fb, 2, "NOPR", "NOCR", "NONC", "- Number of roots (NOCR + NONC) is zero or NA but number of plants with roots (NOPR) is greater than zero:")
  spc02(fb, 1, "NOPR", "CRW", "NCRW", "- Number of plants with roots (NOPR) is zero or NA but root weight (CRW + NCRW) is greater than zero:")
  spc02(fb, 2, "NOPR", "CRW", "NCRW", "- Root weight (CRW + NCRW) is zero or NA but number of plants with roots (NOPR) is greater than zero:")
  spc01(fb, 3, "NOPR", "TRW", "- Number of plants with roots (NOPR) is zero or NA but total root weight (TRW) is greater than zero::")
  spc01(fb, 3, "TRW", "NOPR", "- Total root weight (TRW) is zero or NA but number of plants with roots (NOPR) is greater than zero:")
  spc01(fb, 3, "NOPR", "RYTHA", "- Number of plants with roots (NOPR) is zero or NA but root yield in tons per hectare (rytha) is greater than zero::")
  spc01(fb, 3, "RYTHA", "NOPR", "- Root yield in tons per hectare (RYTHA) is zero or NA but number of plants with roots (NOPR) is greater than zero:")
  
  # Number of roots and root weight
  
  spc01(fb, 3, "NOCR", "CRW", "- Number of commercial roots (NOCR) is zero or NA but the commercial root weight (CRW) is greater than zero:") 
  spc01(fb, 3, "CRW", "NOCR", "- Commercial root weight (CRW) is zero or NA but the number of commercial roots (NOCR) is greater than zero:") 
  spc01(fb, 3, "NONC", "NCRW", "- Number of non commercial roots (NONC) is zero or NA but the non commercial root weight (NCRW) is greater than zero:")
  spc01(fb, 3, "NCRW", "NONC", "- Non commercial root weight (NCRW) is zero or NA but the number of non commercial roots (NONC) is greater than zero:")
  spc02(fb, 3, "TRW", "NOCR", "NONC", "- Total root weight (TRW) is zero or NA but number of roots (NOCR + NONC) is greater than zero:")
  spc02(fb, 2, "TRW", "NOCR", "NONC", "- Number of roots (NOCR + NONC) is zero or NA but total root weight (TRW) is greater than zero:")
  
  # Roots and dependencies

  temp <- array(FALSE, dim(fb)[1])
  
  if (exists("NOPR", fb))
    temp <- temp | (!is.na(fb$NOPR) & fb$NOPR > 0)
  if (exists("NOCR", fb))
    temp <- temp | (!is.na(fb$NOCR) & fb$NOCR > 0)
  if (exists("NONC", fb))
    temp <- temp | (!is.na(fb$NONC) & fb$NONC > 0)
  if (exists("CRW", fb))
    temp <- temp | (!is.na(fb$CRW) & fb$CRW > 0)
  if (exists("NCRW", fb))
    temp <- temp | (!is.na(fb$NCRW) & fb$NCRW > 0)
  if (exists("TRW", fb))
    temp <- temp | (!is.na(fb$TRW) & fb$TRW > 0)
  if (exists("RYTHA", fb))
    temp <- temp | (!is.na(fb$RYTHA) & fb$RYTHA > 0)
  if (exists("CYTHA", fb))
    temp <- temp | (!is.na(fb$CYTHA) & fb$CYTHA > 0)
  
  spc03(fb, temp, "RFCP", "- There are no roots but there is data for root primary flesh color (RFCP):")
  spc03(fb, temp, "RFCS", "- There are no roots but there is data for root secondary flesh color (RFCS):")
  spc03(fb, temp, "SCOL", "- There are no roots but there is data for storage root skin color (SCOL):")
  spc03(fb, temp, "FCOL", "- There are no roots but there is data for storage root flesh color (FCOL):")
  spc03(fb, temp, "RS", "- There are no roots but there is data for root size (RS):")
  spc03(fb, temp, "RF", "- There are no roots but there is data for root form (RF):")
  spc03(fb, temp, "DAMR", "- There are no roots but there is data for root defects (DAMR):")
  spc03(fb, temp, "RSPR", "- There are no roots but there is data for root sprouting (RSPR):")
  spc03(fb, temp, "WED", "- There are no roots but there is data for weevil damage (WED):")
  spc03(fb, temp, "WED1", "- There are no roots but there is data for weevil damage first evaluation (WED1):")
  spc03(fb, temp, "WED2", "- There are no roots but there is data for weevil damage second evaluation (WED2):")
  spc03(fb, temp, "DMF", "- There are no roots but there is data for fresh weight of roots for dry matter assessment (DMF):")
  spc03(fb, temp, "DMD", "- There are no roots but there is data for dry weight of roots for dry matter assessment (DMD):")
  spc01(fb, 1, "DMD", "DMF", "- Dry weight of roots for dry matter assessment (DMD) is greater than fresh weight of roots for dry matter assessment (DMF):")
  spc03(fb, temp, "FRAW", "- There are no roots but there is data for root fiber (FRAW):")
  spc03(fb, temp, "SURAW", "- There are no roots but there is data for root sugar (SURAW):")
  spc03(fb, temp, "STRAW", "- There are no roots but there is data for root starch (STRAW):")
  spc03(fb, temp, "COOF", "- There are no roots but there is data for cooked fiber (COOF):")
  spc03(fb, temp, "COOSU", "- There are no roots but there is data for cooked sugars (COOSU):")
  spc03(fb, temp, "COOST", "- There are no roots but there is data for cooked starch (COOST):")
  spc03(fb, temp, "COOT", "- There are no roots but there is data for cooked taste (COOT):")
  spc03(fb, temp, "COOAP", "- There are no roots but there is data for cooked appearance (COOAP):")
  spc03(fb, temp, "FRAW1", "- There are no roots but there is data for root fiber first determination (FRAW1):")
  spc03(fb, temp, "SURAW1", "- There are no roots but there is data for root sugar first determination (SURAW1):")
  spc03(fb, temp, "STRAW1", "- There are no roots but there is data for root starch first determination (STRAW1):")
  spc03(fb, temp, "COOF1", "- There are no roots but there is data for cooked fiber first evaluation (COOF1):")
  spc03(fb, temp, "COOSU1", "- There are no roots but there is data for cooked sugars first evaluation (COOSU1):")
  spc03(fb, temp, "COOST1", "- There are no roots but there is data for cooked starch first evaluation (COOST1):")
  spc03(fb, temp, "COOT1", "- There are no roots but there is data for cooked taste first evaluation (COOT1):")
  spc03(fb, temp, "COOAP1", "- There are no roots but there is data for cooked appearance first evaluation (COOAP1):")
  spc03(fb, temp, "FRAW2", "- There are no roots but there is data for root fiber second determination (FRAW2):")
  spc03(fb, temp, "SURAW2", "- There are no roots but there is data for root sugar second determination (SURAW2):")
  spc03(fb, temp, "STRAW2", "- There are no roots but there is data for root starch second determination (STRAW2):")
  spc03(fb, temp, "COOF2", "- There are no roots but there is data for cooked fiber second evaluation (COOF2):")
  spc03(fb, temp, "COOSU2", "- There are no roots but there is data for cooked sugars second evaluation (COOSU2):")
  spc03(fb, temp, "COOST2", "- There are no roots but there is data for cooked starch second evaluation (COOST2):")
  spc03(fb, temp, "COOT2", "- There are no roots but there is data for cooked taste second evaluation (COOT2):")
  spc03(fb, temp, "COOAP2", "- There are no roots but there is data for cooked appearance second evaluation (COOAP2):")
  spc03(fb, temp, "PROT", "- There are no roots but there is data for protein (PROT):")
  spc03(fb, temp, "FE", "- There are no roots but there is data for iron in dry weight (FE):")
  spc03(fb, temp, "ZN", "- There are no roots but there is data for zinc in dry weight (ZN):")
  spc03(fb, temp, "CA", "- There are no roots but there is data for calcium in dry weight (CA):")
  spc03(fb, temp, "MG", "- There are no roots but there is data for magnesium in dry weight (MG):")
  spc03(fb, temp, "BC", "- There are no roots but there is data for beta-carotene in dry weight (BC):")
  spc03(fb, temp, "TC", "- There are no roots but there is data for total carotenoids in dry weight (TC):")
  spc03(fb, temp, "STAR", "- There are no roots but there is data for starch (STAR):")
  spc03(fb, temp, "FRUC", "- There are no roots but there is data for fructose (FRUC):")
  spc03(fb, temp, "GLUC", "- There are no roots but there is data for glucose (GLUC):")
  spc03(fb, temp, "SUCR", "- There are no roots but there is data for sucrose (SUCR):")
  spc03(fb, temp, "MALT", "- There are no roots but there is data for maltose (MALT):")
  
  # Calculated variables
  
  spc05(fb, 1, plot.size, "TRW", "CRW", "NCRW", "- Total root weight (TRW) different from CRW + NCRW:")
  spc04(fb, plot.size, "CYTHA", "CRW", "- Commercial root yield in tons per hectare (CYTHA) is different from CRW * 10 / plot.size:")
  spc05(fb, 2, plot.size, "RYTHA", "CRW", "NCRW", "- Total root yield in tons per hectare (RYTHA) is different from (CRW + NCRW) * 10 / plot.size:")
  spc05(fb, 3, plot.size, "ACRW", "CRW", "NOCR", "- Average commercial root weight (ACRW) is different from CRW / NOCR:")
  spc06(fb, 1, plot.size, "NRPP", "NOCR", "NONC", "NOPH", "- Number of roots per plant (NRPP) is different from (NOCR + NONC) / NOPH:")
  spc06(fb, 1, plot.size, "YPP", "CRW", "NCRW", "NOPH", "- Yield per plant (YPP) is different from (CRW + NCRW) / NOPH:")
  spc05(fb, 4, plot.size, "CI", "NOCR", "NONC", "- Commercial index (CI) is different from NOCR / (NOCR + NONC) * 100:")
  spc06(fb, 2, plot.size, "HI", "CRW", "NCRW", "VW", "- Harvest index (HI) is different from (CRW + NCRW) / (VW + CRW + NCRW) * 100:")
  spc05(fb, 5, plot.size, "SHI", "NOPH", "NOPS", "- Harvest sowing index (SHI) is different from NOPH / NOPS * 100:")
  spc06(fb, 3, plot.size, "BIOM", "CRW", "NCRW", "VW", "- Biomass yield (BIOM) is different from (CRW + NCRW + VW) * 10 / plot.size:")
  spc04(fb, plot.size, "FYTHA", "VW", "- Foliage total yield in tons per hectare (FYTHA) is different from VW * 10 / plot.size:")
  spc05(fb, 5, plot.size, "DM", "DMD", "DMF", "- Storage root dry matter content (DM) is different from DMD / DMF * 100:")
  spc05(fb, 5, plot.size, "DMV", "DMVD", "DMVF", "- Vine dry matter content (DMV) is different from DMVD / DMVF * 100:")
  spc06(fb, 4, plot.size, "DMFY", "VW", "DMVD", "DMVF", "- Dry matter foliage yield (DMFY) is different from VW * 10 / plot.size * DMVD / DMVF:")
  spc07(fb, plot.size, "DMRY", "CRW", "NCRW", "DMD", "DMF", "- Dry matter root yield (DMRY) is different from (CRW + NCRW) * 10 / plot.size * DMD / DMF:")
  spc08(fb, "RFR", "CRW", "NCRW", "DMD", "DMF", "VW", "DMVD", "DMVF", "- Root foliage ratio (RFR) is different from (CRW + NCRW) * (DMD / DMF) / (VW * DMVD / DMVF) * 100:")
  
  # Outliers detection and values out of range for field data

  spc09(fb, c(0:100, NA), "NOPE", "- Out of range values for number of plants established (NOPE):")
  spc09(fb, c(1:9, NA), "VIR", "- Out of range values for virus symptoms (VIR):")
  spc09(fb, c(1:9, NA), "VIR1", "- Out of range values for virus symptoms first evaluation (VIR1):")
  spc09(fb, c(1:9, NA), "VIR2", "- Out of range values for virus symptoms second evaluation (VIR2):")
  spc09(fb, c(1:9, NA), "VIR2", "- Out of range values for virus symptoms third evaluation (VIR3):")
  spc09(fb, c(1:9, NA), "ALT", "- Out of range values for alternaria symptoms (ALT):")
  spc09(fb, c(1:9, NA), "ALT1", "- Out of range values for alternaria symptoms first evaluation (ALT1):")
  spc09(fb, c(1:9, NA), "ALT2", "- Out of range values for alternaria symptoms second evaluation (ALT2):")
  spc09(fb, c(1:9, NA), "VV", "- Out of range values for vine vigor (VV):")
  spc09(fb, c(1:9, NA), "VV1", "- Out of range values for vine vigor first evaluation (VV1):")
  spc09(fb, c(1:9, NA), "VV2", "- Out of range values for vine vigor second evaluation (VV2):")
  spc10(fb, "lower", "VW", "- Out of range values for vine weight (VW):")
  spc11(fb, f, "low", "VW", "- Extreme low values for vine weight (VW):")
  spc11(fb, f, "high", "VW", "- Extreme high values for vine weight (VW):")
  spc10(fb, "lower", "NOPH", "- Out of range values for number of plants harvested (NOPH):")
  spc10(fb, "lower", "NOPR", "- Out of range values for number of plants with roots (NOPR):")
  spc10(fb, "lower", "NOCR", "- Out of range values for number of commercial roots (NOCR):")
  spc11(fb, f, "low", "NOCR", "- Extreme low values for number of commercial roots (NOCR):")
  spc11(fb, f, "high", "NOCR", "- Extreme high values for number of commercial roots (NOCR):")
  spc10(fb, "lower", "NONC", "- Out of range values for number of non commercial roots (NONC):")
  spc11(fb, f, "low", "NONC", "- Extreme low values for number of non commercial roots (NONC):")
  spc11(fb, f, "high", "NONC", "- Extreme high values for number of non commercial roots (NONC):")
  spc10(fb, "lower", "CRW", "- Out of range values for commercial root weight (CRW):")
  spc11(fb, f, "low", "CRW", "- Extreme low values for commercial root weight (CRW):")
  spc11(fb, f, "high", "CRW", "- Extreme high values for commercial root weight (CRW):")
  spc10(fb, "lower", "NCRW", "- Out of range values for non commercial root weight (NCRW):")
  spc11(fb, f, "low", "NCRW", "- Extreme low values for non commercial root weight (NCRW):")
  spc11(fb, f, "high", "NCRW", "- Extreme high values for non commercial root weight (NCRW):")
  spc09(fb, c(1:9, NA), "SCOL", "- Out of range values for storage root skin color (SCOL):")
  spc09(fb, c(1:9, NA), "FCOL", "- Out of range values for storage root flesh color (FCOL):")
  spc09(fb, c(1:9, NA), "RS", "- Out of range values for root size (RS):")
  spc09(fb, c(1:9, NA), "RF", "- Out of range values for root form (RF):")
  spc09(fb, c(1:9, NA), "DAMR", "- Out of range values for root defects (DAMR):")
  spc09(fb, c(1:9, NA), "RSPR", "- Out of range values for root sprouting (RSPR):")
  spc09(fb, c(1:9, NA), "WED", "- Out of range values for weevil damage (WED):")
  spc09(fb, c(1:9, NA), "WED1", "- Out of range values for weevil damage first evaluation (WED1):")
  spc09(fb, c(1:9, NA), "WED2", "- Out of range values for weevil damage second evaluation (WED2):")
  
  # Outliers detection and values out of range for DM data
  
  spc10(fb, "lower", "DMF", "- Out of range values for fresh weight of roots for dry matter assessment (DMF):")
  spc11(fb, f, "low", "DMF", "- Extreme low values for fresh weight of roots for dry matter assessment (DMF):")
  spc11(fb, f, "high", "DMF", "- Extreme high values for fresh weight of roots for dry matter assessment (DMF):")
  spc10(fb, "lower", "DMD", "- Out of range values for dry weight of roots for dry matter assessment (DMD):")
  spc11(fb, f, "low", "DMD", "- Extreme low values for dry weight of roots for dry matter assessment (DMD):")
  spc11(fb, f, "high", "DMD", "- Extreme high values for dry weight of roots for dry matter assessment (DMD):")
  spc10(fb, "lower", "DMVF", "- Out of range values for fresh weight vines for dry matter assessment (DMVF):")
  spc11(fb, f, "low", "DMVF", "- Extreme low values for fresh weight of vines for dry matter assessment (DMVF):")
  spc11(fb, f, "high", "DMVF", "- Extreme high values for fresh weight of vines for dry matter assessment (DMVF):")
  spc10(fb, "lower", "DMVD", "- Out of range values for dry weight of vines for dry matter assessment (DMVD):")
  spc11(fb, f, "low", "DMVD", "- Extreme low values for dry weight of vines for dry matter assessment (DMVD):")
  spc11(fb, f, "high", "DMVD", "- Extreme high values for dry weight of vines for dry matter assessment (DMVD):")
  spc10(fb, "lower", "DM", "- Out of range values for storage root dry matter content (DM):")
  spc11(fb, f, "low", "DM", "- Extreme low values for storage root dry matter content (DM):")
  spc11(fb, f, "high", "DM", "- Extreme high values for storage root dry matter content (DM):")
  spc10(fb, "lower", "DMV", "- Out of range values for vine dry matter content (DMV):")
  spc11(fb, f, "low", "DMV", "- Extreme low values for vine dry matter content (DMV):")
  spc11(fb, f, "high", "DMV", "- Extreme high values for vine dry matter content (DMV):")
  spc10(fb, "lower", "DMFY", "- Out of range values for dry matter foliage yield (DMFY):")
  spc11(fb, f, "low", "DMFY", "- Extreme low values for dry matter foliage yield (DMFY):")
  spc11(fb, f, "high", "DMFY", "- Extreme high values for dry matter foliage yield (DMFY):")
  spc10(fb, "lower", "DMRY", "- Out of range values for dry matter root yield (DMRY):")
  spc11(fb, f, "low", "DMRY", "- Extreme low values for dry matter root yield (DMRY):")
  spc11(fb, f, "high", "DMRY", "- Extreme high values for dry matter root yield (DMRY):")
  
  # Outliers detection and values out of range for cooked traits
  
  spc09(fb, c(1:9, NA), "FRAW", "- Out of range values for root fiber (FRAW):")
  spc09(fb, c(1:9, NA), "SURAW", "- Out of range values for root sugar (SURAW):")
  spc09(fb, c(1:9, NA), "STRAW", "- Out of range values for root starch (STRAW):")
  spc09(fb, c(1:9, NA), "COOF", "- Out of range values for cooked fiber (COOF):")
  spc09(fb, c(1:9, NA), "COOSU", "- Out of range values for cooked sugars (COOSU):")
  spc09(fb, c(1:9, NA), "COOST", "- Out of range values for cooked starch (COOST):")
  spc09(fb, c(1:9, NA), "COOT", "- Out of range values for cooked taste (COOT):")
  spc09(fb, c(1:9, NA), "COOAP", "- Out of range values for cooked appearance (COOAP):")
  spc09(fb, c(1:9, NA), "FRAW1", "- Out of range values for root fiber first determination (FRAW1):")
  spc09(fb, c(1:9, NA), "SURAW1", "- Out of range values for root sugar first determination (SURAW1):")
  spc09(fb, c(1:9, NA), "STRAW1", "- Out of range values for root starch first determination (STRAW1):")
  spc09(fb, c(1:9, NA), "COOF1", "- Out of range values for cooked fiber first evaluation (COOF1):")
  spc09(fb, c(1:9, NA), "COOSU1", "- Out of range values for cooked sugars first evaluation (COOSU1):")
  spc09(fb, c(1:9, NA), "COOST1", "- Out of range values for cooked starch first evaluation (COOST1):")
  spc09(fb, c(1:9, NA), "COOT1", "- Out of range values for cooked taste first evaluation (COOT1):")
  spc09(fb, c(1:9, NA), "COOAP1", "- Out of range values for cooked appearance first evaluation (COOAP1):")
  spc09(fb, c(1:9, NA), "FRAW2", "- Out of range values for root fiber second determination (FRAW2):")
  spc09(fb, c(1:9, NA), "SURAW2", "- Out of range values for root sugar second determination (SURAW2):")
  spc09(fb, c(1:9, NA), "STRAW2", "- Out of range values for root starch second determination (STRAW2):")
  spc09(fb, c(1:9, NA), "COOF2", "- Out of range values for cooked fiber second evaluation (COOF2):")
  spc09(fb, c(1:9, NA), "COOSU2", "- Out of range values for cooked sugars second evaluation (COOSU2):")
  spc09(fb, c(1:9, NA), "COOST2", "- Out of range values for cooked starch second evaluation (COOST2):")
  spc09(fb, c(1:9, NA), "COOT2", "- Out of range values for cooked taste second evaluation (COOT2):")
  spc09(fb, c(1:9, NA), "COOAP2", "- Out of range values for cooked appearance second evaluation (COOAP2):")
  
  # Outliers detection and values out of range for lab data
  
  spc10(fb, "lower", "PROT", "- Out of range values for protein (PROT):")
  spc11(fb, f, "low", "PROT", "- Extreme low values for protein (PROT):")
  spc11(fb, f, "high", "PROT", "- Extreme high values for protein (PROT):")
  spc10(fb, "lower", "FE", "- Out of range values for iron in dry weight (FE):")
  spc11(fb, f, "low", "FE", "- Extreme low values for iron in dry weight (FE):")
  spc11(fb, f, "high", "FE", "- Extreme high values for iron in dry weight (FE):")
  spc10(fb, "lower", "ZN", "- Out of range values for zinc in dry weight (ZN):")
  spc11(fb, f, "low", "ZN", "- Extreme low values for zinc in dry weight (ZN):")
  spc11(fb, f, "high", "ZN", "- Extreme high values for zinc in dry weight (ZN):")
  spc10(fb, "lower", "CA", "- Out of range values for calcium in dry weight (CA):")
  spc11(fb, f, "low", "CA", "- Extreme low values for calcium in dry weight (CA):")
  spc11(fb, f, "high", "CA", "- Extreme high values for calcium in dry weight (CA):")
  spc10(fb, "lower", "MG", "- Out of range values for magnesium in dry weight (MG):")
  spc11(fb, f, "low", "MG", "- Extreme low values for magnesium in dry weight (MG):")
  spc11(fb, f, "high", "MG", "- Extreme high values for magnesium in dry weight (MG):")
  spc10(fb, "lower", "BC", "- Out of range values for beta-carotene in dry weight (BC):")
  spc11(fb, f, "low", "BC", "- Extreme low values for beta-carotene in dry weight (BC):")
  spc11(fb, f, "high", "BC", "- Extreme high values for beta-carotene in dry weight (BC):")
  bc.cc.values <- c(0.03, 0, 0.12, 0.02, 0.15, 1.38, 1.65, 1.5, 1.74, 1.76, 0.69, 1.17, 1.32,
                    1.04, 4.41, 4.92, 6.12, 5.46, 3.96, 5.49, 3.03, 3.76, 4.61, 7.23, 7.76,
                    10.5, 11.03, 12.39, 14.37, NA)
  spc09(fb, bc.cc.values, "BC.CC", "- Out of range values for beta-carotene with color chart (BC.CC):")
  spc10(fb, "lower", "TC", "- Out of range values for total carotenoids in dry weight (TC):")
  spc11(fb, f, "low", "TC", "- Extreme low values for total carotenoids in dry weight (TC):")
  spc11(fb, f, "high", "TC", "- Extreme high values for total carotenoids in dry weight (TC):")
  spc10(fb, "lower", "STAR", "- Out of range values for starch (STAR):")
  spc11(fb, f, "low", "STAR", "- Extreme low values for starch (STAR):")
  spc11(fb, f, "high", "STAR", "- Extreme high values for starch (STAR):")
  spc10(fb, "lower", "FRUC", "- Out of range values for fructose (FRUC):")
  spc11(fb, f, "low", "FRUC", "- Extreme low values for fructose (FRUC):")
  spc11(fb, f, "high", "FRUC", "- Extreme high values for fructose (FRUC):")
  spc10(fb, "lower", "GLUC", "- Out of range values for glucose (GLUC):")
  spc11(fb, f, "low", "GLUC", "- Extreme low values for glucose (GLUC):")
  spc11(fb, f, "high", "GLUC", "- Extreme high values for glucose (GLUC):")
  spc10(fb, "lower", "SUCR", "- Out of range values for sucrose (SUCR):")
  spc11(fb, f, "low", "SUCR", "- Extreme low values for sucrose (SUCR):")
  spc11(fb, f, "high", "SUCR", "- Extreme high values for sucrose (SUCR):")
  spc10(fb, "lower", "MALT", "- Out of range values for maltose (MALT):")
  spc11(fb, f, "low", "MALT", "- Extreme low values for maltose (MALT):")
  spc11(fb, f, "high", "MALT", "- Extreme high values for maltose (MALT):")
  
  # Outliers detection and values out of range for derived variables

  spc10(fb, "lower", "TRW", "- Out of range values for total root weight (TRW):")
  spc11(fb, f, "low", "TRW", "- Extreme low values for total root weight (TRW):")
  spc11(fb, f, "high", "TRW", "- Extreme high values for total root weight (TRW):")
  spc10(fb, "lower", "CYTHA", "- Out of range values for commercial root yield in tons per hectare (CYTHA):")
  spc11(fb, f, "low", "CYTHA", "- Extreme low values for commercial root yield in tons per hectare (CYTHA):")
  spc11(fb, f, "high", "CYTHA", "- Extreme high values for commercial root yield in tons per hectare (CYTHA):")
  spc10(fb, "lower", "RYTHA", "- Out of range values for total root yield in tons per hectare (RYTHA):")
  spc11(fb, f, "low", "RYTHA", "- Extreme low values for total root yield in tons per hectare (RYTHA):")
  spc11(fb, f, "high", "RYTHA", "- Extreme high values for total root yield in tons per hectare (RYTHA):")
  spc10(fb, "lower", "ACRW", "- Out of range values for average commercial root weight (ACRW):")
  spc11(fb, f, "low", "ACRW", "- Extreme low values for average commercial root weight (ACRW):")
  spc11(fb, f, "high", "ACRW", "- Extreme high values for average commercial root weight (ACRW):")
  spc10(fb, "lower", "NRPP", "- Out of range values for number of roots per plant (NRPP):")
  spc11(fb, f, "low", "NRPP", "- Extreme low values for number of roots per plant (NRPP):")
  spc11(fb, f, "high", "NRPP", "- Extreme high values for number of roots per plant (NRPP):")
  spc10(fb, "lower", "YPP", "- Out of range values for yield per plant (YPP):")
  spc11(fb, f, "low", "YPP", "- Extreme low values for yield per plant (YPP):")
  spc11(fb, f, "high", "YPP", "- Extreme high values for yield per plant (YPP):")
  spc10(fb, "both", "CI", "- Out of range values for commercial index (CI):")
  spc11(fb, f, "low", "CI", "- Extreme low values for commercial index (CI):")
  spc11(fb, f, "high", "CI", "- Extreme high values for commercial index (CI):")
  spc10(fb, "both", "HI", "- Out of range values for harvest index (HI):")
  spc11(fb, f, "low", "HI", "- Extreme low values for harvest index (HI):")
  spc11(fb, f, "high", "HI", "- Extreme high values for harvest index (HI):")
  spc10(fb, "both", "SHI", "- Out of range values for harvest sowing index (SHI):")
  spc11(fb, f, "low", "SHI", "- Extreme low values for harvest sowing index (SHI):")
  spc11(fb, f, "high", "SHI", "- Extreme high values for harvest sowing index (SHI):")
  spc10(fb, "lower", "BIOM", "- Out of range values for biomass yield (BIOM):")
  spc11(fb, f, "low", "BIOM", "- Extreme low values for biomass yield (BIOM):")
  spc11(fb, f, "high", "BIOM", "- Extreme high values for biomass yield (BIOM):")
  spc10(fb, "lower", "FYTHA", "- Out of range values for foliage total yield in tons per hectare (FYTHA):")
  spc11(fb, f, "low", "FYTHA", "- Extreme low values for foliage total yield in tons per hectare (FYTHA):")
  spc11(fb, f, "high", "FYTHA", "- Extreme high values for foliage total yield in tons per hectare (FYTHA):")
  spc10(fb, "both", "RFR", "- Out of range values for root foliage ratio (RFR):")
  spc11(fb, f, "low", "RFR", "- Extreme low values for root foliage ratio (RFR):")
  spc11(fb, f, "high", "RFR", "- Extreme high values for root foliage ratio (RFR):")
  
  if (file == TRUE) sink()
}

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
spc03 <- function(fb, temp, t1, tx) {
  if (exists(t1, fb)) {
    cond <- (temp == FALSE) & (fb[, t1] > 0)
    output(fb, cond, tx)
  }
}

# Formula for derived traits with two traits
spc04 <- function(fb, plot.size, t1, t2, tx) {
  if (exists(t1, fb) & exists(t2, fb) & !is.null(plot.size)) {
    cond <- abs(fb[, t1] - fb[, t2] * 10 / plot.size) > 1e-10
    output(fb, cond, tx)
  }
}

# Formula for derived traits with three traits
spc05 <- function(fb, dtf, plot.size, t1, t2, t3, tx) {
  if (exists(t1, fb) & exists(t2, fb) & exists(t3, fb) & !is.null(plot.size)) {
    if (dtf == 1)
      cond <- abs(fb[, t1] - suma(fb[, t2], fb[, t3])) > 1e-10
    if (dtf == 2)
      cond <- abs(fb[, t1] - suma(fb[, t2], fb[, t3]) * 10 / plot.size) > 1e-10
    if (dtf == 3)
      cond <- abs(fb[, t1] - fb[, t2] / fb[, t3]) > 1e-10
    if (dtf == 4)
      cond <- abs(fb[, t1] - fb[, t2] / suma(fb[, t2], fb[, t3]) * 100) > 1e-10
    if (dtf == 5)
      cond <- abs(fb[, t1] - fb[, t2] / fb[, t3] * 100) > 1e-10
    output(fb, cond, tx)
  }
}

# Formula for derived traits with four traits
spc06 <- function(fb, dtf, plot.size, t1, t2, t3, t4, tx) {
  if (exists(t1, fb) & exists(t2, fb) & exists(t3, fb) & exists(t4, fb) & !is.null(plot.size)) {
    if (dtf == 1)
      cond <- abs(fb[, t1] - suma(fb[, t2], fb[, t3]) / fb[, t4]) > 1e-10
    if (dtf == 2)
      cond <- abs(fb[, t1] -  suma(fb[, t2], fb[, t3]) / suma(suma(fb[, t2], fb[, t3]), fb[, t4]) * 100) > 1e-10
    if (dtf == 3)
      cond <- abs(fb[, t1] -  suma(suma(fb[, t2], fb[, t3]), fb[, t4]) / plot.size * 10) > 1e-10
    if (dtf == 4)
      cond <- abs(fb[, t1] - fb[, t2] * 10 / plot.size * fb[, t3] / fb[, t4]) > 1e-10
    output(fb, cond, tx)
  }
}

# Formula for derived traits with five traits
spc07 <- function(fb, plot.size, t1, t2, t3, t4, t5, tx) {
  if (exists(t1, fb) & exists(t2, fb) & exists(t3, fb) & exists(t4, fb) & exists(t5, fb) & !is.null(plot.size)) {
    cond <- abs(fb[, t1] - suma(fb[, t2], fb[, t3]) * 10 / plot.size * fb[, t4] / fb[, t5]) > 1e-10
    output(fb, cond, tx)
  }
}

# Formula for derived traits with eight traits
spc08 <- function(fb, t1, t2, t3, t4, t5, t6, t7, t8, tx) {
  if (exists(t1, fb) & exists(t2, fb) & exists(t3, fb) & exists(t4, fb) & exists(t5, fb)
      & exists(t6, fb) & exists(t7, fb) & exists(t8, fb)) {
    cond <- abs(fb[, t1] - suma(fb[, t2], fb[, t3]) * fb[, t4] / fb[, t5]
                / (fb[, t6] * fb[, t7] / fb[, t8])) * 100 > 1e-10
    output(fb, cond, tx)
  }
}

# Detect out of discrete range
spc09 <- function(fb, vv, t1, tx) {
  if (exists(t1, fb)) {
    cond <- !(fb[, t1] %in% vv)
    output(fb, cond, tx)
  }
}

# Detect out of continous range
spc10 <- function(fb, ex, t1, tx) { 
  if (exists(t1, fb)) {
    if (ex == "lower")
      cond <- fb[, t1] < 0
    if (ex == "both")
      cond <- fb[, t1] < 0 | fb[, t1] > 100
    output(fb, cond, tx)
  }
}

# Extreme values
spc11 <- function(fb, f, ex, t1, tx) { 
  if (exists(t1, fb)) {
    if (ex == "low")
      cond <- fb[, t1] < quantile(fb[, t1], 0.25, na.rm = TRUE) - f * IQR(fb[, t1], na.rm = TRUE)
    if (ex == "high")
      cond <- fb[, t1] > quantile(fb[, t1], 0.75, na.rm = TRUE) + f * IQR(fb[, t1], na.rm = TRUE)
    output(fb, cond, tx)
  }
}
