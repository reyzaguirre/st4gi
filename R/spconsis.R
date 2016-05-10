#' Check consistency for sweetpotato experimental data
#'
#' Set of rules to check for consistency of sweetpotato experimental data.
#' Data labels must be defined as specified in the PROCEDURES FOR THE EVALUATION
#' AND ANALYSIS OF SWEETPOTATO TRIALS document.
#' @param fb The name of the fieldbook data frame.
#' @param plot.size Plot size in square meters.
#' @param f Factor for extreme values detection. See details.
#' @param width Number of columns for the output.
#' @param file Logigal, if TRUE the output goes to a file.
#' @details The data frame must use the labels (lower or upper case) listed in function \code{checknames}.
#' See \code{?checknames} for details.
#' Extreme values are detected using the interquartile range.
#' The rule is to detect any value out of the interval 
#' \eqn{[Q_1 - f \times IQR; Q_3 + f \times IQR]}. By default \code{f = 3}.
#' @return If \code{file = TRUE} it returns a file with name checks.txt with a list of
#' all rows with some kind of inconsistency and all rows with outliers. If \code{file = FALSE}
#' the output is shown in the R console.
#' @author Raul Eyzaguirre.
#' @examples
#' spconsis(pjpz09, 4.5)
#' @export

spconsis <- function(fb, plot.size, f = 3, width = 240, file = TRUE) {

  options(width = width)
  
  fb <- checknames(fb)

  if (file == TRUE) sink(paste(getwd(), "/checks.txt", sep = ""))

  # Inconsistencies for NOPS > NOPE > NOPH > NOPR.

  spc01(fb, "NOPE", "NOPS", "- Number of plants established (NOPE) is greater than number of plants sowed (NOPS):")
  spc01(fb, "NOPH", "NOPS", "- Number of plants harvested (NOPH)  is greater than number of plants sowed (NOPS):")
  spc01(fb, "NOPR", "NOPS", "- Number of plants with roots (NOPR) is greater than number of plants sowed (NOPS):")
  spc01(fb, "NOPH", "NOPE", "- Number of plants harvested (NOPH) is greater than number of plants established (NOPE):")
  spc01(fb, "NOPR", "NOPE", "- Number of plants with roots (NOPR) is greater than number of plants established (NOPE):")
  spc01(fb, "NOPR", "NOPH", "- Number of plants with roots (NOPR) is greater than number of plants harvested (NOPH):")

  # Inconsistencies for NOPE and dependencies.

  spc02(fb, "NOPE", "VIR1", "- Number of plants established (NOPE) is zero or NA but there is data for virus symptoms first evaluation (VIR1):")
  spc02(fb, "NOPE", "VIR2", "- Number of plants established (NOPE) is zero or NA but there is data for virus symptoms second evaluation (VIR2):")
  spc02(fb, "NOPE", "ALT1", "- Number of plants established (NOPE) is zero or NA but there is data for alternaria symptoms first evaluation (ALT1):")
  spc02(fb, "NOPE", "ALT2", "- Number of plants established (NOPE) is zero or NA but there is data for alternaria symptoms second evaluation (ALT2):")
  spc02(fb, "NOPE", "VV1", "- Number of plants established (NOPE) is zero or NA but there is data for vine vigor first evaluation (VV1):")
  
  # NOPH and VW
  
  spc03(fb, "NOPH", "VW", "- Number of plants harvested (NOPH) is zero or NA but vine weight (VW) is greater than zero:") 
  spc03(fb, "VW", "NOPH", "- Vine weight (VW) is zero or NA but the number of plants harvested (NOPH) is greater than zero:") 
  spc03(fb, "NOPH", "FYTHA", "- Number of plants harvested (NOPH) is zero or NA but foliage yield in tons per hectare (FYTHA) is greater than zero:") 
  spc03(fb, "FYTHA", "NOPH", "- Foliage yield in tons per hectare (FYTHA) is zero or NA but the number of plants harvested (NOPH) is greater than zero:") 
  
  # VW and dependencies
  
  spc03(fb, "VW", "DMVF", "- Vine weight (VW) is zero or NA but there is fresh weight vines for dry matter assessment (DMVF):") 
  spc03(fb, "VW", "DMVD", "- Vine weight (VW) is zero or NA but there is dry weight vines for dry matter assessment (DMVD):") 
  spc01(fb, "DMVD", "DMVF", "- Dry weight vines for dry matter assessment (DMVD) is greater than fresh weight vines for dry matter assessment (DBVF):")
  spc02(fb, "VW", "VV2", "- Vine weight (VW) is zero or NA but there is data for vine vigor second evaluation (VV2):")
  spc02(fb, "VW", "VIR3", "- Vine weight (VW) is zero or NA but there is data for virus symptoms third evaluation (VIR3):")
  
  # NOPR and roots
  
  spc04(fb, "NOPR", "NOCR", "NONC", "- Number of plants with roots (NOPR) is zero or NA but number of roots (NOCR + NONC) is greater than zero:")
  spc05(fb, "NOPR", "NOCR", "NONC", "- Number of roots (NOCR + NONC) is zero or NA but number of plants with roots (NOPR) is greater than zero:")
  spc04(fb, "NOPR", "CRW", "NCRW", "- Number of plants with roots (NOPR) is zero or NA but root weight (CRW + NCRW) is greater than zero:")
  spc05(fb, "NOPR", "CRW", "NCRW", "- Root weight (CRW + NCRW) is zero or NA but number of plants with roots (NOPR) is greater than zero:")
  spc03(fb, "NOPR", "TRW", "- Number of plants with roots (NOPR) is zero or NA but total root weight (TRW) is greater than zero::")
  spc03(fb, "TRW", "NOPR", "- Total root weight (TRW) is zero or NA but number of plants with roots (NOPR) is greater than zero:")
  spc03(fb, "NOPR", "RYTHA", "- Number of plants with roots (NOPR) is zero or NA but root yield in tons per hectare (rytha) is greater than zero::")
  spc03(fb, "RYTHA", "NOPR", "- Root yield in tons per hectare (RYTHA) is zero or NA but number of plants with roots (NOPR) is greater than zero:")
  
  # Number of roots and root weight
  
  spc03(fb, "NOCR", "CRW", "- Number of commercial roots (NOCR) is zero or NA but the commercial root weight (CRW) is greater than zero:") 
  spc03(fb, "CRW", "NOCR", "- Commercial root weight (CRW) is zero or NA but the number of commercial roots (NOCR) is greater than zero:") 
  spc03(fb, "NONC", "NCRW", "- Number of non commercial roots (NONC) is zero or NA but the non commercial root weight (NCRW) is greater than zero:")
  spc03(fb, "NCRW", "NONC", "- Non commercial root weight (NCRW) is zero or NA but the number of non commercial roots (NONC) is greater than zero:")
  spc06(fb, "TRW", "NOCR", "NONC", "- Total root weight (TRW) is zero or NA but number of roots (NOCR + NONC) is greater than zero:")
  spc05(fb, "TRW", "NOCR", "NONC", "- Number of roots (NOCR + NONC) is zero or NA but total root weight (TRW) is greater than zero:")
  
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
  
  spc07(fb, temp, "RFCP", "- There are no roots but there is data for root primary flesh color (RFCP):")
  spc07(fb, temp, "RFCS", "- There are no roots but there is data for root secondary flesh color (RFCS):")
  spc07(fb, temp, "SCOL", "- There are no roots but there is data for storage root skin color (SCOL):")
  spc07(fb, temp, "FCOL", "- There are no roots but there is data for storage root flesh color (FCOL):")
  spc07(fb, temp, "RS", "- There are no roots but there is data for root size (RS):")
  spc07(fb, temp, "RF", "- There are no roots but there is data for root form (RF):")
  spc07(fb, temp, "DAMR", "- There are no roots but there is data for root defects (DAMR):")
  spc07(fb, temp, "RSPR", "- There are no roots but there is data for root sprouting (RSPR):")
  spc07(fb, temp, "WED1", "- There are no roots but there is data for weevil damage first evaluation (WED1):")
  spc07(fb, temp, "WED2", "- There are no roots but there is data for weevil damage second evaluation (WED2):")
  spc07(fb, temp, "DMF", "- There are no roots but there is data for fresh weight of roots for dry matter assessment (DMF):")
  spc07(fb, temp, "DMD", "- There are no roots but there is data for dry weight of roots for dry matter assessment (DMD):")
  spc01(fb, "DMD", "DMF", "- Dry weight of roots for dry matter assessment (DMD) is greater than fresh weight of roots for dry matter assessment (DMF):")
  spc07(fb, temp, "FRAW1", "- There are no roots but there is data for root fiber first determination (FRAW1):")
  spc07(fb, temp, "SURAW1", "- There are no roots but there is data for root sugar first determination (SURAW1):")
  spc07(fb, temp, "STRAW1", "- There are no roots but there is data for root starch first determination (STRAW1):")
  spc07(fb, temp, "COOF1", "- There are no roots but there is data for cooked fiber first evaluation (COOF1):")
  spc07(fb, temp, "COOSU1", "- There are no roots but there is data for cooked sugars first evaluation (COOSU1):")
  spc07(fb, temp, "COOST1", "- There are no roots but there is data for cooked starch first evaluation (COOST1):")
  spc07(fb, temp, "COOT1", "- There are no roots but there is data for cooked taste first evaluation (COOT1):")
  spc07(fb, temp, "COOAP1", "- There are no roots but there is data for cooked appearance first evaluation (COOAP1):")
  spc07(fb, temp, "FRAW2", "- There are no roots but there is data for root fiber second determination (FRAW2):")
  spc07(fb, temp, "SURAW2", "- There are no roots but there is data for root sugar second determination (SURAW2):")
  spc07(fb, temp, "STRAW2", "- There are no roots but there is data for root starch second determination (STRAW2):")
  spc07(fb, temp, "COOF2", "- There are no roots but there is data for cooked fiber second evaluation (COOF2):")
  spc07(fb, temp, "COOSU2", "- There are no roots but there is data for cooked sugars second evaluation (COOSU2):")
  spc07(fb, temp, "COOST2", "- There are no roots but there is data for cooked starch second evaluation (COOST2):")
  spc07(fb, temp, "COOT2", "- There are no roots but there is data for cooked taste second evaluation (COOT2):")
  spc07(fb, temp, "COOAP2", "- There are no roots but there is data for cooked appearance second evaluation (COOAP2):")
  spc07(fb, temp, "PROT", "- There are no roots but there is data for protein (PROT):")
  spc07(fb, temp, "FE", "- There are no roots but there is data for iron in dry weight (FE):")
  spc07(fb, temp, "ZN", "- There are no roots but there is data for zinc in dry weight (ZN):")
  spc07(fb, temp, "CA", "- There are no roots but there is data for calcium in dry weight (CA):")
  spc07(fb, temp, "MG", "- There are no roots but there is data for magnesium in dry weight (MG):")
  spc07(fb, temp, "BC", "- There are no roots but there is data for beta-carotene in dry weight (BC):")
  spc07(fb, temp, "TC", "- There are no roots but there is data for total carotenoids in dry weight (TC):")
  spc07(fb, temp, "STAR", "- There are no roots but there is data for starch (STAR):")
  spc07(fb, temp, "FRUC", "- There are no roots but there is data for fructose (FRUC):")
  spc07(fb, temp, "GLUC", "- There are no roots but there is data for glucose (GLUC):")
  spc07(fb, temp, "SUCR", "- There are no roots but there is data for sucrose (SUCR):")
  spc07(fb, temp, "MALT", "- There are no roots but there is data for maltose (MALT):")
  
  # Calculated variables
  
  spc08(fb, "TRW", "CRW", "NCRW", "- Total root weight (TRW) different from CRW + NCRW:")
  spc09(fb, plot.size, "CYTHA", "CRW", "- Commercial root yield in tons per hectare (CYTHA) is different from CRW * 10 / plot.size:")
  spc10(fb, plot.size, "RYTHA", "CRW", "NCRW", "- Total root yield in tons per hectare (RYTHA) is different from (CRW + NCRW) * 10 / plot.size:")
  spc11(fb, "ACRW", "CRW", "NOCR", "- Average commercial root weight (ACRW) is different from CRW / NOCR:")
  spc12(fb, "NRPP", "NOCR", "NONC", "NOPH", "- Number of roots per plant (NRPP) is different from (NOCR + NONC) / NOPH:")
  spc12(fb, "YPP", "CRW", "NCRW", "NOPH", "- Yield per plant (YPP) is different from (CRW + NCRW) / NOPH:")
  spc13(fb, "CI", "NOCR", "NONC", "- Commercial index (CI) is different from NOCR / (NOCR + NONC) * 100:")
  spc14(fb, "HI", "CRW", "NCRW", "VW", "- Harvest index (HI) is different from (CRW + NCRW) / (VW + CRW + NCRW) * 100:")
  spc15(fb, "SHI", "NOPH", "NOPS", "- Harvest sowing index (SHI) is different from NOPH / NOPS * 100:")
  spc16(fb, plot.size, "BIOM", "CRW", "NCRW", "VW", "- Biomass yield (BIOM) is different from (CRW + NCRW + VW) * 10 / plot.size:")
  spc09(fb, plot.size, "FYTHA", "VW", "- Foliage total yield in tons per hectare (FYTHA) is different from VW * 10 / plot.size:")
  spc15(fb, "DM", "DMD", "DMF", "- Storage root dry matter content (DM) is different from DMD / DMF * 100:")
  spc15(fb, "DMV", "DMVD", "DMVF", "- Vine dry matter content (DMV) is different from DMVD / DMVF * 100:")
  spc17(fb, plot.size, "DMFY", "VW", "DMVD", "DMVF", "- Dry matter foliage yield (DMFY) is different from VW * 10 / plot.size * DMVD / DMVF:")
  spc18(fb, plot.size, "DMRY", "CRW", "NCRW", "DMD", "DMF", "- Dry matter root yield (DMRY) is different from (CRW + NCRW) * 10 / plot.size * DMD / DMF:")
  spc19(fb, "RFR", "CRW", "NCRW", "DMD", "DMF", "VW", "DMVD", "DMVF", "- Root foliage ratio (RFR) is different from (CRW + NCRW) * (DMD / DMF) / (VW * DMVD / DMVF) * 100:")
  
  # Outliers detection and values out of range for field data

  spc20(fb, "NOPE", c(0:100, NA), "- Out of range values for number of plants established (NOPE):")
  spc20(fb, "VIR1", c(1:9, NA), "- Out of range values for virus symptoms first evaluation (VIR1):")
  spc20(fb, "VIR2", c(1:9, NA), "- Out of range values for virus symptoms second evaluation (VIR2):")
  spc20(fb, "VIR2", c(1:9, NA), "- Out of range values for virus symptoms third evaluation (VIR3):")
  spc20(fb, "ALT1", c(1:9, NA), "- Out of range values for alternaria symptoms first evaluation (ALT1):")
  spc20(fb, "ALT2", c(1:9, NA), "- Out of range values for alternaria symptoms second evaluation (ALT2):")
  spc20(fb, "VV1", c(1:9, NA), "- Out of range values for vine vigor first evaluation (VV1):")
  spc20(fb, "VV2", c(1:9, NA), "- Out of range values for vine vigor second evaluation (VV2):")
  spc21(fb, "VW", "- Out of range values for vine weight (VW):")
  spc23(fb, f, "VW", "- Extreme low values for vine weight (VW):")
  spc24(fb, f, "VW", "- Extreme high values for vine weight (VW):")
  spc21(fb, "NOPH", "- Out of range values for number of plants harvested (NOPH):")
  spc21(fb, "NOPR", "- Out of range values for number of plants with roots (NOPR):")
  spc21(fb, "NOCR", "- Out of range values for number of commercial roots (NOCR):")
  spc23(fb, f, "NOCR", "- Extreme low values for number of commercial roots (NOCR):")
  spc24(fb, f, "NOCR", "- Extreme high values for number of commercial roots (NOCR):")
  spc21(fb, "NONC", "- Out of range values for number of non commercial roots (NONC):")
  spc23(fb, f, "NONC", "- Extreme low values for number of non commercial roots (NONC):")
  spc24(fb, f, "NONC", "- Extreme high values for number of non commercial roots (NONC):")
  spc21(fb, "CRW", "- Out of range values for commercial root weight (CRW):")
  spc23(fb, f, "CRW", "- Extreme low values for commercial root weight (CRW):")
  spc24(fb, f, "CRW", "- Extreme high values for commercial root weight (CRW):")
  spc21(fb, "NCRW", "- Out of range values for non commercial root weight (NCRW):")
  spc23(fb, f, "NCRW", "- Extreme low values for non commercial root weight (NCRW):")
  spc24(fb, f, "NCRW", "- Extreme high values for non commercial root weight (NCRW):")
  spc20(fb, "SCOL", c(1:9, NA), "- Out of range values for storage root skin color (SCOL):")
  spc20(fb, "FCOL", c(1:9, NA), "- Out of range values for storage root flesh color (FCOL):")
  spc20(fb, "RS", c(1:9, NA), "- Out of range values for root size (RS):")
  spc20(fb, "RF", c(1:9, NA), "- Out of range values for root form (RF):")
  spc20(fb, "DAMR", c(1:9, NA), "- Out of range values for root defects (DAMR):")
  spc20(fb, "RSPR", c(1:9, NA), "- Out of range values for root sprouting (RSPR):")
  spc20(fb, "WED1", c(1:9, NA), "- Out of range values for weevil damage first evaluation (WED1):")
  spc20(fb, "WED2", c(1:9, NA), "- Out of range values for weevil damage second evaluation (WED2):")
  
  # Outliers detection and values out of range for DM data
  
  spc21(fb, "DMF", "- Out of range values for fresh weight of roots for dry matter assessment (DMF):")
  spc23(fb, f, "DMF", "- Extreme low values for fresh weight of roots for dry matter assessment (DMF):")
  spc24(fb, f, "DMF", "- Extreme high values for fresh weight of roots for dry matter assessment (DMF):")
  spc21(fb, "DMD", "- Out of range values for dry weight of roots for dry matter assessment (DMD):")
  spc23(fb, f, "DMD", "- Extreme low values for dry weight of roots for dry matter assessment (DMD):")
  spc24(fb, f, "DMD", "- Extreme high values for dry weight of roots for dry matter assessment (DMD):")
  spc21(fb, "DMVF", "- Out of range values for fresh weight vines for dry matter assessment (DMVF):")
  spc23(fb, f, "DMVF", "- Extreme low values for fresh weight of vines for dry matter assessment (DMVF):")
  spc24(fb, f, "DMVF", "- Extreme high values for fresh weight of vines for dry matter assessment (DMVF):")
  spc21(fb, "DMVD", "- Out of range values for dry weight of vines for dry matter assessment (DMVD):")
  spc23(fb, f, "DMVD", "- Extreme low values for dry weight of vines for dry matter assessment (DMVD):")
  spc24(fb, f, "DMVD", "- Extreme high values for dry weight of vines for dry matter assessment (DMVD):")
  spc21(fb, "DM", "- Out of range values for storage root dry matter content (DM):")
  spc23(fb, f, "DM", "- Extreme low values for storage root dry matter content (DM):")
  spc24(fb, f, "DM", "- Extreme high values for storage root dry matter content (DM):")
  spc21(fb, "DMV", "- Out of range values for vine dry matter content (DMV):")
  spc23(fb, f, "DMV", "- Extreme low values for vine dry matter content (DMV):")
  spc24(fb, f, "DMV", "- Extreme high values for vine dry matter content (DMV):")
  spc21(fb, "DMFY", "- Out of range values for dry matter foliage yield (DMFY):")
  spc23(fb, f, "DMFY", "- Extreme low values for dry matter foliage yield (DMFY):")
  spc24(fb, f, "DMFY", "- Extreme high values for dry matter foliage yield (DMFY):")
  spc21(fb, "DMRY", "- Out of range values for dry matter root yield (DMRY):")
  spc23(fb, f, "DMRY", "- Extreme low values for dry matter root yield (DMRY):")
  spc24(fb, f, "DMRY", "- Extreme high values for dry matter root yield (DMRY):")
  
  # Outliers detection and values out of range for cooked traits
  
  spc20(fb, "FRAW1", c(1:9, NA), "- Out of range values for root fiber first determination (FRAW1):")
  spc20(fb, "SURAW1", c(1:9, NA), "- Out of range values for root sugar first determination (SURAW1):")
  spc20(fb, "STRAW1", c(1:9, NA), "- Out of range values for root starch first determination (STRAW1):")
  spc20(fb, "COOF1", c(1:9, NA), "- Out of range values for cooked fiber first evaluation (COOF1):")
  spc20(fb, "COOSU1", c(1:9, NA), "- Out of range values for cooked sugars first evaluation (COOSU1):")
  spc20(fb, "COOST1", c(1:9, NA), "- Out of range values for cooked starch first evaluation (COOST1):")
  spc20(fb, "COOT1", c(1:9, NA), "- Out of range values for cooked taste first evaluation (COOT1):")
  spc20(fb, "COOAP1", c(1:9, NA), "- Out of range values for cooked appearance first evaluation (COOAP1):")
  spc20(fb, "FRAW2", c(1:9, NA), "- Out of range values for root fiber second determination (FRAW2):")
  spc20(fb, "SURAW2", c(1:9, NA), "- Out of range values for root sugar second determination (SURAW2):")
  spc20(fb, "STRAW2", c(1:9, NA), "- Out of range values for root starch second determination (STRAW2):")
  spc20(fb, "COOF2", c(1:9, NA), "- Out of range values for cooked fiber second evaluation (COOF2):")
  spc20(fb, "COOSU2", c(1:9, NA), "- Out of range values for cooked sugars second evaluation (COOSU2):")
  spc20(fb, "COOST2", c(1:9, NA), "- Out of range values for cooked starch second evaluation (COOST2):")
  spc20(fb, "COOT2", c(1:9, NA), "- Out of range values for cooked taste second evaluation (COOT2):")
  spc20(fb, "COOAP2", c(1:9, NA), "- Out of range values for cooked appearance second evaluation (COOAP2):")
  
  # Outliers detection and values out of range for lab data
  
  spc21(fb, "PROT", "- Out of range values for protein (PROT):")
  spc23(fb, f, "PROT", "- Extreme low values for protein (PROT):")
  spc24(fb, f, "PROT", "- Extreme high values for protein (PROT):")
  spc21(fb, "FE", "- Out of range values for iron in dry weight (FE):")
  spc23(fb, f, "FE", "- Extreme low values for iron in dry weight (FE):")
  spc24(fb, f, "FE", "- Extreme high values for iron in dry weight (FE):")
  spc21(fb, "ZN", "- Out of range values for zinc in dry weight (ZN):")
  spc23(fb, f, "ZN", "- Extreme low values for zinc in dry weight (ZN):")
  spc24(fb, f, "ZN", "- Extreme high values for zinc in dry weight (ZN):")
  spc21(fb, "CA", "- Out of range values for calcium in dry weight (CA):")
  spc23(fb, f, "CA", "- Extreme low values for calcium in dry weight (CA):")
  spc24(fb, f, "CA", "- Extreme high values for calcium in dry weight (CA):")
  spc21(fb, "MG", "- Out of range values for magnesium in dry weight (MG):")
  spc23(fb, f, "MG", "- Extreme low values for magnesium in dry weight (MG):")
  spc24(fb, f, "MG", "- Extreme high values for magnesium in dry weight (MG):")
  spc21(fb, "BC", "- Out of range values for beta-carotene in dry weight (BC):")
  spc23(fb, f, "BC", "- Extreme low values for beta-carotene in dry weight (BC):")
  spc24(fb, f, "BC", "- Extreme high values for beta-carotene in dry weight (BC):")
  bc.cc.values <- c(0.03, 0, 0.12, 0.02, 0.15, 1.38, 1.65, 1.5, 1.74, 1.76, 0.69, 1.17, 1.32,
                    1.04, 4.41, 4.92, 6.12, 5.46, 3.96, 5.49, 3.03, 3.76, 4.61, 7.23, 7.76,
                    10.5, 11.03, 12.39, 14.37, NA)
  spc20(fb, "BC.CC", bc.cc.values, "- Out of range values for beta-carotene with color chart (BC.CC):")
  spc21(fb, "TC", "- Out of range values for total carotenoids in dry weight (TC):")
  spc23(fb, f, "TC", "- Extreme low values for total carotenoids in dry weight (TC):")
  spc24(fb, f, "TC", "- Extreme high values for total carotenoids in dry weight (TC):")
  spc21(fb, "STAR", "- Out of range values for starch (STAR):")
  spc23(fb, f, "STAR", "- Extreme low values for starch (STAR):")
  spc24(fb, f, "STAR", "- Extreme high values for starch (STAR):")
  spc21(fb, "FRUC", "- Out of range values for fructose (FRUC):")
  spc23(fb, f, "FRUC", "- Extreme low values for fructose (FRUC):")
  spc24(fb, f, "FRUC", "- Extreme high values for fructose (FRUC):")
  spc21(fb, "GLUC", "- Out of range values for glucose (GLUC):")
  spc23(fb, f, "GLUC", "- Extreme low values for glucose (GLUC):")
  spc24(fb, f, "GLUC", "- Extreme high values for glucose (GLUC):")
  spc21(fb, "SUCR", "- Out of range values for sucrose (SUCR):")
  spc23(fb, f, "SUCR", "- Extreme low values for sucrose (SUCR):")
  spc24(fb, f, "SUCR", "- Extreme high values for sucrose (SUCR):")
  spc21(fb, "MALT", "- Out of range values for maltose (MALT):")
  spc23(fb, f, "MALT", "- Extreme low values for maltose (MALT):")
  spc24(fb, f, "MALT", "- Extreme high values for maltose (MALT):")
  
  # Outliers detection and values out of range for derived variables

  spc21(fb, "TRW", "- Out of range values for total root weight (TRW):")
  spc23(fb, f, "TRW", "- Extreme low values for total root weight (TRW):")
  spc24(fb, f, "TRW", "- Extreme high values for total root weight (TRW):")
  spc21(fb, "CYTHA", "- Out of range values for commercial root yield in tons per hectare (CYTHA):")
  spc23(fb, f, "CYTHA", "- Extreme low values for commercial root yield in tons per hectare (CYTHA):")
  spc24(fb, f, "CYTHA", "- Extreme high values for commercial root yield in tons per hectare (CYTHA):")
  spc21(fb, "RYTHA", "- Out of range values for total root yield in tons per hectare (RYTHA):")
  spc23(fb, f, "RYTHA", "- Extreme low values for total root yield in tons per hectare (RYTHA):")
  spc24(fb, f, "RYTHA", "- Extreme high values for total root yield in tons per hectare (RYTHA):")
  spc21(fb, "ACRW", "- Out of range values for average commercial root weight (ACRW):")
  spc23(fb, f, "ACRW", "- Extreme low values for average commercial root weight (ACRW):")
  spc24(fb, f, "ACRW", "- Extreme high values for average commercial root weight (ACRW):")
  spc21(fb, "NRPP", "- Out of range values for number of roots per plant (NRPP):")
  spc23(fb, f, "NRPP", "- Extreme low values for number of roots per plant (NRPP):")
  spc24(fb, f, "NRPP", "- Extreme high values for number of roots per plant (NRPP):")
  spc21(fb, "YPP", "- Out of range values for yield per plant (YPP):")
  spc23(fb, f, "YPP", "- Extreme low values for yield per plant (YPP):")
  spc24(fb, f, "YPP", "- Extreme high values for yield per plant (YPP):")
  spc22(fb, "CI", "- Out of range values for commercial index (CI):")
  spc23(fb, f, "CI", "- Extreme low values for commercial index (CI):")
  spc24(fb, f, "CI", "- Extreme high values for commercial index (CI):")
  spc22(fb, "HI", "- Out of range values for harvest index (HI):")
  spc23(fb, f, "HI", "- Extreme low values for harvest index (HI):")
  spc24(fb, f, "HI", "- Extreme high values for harvest index (HI):")
  spc22(fb, "SHI", "- Out of range values for harvest sowing index (SHI):")
  spc23(fb, f, "SHI", "- Extreme low values for harvest sowing index (SHI):")
  spc24(fb, f, "SHI", "- Extreme high values for harvest sowing index (SHI):")
  spc21(fb, "BIOM", "- Out of range values for biomass yield (BIOM):")
  spc23(fb, f, "BIOM", "- Extreme low values for biomass yield (BIOM):")
  spc24(fb, f, "BIOM", "- Extreme high values for biomass yield (BIOM):")
  spc21(fb, "FYTHA", "- Out of range values for foliage total yield in tons per hectare (FYTHA):")
  spc23(fb, f, "FYTHA", "- Extreme low values for foliage total yield in tons per hectare (FYTHA):")
  spc24(fb, f, "FYTHA", "- Extreme high values for foliage total yield in tons per hectare (FYTHA):")
  spc22(fb, "RFR", "- Out of range values for root foliage ratio (RFR):")
  spc23(fb, f, "RFR", "- Extreme low values for root foliage ratio (RFR):")
  spc24(fb, f, "RFR", "- Extreme high values for root foliage ratio (RFR):")
  
  if (file == TRUE) sink()
}

output <- function(fb, cond, tx) {
  if (sum(cond, na.rm = TRUE) > 0) {
    cat("\n", tx, "\n", sep = "")
    print(subset(fb, cond))
  }
}

spc01 <- function(fb, t1, t2, tx) {
  if (exists(t1, fb) & exists(t2, fb)) {
    cond <- fb[, t1] > fb[, t2]
    output(fb, cond, tx)
  }
}
  
spc02 <- function(fb, t1, t2, tx) {
  if (exists(t1, fb) & exists(t2, fb)) {
    cond <- (fb[, t1] == 0 | is.na(fb[, t1])) & !is.na(fb[, t2])
    output(fb, cond, tx)
  }
}

spc03 <- function(fb, t1, t2, tx) {
  if (exists(t1, fb) & exists(t2, fb)) {
    cond <- (fb[, t1] == 0 | is.na(fb[, t1])) & fb[, t2] > 0
    output(fb, cond, tx)
  }
}

spc04 <- function(fb, t1, t2, t3, tx) {
  if (exists(t1, fb) & exists(t2, fb) & exists(t3, fb)) {
    cond <- (fb[, t1] == 0 | is.na(fb[, t1])) & (fb[, t2] > 0 | fb[, t3] > 0)
    output(fb, cond, tx)
  }
}

spc05 <- function(fb, t1, t2, t3, tx) {
  if (exists(t1, fb) & exists(t2, fb) & exists(t3, fb)) {
    cond <- fb[, t1] > 0 & (suma(fb[, t2], fb[, t3]) == 0 | is.na(suma(fb[, t2], fb[, t3])))
    output(fb, cond, tx)
  }
}

spc06 <- function(fb, t1, t2, t3, tx) {
  if (exists(t1, fb) & exists(t2, fb) & exists(t3, fb)) {
    cond <- (fb[, t1] == 0 | is.na(fb[, t1])) & suma(fb[, t2], fb[, t3]) > 0
    output(fb, cond, tx)
  }
}

spc07 <- function(fb, temp, t1, tx) {
  if (exists(t1, fb)) {
    cond <- temp == FALSE
    output(fb, cond, tx)
  }
}

spc08 <- function(fb, t1, t2, t3, tx) {
  if (exists(t1, fb) & exists(t2, fb) & exists(t3, fb)) {
    cond <- abs(fb[, t1] - suma(fb[, t2], fb[, t3])) > 1e-10
    output(fb, cond, tx)
  }
}

spc09 <- function(fb, plot.size, t1, t2, tx) {
  if (exists(t1, fb) & exists(t2, fb)) {
    cond <- abs(fb[, t1] - fb[, t2] * 10 / plot.size) > 1e-10
    output(fb, cond, tx)
  }
}

spc10 <- function(fb, plot.size, t1, t2, t3, tx) {
  if (exists(t1, fb) & exists(t2, fb) & exists(t3, fb)) {
    cond <- abs(fb[, t1] - suma(fb[, t2], fb[, t3]) * 10 / plot.size) > 1e-10
    output(fb, cond, tx)
  }
}

spc11 <- function(fb, t1, t2, t3, tx) {
  if (exists(t1, fb) & exists(t2, fb) & exists(t3, fb)) {
    cond <- abs(fb[, t1] - fb[, t2] / fb[, t3]) > 1e-10
    output(fb, cond, tx)
  }
}

spc12 <- function(fb, t1, t2, t3, t4, tx) {
  if (exists(t1, fb) & exists(t2, fb) & exists(t3, fb) & exists(t4, fb)) {
    cond <- abs(fb[, t1] - suma(fb[, t2], fb[, t3]) / fb[, t4]) > 1e-10
    output(fb, cond, tx)
  }
}

spc13 <- function(fb, t1, t2, t3, tx) {
  if (exists(t1, fb) & exists(t2, fb) & exists(t3, fb)) {
    cond <- abs(fb[, t1] - fb[, t2] / suma(fb[, t2], fb[, t3]) * 100) > 1e-10
    output(fb, cond, tx)
  }
}

spc14 <- function(fb, t1, t2, t3, t4, tx) {
  if (exists(t1, fb) & exists(t2, fb) & exists(t3, fb) & exists(t4, fb)) {
    cond <- abs(fb[, t1] -  suma(fb[, t2], fb[, t3]) / suma(suma(fb[, t2], fb[, t3]), fb[, t4]) * 100) > 1e-10
    output(fb, cond, tx)
  }
}

spc15 <- function(fb, t1, t2, t3, tx) {
  if (exists(t1, fb) & exists(t2, fb) & exists(t3, fb)) {
    cond <- abs(fb[, t1] -  fb[, t2] / fb[, t3] * 100) > 1e-10
    output(fb, cond, tx)
  }
}

spc16 <- function(fb, plot.size, t1, t2, t3, t4, tx) {
  if (exists(t1, fb) & exists(t2, fb) & exists(t3, fb) & exists(t4, fb)) {
    cond <- abs(fb[, t1] -  suma(suma(fb[, t2], fb[, t3]), fb[, t4]) / plot.size * 10) > 1e-10
    output(fb, cond, tx)
  }
}

spc17 <- function(fb, plot.size, t1, t2, t3, t4, tx) {
  if (exists(t1, fb) & exists(t2, fb) & exists(t3, fb) & exists(t4, fb)) {
    cond <- abs(fb[, t1] - fb[, t2] * 10 / plot.size * fb[, t3] / fb[, t4]) > 1e-10
    output(fb, cond, tx)
  }
}

spc18 <- function(fb, plot.size, t1, t2, t3, t4, t5, tx) {
  if (exists(t1, fb) & exists(t2, fb) & exists(t3, fb) & exists(t4, fb) & exists(t5, fb)) {
    cond <- abs(fb[, t1] - suma(fb[, t2], fb[, t3]) * 10 / plot.size * fb[, t4] / fb[, t5]) > 1e-10
    output(fb, cond, tx)
  }
}

spc19 <- function(fb, t1, t2, t3, t4, t5, t6, t7, t8, tx) {
  if (exists(t1, fb) & exists(t2, fb) & exists(t3, fb) & exists(t4, fb) & exists(t5, fb)
      & exists(t6, fb) & exists(t7, fb) & exists(t8, fb)) {
    cond <- abs(fb[, t1] - suma(fb[, t2], fb[, t3]) * fb[, t4] / fb[, t5]
                / (fb[, t6] * fb[, t7] / fb[, t8])) * 100 > 1e-10
    output(fb, cond, tx)
  }
}

spc20 <- function(fb, t1, vv, tx) {
  if (exists(t1, fb)) {
    cond <- !(fb[, t1] %in% vv)
    output(fb, cond, tx)
  }
}

spc21 <- function(fb, t1, tx) {
  if (exists(t1, fb)) {
    cond <- fb[, t1] < 0
    output(fb, cond, tx)
  }
}

spc22 <- function(fb, t1, tx) {
  if (exists(t1, fb)) {
    cond <- fb[, t1] < 0 | fb[, t1] > 100
    output(fb, cond, tx)
  }
}

spc23 <- function(fb, f, t1, tx) {
  if (exists(t1, fb)) {
    cond <- fb[, t1] < quantile(fb[, t1], 0.25, na.rm = TRUE) - f * IQR(fb[, t1], na.rm = TRUE)
    output(fb, cond, tx)
  }
}

spc24 <- function(fb, f, t1, tx) {
  if (exists(t1, fb)) {
    cond <- fb[, t1] > quantile(fb[, t1], 0.75, na.rm = TRUE) + f * IQR(fb[, t1], na.rm = TRUE)
    output(fb, cond, tx)
  }
}
