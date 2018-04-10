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
#' @details The data frame must use the labels (lower or upper case) listed in function
#' \code{check.names}. See \code{?check.names} for details.
#' 
#' Extreme low and high values are detected using the interquartile range.
#' The rule is to detect any value out of the interval 
#' \eqn{[Q_1 - f \times IQR; Q_3 + f \times IQR]}. By default \code{f = 3}.
#' 
#' Outliers are detected based on standardized residuals for some statistical models.
#' Options are \code{'rcbd'} and \code{'met'} for a randomized complete block design
#' and a multi environment trial with RCBD in each environment. By default the threshold
#' value is \code{out.max = 4}.
#' @return It returns all rows with some kind of inconsistency or outliers.
#' @author Raul Eyzaguirre.
#' @examples
#' check.data(pjpz09)
#' @importFrom stats IQR quantile rstandard
#' @export

check.data <- function(fb, f = 3, out.mod = c("none", "rcbd", "met"),
                       out.max = 4, aqt = NULL) {
  
  out.mod = match.arg(out.mod)
  
  # Check names
  
  fb <- check.names(fb, aqt)
  if (!is.null(aqt))
    aqt <- tolower(aqt)

  # Inconsistencies for nops > nope > noph > nopr.

  sp1(fb, 1, "nope", "nops", "- Number of plants established (nope) is greater than number of plants sowed (nops):")
  sp1(fb, 1, "noph", "nops", "- Number of plants harvested (noph)  is greater than number of plants sowed (nops):")
  sp1(fb, 1, "nopr", "nops", "- Number of plants with roots (nopr) is greater than number of plants sowed (nops):")
  sp1(fb, 1, "noph", "nope", "- Number of plants harvested (noph) is greater than number of plants established (nope):")
  sp1(fb, 1, "nopr", "nope", "- Number of plants with roots (nopr) is greater than number of plants established (nope):")
  sp1(fb, 1, "nopr", "noph", "- Number of plants with roots (nopr) is greater than number of plants harvested (noph):")

  # Inconsistencies for nope and dependencies.

  sp1(fb, 2, "nope", "vir", "- Number of plants established (nope) is zero or NA but there is data for virus symptoms (vir):")
  sp1(fb, 2, "nope", "vir1", "- Number of plants established (nope) is zero or NA but there is data for virus symptoms first evaluation (vir1):")
  sp1(fb, 2, "nope", "vir2", "- Number of plants established (nope) is zero or NA but there is data for virus symptoms second evaluation (vir2):")
  sp1(fb, 2, "nope", "alt", "- Number of plants established (nope) is zero or NA but there is data for alternaria symptoms (alt):")
  sp1(fb, 2, "nope", "alt1", "- Number of plants established (nope) is zero or NA but there is data for alternaria symptoms first evaluation (alt1):")
  sp1(fb, 2, "nope", "alt2", "- Number of plants established (nope) is zero or NA but there is data for alternaria symptoms second evaluation (alt2):")
  sp1(fb, 2, "nope", "vv", "- Number of plants established (nope) is zero or NA but there is data for vine vigor (vv):")

  # noph and vw
  
  sp1(fb, 3, "noph", "vw", "- Number of plants harvested (noph) is zero or NA but vine weight (vw) is greater than zero:") 
  sp1(fb, 3, "vw", "noph", "- Vine weight (vw) is zero or NA but the number of plants harvested (noph) is greater than zero:") 
  sp1(fb, 3, "noph", "fytha", "- Number of plants harvested (noph) is zero or NA but foliage yield in tons per hectare (fytha) is greater than zero:") 
  sp1(fb, 3, "fytha", "noph", "- Foliage yield in tons per hectare (fytha) is zero or NA but the number of plants harvested (noph) is greater than zero:") 
  sp1(fb, 3, "noph", "fytha.aj", "- Number of plants harvested (noph) is zero or NA but foliage yield in tons per hectare (fytha.aj) is greater than zero:") 
  sp1(fb, 3, "fytha.aj", "noph", "- Foliage yield in tons per hectare (fytha.aj) is zero or NA but the number of plants harvested (noph) is greater than zero:") 
  
  # vw and dependencies
  
  sp1(fb, 3, "vw", "dmvf", "- Vine weight (vw) is zero or NA but there is fresh weight vines for dry matter assessment (dmvf):") 
  sp1(fb, 3, "vw", "dmvd", "- Vine weight (vw) is zero or NA but there is dry weight vines for dry matter assessment (dmvd):") 
  sp1(fb, 1, "dmvd", "dmvf", "- Dry weight vines for dry matter assessment (dmvd) is greater than fresh weight vines for dry matter assessment (dmvf):")

  # nopr and roots
  
  sp1(fb, 3, "nopr", "tnr", "- Number of plants with roots (nopr) is zero or NA but total number of roots (tnr) is greater than zero:")
  sp1(fb, 3, "tnr", "nopr", "- Number of roots (tnr) is zero or NA but number of plants with roots (nopr) is greater than zero:")
  sp1(fb, 3, "nopr", "trw", "- Number of plants with roots (nopr) is zero or NA but total root weight (trw) is greater than zero::")
  sp1(fb, 3, "trw", "nopr", "- Total root weight (trw) is zero or NA but number of plants with roots (nopr) is greater than zero:")
  sp1(fb, 3, "nopr", "rytha", "- Number of plants with roots (nopr) is zero or NA but root yield in tons per hectare (rytha) is greater than zero::")
  sp1(fb, 3, "rytha", "nopr", "- Root yield in tons per hectare (rytha) is zero or NA but number of plants with roots (nopr) is greater than zero:")
  sp1(fb, 3, "nopr", "rytha.aj", "- Number of plants with roots (nopr) is zero or NA but root yield in tons per hectare (rytha.aj) is greater than zero::")
  sp1(fb, 3, "rytha.aj", "nopr", "- Root yield in tons per hectare (rytha.aj) is zero or NA but number of plants with roots (nopr) is greater than zero:")
  sp1(fb, 3, "nopr", "cytha", "- Number of plants with roots (nopr) is zero or NA but commercial root yield in tons per hectare (cytha) is greater than zero::")
  sp1(fb, 3, "nopr", "cytha.aj", "- Number of plants with roots (nopr) is zero or NA but commercial root yield in tons per hectare (cytha.aj) is greater than zero::")
  sp1(fb, 2, "nopr", "wed", "- Number of plants with roots (nopr) is zero or NA but there is data for weevil damage (wed):")
  
  # Number of roots and root weight
  
  sp1(fb, 3, "nocr", "crw", "- Number of commercial roots (nocr) is zero or NA but the commercial root weight (crw) is greater than zero:") 
  sp1(fb, 3, "crw", "nocr", "- Commercial root weight (crw) is zero or NA but the number of commercial roots (nocr) is greater than zero:") 
  sp1(fb, 3, "nocr", "cytha", "- Number of commercial roots (nocr) is zero or NA but the commercial root yield in tons per hectare (cytha) is greater than zero:") 
  sp1(fb, 3, "cytha", "nocr", "- Commercial root yield in tons per hectare (cytha) is zero or NA but the number of commercial roots (nocr) is greater than zero:") 
  sp1(fb, 3, "nocr", "cytha.aj", "- Number of commercial roots (nocr) is zero or NA but the commercial root yield in tons per hectare (cytha.aj) is greater than zero:") 
  sp1(fb, 3, "cytha.aj", "nocr", "- Commercial root yield in tons per hectare (cytha.aj) is zero or NA but the number of commercial roots (nocr) is greater than zero:") 
  sp1(fb, 3, "nonc", "ncrw", "- Number of non commercial roots (nonc) is zero or NA but the non commercial root weight (ncrw) is greater than zero:")
  sp1(fb, 3, "ncrw", "nonc", "- Non commercial root weight (ncrw) is zero or NA but the number of non commercial roots (nonc) is greater than zero:")
  sp1(fb, 3, "trw", "tnr", "- Total root weight (trw) is zero or NA but total number of roots (tnr) is greater than zero:")
  sp1(fb, 3, "tnr", "trw", "- Total number of roots (tnr) is zero or NA but total root weight (trw) is greater than zero:")
  
  # Roots and dependencies

  temp <- array(FALSE, dim(fb)[1])
  do <- FALSE
  
  if (exists("nopr", fb)) {
    temp <- temp | (!is.na(fb$nopr) & fb$nopr > 0)
    do <- TRUE
  }
  if (exists("nocr", fb)) {
    temp <- temp | (!is.na(fb$nocr) & fb$nocr > 0)
    do <- TRUE
  }
  if (exists("nonc", fb)) {
    temp <- temp | (!is.na(fb$nonc) & fb$nonc > 0)
    do <- TRUE
  }
  if (exists("crw", fb)) {
    temp <- temp | (!is.na(fb$crw) & fb$crw > 0)
    do <- TRUE
  }
  if (exists("ncrw", fb)) {
    temp <- temp | (!is.na(fb$ncrw) & fb$ncrw > 0)
    do <- TRUE
  }
  if (exists("trw", fb)) {
    temp <- temp | (!is.na(fb$trw) & fb$trw > 0)
    do <- TRUE
  }
  if (exists("rytha", fb)) {
    temp <- temp | (!is.na(fb$rytha) & fb$rytha > 0)
    do <- TRUE
  }
  if (exists("rytha.aj", fb)) {
    temp <- temp | (!is.na(fb$rytha.aj) & fb$rytha.aj > 0)
    do <- TRUE
  }

  sp2(fb, temp, do, "fcol.cc", "- There are no roots but there is data for root flesh color using RHS color charts (fcol.cc):")
  sp2(fb, temp, do, "scol", "- There are no roots but there is data for storage root skin color (scol):")
  sp2(fb, temp, do, "fcol", "- There are no roots but there is data for storage root flesh color (fcol):")
  sp2(fb, temp, do, "rs", "- There are no roots but there is data for root size (rs):")
  sp2(fb, temp, do, "rf", "- There are no roots but there is data for root form (rf):")
  sp2(fb, temp, do, "damr", "- There are no roots but there is data for root defects (damr):")
  sp2(fb, temp, do, "rspr", "- There are no roots but there is data for root sprouting (rspr):")
  sp2(fb, temp, do, "wed", "- There are no roots but there is data for weevil damage (wed):")
  sp2(fb, temp, do, "dmf", "- There are no roots but there is data for fresh weight of roots for dry matter assessment (dmf):")
  sp2(fb, temp, do, "dmd", "- There are no roots but there is data for dry weight of roots for dry matter assessment (dmd):")
  sp1(fb, 1, "dmd", "dmf", "- Dry weight of roots for dry matter assessment (dmd) is greater than fresh weight of roots for dry matter assessment (dmf):")
  sp2(fb, temp, do, "fraw", "- There are no roots but there is data for root fiber (fraw):")
  sp2(fb, temp, do, "suraw", "- There are no roots but there is data for root sugar (suraw):")
  sp2(fb, temp, do, "straw", "- There are no roots but there is data for root starch (straw):")
  sp2(fb, temp, do, "coof", "- There are no roots but there is data for cooked fiber (coof):")
  sp2(fb, temp, do, "coosu", "- There are no roots but there is data for cooked sugars (coosu):")
  sp2(fb, temp, do, "coost", "- There are no roots but there is data for cooked starch (coost):")
  sp2(fb, temp, do, "coot", "- There are no roots but there is data for cooked taste (coot):")
  sp2(fb, temp, do, "cooap", "- There are no roots but there is data for cooked appearance (cooap):")
  sp2(fb, temp, do, "fraw1", "- There are no roots but there is data for root fiber first determination (fraw1):")
  sp2(fb, temp, do, "suraw1", "- There are no roots but there is data for root sugar first determination (suraw1):")
  sp2(fb, temp, do, "straw1", "- There are no roots but there is data for root starch first determination (straw1):")
  sp2(fb, temp, do, "coof1", "- There are no roots but there is data for cooked fiber first evaluation (coof1):")
  sp2(fb, temp, do, "coosu1", "- There are no roots but there is data for cooked sugars first evaluation (coosu1):")
  sp2(fb, temp, do, "coost1", "- There are no roots but there is data for cooked starch first evaluation (coost1):")
  sp2(fb, temp, do, "coot1", "- There are no roots but there is data for cooked taste first evaluation (coot1):")
  sp2(fb, temp, do, "cooap1", "- There are no roots but there is data for cooked appearance first evaluation (cooap1):")
  sp2(fb, temp, do, "fraw2", "- There are no roots but there is data for root fiber second determination (fraw2):")
  sp2(fb, temp, do, "suraw2", "- There are no roots but there is data for root sugar second determination (suraw2):")
  sp2(fb, temp, do, "straw2", "- There are no roots but there is data for root starch second determination (straw2):")
  sp2(fb, temp, do, "coof2", "- There are no roots but there is data for cooked fiber second evaluation (coof2):")
  sp2(fb, temp, do, "coosu2", "- There are no roots but there is data for cooked sugars second evaluation (coosu2):")
  sp2(fb, temp, do, "coost2", "- There are no roots but there is data for cooked starch second evaluation (coost2):")
  sp2(fb, temp, do, "coot2", "- There are no roots but there is data for cooked taste second evaluation (coot2):")
  sp2(fb, temp, do, "cooap2", "- There are no roots but there is data for cooked appearance second evaluation (cooap2):")
  sp2(fb, temp, do, "prot", "- There are no roots but there is data for protein (prot):")
  sp2(fb, temp, do, "fe", "- There are no roots but there is data for iron (fe):")
  sp2(fb, temp, do, "zn", "- There are no roots but there is data for zinc (zn):")
  sp2(fb, temp, do, "ca", "- There are no roots but there is data for calcium (ca):")
  sp2(fb, temp, do, "mg", "- There are no roots but there is data for magnesium (mg):")
  sp2(fb, temp, do, "bc", "- There are no roots but there is data for beta-carotene (bc):")
  sp2(fb, temp, do, "bc.cc", "- There are no roots but there is data for beta-carotene (bc.cc):")
  sp2(fb, temp, do, "tc", "- There are no roots but there is data for total carotenoids (tc):")
  sp2(fb, temp, do, "star", "- There are no roots but there is data for starch (star):")
  sp2(fb, temp, do, "fruc", "- There are no roots but there is data for fructose (fruc):")
  sp2(fb, temp, do, "gluc", "- There are no roots but there is data for glucose (gluc):")
  sp2(fb, temp, do, "sucr", "- There are no roots but there is data for sucrose (sucr):")
  sp2(fb, temp, do, "malt", "- There are no roots but there is data for maltose (malt):")
  
  # Extreme values detection and values out of range for field data

  sp3(fb, c(0:100, NA), "nope", "- Out of range values for number of plants established (nope):")
  sp3(fb, c(1:9, NA), "vir", "- Out of range values for virus symptoms (vir):")
  sp3(fb, c(1:9, NA), "vir1", "- Out of range values for virus symptoms first evaluation (vir1):")
  sp3(fb, c(1:9, NA), "vir2", "- Out of range values for virus symptoms second evaluation (vir2):")
  sp3(fb, c(1:9, NA), "alt", "- Out of range values for alternaria symptoms (alt):")
  sp3(fb, c(1:9, NA), "alt1", "- Out of range values for alternaria symptoms first evaluation (alt1):")
  sp3(fb, c(1:9, NA), "alt2", "- Out of range values for alternaria symptoms second evaluation (alt2):")
  sp3(fb, c(1:9, NA), "vv", "- Out of range values for vine vigor (vv):")
  sp4(fb, "lower", "vw", "- Out of range values for vine weight (vw):")
  sp5(fb, f, "low", "vw", "- Extreme low values for vine weight (vw):")
  sp5(fb, f, "high", "vw", "- Extreme high values for vine weight (vw):")
  sp4(fb, "lower", "noph", "- Out of range values for number of plants harvested (noph):")
  sp4(fb, "lower", "nopr", "- Out of range values for number of plants with roots (nopr):")
  sp4(fb, "lower", "nocr", "- Out of range values for number of commercial roots (nocr):")
  sp5(fb, f, "low", "nocr", "- Extreme low values for number of commercial roots (nocr):")
  sp5(fb, f, "high", "nocr", "- Extreme high values for number of commercial roots (nocr):")
  sp4(fb, "lower", "nonc", "- Out of range values for number of non commercial roots (nonc):")
  sp5(fb, f, "low", "nonc", "- Extreme low values for number of non commercial roots (nonc):")
  sp5(fb, f, "high", "nonc", "- Extreme high values for number of non commercial roots (nonc):")
  sp4(fb, "lower", "crw", "- Out of range values for commercial root weight (crw):")
  sp5(fb, f, "low", "crw", "- Extreme low values for commercial root weight (crw):")
  sp5(fb, f, "high", "crw", "- Extreme high values for commercial root weight (crw):")
  sp4(fb, "lower", "ncrw", "- Out of range values for non commercial root weight (ncrw):")
  sp5(fb, f, "low", "ncrw", "- Extreme low values for non commercial root weight (ncrw):")
  sp5(fb, f, "high", "ncrw", "- Extreme high values for non commercial root weight (ncrw):")
  sp3(fb, c(1:30, NA), "fcol.cc", "- Out of range values for root flesh color using RHS color charts (fcol.cc):")
  sp3(fb, c(1:9, NA), "scol", "- Out of range values for storage root skin color (scol):")
  sp3(fb, c(1:9, NA), "fcol", "- Out of range values for storage root flesh color (fcol):")
  sp3(fb, c(1:9, NA), "rs", "- Out of range values for root size (rs):")
  sp3(fb, c(1:9, NA), "rf", "- Out of range values for root form (rf):")
  sp3(fb, c(1:9, NA), "damr", "- Out of range values for root defects (damr):")
  sp3(fb, c(1:9, NA), "rspr", "- Out of range values for root sprouting (rspr):")
  sp3(fb, c(1:9, NA), "wed", "- Out of range values for weevil damage (wed):")

  # Extreme values detection and values out of range for dm data
  
  sp4(fb, "lower", "dmf", "- Out of range values for fresh weight of roots for dry matter assessment (dmf):")
  sp5(fb, f, "low", "dmf", "- Extreme low values for fresh weight of roots for dry matter assessment (dmf):")
  sp5(fb, f, "high", "dmf", "- Extreme high values for fresh weight of roots for dry matter assessment (dmf):")
  sp4(fb, "lower", "dmd", "- Out of range values for dry weight of roots for dry matter assessment (dmd):")
  sp5(fb, f, "low", "dmd", "- Extreme low values for dry weight of roots for dry matter assessment (dmd):")
  sp5(fb, f, "high", "dmd", "- Extreme high values for dry weight of roots for dry matter assessment (dmd):")
  sp4(fb, "lower", "dmvf", "- Out of range values for fresh weight vines for dry matter assessment (dmvf):")
  sp5(fb, f, "low", "dmvf", "- Extreme low values for fresh weight of vines for dry matter assessment (dmvf):")
  sp5(fb, f, "high", "dmvf", "- Extreme high values for fresh weight of vines for dry matter assessment (dmvf):")
  sp4(fb, "lower", "dmvd", "- Out of range values for dry weight of vines for dry matter assessment (dmvd):")
  sp5(fb, f, "low", "dmvd", "- Extreme low values for dry weight of vines for dry matter assessment (dmvd):")
  sp5(fb, f, "high", "dmvd", "- Extreme high values for dry weight of vines for dry matter assessment (dmvd):")
  sp4(fb, "lower", "dm", "- Out of range values for storage root dry matter content (dm):")
  sp5(fb, f, "low", "dm", "- Extreme low values for storage root dry matter content (dm):")
  sp5(fb, f, "high", "dm", "- Extreme high values for storage root dry matter content (dm):")
  sp4(fb, "lower", "dmv", "- Out of range values for vine dry matter content (dmv):")
  sp5(fb, f, "low", "dmv", "- Extreme low values for vine dry matter content (dmv):")
  sp5(fb, f, "high", "dmv", "- Extreme high values for vine dry matter content (dmv):")
  sp4(fb, "lower", "dmry", "- Out of range values for dry matter root yield (dmry):")
  sp5(fb, f, "low", "dmry", "- Extreme low values for dry matter root yield (dmry):")
  sp5(fb, f, "high", "dmry", "- Extreme high values for dry matter root yield (dmry):")
  sp4(fb, "lower", "dmvy", "- Out of range values for dry matter vine yield (dmvy):")
  sp5(fb, f, "low", "dmvy", "- Extreme low values for dry matter vine yield (dmvy):")
  sp5(fb, f, "high", "dmvy", "- Extreme high values for dry matter vine yield (dmvy):")
  sp4(fb, "lower", "dmbiom", "- Out of range values for dry matter biomass (dmbiom):")
  sp5(fb, f, "low", "dmbiom", "- Extreme low values for dry matter biomass (dmbiom):")
  sp5(fb, f, "high", "dmbiom", "- Extreme high values for dry matter biomass (dmbiom):")
  
  # Extreme values detection and values out of range for cooked traits
  
  sp3(fb, c(1:9, NA), "fraw", "- Out of range values for root fiber (fraw):")
  sp3(fb, c(1:9, NA), "suraw", "- Out of range values for root sugar (suraw):")
  sp3(fb, c(1:9, NA), "straw", "- Out of range values for root starch (straw):")
  sp3(fb, c(1:9, NA), "coof", "- Out of range values for cooked fiber (coof):")
  sp3(fb, c(1:9, NA), "coosu", "- Out of range values for cooked sugars (coosu):")
  sp3(fb, c(1:9, NA), "coost", "- Out of range values for cooked starch (coost):")
  sp3(fb, c(1:9, NA), "coot", "- Out of range values for cooked taste (coot):")
  sp3(fb, c(1:9, NA), "cooap", "- Out of range values for cooked appearance (cooap):")
  sp3(fb, c(1:9, NA), "fraw1", "- Out of range values for root fiber first determination (fraw1):")
  sp3(fb, c(1:9, NA), "suraw1", "- Out of range values for root sugar first determination (suraw1):")
  sp3(fb, c(1:9, NA), "straw1", "- Out of range values for root starch first determination (straw1):")
  sp3(fb, c(1:9, NA), "coof1", "- Out of range values for cooked fiber first evaluation (coof1):")
  sp3(fb, c(1:9, NA), "coosu1", "- Out of range values for cooked sugars first evaluation (coosu1):")
  sp3(fb, c(1:9, NA), "coost1", "- Out of range values for cooked starch first evaluation (coost1):")
  sp3(fb, c(1:9, NA), "coot1", "- Out of range values for cooked taste first evaluation (coot1):")
  sp3(fb, c(1:9, NA), "cooap1", "- Out of range values for cooked appearance first evaluation (cooap1):")
  sp3(fb, c(1:9, NA), "fraw2", "- Out of range values for root fiber second determination (fraw2):")
  sp3(fb, c(1:9, NA), "suraw2", "- Out of range values for root sugar second determination (suraw2):")
  sp3(fb, c(1:9, NA), "straw2", "- Out of range values for root starch second determination (straw2):")
  sp3(fb, c(1:9, NA), "coof2", "- Out of range values for cooked fiber second evaluation (coof2):")
  sp3(fb, c(1:9, NA), "coosu2", "- Out of range values for cooked sugars second evaluation (coosu2):")
  sp3(fb, c(1:9, NA), "coost2", "- Out of range values for cooked starch second evaluation (coost2):")
  sp3(fb, c(1:9, NA), "coot2", "- Out of range values for cooked taste second evaluation (coot2):")
  sp3(fb, c(1:9, NA), "cooap2", "- Out of range values for cooked appearance second evaluation (cooap2):")
  
  # Extreme values detection and values out of range for lab data
  
  sp4(fb, "lower", "prot", "- Out of range values for protein (prot):")
  sp5(fb, f, "low", "prot", "- Extreme low values for protein (prot):")
  sp5(fb, f, "high", "prot", "- Extreme high values for protein (prot):")
  sp4(fb, "lower", "fe", "- Out of range values for iron (fe):")
  sp5(fb, f, "low", "fe", "- Extreme low values for iron (fe):")
  sp5(fb, f, "high", "fe", "- Extreme high values for iron (fe):")
  sp4(fb, "lower", "zn", "- Out of range values for zinc (zn):")
  sp5(fb, f, "low", "zn", "- Extreme low values for zinc (zn):")
  sp5(fb, f, "high", "zn", "- Extreme high values for zinc (zn):")
  sp4(fb, "lower", "ca", "- Out of range values for calcium (ca):")
  sp5(fb, f, "low", "ca", "- Extreme low values for calcium (ca):")
  sp5(fb, f, "high", "ca", "- Extreme high values for calcium (ca):")
  sp4(fb, "lower", "mg", "- Out of range values for magnesium (mg):")
  sp5(fb, f, "low", "mg", "- Extreme low values for magnesium (mg):")
  sp5(fb, f, "high", "mg", "- Extreme high values for magnesium (mg):")
  sp4(fb, "lower", "bc", "- Out of range values for beta-carotene (bc):")
  sp5(fb, f, "low", "bc", "- Extreme low values for beta-carotene (bc):")
  sp5(fb, f, "high", "bc", "- Extreme high values for beta-carotene (bc):")
  bc.cc.values <- c(0.03, 0, 0.12, 0.02, 0.15, 1.38, 1.65, 1.5, 1.74, 1.76, 0.69, 1.17, 1.32,
                    1.04, 4.41, 4.92, 6.12, 5.46, 3.96, 5.49, 3.03, 3.76, 4.61, 7.23, 7.76,
                    10.5, 11.03, 12.39, 14.37, NA)
  sp3(fb, bc.cc.values, "bc.cc", "- Out of range values for beta-carotene (bc.cc):")
  sp4(fb, "lower", "tc", "- Out of range values for total carotenoids (tc):")
  sp5(fb, f, "low", "tc", "- Extreme low values for total carotenoids (tc):")
  sp5(fb, f, "high", "tc", "- Extreme high values for total carotenoids (tc):")
  sp4(fb, "lower", "star", "- Out of range values for starch (star):")
  sp5(fb, f, "low", "star", "- Extreme low values for starch (star):")
  sp5(fb, f, "high", "star", "- Extreme high values for starch (star):")
  sp4(fb, "lower", "fruc", "- Out of range values for fructose (fruc):")
  sp5(fb, f, "low", "fruc", "- Extreme low values for fructose (fruc):")
  sp5(fb, f, "high", "fruc", "- Extreme high values for fructose (fruc):")
  sp4(fb, "lower", "gluc", "- Out of range values for glucose (gluc):")
  sp5(fb, f, "low", "gluc", "- Extreme low values for glucose (gluc):")
  sp5(fb, f, "high", "gluc", "- Extreme high values for glucose (gluc):")
  sp4(fb, "lower", "sucr", "- Out of range values for sucrose (sucr):")
  sp5(fb, f, "low", "sucr", "- Extreme low values for sucrose (sucr):")
  sp5(fb, f, "high", "sucr", "- Extreme high values for sucrose (sucr):")
  sp4(fb, "lower", "malt", "- Out of range values for maltose (malt):")
  sp5(fb, f, "low", "malt", "- Extreme low values for maltose (malt):")
  sp5(fb, f, "high", "malt", "- Extreme high values for maltose (malt):")
  
  # Extreme values detection and values out of range for derived variables

  sp4(fb, "lower", "trw", "- Out of range values for total root weight (trw):")
  sp5(fb, f, "low", "trw", "- Extreme low values for total root weight (trw):")
  sp5(fb, f, "high", "trw", "- Extreme high values for total root weight (trw):")
  sp4(fb, "lower", "cytha", "- Out of range values for commercial root yield in tons per hectare (cytha):")
  sp5(fb, f, "low", "cytha", "- Extreme low values for commercial root yield in tons per hectare (cytha):")
  sp5(fb, f, "high", "cytha", "- Extreme high values for commercial root yield in tons per hectare (cytha):")
  sp4(fb, "lower", "cytha.aj", "- Out of range values for commercial root yield in tons per hectare (cytha.aj):")
  sp5(fb, f, "low", "cytha.aj", "- Extreme low values for commercial root yield in tons per hectare (cytha.aj):")
  sp5(fb, f, "high", "cytha.aj", "- Extreme high values for commercial root yield in tons per hectare (cytha.aj):")
  sp4(fb, "lower", "rytha", "- Out of range values for total root yield in tons per hectare (rytha):")
  sp5(fb, f, "low", "rytha", "- Extreme low values for total root yield in tons per hectare (rytha):")
  sp5(fb, f, "high", "rytha", "- Extreme high values for total root yield in tons per hectare (rytha):")
  sp4(fb, "lower", "rytha.aj", "- Out of range values for total root yield in tons per hectare (rytha.aj):")
  sp5(fb, f, "low", "rytha.aj", "- Extreme low values for total root yield in tons per hectare (rytha.aj):")
  sp5(fb, f, "high", "rytha.aj", "- Extreme high values for total root yield in tons per hectare (rytha.aj):")
  sp4(fb, "lower", "acrw", "- Out of range values for average commercial root weight (acrw):")
  sp5(fb, f, "low", "acrw", "- Extreme low values for average commercial root weight (acrw):")
  sp5(fb, f, "high", "acrw", "- Extreme high values for average commercial root weight (acrw):")
  sp4(fb, "lower", "nrpp", "- Out of range values for number of roots per harvested plant (nrpp):")
  sp5(fb, f, "low", "nrpp", "- Extreme low values for number of roots per harvested plant (nrpp):")
  sp5(fb, f, "high", "nrpp", "- Extreme high values for number of roots per harvested plant (nrpp):")
  sp4(fb, "lower", "nrpsp", "- Out of range values for number of roots per sowed plant (nrpsp):")
  sp5(fb, f, "low", "nrpsp", "- Extreme low values for number of roots per sowed plant (nrpsp):")
  sp5(fb, f, "high", "nrpsp", "- Extreme high values for number of roots per sowed plant (nrpsp):")
  sp4(fb, "lower", "ncrpp", "- Out of range values for number of commercial roots per harvested plant (ncrpp):")
  sp5(fb, f, "low", "ncrpp", "- Extreme low values for number of commercial roots per harvested plant (ncrpp):")
  sp5(fb, f, "high", "ncrpp", "- Extreme high values for number of commercial roots per harvested plant (ncrpp):")
  sp4(fb, "lower", "ncrpsp", "- Out of range values for number of commercial roots per sowed plant (ncrpsp):")
  sp5(fb, f, "low", "ncrpsp", "- Extreme low values for number of commercial roots per sowed plant (ncrpsp):")
  sp5(fb, f, "high", "ncrpsp", "- Extreme high values for number of commercial roots per sowed plant (ncrpsp):")
  sp4(fb, "lower", "ypp", "- Out of range values for yield per harvested plant (ypp):")
  sp5(fb, f, "low", "ypp", "- Extreme low values for yield per harvested plant (ypp):")
  sp5(fb, f, "high", "ypp", "- Extreme high values for yield per harvested plant (ypp):")
  sp4(fb, "lower", "ypsp", "- Out of range values for yield per sowed plant (ypsp):")
  sp5(fb, f, "low", "ypsp", "- Extreme low values for yield per sowed plant (ypsp):")
  sp5(fb, f, "high", "ypsp", "- Extreme high values for yield per sowed plant (ypsp):")
  sp4(fb, "lower", "vpp", "- Out of range values for vine weight per harvested plant (vpp):")
  sp5(fb, f, "low", "vpp", "- Extreme low values for vine weight per harvested plant (vpp):")
  sp5(fb, f, "high", "vpp", "- Extreme high values for vine weight per harvested plant (vpp):")
  sp4(fb, "lower", "vpsp", "- Out of range values for vine weight per sowed plant (vpsp):")
  sp5(fb, f, "low", "vpsp", "- Extreme low values for vine weight per sowed plant (vpsp):")
  sp5(fb, f, "high", "vpsp", "- Extreme high values for vine weight per sowed plant (vpsp):")
  sp4(fb, "both", "ci", "- Out of range values for commercial index (ci):")
  sp5(fb, f, "low", "ci", "- Extreme low values for commercial index (ci):")
  sp5(fb, f, "high", "ci", "- Extreme high values for commercial index (ci):")
  sp4(fb, "both", "hi", "- Out of range values for harvest index (hi):")
  sp5(fb, f, "low", "hi", "- Extreme low values for harvest index (hi):")
  sp5(fb, f, "high", "hi", "- Extreme high values for harvest index (hi):")
  sp4(fb, "both", "shi", "- Out of range values for harvest sowing index (shi):")
  sp5(fb, f, "low", "shi", "- Extreme low values for harvest sowing index (shi):")
  sp5(fb, f, "high", "shi", "- Extreme high values for harvest sowing index (shi):")
  sp4(fb, "lower", "biom", "- Out of range values for biomass yield (biom):")
  sp5(fb, f, "low", "biom", "- Extreme low values for biomass yield (biom):")
  sp5(fb, f, "high", "biom", "- Extreme high values for biomass yield (biom):")
  sp4(fb, "lower", "biom.aj", "- Out of range values for biomass yield (biom.aj):")
  sp5(fb, f, "low", "biom.aj", "- Extreme low values for biomass yield (biom.aj):")
  sp5(fb, f, "high", "biom.aj", "- Extreme high values for biomass yield (biom.aj):")
  sp4(fb, "lower", "fytha", "- Out of range values for foliage total yield in tons per hectare (fytha):")
  sp5(fb, f, "low", "fytha", "- Extreme low values for foliage total yield in tons per hectare (fytha):")
  sp5(fb, f, "high", "fytha", "- Extreme high values for foliage total yield in tons per hectare (fytha):")
  sp4(fb, "lower", "fytha.aj", "- Out of range values for foliage total yield in tons per hectare (fytha.aj):")
  sp5(fb, f, "low", "fytha.aj", "- Extreme low values for foliage total yield in tons per hectare (fytha.aj):")
  sp5(fb, f, "high", "fytha.aj", "- Extreme high values for foliage total yield in tons per hectare (fytha.aj):")
  sp4(fb, "both", "rfr", "- Out of range values for root foliage ratio (rfr):")
  sp5(fb, f, "low", "rfr", "- Extreme low values for root foliage ratio (rfr):")
  sp5(fb, f, "high", "rfr", "- Extreme high values for root foliage ratio (rfr):")

  # Extreme values detection and values out of range for additional traits
  
  if (!is.null(aqt)) {
    for (i in 1:length(aqt)) {
      sp5(fb, f, "low", aqt[i], paste("- Extreme low values for (", aqt[i], "):", sep = ""))
      sp5(fb, f, "high", aqt[i], paste("- Extreme high values for (", aqt[i], "):", sep = ""))
    }
  }
  
  # Outliers' detection
  
  oc <- 0
  
  if (out.mod == "rcbd") {
    oc <- 1
    if (exists('g', fb)) {
      geno <- fb[, 'g']
    } else {
      if (exists('geno', fb)) {
        geno <- fb[, 'geno']
      } else {
        if (exists('name', fb)) {
          geno <- fb[, 'name']
        } else {
          oc <- 0
          warning("Genotypes are not defined. Use g or geno as labels.", call. = FALSE)  
        }
      }
    }
    if (exists('r', fb)) {
      rep <- fb[, 'r']
    } else {
      if (exists('rep', fb)) {
        rep <- fb[, 'rep']
      } else {
        oc <- 0
        warning('Blocks are not defined. Use r or rep as labels.', call. = FALSE)
      }
    }
    env <- NULL
  }
  
  if (out.mod == "met") {
    oc <- 1
    if (exists('g', fb)) {
      geno <- fb[, 'g']
    } else {
      if (exists('geno', fb)) {
        geno <- fb[, 'geno']
      } else {
        if (exists('name', fb)) {
          geno <- fb[, 'name']
        } else {
          oc <- 0
          warning("Genotypes are not defined. Use g or geno as labels.", call. = FALSE)  
        }
      }
    }
    if (exists('e', fb)) {
      env <- fb[, 'e']
    } else {
      if (exists('env', fb)) {
        env <- fb[, 'env']
      } else {
        oc <- 0
        warning('Environments are not defined. Use e or env as labels.', call. = FALSE)
      }
    }
    if (exists('r', fb)) {
      rep <- fb[, 'r']
    } else {
      if (exists('rep', fb)) {
        rep <- fb[, 'rep']
      } else {
        oc <- 0
        warning('Blocks are not defined. Use r or rep as labels.', call. = FALSE)
      }
    }
  }
  
  if (oc == 1) {
    sp6(fb, geno, env, rep, "vw", out.mod, out.max, "- Outliers for vine weight (vw):")
    sp6(fb, geno, env, rep, "nocr", out.mod, out.max, "- Outliers for number of commercial roots (nocr):")
    sp6(fb, geno, env, rep, "nonc", out.mod, out.max, "- Outliers for number of non commercial roots (nonc):")
    sp6(fb, geno, env, rep, "crw", out.mod, out.max, "- Outliers for commercial root weight (crw):")
    sp6(fb, geno, env, rep, "ncrw", out.mod, out.max, "- Outliers for non commercial root weight (ncrw):")
    sp6(fb, geno, env, rep, "dm", out.mod, out.max, "- Outliers for storage root dry matter content (dm):")
    sp6(fb, geno, env, rep, "dmry", out.mod, out.max, "- Outliers for dry matter root yield (dmry):")
    sp6(fb, geno, env, rep, "dmv", out.mod, out.max, "- Outliers for vines dry matter content (dmv):")
    sp6(fb, geno, env, rep, "dmvy", out.mod, out.max, "- Outliers for dry matter vine yield (dmvy):")
    sp6(fb, geno, env, rep, "dmbiom", out.mod, out.max, "- Outliers for dry matter biomass (dmbiom):")
    sp6(fb, geno, env, rep, "prot", out.mod, out.max, "- Outliers for protein (prot):")
    sp6(fb, geno, env, rep, "fe", out.mod, out.max, "- Outliers for iron (fe):")
    sp6(fb, geno, env, rep, "zn", out.mod, out.max, "- Outliers for zinc (zn):")
    sp6(fb, geno, env, rep, "ca", out.mod, out.max, "- Outliers for calcium (ca):")
    sp6(fb, geno, env, rep, "mg", out.mod, out.max, "- Outliers for magnesium (mg):")
    sp6(fb, geno, env, rep, "bc", out.mod, out.max, "- Outliers for beta-carotene (bc):")
    sp6(fb, geno, env, rep, "bc.cc", out.mod, out.max, "- Outliers for beta-carotene (bc.cc):")
    sp6(fb, geno, env, rep, "tc", out.mod, out.max, "- Outliers for total carotenoids (tc):")
    sp6(fb, geno, env, rep, "star", out.mod, out.max, "- Outliers for starch (star):")
    sp6(fb, geno, env, rep, "fruc", out.mod, out.max, "- Outliers for fructose (fruc):")
    sp6(fb, geno, env, rep, "gluc", out.mod, out.max, "- Outliers for glucose (gluc):")
    sp6(fb, geno, env, rep, "sucr", out.mod, out.max, "- Outliers for sucrose (sucr):")
    sp6(fb, geno, env, rep, "malt", out.mod, out.max, "- Outliers for maltose (malt):")
    sp6(fb, geno, env, rep, "trw", out.mod, out.max, "- Outliers for total root weight (trw):")
    sp6(fb, geno, env, rep, "cytha", out.mod, out.max, "- Outliers for commercial root yield in tons per hectare (cytha):")
    sp6(fb, geno, env, rep, "cytha.aj", out.mod, out.max, "- Outliers for commercial root yield in tons per hectare (cytha.aj):")
    sp6(fb, geno, env, rep, "rytha", out.mod, out.max, "- Outliers for total root yield in tons per hectare (rytha):")
    sp6(fb, geno, env, rep, "rytha.aj", out.mod, out.max, "- Outliers for total root yield in tons per hectare (rytha.aj):")
    sp6(fb, geno, env, rep, "acrw", out.mod, out.max, "- Outliers for average commercial root weight (acrw):")
    sp6(fb, geno, env, rep, "nrpp", out.mod, out.max, "- Outliers for number of roots per harvested plant (nrpp):")
    sp6(fb, geno, env, rep, "nrpsp", out.mod, out.max, "- Outliers for number of roots per sowed plant (nrpsp):")
    sp6(fb, geno, env, rep, "ncrpp", out.mod, out.max, "- Outliers for number of commercial roots per harvested plant (ncrpp):")
    sp6(fb, geno, env, rep, "ncrpsp", out.mod, out.max, "- Outliers for number of commercial roots per sowed plant (ncrpsp):")
    sp6(fb, geno, env, rep, "ypp", out.mod, out.max, "- Outliers for yield per harvested plant (ypp):")
    sp6(fb, geno, env, rep, "ypsp", out.mod, out.max, "- Outliers for yield per sowed plant (ypsp):")
    sp6(fb, geno, env, rep, "vpp", out.mod, out.max, "- Outliers for vine weight per harvested plant (vpp):")
    sp6(fb, geno, env, rep, "vpsp", out.mod, out.max, "- Outliers for vine weight per sowed plant (vpsp):")
    sp6(fb, geno, env, rep, "ci", out.mod, out.max, "- Outliers for commercial index (ci):")
    sp6(fb, geno, env, rep, "hi", out.mod, out.max, "- Outliers for harvest index (hi):")
    sp6(fb, geno, env, rep, "shi", out.mod, out.max, "- Outliers for harvest sowing index (shi):")
    sp6(fb, geno, env, rep, "biom", out.mod, out.max, "- Outliers for biomass yield (biom):")
    sp6(fb, geno, env, rep, "biom.aj", out.mod, out.max, "- Outliers for biomass yield (biom.aj):")
    sp6(fb, geno, env, rep, "fytha", out.mod, out.max, "- Outliers for foliage total yield in tons per hectare (fytha):")
    sp6(fb, geno, env, rep, "fytha.aj", out.mod, out.max, "- Outliers for foliage total yield in tons per hectare (fytha.aj):")
    sp6(fb, geno, env, rep, "rfr", out.mod, out.max, "- Outliers for root foliage ratio (rfr):")
    
    # Outliers' detection for additional traits
    
    if (!is.null(aqt))
      for (i in 1:length(aqt))
        sp6(fb, geno, env, rep, aqt[i], out.mod, out.max, paste("- Outlilers for (", aqt[i], "):", sep = ""))
    
  }
}

# Print results function

output <- function(fb, cond, tx) {
  if (sum(cond, na.rm = TRUE) > 0) {
    cat("\n", tx, "\n", sep = "")
    print(subset(fb, cond))
  }
}

# Two traits conditions

sp1 <- function(fb, type, t1, t2, tx) {
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

# Roots and dependencies

sp2 <- function(fb, temp, do, t1, tx) {
  if (exists(t1, fb) & do == TRUE) {
    cond <- (temp == FALSE) & (fb[, t1] > 0)
    output(fb, cond, tx)
  }
}

# Detect out of discrete range

sp3 <- function(fb, vv, t1, tx) {
  if (exists(t1, fb)) {
    cond <- !(fb[, t1] %in% vv)
    output(fb, cond, tx)
  }
}

# Detect out of continous range

sp4 <- function(fb, ex, t1, tx) { 
  if (exists(t1, fb)) {
    if (ex == "lower")
      cond <- fb[, t1] < 0
    if (ex == "both")
      cond <- fb[, t1] < 0 | fb[, t1] > 100
    output(fb, cond, tx)
  }
}

# Extreme values

sp5 <- function(fb, f, ex, t1, tx) { 
  if (exists(t1, fb)) {
    if (ex == "low")
      cond <- fb[, t1] < quantile(fb[, t1], 0.25, na.rm = TRUE) - f * IQR(fb[, t1], na.rm = TRUE)
    if (ex == "high")
      cond <- fb[, t1] > quantile(fb[, t1], 0.75, na.rm = TRUE) + f * IQR(fb[, t1], na.rm = TRUE)
    output(fb, cond, tx)
  } 
}

# Outliers' detection

sp6 <- function(fb, geno, env, rep, t1, out.mod, out.max, tx) {
  if (exists(t1, fb)) {
    fb$id <- 1:dim(fb)[1]
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
