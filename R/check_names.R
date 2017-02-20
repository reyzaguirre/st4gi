#' Check fieldbook traits names
#'
#' Check that fieldbook traits names correspond with the names defined in the document
#' "PROCEDURES  FOR THE EVALUATION AND ANALYSIS OF SWEETPOTATO TRIALS".
#' @param fb The name of the fieldbook data frame.
#' @param aqt Additional quantitative traits.
#' @details The data frame must use the following labels (lower or upper case):
#' \itemize{
#'  \item \code{l}       : Locations (\code{loc} is also valid)
#'  \item \code{y}       : Years
#'  \item \code{s}       : Seasons
#'  \item \code{g}       : Genotypes (\code{geno} is also valid)
#'  \item \code{name}    : Names for genotypes
#'  \item \code{e}       : Environments (\code{env} is also valid)
#'  \item \code{r}       : Replications (\code{rep} is also valid)
#'  \item \code{nops}    : Number of plants sowed
#'  \item \code{nope}    : Number of plants established
#'  \item \code{vir1}    : Virus symptoms (1-9), first evaluation, 4-8 weeks after planting (\code{vir} is also valid)
#'  \item \code{vir2}    : Virus symptoms (1-9), second evaluation, 1 month before harvest
#'  \item \code{alt1}    : Alternaria symptoms (1-9), first evaluation, 4-8 weeks after planting (\code{alt} is also valid)
#'  \item \code{alt2}    : Alternaria symptoms (1-9), second evaluation, 1 month before harvest
#'  \item \code{vv}      : Vine vigor (1-9), 1 month before harvest
#'  \item \code{vw}      : Vine weight (kg/plot)
#'  \item \code{noph}    : Number of plants harvested
#'  \item \code{nopr}    : Number of plants with roots
#'  \item \code{nocr}    : Number of commercial roots
#'  \item \code{nonc}    : Number of non commercial roots
#'  \item \code{crw}     : Commercial root weight (kg/plot)
#'  \item \code{ncrw}    : Non commercial root weight (kg/plot)
#'  \item \code{rfcp.cc} : Root primary flesh color using CIP color charts
#'  \item \code{rfcs.cc} : Root secondary flesh color using CIP color charts
#'  \item \code{scol}    : Storage root skin color (1-9)
#'  \item \code{fcol}    : Storage root flesh color (1-9)
#'  \item \code{rfcp}    : Storage root primary flesh color (1-9)
#'  \item \code{rfcs}    : Storage root secondary flesh color (1-9)
#'  \item \code{rs}      : Root size (1-9)
#'  \item \code{rf}      : Root form (1-9)
#'  \item \code{damr}    : Root defects (1-9)
#'  \item \code{rspr}    : Root sprouting (1-9)
#'  \item \code{wed}     : Weevil damage (1-9)
#'  \item \code{dmf}     : Fresh weight of roots for dry matter assessment (kg)
#'  \item \code{dmd}     : Dry weight of dmf samples (kg)
#'  \item \code{dm}      : Storage root dry matter content (\%)
#'  \item \code{dmvf}    : Fresh weight vines for dry matter assessment (kg)
#'  \item \code{dmvd}    : Dry weight of dmvf samples (kg)
#'  \item \code{dmv}     : Vines dry matter content (\%)
#'  \item \code{fraw1}   : Root fiber (1-9), first determination (\code{fraw} is also valid)
#'  \item \code{suraw1}  : Root sugar (1-9), first determination (\code{suraw} is also valid)
#'  \item \code{straw1}  : Root starch (1-9), first determination (\code{straw} is also valid)
#'  \item \code{coof1}   : Cooked fiber (1-9), first evaluation (\code{coof} is also valid)
#'  \item \code{coosu1}  : Cooked sugars (1-9), first evaluation (\code{coosu} is also valid)
#'  \item \code{coost1}  : Cooked starch (1-9), first evaluation (\code{coost} is also valid)
#'  \item \code{coot1}   : Cooked taste (1-9), first evaluation (\code{coot} is also valid)
#'  \item \code{cooap1}  : Cooked appearance (1-9), first evaluation (\code{cooap} is also valid)
#'  \item \code{fraw2}   : Root fiber (1-9), second determination
#'  \item \code{suraw2}  : Root sugar (1-9), second determination
#'  \item \code{straw2}  : Root starch (1-9), second determination
#'  \item \code{coof2}   : Cooked fiber (1-9), second evaluation
#'  \item \code{coosu2}  : Cooked sugars (1-9), second evaluation
#'  \item \code{coost2}  : Cooked starch (1-9), second evaluation
#'  \item \code{coot2}   : Cooked taste (1-9), second evaluation
#'  \item \code{cooap2}  : Cooked appearance (1-9), second evaluation
#'  \item \code{prot}    : Protein (\%)
#'  \item \code{fe}      : Iron (mg/100 g dry weight)
#'  \item \code{zn}      : Zinc (mg/100 g dry weight)
#'  \item \code{ca}      : Calcium (mg/100 g dry weight)
#'  \item \code{mg}      : Magnesium (mg/100 g dry weight)
#'  \item \code{bc}      : Beta-carotene (mg/100 g dry weight)
#'  \item \code{bc.cc}   : Beta-carotene with color charts (mg/100 g fresh weight)
#'  \item \code{tc}      : Total carotenoids (mg/100 g dry weight)
#'  \item \code{star}    : Starch (\%)
#'  \item \code{fruc}    : Fructose (\%)
#'  \item \code{gluc}    : Glucose (\%)
#'  \item \code{sucr}    : Sucrose (\%)
#'  \item \code{malt}    : Maltose (\%)
#'  \item \code{trw}     : Total root weight (kg/plot)
#'  \item \code{trw.d}   : Total root dry weight (kg/plot)
#'  \item \code{cytha}   : Commercial root yield (t/ha)
#'  \item \code{cytha.aj}: Commercial root yield (t/ha) adjusted by number of harvested plants
#'  \item \code{rytha}   : Total root yield (t/ha)
#'  \item \code{rytha.aj}: Total root yield (t/ha) adjusted by number of harvested plants
#'  \item \code{dmry}    : Dry matter root yield (t/ha)
#'  \item \code{dmry.aj} : Dry matter root yield (t/ha) adjusted by number of harvested plants
#'  \item \code{vw.d}    : Vine dry weight (kg/plot)
#'  \item \code{fytha}   : Foliage total yield (t/ha)
#'  \item \code{fytha.aj}: Foliage total yield (t/ha) adjusted by number of harvested plants
#'  \item \code{dmvy}    : Dry matter vine yield (t/ha)
#'  \item \code{dmvy.aj} : Dry matter vine yield (t/ha) adjusted by number of harvested plants
#'  \item \code{biom}    : Biomass yield (t/ha)
#'  \item \code{biom.aj} : Biomass yield (t/ha) adjusted by number of harvested plants
#'  \item \code{acrw}    : Average commercial root weight (kg/root)
#'  \item \code{nrpp}    : Number of roots per plant
#'  \item \code{ncrpp}   : Number of commercial roots per plant
#'  \item \code{ypp}     : Yield per plant (kg/plant)
#'  \item \code{ci}      : Percent marketable roots (commercial index)
#'  \item \code{hi}      : Harvest index
#'  \item \code{shi}     : Harvest sowing index (survival)
#'  \item \code{rfr}     : Root foliage ratio
#'  }
#' @return It returns a data frame with all traits names in lower case, and a list of the
#' traits with names not included in the list shown above.
#' @author Raul Eyzaguirre.
#' @examples
#' checknames(pjpz09)
#' @export

checknames <- function(fb, aqt = NULL) {

  colnames.valid <- c("l", "loc", "y", "s", "g", "geno", "name", "e", "env", "r", "rep",
                      "nops", "nope", "vir", "vir1", "vir2", "alt", "alt1", "alt2", "vv",
                      "vw", "noph", "nopr", "nocr", "nonc", "crw", "ncrw", "rfcp.cc",
                      "rfcs.cc", "scol", "fcol", "rfcp", "rfcs", "rs", "rf", "damr",
                      "rspr", "wed", "dmf", "dmd", "dm", "dmvf", "dmvd", "dmv", "fraw",
                      "fraw1", "suraw", "suraw1", "straw", "straw1", "coof", "coof1",
                      "coosu", "coosu1", "coost", "coost1", "coot", "coot1", "cooap",
                      "cooap1", "fraw2", "suraw2", "straw2", "coof2", "coosu2", "coost2",
                      "coot2", "cooap2", "prot", "fe", "zn", "ca", "mg", "bc", "bc.cc",
                      "tc", "star", "fruc", "gluc", "sucr", "malt", "trw", "trw.d", "cytha",
                      "cytha.aj", "rytha", "rytha.aj", "dmry", "dmry.aj", "vw.d", "fytha",
                      "fytha.aj", "dmvy", "dmvy.aj", "biom", "biom.aj", "acrw", "nrpp",
                      "ncrpp", "ypp", "ci", "hi", "shi", "rfr", tolower(aqt))
    
  colnames.list <- colnames(fb)
  
  check.list.1 <- !(tolower(colnames.list) %in% colnames.valid) # which are not valid
  temp <- colnames.list[!check.list.1]                          # list of valid names
  check.list.2 <- !(temp %in% colnames.valid)                   # which are valid but lower case
  
  colnames(fb) <- tolower(colnames(fb))
    
  # Warnings
  
  if (max(check.list.1) == 1)
    warning("Some columns with invalid names: ", list(colnames.list[check.list.1]), call. = FALSE)
  
  if (max(check.list.2) == 1)
    warning("Some labels converted to lower case: ", list(temp[check.list.2]), call. = FALSE)
  
  fb
}
