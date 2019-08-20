#' Check fieldbook traits names for sweetpotato
#'
#' Check that fieldbook traits names correspond with the names defined in 
#' Procedures for the evaluation and analysis of sweetpotato trials, ISBN 978-92-9060-522-5.
#' @param dfr The name of the data frame.
#' @param add Additional traits.
#' @details The data frame must use the following labels (lower or upper case):
#' 
#'  -------------------- Plot identifiers --------------------
#'  \itemize{
#'  \item \code{plot}      : Plot number
#'  \item \code{row}       : Row number
#'  \item \code{col}       : Column number
#'  }
#' -------------------- Classification factors --------------------
#'  \itemize{
#'  \item \code{l}         : Locations (\code{loc} is also valid)
#'  \item \code{y}         : Years (\code{year} is also valid)
#'  \item \code{s}         : Seasons (\code{season} is also valid) 
#'  \item \code{e}         : Environments (\code{env} is also valid)
#'  \item \code{g}         : Genotypes (\code{geno} is also valid)
#'  \item \code{cipno}     : Institutional CIP number
#'  \item \code{r}         : Replications (\code{rep} is also valid)
#'  \item \code{b}         : Blocks (\code{block} is also valid)
#'  \item \code{treat}     : Treatment
#'  \item \code{harvest}   : Harvest time
#'  }
#' -------------------- On field evaluated traits --------------------
#'  \itemize{
#'  \item \code{nops}      : Number of plants sowed
#'  \item \code{nope}      : Number of plants established
#'  \item \code{noph}      : Number of plants harvested
#'  \item \code{vir1}      : Virus symptoms (1-9), first evaluation, 4-8 weeks after planting (\code{vir} is also valid)
#'  \item \code{vir2}      : Virus symptoms (1-9), second evaluation, 1 month before harvest
#'  \item \code{alt1}      : Alternaria symptoms (1-9), first evaluation, 4-8 weeks after planting (\code{alt} is also valid)
#'  \item \code{alt2}      : Alternaria symptoms (1-9), second evaluation, 1 month before harvest
#'  \item \code{vv}        : Vine vigor (1-9), first evaluation, 1 month before harvest
#'  \item \code{vw}        : Vine weight (kg/plot)
#'  \item \code{nopr}      : Number of plants with roots
#'  \item \code{nocr}      : Number of commercial roots
#'  \item \code{nonc}      : Number of non commercial roots
#'  \item \code{crw}       : Commercial root weight (kg/plot)
#'  \item \code{ncrw}      : Non commercial root weight (kg/plot)
#'  \item \code{scol}      : Storage root skin color (1-9)
#'  \item \code{fcol}      : Storage root flesh color (1-9)
#'  \item \code{fcol.cc}   : Root flesh color using RHS color charts (1-30)
#'  \item \code{rs}        : Root size (1-9)
#'  \item \code{rf}        : Root form (1-9)
#'  \item \code{damr}      : Root defects (1-9)
#'  \item \code{rspr}      : Root sprouting (1-9)
#'  \item \code{wed}       : Weevil damage (1-9)
#'  }
#' -------------------- Dry matter assesment --------------------
#'  \itemize{
#'  \item \code{dmf}       : Fresh weight of roots for dry matter assessment (g)
#'  \item \code{dmd}       : Dry weight of roots for dry matter assessment (g)
#'  \item \code{dm}        : Storage root dry matter content (\%)
#'  \item \code{dmvf}      : Fresh weight vines for dry matter assessment (g)
#'  \item \code{dmvd}      : Dry weight vines for dry matter assessment (g)
#'  \item \code{dmv}       : Vines dry matter content (\%)
#'  }
#' -------------------- Raw and cooked roots attributes --------------------
#'  \itemize{
#'  \item \code{fraw1}     : Raw root fiber (1-9), first determination (\code{fraw} is also valid)
#'  \item \code{suraw1}    : Raw root sugar (1-9), first determination (\code{suraw} is also valid)
#'  \item \code{straw1}    : Raw root starch (1-9), first determination (\code{straw} is also valid)
#'  \item \code{coof1}     : Cooked root fiber (1-9), first evaluation (\code{coof} is also valid)
#'  \item \code{coosu1}    : Cooked root sugars (1-9), first evaluation (\code{coosu} is also valid)
#'  \item \code{coost1}    : Cooked root starch (1-9), first evaluation (\code{coost} is also valid)
#'  \item \code{coot1}     : Cooked root taste (1-9), first evaluation (\code{coot} is also valid)
#'  \item \code{cooap1}    : Cooked root appearance (1-9), first evaluation (\code{cooap} is also valid)
#'  \item \code{fraw2}     : Raw root fiber (1-9), second determination
#'  \item \code{suraw2}    : Raw root sugar (1-9), second determination
#'  \item \code{straw2}    : Raw root starch (1-9), second determination
#'  \item \code{coof2}     : Cooked root fiber (1-9), second evaluation
#'  \item \code{coosu2}    : Cooked root sugars (1-9), second evaluation
#'  \item \code{coost2}    : Cooked root starch (1-9), second evaluation
#'  \item \code{coot2}     : Cooked root taste (1-9), second evaluation
#'  \item \code{cooap2}    : Cooked root appearance (1-9), second evaluation
#'  }
#' -------------------- Nutrients evaluations --------------------
#'  \itemize{
#'  \item \code{prot}      : Protein (\% raw fresh)
#'  \item \code{fe}        : Iron (mg/100g raw dry weight)
#'  \item \code{zn}        : Zinc (mg/100g raw dry weight)
#'  \item \code{ca}        : Calcium (mg/100g raw dry weight)
#'  \item \code{mg}        : Magnesium (mg/100g raw dry weight)
#'  \item \code{bc}        : Beta-carotene (mg/100g raw dry weight)
#'  \item \code{bc.cc}     : Beta-carotene with RHS color charts (mg/100g raw fresh weight)
#'  \item \code{tc}        : Total carotenoids (mg/100g raw dry weight)
#'  \item \code{star}      : Starch (\% raw fresh)
#'  \item \code{fruc}      : Fructose (\% raw fresh)
#'  \item \code{gluc}      : Glucose (\% raw fresh)
#'  \item \code{sucr}      : Sucrose (\% raw fresh)
#'  \item \code{malt}      : Maltose (\% raw fresh)
#'  }
#' -------------------- Calculated traits --------------------
#'  \itemize{
#'  \item \code{tnr}       : Total number of roots per plot
#'  \item \code{trw}       : Total root weight (kg/plot)
#'  \item \code{trw.d}     : Total root dry weight (kg/plot)
#'  \item \code{biom}      : Biomass yield (kg/plot)
#'  \item \code{biom.d}    : Biomass dry yield (kg/plot)
#'  \item \code{cytha}     : Commercial root yield (t/ha)
#'  \item \code{cytha.aj}  : Commercial root yield (t/ha) adjusted by number of harvested plants
#'  \item \code{rytha}     : Total root yield (t/ha)
#'  \item \code{rytha.aj}  : Total root yield (t/ha) adjusted by number of harvested plants
#'  \item \code{dmry}      : Dry matter root yield (t/ha)
#'  \item \code{dmry.aj}   : Dry matter root yield (t/ha) adjusted by number of harvested plants
#'  \item \code{vw.d}      : Vine dry weight (kg/plot)
#'  \item \code{fytha}     : Foliage total yield (t/ha)
#'  \item \code{fytha.aj}  : Foliage total yield (t/ha) adjusted by number of harvested plants
#'  \item \code{dmvy}      : Dry matter vine yield (t/ha)
#'  \item \code{dmvy.aj}   : Dry matter vine yield (t/ha) adjusted by number of harvested plants
#'  \item \code{bytha}     : Biomass yield (t/ha)
#'  \item \code{bytha.aj}  : Biomass yield (t/ha) adjusted by number of harvested plants
#'  \item \code{dmby}      : Dry matter biomass (t/ha)
#'  \item \code{dmby.aj}   : Dry matter biomass (t/ha) adjusted by number of harvested plants
#'  \item \code{acrw}      : Average commercial root weight (kg/root)
#'  \item \code{ancrw}     : Average non commercial root weight (kg/root)
#'  \item \code{atrw}      : Average total root weight (kg/root)
#'  \item \code{nrpp}      : Number of roots per harvested plant
#'  \item \code{nrpsp}     : Number of roots per sowed plant
#'  \item \code{ncrpp}     : Number of commercial roots per harvested plant
#'  \item \code{ncrpsp}    : Number of commercial roots per sowed plant
#'  \item \code{ypp}       : Yield per harvested plant (kg/plant)
#'  \item \code{ypsp}      : Yield per sowed plant (kg/plant)
#'  \item \code{vpp}       : Vine weight per harvested plant (kg/plant)
#'  \item \code{vpsp}      : Vine weight per sowed plant (kg/plant)
#'  \item \code{ci}        : Commercial index (\%)
#'  \item \code{hi}        : Harvest index (\%)
#'  \item \code{shi}       : Harvest sowing index (\%)
#'  \item \code{rfr}       : Root foliage ratio (\%)
#'  }
#' @return It returns a data frame with all traits names in lower case, and a list of the
#' traits with names not included in the list shown above.
#' @author Raul Eyzaguirre.
#' @examples
#' check.names.sp(pjpz09)
#' @export

check.names.sp <- function(dfr, add = NULL) {
  
  plot.id <- c("plot", "row", "col")
  
  factors <- c("l", "loc", "y", "year", "s", "season", "e", "env", "g", "geno", "cipno",
               "r", "rep", "b", "block", "treat", "harvest")
  
  traits <- c("nops", "nope", "noph", "vir", "vir1", "vir2", "alt", "alt1", "alt2", "vv",
              "vw", "nopr", "nocr", "nonc", "crw", "ncrw", "scol", "fcol", "fcol.cc",
              "rs", "rf", "damr", "rspr", "wed", "dmf", "dmd", "dm", "dmvf", "dmvd",
              "dmv", "fraw", "fraw1", "suraw", "suraw1", "straw", "straw1", "coof",
              "coof1", "coosu", "coosu1", "coost", "coost1", "coot", "coot1", "cooap",
              "cooap1", "fraw2", "suraw2", "straw2", "coof2", "coosu2", "coost2", "coot2",
              "cooap2", "prot", "fe", "zn", "ca", "mg", "bc", "bc.cc", "tc", "star",
              "fruc", "gluc", "sucr", "malt", "tnr", "trw", "trw.d", "biom", "biom.d",
              "cytha", "cytha.aj", "rytha", "rytha.aj", "dmry", "dmry.aj", "vw.d", "fytha",
              "fytha.aj", "dmvy", "dmvy.aj", "bytha", "bytha.aj", "dmby", "dmby.aj", "acrw",
              "ancrw", "atrw", "nrpp", "nrpsp", "ncrpp", "ncrpsp", "ypp", "ypsp", "vpp",
              "vpsp", "ci", "hi", "shi", "rfr")
  
  colnames.valid <- c(plot.id, factors, traits, tolower(add))
    
  colnames.list <- colnames(dfr)
  
  check.list.1 <- !(tolower(colnames.list) %in% colnames.valid) # which are not valid
  temp <- colnames.list[!check.list.1]                          # list of valid names
  check.list.2 <- !(temp %in% colnames.valid)                   # which are valid but lower case
  
  colnames(dfr) <- tolower(colnames(dfr))
    
  # Warnings
  
  if (max(check.list.1) == 1)
    warning("Some columns with invalid names: ", list(colnames.list[check.list.1]), call. = FALSE)
  
  if (max(check.list.2) == 1)
    warning("Some labels converted to lower case: ", list(temp[check.list.2]), call. = FALSE)
  
  # Return
  
  dfr
  
}
