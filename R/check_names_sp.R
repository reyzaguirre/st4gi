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
#'  \item \code{plot}      : Plot number.
#'  \item \code{row}       : Row number.
#'  \item \code{col}       : Column number.
#'  }
#' -------------------- Classification factors --------------------
#'  \itemize{
#'  \item \code{l}         : Locations (\code{loc} is also valid).
#'  \item \code{y}         : Years (\code{year} is also valid).
#'  \item \code{s}         : Seasons (\code{season} is also valid). 
#'  \item \code{e}         : Environments (\code{env} is also valid).
#'  \item \code{g}         : Genotypes (\code{geno} is also valid).
#'  \item \code{cipno}     : Institutional CIP number.
#'  \item \code{r}         : Replications (\code{rep} is also valid).
#'  \item \code{b}         : Blocks (\code{block} is also valid).
#'  \item \code{treat}     : Treatment.
#'  \item \code{harvest}   : Harvest time.
#'  }
#' -------------------- On field evaluated traits --------------------
#'  \itemize{
#'  \item \code{nops}      : Number of plants sowed (CO_331:0000678).
#'  \item \code{nope}      : Number of plants established (CO_331:0000192).
#'  \item \code{noph}      : Number of plants harvested (CO_331:0000679).
#'  \item \code{vir}       : Virus symptoms, 1-9 scale.
#'  \item \code{vir1}      : Virus symptoms, 1-9 scale, first evaluation, 4-8 weeks after planting (CO_331:0000193).
#'  \item \code{vir2}      : Virus symptoms, 1-9 scale, second evaluation, 1 month before harvest.
#'  \item \code{alt}       : Alternaria symptoms, 1-9 scale.
#'  \item \code{alt1}      : Alternaria symptoms, 1-9 scale, first evaluation, 4-8 weeks after planting (CO_331:0000198).
#'  \item \code{alt2}      : Alternaria symptoms, 1-9 scale, second evaluation, 1 month before harvest.
#'  \item \code{vv}        : Vine vigor, 1-9 scale (CO_331:0000197).
#'  \item \code{vw}        : Vine weight in kg/plot (CO_331:0000227).
#'  \item \code{nopr}      : Number of plants with roots (CO_331:0000211).
#'  \item \code{nocr}      : Number of commercial roots (CO_331:0000214).
#'  \item \code{nonc}      : Number of non commercial roots (CO_331:0000217).
#'  \item \code{crw}       : Commercial root weight in kg/plot (CO_331:0000220).
#'  \item \code{ncrw}      : Non commercial root weight in kg/plot (CO_331:0000223).
#'  \item \code{scol}      : Storage root skin color, 1-9 scale, (CO_331:0000175).
#'  \item \code{fcol}      : Storage root flesh color, 1-9 scale (CO_331:0000178).
#'  \item \code{fcol.cc}   : Root flesh color using RHS color charts, 1-30 scale.
#'  \item \code{rs}        : Root size, 1-9 scale, (CO_331:0000184).
#'  \item \code{rf}        : Root form, 1-9 scale, (CO_331:0000202).
#'  \item \code{damr}      : Root defects, 1-9 scale, (CO_331:0000206).
#'  \item \code{rspr}      : Root sprouting, 1-9 scale, (CO_331:0000277).
#'  \item \code{alcdam}    : Alcidodes sp. damage, 1-9 scale (CO_331:0000806).
#'  \item \code{wed}       : Weevil damage, 1-9 scale (CO_331:0000207).
#'  }
#' -------------------- Dry matter assesment --------------------
#'  \itemize{
#'  \item \code{dmf}       : Fresh weight of roots for dry matter assessment in g (CO_331:0000243). 
#'  \item \code{dmd}       : Dry weight of roots for dry matter assessment in g (CO_331:0000247).
#'  \item \code{dm}        : Storage root dry matter content \% (CO_331:0000297).
#'  \item \code{dmvf}      : Fresh weight vines for dry matter assessment in g (CO_331:0000251).
#'  \item \code{dmvd}      : Dry weight vines for dry matter assessment in g (CO_331:0000255).
#'  \item \code{dmv}       : Vines dry matter content \% (CO_331:0001014).
#'  }
#' -------------------- Raw and cooked roots attributes --------------------
#'  \itemize{
#'  \item \code{fraw}      : Raw root fiber, 1-9 scale (CO_331:0001038).
#'  \item \code{suraw}     : Raw root sugar, 1-9 scale (CO_331:0000766).
#'  \item \code{straw}     : Raw root starch, 1-9 scale (CO_331:0001055).
#'  \item \code{coof}      : Cooked root fiber, 1-9 (CO_331:0000259).
#'  \item \code{coosu}     : Cooked root sugars, 1-9 scale (CO_331:0000261).
#'  \item \code{coost}     : Cooked root starch, 1-9 scale (CO_331:0000265).
#'  \item \code{coot}      : Cooked root taste, 1-9 scale (CO_331:0000269).
#'  \item \code{cooap}     : Cooked root appearance, 1-9 scale (CO_331:0000273).
#'  }
#' -------------------- Nutrients evaluations --------------------
#'  \itemize{
#'  \item \code{prot}      : Protein, g/100g raw fresh (CO_331:0001010).
#'  \item \code{fe}        : Iron, mg/100g raw dry weight measured by NIRS (CO_331:0001016).
#'  \item \code{zn}        : Zinc, mg/100g raw dry weight measured by NIRS (CO_331:0001017).
#'  \item \code{ca}        : Calcium, mg/100g raw dry weight measured by NIRS (CO_331:0001029).
#'  \item \code{mg}        : Magnesium, mg/100g raw dry weight measured by NIRS (CO_331:0001030).
#'  \item \code{bc}        : Beta-carotene, mg/100g raw dry weight measured by NIRS (CO_331:0000289).
#'  \item \code{bc.cc}     : Beta-carotene with RHS color charts, mg/100g raw fresh weight (CO_331:0001023).
#'  \item \code{tc}        : Total carotenoids, mg/100g raw dry weight (CO_331:0000290).
#'  \item \code{star}      : Starch, g/100g raw fresh (CO_331:0001012).
#'  \item \code{fruc}      : Fructose, g/100g raw fresh (CO_331:0000292 - CO_331:0000979).
#'  \item \code{gluc}      : Glucose, g/100g raw fresh (CO_331:0000293 - CO_331:0000981).
#'  \item \code{sucr}      : Sucrose, g/100g raw fresh (CO_331:0000294 - CO_331:0000983).
#'  \item \code{malt}      : Maltose, g/100g raw fresh (CO_331:0000295 - CO_331:0000985).
#'  }
#' -------------------- Calculated traits --------------------
#'  \itemize{
#'  \item \code{tnr}       : Total number of roots per plot (CO_331:0000233).
#'  \item \code{trw}       : Total root weight in kg/plot (CO_331:0000237).
#'  \item \code{trw.d}     : Total root dry weight in kg/plot (CO_331:0000995).
#'  \item \code{biom}      : Biomass yield in kg/plot (CO_331:0000683).
#'  \item \code{biom.d}    : Biomass dry yield in kg/plot (CO_331:0000997).
#'  \item \code{cytha}     : Commercial root yield in t/ha (CO_331:0000809).
#'  \item \code{cytha.aj}  : Commercial root yield in t/ha adjusted by number of harvested plants (CO_331:0000998).
#'  \item \code{rytha}     : Total root yield in t/ha (CO_331:0000296).
#'  \item \code{rytha.aj}  : Total root yield in t/ha adjusted by number of harvested plants (CO_331:0000999).
#'  \item \code{dmry}      : Dry matter root yield in t/ha (CO_331:0001000).
#'  \item \code{dmry.aj}   : Dry matter root yield in t/ha adjusted by number of harvested plants (CO_331:0001001).
#'  \item \code{vw.d}      : Vine dry weight in kg/plot (CO_331:0000996).
#'  \item \code{fytha}     : Foliage total yield in t/ha (CO_331:0000684).
#'  \item \code{fytha.aj}  : Foliage total yield in t/ha adjusted by number of harvested plants (CO_331:0001002).
#'  \item \code{dmvy}      : Dry matter vine yield in t/ha (CO_331:0001005).
#'  \item \code{dmvy.aj}   : Dry matter vine yield (t/ha) adjusted by number of harvested plants (CO_331:0001006).
#'  \item \code{bytha}     : Biomass yield in t/ha (CO_331:0001003).
#'  \item \code{bytha.aj}  : Biomass yield in t/ha adjusted by number of harvested plants (CO_331:0001004).
#'  \item \code{dmby}      : Dry matter biomass in t/ha (CO_331:0001007).
#'  \item \code{dmby.aj}   : Dry matter biomass in t/ha adjusted by number of harvested plants (CO_331:0001008).
#'  \item \code{acrw}      : Average commercial root weight in kg/root (CO_331:0000680).
#'  \item \code{ancrw}     : Average non commercial root weight in kg/root (CO_331:0002039).
#'  \item \code{atrw}      : Average total root weight in kg/root (CO_331:0002040).
#'  \item \code{nrpp}      : Number of roots per harvested plant (CO_331:0000230).
#'  \item \code{nrpsp}     : Number of roots per sowed plant (CO_331:0000994).
#'  \item \code{ncrpp}     : Number of commercial roots per harvested plant (CO_331:2000036).
#'  \item \code{ncrpsp}    : Number of commercial roots per sowed plant (CO_331:0000993).
#'  \item \code{ypp}       : Yield per harvested plant in kg (CO_331:0000681).
#'  \item \code{ypsp}      : Yield per sowed plant in kg (CO_331:0000989).
#'  \item \code{vpp}       : Vine weight per harvested plant in kg (CO_331:0000990).
#'  \item \code{vpsp}      : Vine weight per sowed plant in kg (CO_331:0000991).
#'  \item \code{ci}        : Commercial index \% (CO_331:0000682).
#'  \item \code{hi}        : Harvest index \% (CO_331:0000302).
#'  \item \code{shi}       : Harvest sowing index \% (CO_331:0000301).
#'  \item \code{rfr}       : Root foliage ratio \% (CO_331:0001015).
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
  
  traits <- c("nops", "nope", "noph", "vir", "vir1", "vir2", "alt", "alt1", "alt2",
              "vv", "vw", "nopr", "nocr", "nonc", "crw", "ncrw", "scol", "fcol",
              "fcol.cc", "rs", "rf", "damr", "rspr", "alcdam", "wed", "dmf", "dmd",
              "dm", "dmvf", "dmvd", "dmv", "fraw", "suraw", "straw", "coof", "coosu",
              "coost", "coot", "cooap", "prot", "fe", "zn", "ca", "mg", "bc", "bc.cc",
              "tc", "star", "fruc", "gluc", "sucr", "malt", "tnr", "trw", "trw.d",
              "biom", "biom.d", "cytha", "cytha.aj", "rytha", "rytha.aj", "dmry",
              "dmry.aj", "vw.d", "fytha", "fytha.aj", "dmvy", "dmvy.aj", "bytha",
              "bytha.aj", "dmby", "dmby.aj", "acrw", "ancrw", "atrw", "nrpp", "nrpsp",
              "ncrpp", "ncrpsp", "ypp", "ypsp", "vpp", "vpsp", "ci", "hi", "shi", "rfr")
  
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
