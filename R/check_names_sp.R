#' Check fieldbook traits names for sweetpotato
#'
#' Check that fieldbook traits names correspond with the names defined in 
#' Procedures for the evaluation and analysis of sweetpotato trials, ISBN 978-92-9060-522-5.
#' @param dfr The name of the data frame.
#' @param add Additional traits.
#' @details The data frame must use the labels (lower or upper case) listed below.
#' Between parentheses the CO numbers for variable/trait according to
#' https://github.com/Planteome/CO_331-sweetpotato-traits/blob/master/sweetpotato-trait-ontology_withcropname.obo
#' or COMP (compound trait) number according to https://sweetpotatobase.org/
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
#' -------------------- Roots and vines --------------------
#'  \itemize{
#'  \item \code{nops}      : Number of plants sowed per plot (CO_331:0000678/CO_331:0000303).
#'  \item \code{nope}      : Number of plants established per plot (CO_331:0000192/CO_331:0000189).
#'  \item \code{noph}      : Number of plants harvested per plot (CO_331:0000679/CO_331:0000304).
#'  \item \code{vir}       : Virus symptoms, 1-9 scale (CO_331:0000193/CO_331:0000094).
#'  \item \code{vir1}      : Virus symptoms, 1-9 scale, first evaluation, 1 month after planting (COMP:0000023).
#'  \item \code{vir2}      : Virus symptoms, 1-9 scale, second evaluation, 1 month before harvest (COMP:0000024).
#'  \item \code{alt}       : Alternaria symptoms, 1-9 scale (CO_331:0000198/CO_331:0000091).
#'  \item \code{alt1}      : Alternaria symptoms, 1-9 scale, first evaluation, 1 month after planting (COMP:0000026).
#'  \item \code{alt2}      : Alternaria symptoms, 1-9 scale, second evaluation, 1 month before harvest (COMP:0000022).
#'  \item \code{vv}        : Vine vigor, 1-9 scale (CO_331:0000197/CO_331:0000194).
#'  \item \code{vw}        : Vine weight in kg/plot (CO_331:0000227/CO_331:0000224).
#'  \item \code{nopr}      : Number of plants with roots per plot (CO_331:0000211/CO_331:0000208).
#'  \item \code{nocr}      : Number of commercial roots per plot (CO_331:0000214/CO_331:0000212).
#'  \item \code{nonc}      : Number of non commercial roots per plot (CO_331:0000217/CO_331:0000215).
#'  \item \code{crw}       : Commercial root weight in kg/plot (CO_331:0000220/CO_331:0000218).
#'  \item \code{ncrw}      : Non commercial root weight in kg/plot (CO_331:0000223/CO_331:0000221).
#'  \item \code{scol}      : Storage root skin color, 1-9 scale, (CO_331:0000175/CO_331:0000049).
#'  \item \code{fcol}      : Storage root predominant flesh color, 1-9 scale (CO_331:0000178/CO_331:0000058).
#'  \item \code{fcol2}     : Storage root secondary flesh color, 1-9 scale (CO_331:0000179/CO_331:0000061).
#'  \item \code{fcol.cc}   : Root flesh color using RHS color charts, 1-30 scale.
#'  \item \code{rs}        : Root size, 1-9 scale, (CO_331:0000184/CO_331:0000076).
#'  \item \code{rf}        : Root form, 1-9 scale, (CO_331:0000202/CO_331:0000199).
#'  \item \code{rtshp}     : Root shape, 1-9 scale, (CO_331:0000181/CO_331:0000067).
#'  \item \code{damr}      : Root damage, 1-9 scale, (CO_331:0000206/CO_331:0000203).
#'  \item \code{rspr}      : Root sprouting, 1-9 scale, (CO_331:0000277/CO_331:0000274).
#'  \item \code{alcdam}    : Alcidodes sp. damage, 1-9 scale (CO_331:0000806/CO_331:0000439).
#'  \item \code{wed}       : Weevil damage, 1-9 scale (CO_331:0000207/CO_331:0000088).
#'  \item \code{stspwv}    : Reaction to striped weevil, 1-9 scale (CO_331:0000720/CO_331:0000356).
#'  \item \code{milldam}   : Millipede damage, 1-9 scale (CO_331:0000805/CO_331:0000438).
#'  }
#' -------------------- Dry matter assesment --------------------
#'  \itemize{
#'  \item \code{dmf}       : Fresh weight of roots for dry matter assessment in g (CO_331:0000243/CO_331:0000240). 
#'  \item \code{dmd}       : Dry weight of roots for dry matter assessment in g (CO_331:0000247/CO_331:0000244).
#'  \item \code{dm}        : Storage root dry matter content \% (CO_331:0000297/CO_331:0000142).
#'  \item \code{dmvf}      : Fresh weight vines for dry matter assessment in g (CO_331:0000251/CO_331:0000248).
#'  \item \code{dmvd}      : Dry weight vines for dry matter assessment in g (CO_331:0000255/CO_331:0000252).
#'  \item \code{dmv}       : Vines dry matter content \%.
#'  }
#' -------------------- Raw and cooked roots attributes --------------------
#'  \itemize{
#'  \item \code{fraw}      : Raw root fiber, 1-9 scale.
#'  \item \code{suraw}     : Raw root sugar, 1-9 scale (CO_331:0000766/CO_331:0000402).
#'  \item \code{straw}     : Raw root starch, 1-9 scale.
#'  \item \code{coof}      : Cooked root fiber, 1-9 scale (CO_331:0000259/CO_331:0000256).
#'  \item \code{coosu}     : Cooked root sugars, 1-9 scale (CO_331:0000261/CO_331:0000139).
#'  \item \code{coost}     : Cooked root starch, 1-9 scale (CO_331:0000265/CO_331:0000136).
#'  \item \code{coot}      : Cooked root taste, 1-9 scale (CO_331:0000269/CO_331:0000266).
#'  \item \code{cooap}     : Cooked root appearance, 1-9 scale (CO_331:0000273/CO_331:0000270).
#'  }
#' -------------------- Nutrients evaluations --------------------
#'  \itemize{
#'  \item \code{prot}      : Protein, g/100g raw dry weight (CO_331:0000278/CO_331:0000100).
#'  \item \code{fe}        : Iron, mg/100g raw dry weight measured by NIRS (CO_331:0000279/CO_331:0000103).
#'  \item \code{zn}        : Zinc, mg/100g raw dry weight measured by NIRS (CO_331:0000280/CO_331:0000106).
#'  \item \code{ca}        : Calcium, mg/100g raw dry weight measured by NIRS (CO_331:0000284/CO_331:0000281).
#'  \item \code{mg}        : Magnesium, mg/100g raw dry weight measured by NIRS (CO_331:0000288/CO_331:0000285).
#'  \item \code{bc}        : Beta-carotene, mg/100g raw dry weight measured by NIRS (CO_331:0000289/CO_331:0000109).
#'  \item \code{bc.cc}     : Beta-carotene with RHS color charts, mg/100g raw fresh weight.
#'  \item \code{tc}        : Total carotenoids, mg/100g raw dry weight (CO_331:0000290/CO_331:0000112).
#'  \item \code{star}      : Starch, g/100g raw dry weight (CO_331:0000291/CO_331:0000115).
#'  \item \code{fruc}      : Fructose, g/100g raw dry weight (CO_331:0000292/CO_331:0000118).
#'  \item \code{gluc}      : Glucose, g/100g raw dry weight (CO_331:0000293/CO_331:0000121).
#'  \item \code{sucr}      : Sucrose, g/100g raw dry weight (CO_331:0000294/CO_331:0000124).
#'  \item \code{malt}      : Maltose, g/100g raw dry weight (CO_331:0000295/CO_331:0000127).
#'  }
#' -------------------- Calculated traits --------------------
#'  \itemize{
#'  \item \code{tnr}       : Total number of roots per plot (CO_331:0000233/CO_331:0000079).
#'  \item \code{trw}       : Total root weight in kg/plot (CO_331:0000237/O_331:0000234).
#'  \item \code{trw.d}     : Total root dry weight in kg/plot.
#'  \item \code{biom}      : Biomass yield in kg/plot (CO_331:0000683/CO_331:0000311).
#'  \item \code{biom.d}    : Biomass dry yield in kg/plot.
#'  \item \code{cytha}     : Commercial root yield in t/ha.
#'  \item \code{cytha.aj}  : Commercial root yield in t/ha adjusted by number of harvested plants.
#'  \item \code{rytha}     : Total root yield in t/ha (CO_331:0000296/CO_331:0000082).
#'  \item \code{rytha.aj}  : Total root yield in t/ha adjusted by number of harvested plants.
#'  \item \code{dmry}      : Dry matter root yield in t/ha.
#'  \item \code{dmry.aj}   : Dry matter root yield in t/ha adjusted by number of harvested plants.
#'  \item \code{vw.d}      : Vine dry weight in kg/plot.
#'  \item \code{fytha}     : Foliage total yield in t/ha (CO_331:0000684/CO_331:0000312).
#'  \item \code{fytha.aj}  : Foliage total yield in t/ha adjusted by number of harvested plants.
#'  \item \code{dmvy}      : Dry matter vine yield in t/ha.
#'  \item \code{dmvy.aj}   : Dry matter vine yield in t/ha adjusted by number of harvested plants.
#'  \item \code{bytha}     : Biomass yield in t/ha.
#'  \item \code{bytha.aj}  : Biomass yield in t/ha adjusted by number of harvested plants.
#'  \item \code{dmby}      : Dry matter biomass in t/ha.
#'  \item \code{dmby.aj}   : Dry matter biomass in t/ha adjusted by number of harvested plants.
#'  \item \code{acrw}      : Average commercial root weight in kg/root (CO_331:0000680/CO_331:0000308).
#'  \item \code{ancrw}     : Average non commercial root weight in kg/root.
#'  \item \code{atrw}      : Average total root weight in kg/root.
#'  \item \code{nrpp}      : Number of roots per harvested plant (CO_331:0000230/CO_331:0000079).
#'  \item \code{nrpsp}     : Number of roots per sowed plant.
#'  \item \code{ncrpp}     : Number of commercial roots per harvested plant.
#'  \item \code{ncrpsp}    : Number of commercial roots per sowed plant.
#'  \item \code{ypp}       : Yield per harvested plant in kg (CO_331:0000681/CO_331:0000309).
#'  \item \code{ypsp}      : Yield per sowed plant in kg.
#'  \item \code{vpp}       : Vine weight per harvested plant in kg.
#'  \item \code{vpsp}      : Vine weight per sowed plant in kg.
#'  \item \code{rtyldpct}  : Yield as percentage of check (CO_331:0000792/CO_331:0000425).
#'  \item \code{ci}        : Commercial index \% (CO_331:0000682/CO_331:0000310).
#'  \item \code{hi}        : Harvest index \% (CO_331:0000302/CO_331:0000085).
#'  \item \code{shi}       : Harvest sowing index \% (CO_331:0000301/CO_331:0000298).
#'  \item \code{rfr}       : Root foliage ratio \%.
#'  }
#' @return It returns a data frame with all traits names in lower case, and a list of the
#' traits with names not included in the list shown above.
#' @author Raul Eyzaguirre.
#' @examples
#' check.names.sp(pjpz09)
#' @export

check.names.sp <- function(dfr, add = NULL) {
  
  plot.id <- c("plot", "row", "col")
  
  factors <- c("l", "loc", "y", "year", "s", "season", "e", "env", "g", "geno",
               "cipno", "r", "rep", "b", "block", "treat", "harvest")
  
  traits <- c("nops", "nope", "noph", "vir", "vir1", "vir2", "alt", "alt1", "alt2",
              "vv", "vw", "nopr", "nocr", "nonc", "crw", "ncrw", "scol", "fcol",
              "fcol2", "fcol.cc", "rs", "rf", "rtshp", "damr", "rspr", "alcdam",
              "wed", "stspwv", "milldam", "dmf", "dmd", "dm", "dmvf", "dmvd",
              "dmv", "fraw", "suraw", "straw", "coof", "coosu", "coost", "coot",
              "cooap", "prot", "fe", "zn", "ca", "mg", "bc", "bc.cc", "tc", "star",
              "fruc", "gluc", "sucr", "malt", "tnr", "trw", "trw.d", "biom",
              "biom.d", "cytha", "cytha.aj", "rytha", "rytha.aj", "dmry",
              "dmry.aj", "vw.d", "fytha", "fytha.aj", "dmvy", "dmvy.aj", "bytha",
              "bytha.aj", "dmby", "dmby.aj", "acrw", "ancrw", "atrw", "nrpp",
              "nrpsp", "ncrpp", "ncrpsp", "ypp", "ypsp", "vpp", "vpsp", "rtyldpct",
              "ci", "hi", "shi", "rfr")
  
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
