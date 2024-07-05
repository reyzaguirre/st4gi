#' Check fieldbook traits names for sweetpotato
#'
#' Check that fieldbook factors and traits names correspond with the names defined in 
#' \url{https://cropontology.org/term/CO_331:ROOT} and in the 
#' Procedures for the evaluation and analysis of sweetpotato trials, ISBN 978-92-9060-522-5.
#' @param dfr The name of the data frame.
#' @param add Additional traits.
#' @details The data frame must use the labels (lower or upper case) listed below.
#' Between parentheses are the CO numbers for variables defined in
#' \url{https://cropontology.org/term/CO_331:ROOT}
#' or COMP numbers defined in \url{https://sweetpotatobase.org/}.
#' 
#' ----------------------------- Plot identifiers -----------------------------
#'  \itemize{
#'  \item \code{plot}      : Plot number.
#'  \item \code{row}       : Row number position of the plot (\code{row_number} is also valid).
#'  \item \code{col}       : Column number position of the plot (\code{col_number} is also valid).
#'  }
#' --------------------------- Classification factors--------------------------
#'  \itemize{
#'  \item \code{loc}       : Locations.
#'  \item \code{year}      : Years.
#'  \item \code{season}    : Seasons. 
#'  \item \code{env}       : Environments.
#'  \item \code{geno}      : Genotypes (\code{accession_name} and \code{cipno} are also valid).
#'  \item \code{type}      : Entry type (\code{clon}, \code{check}, \code{progeny}, \code{parent}).
#'  \item \code{rep}       : Replications (\code{rep_number} is also valid).
#'  \item \code{block}     : Blocks (\code{block_number} is also valid).
#'  \item \code{treat}     : Treatment.
#'  \item \code{harvest}   : Harvest time.
#'  }
#' ------------------------------ Roots and vines -----------------------------
#'  \itemize{
#'  \item \code{nops}      : Number of plants planted per plot (CO_331:0000678).
#'  \item \code{nope}      : Number of plants established per plot (CO_331:0000192).
#'  \item \code{noph}      : Number of plants harvested per plot (CO_331:0000679).
#'  \item \code{vir}       : Virus symptoms estimating 1-9 (CO_331:0000193).
#'  \item \code{vir1}      : Virus symptoms estimating 1-9, first evaluation, 1 month after planting (COMP:0000023).
#'  \item \code{vir2}      : Virus symptoms estimating 1-9, second evaluation, 1 month before harvest (COMP:0000024).
#'  \item \code{alt}       : Alternaria symptoms estimating 1-9 (CO_331:0000198).
#'  \item \code{alt1}      : Alternaria symptoms estimating 1-9, first evaluation, 1 month after planting (COMP:0000026).
#'  \item \code{alt2}      : Alternaria symptoms estimating 1-9, second evaluation, 1 month before harvest (COMP:0000022).
#'  \item \code{vv}        : Vine vigor estimating 1-9 (CO_331:0000197).
#'  \item \code{vw}        : Weight of vines per plot in kg (CO_331:0000227).
#'  \item \code{nopr}      : Number of plants with storage roots per plot (CO_331:0000211).
#'  \item \code{nocr}      : Number of commercial storage roots per plot (CO_331:0000214).
#'  \item \code{nonc}      : Number of non-commercial storage roots per plot (CO_331:0000217).
#'  \item \code{crw}       : Weight of commercial storage roots per plot in kg (CO_331:0000220).
#'  \item \code{ncrw}      : Weight of non-commercial storage roots per plot in kg (CO_331:0000223).
#'  \item \code{scol}      : Storage root predominant skin color estimating 1-9 (CO_331:0000175).
#'  \item \code{fcol}      : Storage root predominant flesh color estimating 1-9 (CO_331:0000178).
#'  \item \code{fcol2}     : Storage root secondary flesh color estimating 0-9 (CO_331:0000179).
#'  \item \code{fcol.cc}   : Root flesh color from RHS color chart, 1-30.
#'  \item \code{rs}        : Storage root size estimating 1-9 (CO_331:0000184).
#'  \item \code{rf}        : Storage root appearance estimating 1-9 (CO_331:0000202).
#'  \item \code{rtshp}     : Storage root shape estimating 1-9 (CO_331:0000181).
#'  \item \code{damr}      : Storage root damage estimating 1-9 (CO_331:0000206).
#'  \item \code{rspr}      : Sprouting ability estimating 1-9 (CO_331:0000277).
#'  \item \code{alcdam}    : Alcidodes sp. damage estimating 1-9 (CO_331:0000806).
#'  \item \code{wed}       : Reaction to sweet potato weevil estimating 1-9 (CO_331:0000207).
#'  \item \code{stspwv}    : Reaction to striped sweet potato weevil estimating 1-9 (CO_331:0000720).
#'  \item \code{milldam}   : Millipede damage estimating 1-9 (CO_331:0000805).
#'  }
#' --------------------------- Dry matter assesment ---------------------------
#'  \itemize{
#'  \item \code{dmf}       : Fresh weight of storage root samples for dry matter assessment measuring g (CO_331:0000243). 
#'  \item \code{dmd}       : Dry weight of storage root samples for dry matter assessment measuring g (CO_331:0000247).
#'  \item \code{dm}        : Storage root dry matter content computing percent (CO_331:0000297).
#'  \item \code{dmvf}      : Fresh weight of vines samples for dry matter assessment measuring g (CO_331:0000251).
#'  \item \code{dmvd}      : Dry weight of vines samples for dry matter assessment measuring g (CO_331:0000255).
#'  \item \code{dmv}       : Vine dry matter content computing percent (CO_331:0001014).
#'  }
#' ---------------------- Raw and cooked roots attributes ---------------------
#'  \itemize{
#'  \item \code{fraw}      : Fiber content in raw storage roots estimating 1-9 (CO_331:0001038).
#'  \item \code{suraw}     : Content of total sugars in raw storage roots estimating 1-9 (CO_331:0000766).
#'  \item \code{straw}     : Content of starch in raw storage roots estimating 1-9 (CO_331:0001055).
#'  \item \code{coof}      : Fiber content in boiled storage roots estimating 1-9 (CO_331:0000259).
#'  \item \code{coosu}     : Storage root sweetness in cooked samples estimating 1-9 (CO_331:0000261).
#'  \item \code{coost}     : Storage root texture in boiled samples estimating 1-9 (CO_331:0000265).
#'  \item \code{coot}      : Boiled storage root flesh taste estimating 1-9 (CO_331:0000269).
#'  \item \code{cooap}     : Boiled storage root overall appearance estimating 1-9 (CO_331:0000273).
#'  }
#' --------------------------- Nutrients evaluations --------------------------
#'  \itemize{
#'  \item \code{prot}      : Content of protein in dry weight basis in raw storage roots measuring percentage (CO_331:0000278).
#'  \item \code{fe}        : Content of iron in dry weight basis measuring mg per 100 g (CO_331:0001016).
#'  \item \code{zn}        : Content of zinc in dry weight basis measuring mg per 100 g (CO_331:0001017).
#'  \item \code{ca}        : Content of calcium in dry weight basis measuring mg per 100 g (CO_331:0001029).
#'  \item \code{mg}        : Content of magnesium in dry weight basis measuring mg per 100 g (CO_331:0001030).
#'  \item \code{bc}        : Content of beta-carotene in dry weight basis in raw storage roots measuring mg per 100 g (CO_331:0000289).
#'  \item \code{bc.cc}     : Content of beta-carotene in fresh weight basis in raw storage roots estimating from RHS color chart in mg per 100 g	(CO_331:0001023).
#'  \item \code{tc}        : Content of total carotenoids in dry weight basis in raw storage roots measuring mg per 100 g (CO_331:0000290).
#'  \item \code{star}      : Content of starch in dry weight basis in raw storage roots measuring percentage (CO_331:0000291).
#'  \item \code{star.b}    : Content of starch in dry weight basis in boiled storage roots measuring percentage	(CO_331:2000027).
#'  \item \code{fruc}      : Content of fructose in dry weight basis in raw storage roots measuring g per 100 g (CO_331:0001045).
#'  \item \code{fruc.b}    : Content of fructose in dry weight basis in boiled storage roots measuring g per 100 g (CO_331:0001046).
#'  \item \code{gluc}      : Content of glucose in dry weight basis in raw storage roots measuring g per 100 g (CO_331:0001047).
#'  \item \code{gluc.b}    : Content of glucose in dry weight basis in boiled storage roots measuring g per 100 g (CO_331:0001048).
#'  \item \code{sucr}      : Content of sucrose in dry weight basis in raw storage roots measuring g per 100 g (CO_331:0001049).
#'  \item \code{sucr.b}    : Content of sucrose in dry weight basis in boiled storage roots measuring g per 100 g (CO_331:0001050).
#'  \item \code{malt}      : Content of maltose in dry weight basis in raw storage roots measuring g per 100 g (CO_331:0001051).
#'  \item \code{malt.b}    : Content of maltose in dry weight basis in boiled storage roots measuring g per 100 g (CO_331:0001052).
#'  }
#' ----------------------------- Calculated traits ----------------------------
#'  \itemize{
#'  \item \code{tnr}       : Total number of storage roots per plot (CO_331:0000233).
#'  \item \code{trw}       : Total storage root weight per plot in kg (CO_331:0000237).
#'  \item \code{trw.d}     : Total storage root weight in dry weight basis per plot in kg	(CO_331:0000995).
#'  \item \code{biom}      : Biomass yield per plot in kg (CO_331:0000683).
#'  \item \code{biom.d}    : Biomass yield in dry weight basis per plot in kg (CO_331:0000997).
#'  \item \code{cytha}     : Commercial storage root yield in tons per ha (CO_331:0000809).
#'  \item \code{cytha.aj}  : Commercial storage root yield adjusted by number of harvested plants per plot in tons per ha (CO_331:0000998).
#'  \item \code{rytha}     : Total storage roots yield in tons per ha (CO_331:0000296).
#'  \item \code{rytha.aj}  : Total storage roots yield adjusted by number of harvested plants per plot in tons per ha	(CO_331:0000999).
#'  \item \code{dmry}      : Storage root dry matter yield in tons per ha (CO_331:0001000).
#'  \item \code{dmry.aj}   : Storage root dry matter yield adjusted by number of harvested plants per plot in tons per ha	(CO_331:0001001).
#'  \item \code{vw.d}      : Dry weight of vines per plot in kg	(CO_331:0000996).
#'  \item \code{fytha}     : Total foliage yield in tons per ha (CO_331:0000684).
#'  \item \code{fytha.aj}  : Total foliage yield adjusted by number of harvested plants per plot in tons per ha (CO_331:0001002).
#'  \item \code{dmvy}      : Vine dry matter yield in tons per ha (CO_331:0001005).
#'  \item \code{dmvy.aj}   : Vine dry matter yield adjusted by number of harvested plants per plot in tons per ha	(CO_331:0001006).
#'  \item \code{bytha}     : Biomass yield in tons per ha	(CO_331:0001003).
#'  \item \code{bytha.aj}  : Biomass yield adjusted by number of harvested plants per plot in tons per ha (CO_331:0001004).
#'  \item \code{dmby}      : Biomass dry matter yield in tons per ha (CO_331:0001007).
#'  \item \code{dmby.aj}   : Biomass dry matter yield adjusted by number of harvested plants per plot in tons per ha (CO_331:0001008).
#'  \item \code{acrw}      : Average commercial storage root weight in kg (CO_331:0000680).
#'  \item \code{ancrw}     : Average non-commercial storage root weight in kg (CO_331:0002039). 
#'  \item \code{atrw}      : Average total storage root weight in kg (CO_331:0002040).
#'  \item \code{nrpp}      : Number of storage roots per plant harvested (CO_331:0000230).
#'  \item \code{nrpsp}     : Number of storage roots per plant planted (CO_331:0000994).
#'  \item \code{ncrpp}     : Number of commercial storage roots per plant harvested	(CO_331:0006024).
#'  \item \code{ncrpsp}    : Number of commercial storage roots per plant planted	(CO_331:0000993).
#'  \item \code{ypp}       : Storage roots yield per plant harvested in kg (CO_331:0000681).
#'  \item \code{ypsp}      : Storage roots yield per plant planted in kg (CO_331:0000989).
#'  \item \code{vpp}       : Vine weight per plant harvested in kg (CO_331:0000990).
#'  \item \code{vpsp}      : Vine weight per sowed plant in kg (CO_331:0000991).
#'  \item \code{rtyldpct}  : Yield as percent of check (CO_331:0000792).
#'  \item \code{ci}        : Commercial index in percentage (CO_331:0000682).
#'  \item \code{hi}        : Harvest index in percentage (CO_331:0000302).
#'  \item \code{shi}       : Survival index in percentage (CO_331:0000301).
#'  \item \code{rfr}       : Root foliage ratio in percentage	(CO_331:0001015).
#'  }
#' @return It returns:
#' \itemize{
#' \item The fieldbook data frame with all column names in lowercase and
#' with some possible modifications in the names.
#' \item A list of warnings for all the column names that have been changed.
#' \item A list of warnings for all the column names not recognized.
#' }
#' @author Raul Eyzaguirre.
#' @examples
#' check.names.sp(pjpz09)
#' @export

check.names.sp <- function(dfr, add = NULL) {
  
  # Valid names for factors
  
  plot.id <- c('observationunit_name', "plot", "row", 'row_number', "col", 'col_number')
  
  factors <- c("loc", "year", "season", "env",
               "geno", 'accession_name', "cipno", 'type',
               "rep", 'rep_number', "block", 'block_number',
               "treat", "harvest")
  
  # Valid names for traits
  
  traits <- c(spont$Label, 'fcol.cc')
  
  # Valid names for factors and traits
  
  colnames.valid <- c(plot.id, factors, traits, tolower(add))
  
  # Factors and traits in field book
    
  colnames.fb <- colnames(dfr)
  
  # check.list.1: mark names not valid
  
  check.list.1 <- !(tolower(colnames.fb) %in% colnames.valid)
  
  # check.list.2: mark valid names but upper case to be converted to lower case
  
  colnames.fb.valid <- colnames.fb[!check.list.1]
  
  check.list.2 <- !(colnames.fb.valid %in% colnames.valid)
  
  # Convert all fieldbook names to lower case
  
  colnames(dfr) <- tolower(colnames(dfr))
    
  # Solve synonyms for factors
  
  change.names.f <- NULL
  
  old.names.f <- c("rep_number", "block_number", "row_number", "col_number", "accession_name", "cipno")
  new.names.f <- c("rep",        "block",        "row",        "col",        "geno",           "geno")
  
  for (i in 1:length(old.names.f)) {
    
    if (exists(old.names.f[i], dfr) & !exists(new.names.f[i], dfr)) {
      change.names.f <- c(change.names.f, old.names.f[i])
      colnames(dfr)[colnames(dfr) == old.names.f[i]] <- new.names.f[i]
    }
    
  }  
  
  # Warnings
  
  if (!is.null(change.names.f)) {
    change.names.list <- old.names.f %in% change.names.f
    warning("Factors' names ", list(old.names.f[change.names.list]), " changed to ", list(new.names.f[change.names.list]), call. = FALSE)
  }
  
  if (max(check.list.1) == 1)
    warning("Some columns with invalid names: ", list(colnames.fb[check.list.1]), call. = FALSE)
  
  if (max(check.list.2) == 1)
    warning("Some labels converted to lower case: ", list(colnames.fb.valid[check.list.2]), call. = FALSE)
  
  # Return
  
  dfr
  
}
