#' Check fieldbook traits names for potato
#'
#' Check that fieldbook traits names correspond with the names defined in 
#' \url{http://www.cropontology.org/ontology/CO_330/Potato}
#' @param dfr The name of the data frame.
#' @param add Additional traits.
#' @details The data frame must use the following labels (lower or upper case):
#' 
#'  -------------------- Plot identifiers --------------------
#'  \itemize{
#'  \item \code{plot}        : Plot number
#'  \item \code{row}         : Row number
#'  \item \code{col}         : Column number
#'  }
#' -------------------- Classification factors --------------------
#'  \itemize{
#'  \item \code{l}           : Locations (\code{loc} is also valid)
#'  \item \code{y}           : Years
#'  \item \code{s}           : Seasons
#'  \item \code{e}           : Environments (\code{env} is also valid)
#'  \item \code{g}           : Genotypes (\code{geno} is also valid)
#'  \item \code{name}        : Names for genotypes
#'  \item \code{r}           : Replications (\code{rep} is also valid)
#'  \item \code{b}           : Blocks (\code{block} is also valid)
#'  \item \code{treatment}   : Irrigation treatments
#'  \item \code{type}        : Type: early or late
#'  \item \code{instn}       : Code number
#'  }
#' -------------------- Traits N1 group --------------------
#'  \itemize{
#'  \item \code{ntp}         : Number of tubers planted
#'  \item \code{npe}         : Number of plants emerged
#'  \item \code{nph}         : Number of plants harvested
#'  \item \code{ppe}         : Proportion of plants emerged
#'  \item \code{pph}         : Proportion of plants harvested
#'  }
#' -------------------- Traits N2 group --------------------
#'  \itemize{
#'  \item \code{snpp}        : Stem number per plant
#'  \item \code{nipp}        : Number of inflorescences per plant
#'  \item \code{nfwp}        : Number of flowers per main inflorescence
#'  \item \code{nlpp}        : Number of leaves per plant
#'  \item \code{Num_Stolon}  : Number of stolons
#'  \item \code{Leng_Stolon} : Lenght of stolons
#'  }
#' -------------------- Traits N3 group --------------------
#'  \itemize{
#'  \item \code{tntp}        : Total number of tuber per plot
#'  \item \code{tntpl}       : Total number of tuber per plant
#'  \item \code{nmtp}        : Number of marketable tubers per plot
#'  \item \code{nmtpl}       : Number of marketable tubers per plant
#'  \item \code{nnomtp}      : Number of non-marketable tubers per plot
#'  \item \code{nmtci}       : Number of marketable tubers category I per plot
#'  \item \code{nmtcii}      : Number of marketable tubers category II per plot
#'  }
#' -------------------- Traits N4 group --------------------
#'  \itemize{
#'  \item \code{ttwp}        : Total tuber weight per plot
#'  \item \code{ttwpl}       : Total tuber weight per plant
#'  \item \code{mtwp}        : Marketable tuber weight per plot
#'  \item \code{mtwpl}       : Marketable tuber weight per plant
#'  \item \code{nomtwp}      : Non-marketable tuber weight per plot
#'  \item \code{mtwci}       : Marketable tuber weight category I per plot
#'  \item \code{mtwcii}      : Marketable tuber weight category II per plot
#'  }
#' -------------------- Traits N5 group--------------------
#'  \itemize{
#'  \item \code{ttya}        : Total tuber yield adjusted
#'  \item \code{ttyna}       : Total tuber yield no adjusted
#'  \item \code{mtya}        : Marketable tuber yield adjusted
#'  \item \code{mtyna}       : Marketable tuber yiedl no adjusted
#'  \item \code{atw}         : Average of tuber weight in grams (\code{mwt} is also valid)
#'  \item \code{atmw}        : Average of marketable tuber weight in grams (\code{mwmt} is also valid)
#'  }
#' -------------------- Traits N6 group--------------------
#'  \itemize{
#'  \item \code{stlfw}       : Stolon fresh weight per plant
#'  \item \code{sfw}         : Stem fresh weight per plant (\code{stfw} is also valid)
#'  \item \code{lfw}         : Leaf fresh weight per plant
#'  \item \code{rfw}         : Root fresh weight per plant
#'  \item \code{tfw}         : Tuber fresh weight per plant
#'  \item \code{tbfw}        : Total biomass fresh weight per plant
#'  \item \code{hi_fw}       : Harvest index fresh weight
#'  \item \code{fwts}        : Fresh weight of tuber sample
#'  }
#' -------------------- Traits N7 group--------------------
#'  \itemize{
#'  \item \code{stldw}       : Stolon dry weight per plant
#'  \item \code{sdw}         : Stem dry weight per plant (\code{stdw} is also valid)
#'  \item \code{ldw}         : Leaf dry weight per plant
#'  \item \code{rdw}         : Root dry weight per plant
#'  \item \code{tdw}         : Tuber dry weight per plant
#'  \item \code{tbdw}        : Total biomass dry weight per plant
#'  \item \code{hi_dw}       : Harvest index dry weight
#'  \item \code{dwts}        : Dry weight of tuber sample
#'  }
#' -------------------- Traits N8 group--------------------
#'  \itemize{
#'  \item \code{ldmcp}       : Leaf dry matter content per plot
#'  \item \code{sdmcp}       : Stem dry matter content per plot
#'  \item \code{rdmcp}       : Root dry matter content per plot
#'  \item \code{tdmcp}       : Tuber dry matter content per plot
#'  \item \code{pdm}         : Tuber dry matter content (\code{dm} is also valid)
#'  }
#' -------------------- Traits N9 group--------------------
#'  \itemize{
#'  \item \code{twa}         : Tuber weight in air
#'  \item \code{tww}         : Tuber weight in water
#'  \item \code{rsdw}        : Root system dry weight per plant
#'  \item \code{rd}          : Root density
#'  \item \code{rl}          : Root length
#'  \item \code{sg}          : Tuber specific gravity
#'  \item \code{dsi}         : Drought susceptibility index
#'  \item \code{dti}         : Drought tolerance index
#'  }
#' -------------------- Traits N10 group--------------------
#'  \itemize{
#'  \item \code{fedw}        : Tuber iron concentration in dry weight basis
#'  \item \code{fefw}        : Tuber iron concentration in fresh weight basis
#'  \item \code{zndw}        : Tuber zinc concentration in dry weight basis
#'  \item \code{znfw}        : Tuber zinc concentration in fresh weight basis
#'  \item \code{antho_dw}    : Tuber anthocyanin concentration in dry weight basis
#'  \item \code{antho_fw}    : Tuber anthocyanin concentration in fresh weight basis
#'  \item \code{aah_dw}      : Tuber hydrophilic antioxidant activity in dry weight basis
#'  \item \code{aah_fw}      : Tuber hydrophilic antioxidant activity in fresh weight basis
#'  \item \code{aal_dw}      : Tuber lipophilic antioxidant activity in dry weight basis
#'  \item \code{aal_fw}      : Tuber lipophilic antioxidant activity in fresh weight basis
#'  \item \code{asc_dw}      : Tuber ascorbic acid concentration in dry weight basis
#'  \item \code{asc_fw}      : Tuber ascorbic acid concentration in fresh weight basis
#'  \item \code{pro}         : Tuber protein content (\code{protein} is also valid)
#'  \item \code{star}        : Tuber starch content
#'  \item \code{fruc}        : Tuber fructose content
#'  \item \code{gluc}        : Tuber glucose content
#'  \item \code{sucr}        : Tuber sucrose content
#'  \item \code{malt}        : Tuber maltose content
#'  \item \code{fiber}       : Tuber fiber content
#'  }
#' -------------------- Traits N11 group--------------------
#'  \itemize{
#'  \item \code{plant_unif}  : Plant uniformity
#'  \item \code{plant_vigor} : Plant vigor
#'  \item \code{se}          : Senescence
#'  \item \code{tuber_apper} : Tuber appearance
#'  \item \code{tub_unif}    : Tuber uniformity
#'  \item \code{tub_size}    : Tuber size
#'  \item \code{op}          : Osmotic potential
#'  \item \code{tlwp}        : Total leaf water potential
#'  \item \code{ltp}         : Leaf turgor potential
#'  \item \code{dap_chl}     : Chlorophyll relative content - Days after planting # removed
#'  \item \code{dap_tc}      : Canopy temperature - Days after planting           # removed
#'  }
#' -------------------- Traits N12 group--------------------
#'  \itemize{
#'  \item \code{leaflet_twi} : Leaflet turgid weight Evaluation i (i = 1,..., 5)
#'  \item \code{insnppi}     : Increase stem number per plant i (i = 1,..., 5)
#'  \item \code{insdi}       : Increase stem diameter i (i = 1,...,5)
#'  \item \code{inrwci}      : Increase relative water content i (i = 1,..., 5)
#'  }
#' -------------------- Traits N13 group--------------------
#'  \itemize{
#'  \item \code{pw_evi}      : Plant wilting Evaluation i (i = 1,..., 5)
#'  \item \code{inplahei}    : Increase plant height i (i = 1,..., 5)
#'  \item \code{snppi}       : Stem number per plant Evaluation i (i = 1,..., 5)
#'  \item \code{leaflet_fwi} : Leaflet fresh weight Evaluation i (i = 1,..., 5)
#'  \item \code{leaflet_dwi} : Leaflet dry weight Evaluation i (i = 1,..., 5)
#'  }
#' -------------------- Traits N14 group--------------------
#'  \itemize{
#'  \item \code{chci}        : Chlorophyll content sample i (i = 1,..., 5)
#'  \item \code{inchci}      : Increase chlorophyll content i (i = 1,..., 5)
#'  \item \code{leafsdi}     : Leaf stomata density sample i (i = 1,..., 5)
#'  \item \code{av_leafsd}   : Average leaf stomata density
#'  }
#' -------------------- Traits N15 group--------------------
#'  \itemize{
#'  \item \code{plahe_evi}   : Plant height Evaluation i (i = 1,..., 5)
#'  \item \code{sd_evi}      : Stem diameter Evaluation i (i = 1,..., 5)
#'  \item \code{cc_evi}      : Canopy cover Evaluation i (i = 1,..., 5)
#'  \item \code{chlspad_evi} : Chlorophyll content index Evaluation i (i = 1,..., 5)
#'  \item \code{cr_evi}      : Canopy reflectance Evaluation i (i = 1,..., 5)
#'  \item \code{lfa_evi}     : Leaflet area Evaluation i (i = 1,..., 5)
#'  \item \code{rwc_evi}     : Relative water content Evaluation i (i = 1,..., 5)
#'  \item \code{sla_evi}     : Specific leaf area Evaluation i (i = 1,..., 5)
#'  }
#' -------------------- Traits N16 group--------------------
#'  \itemize{
#'  \item \code{plahe_slp}   : Plant height slope
#'  \item \code{sd_slp}      : Stem diameter slope
#'  \item \code{cc_slp}      : Canopy cover slope
#'  \item \code{chlspad_slp} : Chlorophyll content index slope
#'  \item \code{cr_slp}      : Canopy reflectance slope
#'  \item \code{lfa_slp}     : Leaflet area slope
#'  \item \code{sla_slp}     : Specific leaf area slope         
#'  }
#' -------------------- Traits N17 group--------------------
#'  \itemize{
#'  \item \code{lbi}         : Late Blight severity evaluation i (i = 1, 2,... 8)
#'  \item \code{audpc}       : Late blight AUDPC
#'  \item \code{raudpc}      : Late blight relative AUDPC
#'  \item \code{sraudpc}     : Late blight susceptibility
#'  }
#' @return It returns a data frame with all traits names in lower case, and a list of the
#' traits with names not included in the list shown above.
#' @author Raul Eyzaguirre, Johan Ninanya. 
#' @examples
#' check.names.pt(potatoyield)
#' @export

check.names.pt <- function(dfr, add = NULL) {
  
  plot.id <- c("plot", "row", "col")
  
  factors <- c("l", "loc", "y", "s", "g", "e", "env", "geno", "name", "r", "rep",
               "b", "block", "treatment", "type", "instn")
  
  # Repeated measures traits (i = 1...5)
  
  rm.traits <- c("leaflet_tw", "insnpp", "insd", "inrwc", "pw_ev", "inplahe",
                 "snpp", "leaflet_fw", "leaflet_dw", "chc", "inchc", "leafsd",
                 "plahe_ev", "sd_ev", "cc_ev", "chlspad_ev", "cr_ev", "lfa_ev",
                 "rwc_ev", "sla_ev")
  
  rm.traits <- c(rm.traits, sapply(rm.traits, function(x) paste0(x, 1:5)))
  
  # Late blight evaluations
  
  lb.ev <- paste0("lb", 1:8)
  
  # List of traits

  traits <- c("ntp", "npe", "nph", "ppe", "pph",
              "snpp", "nipp", "nfwp", "nlpp", "num_stolon", "leng_stolon",
              "tntp", "tntpl", "nmtp", "nmtpl", "nnomtp", "nmtci", "nmtcii",
              "ttwp", "ttwpl", "mtwp", "mtwpl", "nomtwp", "mtwci", "mtwcii",
              "ttya", "ttyna", "mtya", "mtyna", "atw", "mwt", "atmw", "mwmt",
              "stlfw", "sfw", "stfw", "lfw", "rfw", "tfw", "tbfw", "hi_fw", "fwts",
              "stldw", "sdw", "stdw", "ldw", "rdw", "tdw", "tbdw", "hi_dw", "dwts",
              "ldmcp", "sdmcp", "rdmcp", "tdmcp", "pdm", "dm",
              "twa", "tww", "rsdw", "rd", "rl", "sg", "dsi", "dti",
              "fedw", "fefw", "zndw", "znfw", "antho_dw", "antho_fw",
              "aah_dw", "aah_fw", "aal_dw", "aal_fw", "asc_dw", "asc_fw",
              "pro", "protein", "star", "fruc", "gluc", "sucr", "malt", "fiber",
              "plant_unif", "plant_vigor", "se", "op", "tlwp", "ltp", "dap_chl", "dap_tc",
              "tuber_apper", "tub_unif", "tub_size", "av_leafsd", "plahe_slp",
              "sd_slp", "cc_slp", "chlspad_slp", "cr_slp", "lfa_slp", "rwc_slp", "sla_slp",
              rm.traits, lb.ev, "audpc", "raudpc", "saudpc")
  
  colnames.valid <- c(plot.id, factors, traits, tolower(add))
  
  colnames.list <- colnames(dfr)
  
  check.list.1 <- !(tolower(colnames.list) %in% colnames.valid) # which are not valid
  temp <- colnames.list[!check.list.1]                          # list of valid names
  check.list.2 <- !(temp %in% colnames.valid)                   # which are valid but lower case
  
  colnames(dfr) <- tolower(colnames(dfr))
  
  # Solve synonyms
  
  ch.names <- NULL 
  
  if (exists("mwt", dfr)) {
    ch.names <- c(ch.names, "mwt")
    colnames(dfr)[colnames(dfr) == "mwt"] <- "atw"
  }
  
  if (exists("mwmt", dfr)) {
    ch.names <- c(ch.names, "mwmt")
    colnames(dfr)[colnames(dfr) == "mwmt"] <- "atmw"
  }
  
  if (exists("stfw", dfr)) {
    ch.names <- c(ch.names, "stfw")
    colnames(dfr)[colnames(dfr) == "stfw"] <- "sfw"
  }
  
  if (exists("stdw", dfr)) {
    ch.names <- c(ch.names, "stdw")
    colnames(dfr)[colnames(dfr) == "stdw"] <- "sdw"
  }
  
  if (exists("dm", dfr)) {
    ch.names <- c(ch.names, "dm")
    colnames(dfr)[colnames(dfr) == "dm"] <- "pdm"
  }
  
  if (exists("protein", dfr)) {
    ch.names <- c(ch.names, "protein")
    colnames(dfr)[colnames(dfr) == "protein"] <- "pro"
  }
  
  old.names <- c("mwt", "mwmt", "stfw", "stdw", "dm", "protein")
  new.names <- c("atw", "atmw", "sfw", "sdw", "pdm", "pro")
  
  if (!is.null(ch.names)) {
    ch.names.list <- old.names %in% ch.names
    warning("Traits ", list(old.names[ch.names.list]), " changed to ", list(new.names[ch.names.list]), call. = FALSE)
  }
  
  # Warnings
  
  if (max(check.list.1) == 1)
    warning("Some columns with invalid names: ", list(colnames.list[check.list.1]), call. = FALSE)
  
  if (max(check.list.2) == 1)
    warning("Some labels converted to lower case: ", list(temp[check.list.2]), call. = FALSE)
  
  # Return
  
  dfr
  
}
