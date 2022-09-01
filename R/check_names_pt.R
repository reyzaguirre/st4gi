#' Check fieldbook traits names for potato
#'
#' Check that fieldbook traits names correspond with the names defined in 
#' \url{https://cropontology.org/term/CO_330:ROOT}
#' @param dfr The name of the data frame.
#' @param add Additional traits.
#' @details The data frame must use the labels (lower or upper case) listed below.
#' Between parentheses the CO numbers for variable/trait according to
#' https://github.com/Planteome/CO_330-potato-traits/blob/master/potato_trait.obo
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
#'  \item \code{ntp}         : Number of tubers planted (CO_330:0000265/CO_330:0000263)
#'  \item \code{npe}         : Number of plants emerged (CO_330:0000268/CO_330:0000266)
#'  \item \code{nph}         : Number of plants harvested (CO_330:0000287/CO_330:0000284)
#'  \item \code{ppe}         : Proportion of plants emerged (CO_330:0000283/CO_330:0000281)
#'  \item \code{pph}         : Proportion of plants harvested (CO_330:0000290/CO_330:0000288)
#'  }
#' -------------------- Traits N2 group --------------------
#'  \itemize{
#'  \item \code{snpp}        : Stem number per plant (CO_330:0000851/CO_330:0000848)
#'  \item \code{nipp}        : Number of inflorescences per plant (CO_330:0000226/CO_330:0000009)
#'  \item \code{nfwp}        : Number of flowers per main inflorescence (CO_330:0000227/CO_330:0000012)
#'  \item \code{nlpp}        : Number of leaves per plant
#'  \item \code{Num_Stolon}  : Number of stolons (CO_330:0000340/CO_330:0000337)
#'  \item \code{Leng_Stolon} : Lenght of stolons (CO_330:0000344/CO_330:0000341)
#'  }
#' -------------------- Traits N3 group --------------------
#'  \itemize{
#'  \item \code{tntp}        : Total number of tuber per plot (CO_330:0000304/CO_330:0000301)
#'  \item \code{tntpl}       : Total number of tuber per plant (CO_330:0000305/CO_330:0000108)
#'  \item \code{nmtp}        : Number of marketable tubers per plot (CO_330:0000293/CO_330:0000291)
#'  \item \code{nmtpl}       : Number of marketable tubers per plant (CO_330:0000297/CO_330:0000294)
#'  \item \code{nnomtp}      : Number of non-marketable tubers per plot (CO_330:0000300/CO_330:0000298)
#'  \item \code{nmtci}       : Number of marketable tubers category I per plot
#'  \item \code{nmtcii}      : Number of marketable tubers category II per plot
#'  }
#' -------------------- Traits N4 group --------------------
#'  \itemize{
#'  \item \code{ttwp}        : Total tuber weight per plot (CO_330:0000317/CO_330:0000315)
#'  \item \code{ttwpl}       : Total tuber weight per plant (CO_330:0000321/CO_330:0000318)
#'  \item \code{mtwp}        : Marketable tuber weight per plot (CO_330:0000308/CO_330:0000306)
#'  \item \code{mtwpl}       : Marketable tuber weight per plant (CO_330:0000311/CO_330:0000309)
#'  \item \code{nomtwp}      : Non-marketable tuber weight per plot (CO_330:0000314/CO_330:0000312)
#'  \item \code{mtwci}       : Marketable tuber weight category I per plot
#'  \item \code{mtwcii}      : Marketable tuber weight category II per plot
#'  }
#' -------------------- Traits N5 group--------------------
#'  \itemize{
#'  \item \code{ttya}        : Total tuber yield adjusted (CO_330:0000323/CO_330:0000144)
#'  \item \code{ttyna}       : Total tuber yield no adjusted (CO_330:0000324/CO_330:0000144)
#'  \item \code{mtya}        : Marketable tuber yield adjusted (CO_330:0000327/CO_330:0000325)
#'  \item \code{mtyna}       : Marketable tuber yiedl no adjusted (CO_330:0000330/CO_330:0000328)
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
#'  \item \code{hi_fw}       : Harvest index fresh weight (CO_330:0000738/CO_330:0000733)
#'  \item \code{stldw}       : Stolon dry weight per plant
#'  \item \code{sdw}         : Stem dry weight per plant (\code{stdw} is also valid)
#'  \item \code{ldw}         : Leaf dry weight per plant
#'  \item \code{rdw}         : Root dry weight per plant
#'  \item \code{tdw}         : Tuber dry weight per plant
#'  \item \code{tbdw}        : Total biomass dry weight per plant
#'  \item \code{hi_dw}       : Harvest index dry weight (CO_330:0000735/CO_330:0000733)
#'  \item \code{ldmcp}       : Leaf dry matter content per plot
#'  \item \code{sdmcp}       : Stem dry matter content per plot
#'  \item \code{rdmcp}       : Root dry matter content per plot
#'  \item \code{tdmcp}       : Tuber dry matter content per plot
#'  }
#' -------------------- Traits N7 group--------------------
#'  \itemize{
#'  \item \code{fwts1}        : Fresh weight of tuber sample 1
#'  \item \code{fwts2}        : Fresh weight of tuber sample 2
#'  \item \code{dwts1}        : Dry weight of tuber sample 1
#'  \item \code{dwts2}        : Dry weight of tuber sample 2
#'  \item \code{dm1}          : Tuber dry matter content sample 1
#'  \item \code{dm2}          : Tuber dry matter content sample 2
#'  \item \code{avdm}         : Average tuber dry matter content (CO_330:0000390/CO_330:0000192)
#'  }
#' -------------------- Traits N8 group--------------------
#'  \itemize{
#'  \item \code{twa}         : Tuber weight in air
#'  \item \code{tww}         : Tuber weight in water
#'  \item \code{rsdw}        : Root system dry weight per plant
#'  \item \code{rd}          : Root density (CO_330:0000834/CO_330:0000831)
#'  \item \code{rl}          : Root length (CO_330:0000764/CO_330:0000762)
#'  \item \code{sg}          : Tuber specific gravity (CO_330:0000394/CO_330:0000391)
#'  \item \code{dsi}         : Drought susceptibility index (CO_330:0000776/CO_330:0000774)
#'  \item \code{dti}         : Drought tolerance index (CO_330:0000779/CO_330:0000777)
#'  }
#' -------------------- Traits N9 group--------------------
#'  \itemize{
#'  \item \code{pdm}         : Tuber dry matter content (\code{dm} is also valid)
#'  \item \code{fedw}        : Tuber iron concentration in dry weight basis (CO_330:0000402/CO_330:0000201)
#'  \item \code{fefw}        : Tuber iron concentration in fresh weight basis (CO_330:0000407/CO_330:0000405)
#'  \item \code{zndw}        : Tuber zinc concentration in dry weight basis (CO_330:0000408/CO_330:0000204)
#'  \item \code{znfw}        : Tuber zinc concentration in fresh weight basis (CO_330:0000412/CO_330:0000410)
#'  \item \code{antho_dw}    : Tuber anthocyanin concentration in dry weight basis (CO_330:0000414/CO_330:0000413)
#'  \item \code{antho_fw}    : Tuber anthocyanin concentration in fresh weight basis (CO_330:0000417/CO_330:0000415)
#'  \item \code{aah_dw}      : Tuber hydrophilic antioxidant activity in dry weight basis (CO_330:0000418/CO_330:0000210)
#'  \item \code{aah_fw}      : Tuber hydrophilic antioxidant activity in fresh weight basis (CO_330:0000421/CO_330:0000419)
#'  \item \code{aal_dw}      : Tuber lipophilic antioxidant activity in dry weight basis (CO_330:0000422/CO_330:0000213)
#'  \item \code{aal_fw}      : Tuber lipophilic antioxidant activity in fresh weight basis (CO_330:0000425/CO_330:0000423)
#'  \item \code{asc_dw}      : Tuber ascorbic acid concentration in dry weight basis (CO_330:0000426/CO_330:0000198)
#'  \item \code{asc_fw}      : Tuber ascorbic acid concentration in fresh weight basis (CO_330:0000429/CO_330:0000427)
#'  \item \code{pro}         : Tuber protein content (\code{protein} is also valid) (CO_330:0000432/CO_330:0000430)
#'  \item \code{star}        : Tuber starch content (CO_330:0000435/CO_330:0000433)
#'  \item \code{fruc}        : Tuber fructose content (CO_330:0000438/CO_330:0000436)
#'  \item \code{gluc}        : Tuber glucose content (CO_330:0000441/CO_330:0000439)
#'  \item \code{sucr}        : Tuber sucrose content (CO_330:0000444/CO_330:0000442)
#'  \item \code{malt}        : Tuber maltose content (CO_330:0000447/CO_330:0000445)
#'  \item \code{fiber}       : Tuber fiber content (CO_330:0000450/CO_330:0000448)
#'  }
#' -------------------- Traits N10 group--------------------
#'  \itemize{
#'  \item \code{plant_unif}  : Plant uniformity (CO_330:0000272/CO_330:0000269)
#'  \item \code{plant_vigor} : Plant vigor (CO_330:0000276/CO_330:0000273)
#'  \item \code{se}          : Senescence (CO_330:0000280/CO_330:0000277)
#'  \item \code{tuber_apper} : Tuber appearance (CO_330:0000348/CO_330:0000345)
#'  \item \code{tub_unif}    : Tuber uniformity (CO_330:0000352/CO_330:0000349)
#'  \item \code{tub_size}    : Tuber size (CO_330:0000356/CO_330:0000353)
#'  \item \code{op}          : Osmotic potential
#'  \item \code{tlwp}        : Total leaf water potential
#'  \item \code{ltp}         : Leaf turgor potential
#'  }
#' -------------------- Traits N11 group--------------------
#'  \itemize{
#'  \item \code{leaflet_twi} : Leaflet turgid weight Evaluation i (i = 1,..., 5)
#'  \item \code{insnppi}     : Increase stem number per plant i (i = 1,..., 5)
#'  \item \code{insdi}       : Increase stem diameter i (i = 1,...,5)
#'  \item \code{inrwci}      : Increase relative water content i (i = 1,..., 5)
#'  }
#' -------------------- Traits N12 group--------------------
#'  \itemize{
#'  \item \code{pw_evi}      : Plant wilting Evaluation i (i = 1,..., 5)
#'  \item \code{inplahei}    : Increase plant height i (i = 1,..., 5)
#'  \item \code{snppi}       : Stem number per plant Evaluation i (i = 1,..., 5)
#'  \item \code{leaflet_fwi} : Leaflet fresh weight Evaluation i (i = 1,..., 5)
#'  \item \code{leaflet_dwi} : Leaflet dry weight Evaluation i (i = 1,..., 5)
#'  }
#' -------------------- Traits N13 group--------------------
#'  \itemize{
#'  \item \code{chci}        : Chlorophyll content sample i (i = 1,..., 5)
#'  \item \code{inchci}      : Increase chlorophyll content i (i = 1,..., 5)
#'  \item \code{leafsdi}     : Leaf stomata density sample i (i = 1,..., 5)
#'  \item \code{av_leafsd}   : Average leaf stomata density
#'  }
#' -------------------- Traits N14 group--------------------
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
#' -------------------- Traits N15 group--------------------
#'  \itemize{
#'  \item \code{plahe_slp}   : Plant height slope
#'  \item \code{sd_slp}      : Stem diameter slope
#'  \item \code{cc_slp}      : Canopy cover slope
#'  \item \code{chlspad_slp} : Chlorophyll content index slope
#'  \item \code{cr_slp}      : Canopy reflectance slope
#'  \item \code{lfa_slp}     : Leaflet area slope
#'  \item \code{sla_slp}     : Specific leaf area slope         
#'  }
#' -------------------- Traits N16 group--------------------
#'  \itemize{
#'  \item \code{lbi}         : Late Blight severity evaluation i (i = 1, 2,... 8)
#'  \item \code{audpc}       : Late blight AUDPC (CO_330:0000363/CO_330:0000360)
#'  \item \code{raudpc}      : Late blight relative AUDPC (CO_330:0000367/CO_330:0000364)
#'  \item \code{saudpc}     : Late blight susceptibility (CO_330:0000371/CO_330:0000368)
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
              "stlfw", "sfw", "stfw", "lfw", "rfw", "tfw", "tbfw", "hi_fw",
              "stldw", "sdw", "stdw", "ldw", "rdw", "tdw", "tbdw", "hi_dw",
              "ldmcp", "sdmcp", "rdmcp", "tdmcp", "fwts1", "fwts2", "dwts1",
              "dwts2", "dm1", "dm2", "avdm", "pdm", "dm",
              "twa", "tww", "rsdw", "rd", "rl", "sg", "dsi", "dti",
              "fedw", "fefw", "zndw", "znfw", "antho_dw", "antho_fw",
              "aah_dw", "aah_fw", "aal_dw", "aal_fw", "asc_dw", "asc_fw",
              "pro", "protein", "star", "fruc", "gluc", "sucr", "malt", "fiber",
              "plant_unif", "plant_vigor", "se", "op", "tlwp", "ltp",
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
