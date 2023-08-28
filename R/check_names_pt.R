#' Check fieldbook traits names for potato
#'
#' Check that fieldbook traits names correspond with the names defined in 
#' \url{https://cropontology.org/term/CO_330:ROOT}
#' @param dfr The name of the data frame.
#' @param add Additional traits.
#' @details The data frame must use the labels (lower or upper case) listed below.
#' Between parentheses are the CO numbers for variables defined in
#' \url{https://cropontology.org/term/CO_330:ROOT}.
#' 
#'  -------------------- Plot identifiers --------------------
#'  \itemize{
#'  \item \code{plot}        : Plot number.
#'  \item \code{row}         : Row number.
#'  \item \code{col}         : Column number.
#'  }
#' -------------------- Classification factors --------------------
#'  \itemize{
#'  \item \code{loc}       : Locations.
#'  \item \code{year}      : Years.
#'  \item \code{season}    : Seasons. 
#'  \item \code{env}       : Environments.
#'  \item \code{geno}      : Genotypes (\code{accession_name} is also valid).
#'  \item \code{instn}     : Institutional CIP number.
#'  \item \code{rep}       : Replications (\code{rep_number} is also valid).
#'  \item \code{block}     : Blocks (\code{block_number} is also valid).
#'  \item \code{row_number}: Row number position of the plot.
#'  \item \code{col_number}: Column number position of the plot.
#'  \item \code{treatment} : Irrigation treatments.
#'  }
#' -------------------- Traits N1 group --------------------
#'  \itemize{
#'  \item \code{ntp}         : Number of tubers planted per plot (CO_330:0000265).
#'  \item \code{npe}         : Number of plants emerged per plot (CO_330:0000268).
#'  \item \code{nph}         : Number of plants harvested per plot (CO_330:0000287).
#'  \item \code{ppe}         : Proportion of plants emerged (CO_330:0000283).
#'  \item \code{pph}         : Proportion of plants harvested (CO_330:0000290).
#'  }
#' -------------------- Traits N2 group --------------------
#'  \itemize{
#'  \item \code{snpp}        : Stem number per plant (CO_330:0000851).
#'  \item \code{nipp}        : Number of inflorescences per plant estimating 0-7 (CO_330:0000226).
#'  \item \code{nfwp}        : Number of flowers per main inflorescence estimating 0-5 (CO_330:0000227).
#'  \item \code{Num_Stolon}  : Number of stolons estimating 1-9 (CO_330:0000340).
#'  \item \code{Leng_Stolon} : Lenght of stolons estimating 1-9 (CO_330:0000344).
#'  }
#' -------------------- Traits N3 group --------------------
#'  \itemize{
#'  \item \code{tntp}        : Total number of tuber per plot (CO_330:0000304).
#'  \item \code{tntpl}       : Total number of tuber per plant (CO_330:0000305).
#'  \item \code{nmtp}        : Number of marketable tubers per plot (CO_330:0000293).
#'  \item \code{nmtpl}       : Number of marketable tubers per plant (CO_330:0000297).
#'  \item \code{nnomtp}      : Number of non-marketable tubers per plot (CO_330:0000300).
#'  \item \code{nmtci}       : Number of marketable tubers category I per plot.
#'  \item \code{nmtcii}      : Number of marketable tubers category II per plot.
#'  }
#' -------------------- Traits N4 group --------------------
#'  \itemize{
#'  \item \code{ttwp}        : Total tuber weight per plot in kg (CO_330:0000317).
#'  \item \code{ttwpl}       : Total tuber weight per plant in kg (CO_330:0000321).
#'  \item \code{mtwp}        : Marketable tuber weight per plot in kg (CO_330:0000308).
#'  \item \code{mtwpl}       : Marketable tuber weight per plant in kg (CO_330:0000311).
#'  \item \code{nomtwp}      : Non-marketable tuber weight per plot in kg (CO_330:0000314).
#'  \item \code{mtwci}       : Marketable tuber weight category I per plot in kg.
#'  \item \code{mtwcii}      : Marketable tuber weight category II per plot in kg.
#'  }
#' -------------------- Traits N5 group--------------------
#'  \itemize{
#'  \item \code{ttya}        : Total tuber yield adjusted by number of plants harvested per plot in tons per ha (CO_330:0000323).
#'  \item \code{ttyna}       : Total tuber yield no adjusted in tons per ha (CO_330:0000324).
#'  \item \code{mtya}        : Marketable tuber yield adjusted by number of plants harvested per plot in tons per ha (CO_330:0000327).
#'  \item \code{mtyna}       : Marketable tuber yield no adjusted in tons per ha (CO_330:0000330).
#'  \item \code{atw}         : Average of tuber weight in g (CO_330:0000333).
#'  \item \code{atmw}        : Average of marketable tuber weight in g (CO_330:0000336).
#'  }
#' -------------------- Traits N6 group--------------------
#'  \itemize{
#'  \item \code{tafw}        : Total aerial fresh weight per plot in g.
#'  \item \code{tbfw}        : Total biomass fresh weight per plot in g.
#'  \item \code{tbfwp}       : Total biomass fresh weight per plant in g (CO_330:0000799).
#'  \item \code{hi_fw}       : Harvest index fresh weight in percentage (CO_330:0000738).
#'  \item \code{tadw}        : Total aerial dry weight per plot in g.
#'  \item \code{tbdw}        : Total biomass dry weight per plot in g.
#'  \item \code{tbdwp}       : Total biomass dry weight per plant in g (CO_330:0000799).
#'  \item \code{hi_dw}       : Harvest index dry weight in percentage (CO_330:0000735).
#'  }
#' -------------------- Traits N7 group--------------------
#'  \itemize{
#'  \item \code{fwts}         : Fresh weight of tuber for dry matter determination.
#'  \item \code{dwts}         : Dry weight of tuber for dry matter determination.
#'  \item \code{fwts1}        : Fresh weight of tuber sample 1.
#'  \item \code{fwts2}        : Fresh weight of tuber sample 2.
#'  \item \code{dwts1}        : Dry weight of tuber sample 1.
#'  \item \code{dwts2}        : Dry weight of tuber sample 2.
#'  \item \code{dm1}          : Tuber dry matter content sample 1.
#'  \item \code{dm2}          : Tuber dry matter content sample 2.
#'  \item \code{dm}           : Tuber dry matter content in percentage (CO_330:0000390).
#'  }
#' -------------------- Traits N8 group--------------------
#'  \itemize{
#'  \item \code{rd}          : Root density estimating 1-3 (CO_330:0000834).
#'  \item \code{rl}          : Root length (CO_330:0000764).
#'  \item \code{sg}          : Tuber specific gravity (CO_330:0000394).
#'  \item \code{dsi}         : Drought susceptibility index (CO_330:0000776).
#'  \item \code{dti}         : Drought tolerance index (CO_330:0000779).
#'  }
#' -------------------- Traits N9 group--------------------
#'  \itemize{
#'  \item \code{fedw}        : Tuber iron concentration in dry weight basis in mg per kg (CO_330:0000402).
#'  \item \code{fefw}        : Tuber iron concentration in fresh weight basis in mg per kg (CO_330:0000407).
#'  \item \code{zndw}        : Tuber zinc concentration in dry weight basis in mg per kg (CO_330:0000408).
#'  \item \code{znfw}        : Tuber zinc concentration in fresh weight basis in mg per kg (CO_330:0000412).
#'  \item \code{antho_dw}    : Tuber anthocyanin concentration in dry weight basis in mg per 100 g (CO_330:0000414).
#'  \item \code{antho_fw}    : Tuber anthocyanin concentration in fresh weight basis in mg per 100 g (CO_330:0000417).
#'  \item \code{aah_dw}      : Tuber hydrophilic antioxidant activity in dry weight basis in ugTroloxequiv/g (CO_330:0000418).
#'  \item \code{aah_fw}      : Tuber hydrophilic antioxidant activity in fresh weight basis in ugTroloxequiv/g (CO_330:0000421).
#'  \item \code{aal_dw}      : Tuber lipophilic antioxidant activity in dry weight basis in ugTroloxequiv/g (CO_330:0000422).
#'  \item \code{aal_fw}      : Tuber lipophilic antioxidant activity in fresh weight basis in ugTroloxequiv/g (CO_330:0000425).
#'  \item \code{asc_dw}      : Tuber ascorbic acid concentration in dry weight basis in mg per 100 g (CO_330:0000426).
#'  \item \code{asc_fw}      : Tuber ascorbic acid concentration in fresh weight basis in mg per 100 g (CO_330:0000429).
#'  \item \code{pro}         : Tuber protein content in percentage (CO_330:0000432).
#'  \item \code{star}        : Tuber starch content in percentage (CO_330:0000435).
#'  \item \code{fruc}        : Tuber fructose content in percentage (CO_330:0000438).
#'  \item \code{gluc}        : Tuber glucose content in percentage (CO_330:0000441).
#'  \item \code{sucr}        : Tuber sucrose content in percentage (CO_330:0000444).
#'  \item \code{malt}        : Tuber maltose content in percentage (CO_330:0000447).
#'  \item \code{fiber}       : Tuber fiber content in percentage (CO_330:0000450).
#'  }
#' -------------------- Traits N10 group--------------------
#'  \itemize{
#'  \item \code{plant_unif}  : Plant uniformity estimating 1-9 (CO_330:0000272).
#'  \item \code{plant_vigor} : Plant vigor estimating 1-9 (CO_330:0000276).
#'  \item \code{se}          : Senescence estimating 1-9 (CO_330:0000280).
#'  \item \code{tuber_apper} : Tuber appearance estimating 1-9 (CO_330:0000348).
#'  \item \code{tub_unif}    : Tuber uniformity estimating 1-9 (CO_330:0000352).
#'  \item \code{tub_size}    : Tuber size estimating 1-9 (CO_330:0000356).
#'  }
#' -------------------- Traits N11 group--------------------
#'  \itemize{
#'  \item \code{leaflet_tw}  : Leaflet turgid weight in mg (CO_330:0000818).
#'  \item \code{leaflet_twi} : Leaflet turgid weight evaluation i (i = 1,..., 5).
#'  }
#' -------------------- Traits N12 group--------------------
#'  \itemize{
#'  \item \code{pw_ev}       : Plant wilting estimating 1-9 (CO_330:0000773).
#'  \item \code{pw_evi}      : Plant wilting evaluation i (i = 1,..., 5)
#'  \item \code{snpp}        : Stem number per plant (CO_330:0000851).
#'  \item \code{snppi}       : Stem number per plant evaluation i (i = 1,..., 5).
#'  }
#' -------------------- Traits N13 group--------------------
#'  \itemize{
#'  \item \code{plahe_ev}    : Plant height in cm (CO_330:0000758).
#'  \item \code{plahe_evi}   : Plant height evaluation i (i = 1,..., 5).
#'  \item \code{sd_ev}       : Stem diameter in mm (CO_330:0000847).
#'  \item \code{sd_evi}      : Stem diameter evaluation i (i = 1,..., 5)
#'  \item \code{chlspad_ev}  : Chlorophyll content SPAD units (CO_330:0000805).
#'  \item \code{chlspad_evi} : Chlorophyll content index evaluation i (i = 1,..., 5).
#'  \item \code{cr_ev}       : Canopy reflectance NDVI (CO_330:0000860).
#'  \item \code{cr_evi}      : Canopy reflectance evaluation i (i = 1,..., 5).
#'  \item \code{lfa_ev}      : Leaflet area en cm square (CO_330:0000766).
#'  \item \code{lfa_evi}     : Leaflet area evaluation i (i = 1,..., 5).
#'  \item \code{rwc_ev}      : Relative water content in percentage (CO_330:0000765).
#'  \item \code{rwc_evi}     : Relative water content evaluation i (i = 1,..., 5).
#'  \item \code{sla_ev}      : Leaflet area in cm square per g (CO_330:0000812).
#'  \item \code{sla_evi}     : Leaflet area evaluation i (i = 1,..., 5).
#'  }
#' -------------------- Traits N14 group--------------------
#'  \itemize{
#'  \item \code{lb}          : Late Blight severity in percentage (CO_330:0000359).
#'  \item \code{lbi}         : Late Blight severity evaluation i (i = 1, 2,... 8).
#'  \item \code{audpc}       : Late blight AUDPC (CO_330:0000363).
#'  \item \code{raudpc}      : Late blight relative AUDPC (CO_330:0000367).
#'  \item \code{saudpc}      : Late blight susceptibility (CO_330:0000371).
#'  }
#' @return It returns a data frame with all traits names in lower case, and a list of the
#' traits with names not included in the list shown above.
#' @author Raul Eyzaguirre, Johan Ninanya. 
#' @examples
#' check.names.pt(potatoyield)
#' @export

check.names.pt <- function(dfr, add = NULL) {
  
  # Valid names for factors
  
  plot.id <- c("plot", "row", "col")
  
  factors <- c("loc", "year", "season", "env", "geno", "accession_name", "instn",
               "rep", 'rep_number', "block", "block_number", 'row_number',
               'col_number', "treatment")
  
  # Repeated measures traits (i = 1...5)
  
  rm.traits <- c("leaflet_tw", "pw_ev", "snpp", "plahe_ev", "sd_ev", "chlspad_ev",
                 "cr_ev", "lfa_ev", "rwc_ev", "sla_ev")
  
  rm.traitsi <- c(rm.traits, sapply(rm.traits, function(x) paste0(x, 1:5)))
  
  # Late blight evaluations
  
  lb.ev <- paste0("lb", 1:8)
  
  # Valid names for traits

  traits <- c("ntp", "npe", "nph", "ppe", "pph",
              "snpp", "nipp", "nfwp", "num_stolon", "leng_stolon",
              "tntp", "tntpl", "nmtp", "nmtpl", "nnomtp", "nmtci", "nmtcii",
              "ttwp", "ttwpl", "mtwp", "mtwpl", "nomtwp", "mtwci", "mtwcii",
              "ttya", "ttyna", "mtya", "mtyna", "atw", "atmw",
              'tafw', 'tbfw', "tbfwp", "hi_fw",
              'tadw', 'tbdw', "tbdwp", "hi_dw",
              'fwts', 'dwts', "fwts1", "fwts2", "dwts1", "dwts2", "dm1", "dm2", "dm",
              "rd", "rl", "sg", "dsi", "dti",
              "fedw", "fefw", "zndw", "znfw", "antho_dw", "antho_fw",
              "aah_dw", "aah_fw", "aal_dw", "aal_fw", "asc_dw", "asc_fw",
              "pro", "star", "fruc", "gluc", "sucr", "malt", "fiber",
              "plant_unif", "plant_vigor", "se",
              "tuber_apper", "tub_unif", "tub_size",
              rm.traits, rm.traitsi, "lb", lb.ev, "audpc", "raudpc", "saudpc")
  
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
  
  # Solve synonyms
  
  change.names <- NULL 
  
  if (exists("mwt", dfr)) {
    change.names <- c(change.names, "mwt")
    colnames(dfr)[colnames(dfr) == "mwt"] <- "atw"
  }
  
  if (exists("mwmt", dfr)) {
    change.names <- c(change.names, "mwmt")
    colnames(dfr)[colnames(dfr) == "mwmt"] <- "atmw"
  }
  
  if (exists("stfw", dfr)) {
    change.names <- c(change.names, "stfw")
    colnames(dfr)[colnames(dfr) == "stfw"] <- "sfw"
  }
  
  if (exists("stdw", dfr)) {
    change.names <- c(change.names, "stdw")
    colnames(dfr)[colnames(dfr) == "stdw"] <- "sdw"
  }
  
  if (exists("pdm", dfr)) {
    change.names <- c(change.names, "pdm")
    colnames(dfr)[colnames(dfr) == "pdm"] <- "dm"
  }
  
  if (exists("avdm", dfr)) {
    change.names <- c(change.names, "avdm")
    colnames(dfr)[colnames(dfr) == "avdm"] <- "dm"
  }
  
  if (exists("protein", dfr)) {
    change.names <- c(change.names, "protein")
    colnames(dfr)[colnames(dfr) == "protein"] <- "pro"
  }
  
  old.names <- c("mwt", "mwmt", "stfw", "stdw", "pdm", 'avdm', "protein")
  new.names <- c("atw", "atmw", "sfw", "sdw", "dm", 'dm', "pro")
  
  # Warnings
  
  if (!is.null(change.names)) {
    change.names.list <- old.names %in% change.names
    warning("Trait's names ", list(old.names[change.names.list]), " changed to ", list(new.names[change.names.list]), call. = FALSE)
  }
  
  if (max(check.list.1) == 1)
    warning("Some columns with invalid names: ", list(colnames.fb[check.list.1]), call. = FALSE)
  
  if (max(check.list.2) == 1)
    warning("Some labels converted to lower case: ", list(colnames.fb.valid[check.list.2]), call. = FALSE)
  
  # Return
  
  dfr
  
}
