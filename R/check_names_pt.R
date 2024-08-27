#' Check fieldbook traits names for potato
#'
#' Check that fieldbook traits names correspond with the names defined in 
#' \url{https://cropontology.org/term/CO_330:ROOT},
#' and that all traits are stored as numeric.
#' @param dfr The name of the data frame.
#' @param add Additional traits.
#' @details The data frame must use the labels (lower or upper case) listed below.
#' Between parentheses are the CO numbers for variables defined in
#' \url{https://cropontology.org/term/CO_330:ROOT}.
#' 
#' ----------------------------- Plot identifiers -----------------------------
#'  \itemize{
#'  \item \code{plot}        : Plot number.
#'  \item \code{row}         : Row number position of the plot (\code{row_number} is also valid).
#'  \item \code{col}         : Column number position of the plot (\code{col_number} is also valid).
#'  }
#' --------------------------- Classification factors--------------------------
#'  \itemize{
#'  \item \code{loc}         : Locations.
#'  \item \code{year}        : Years.
#'  \item \code{season}      : Seasons. 
#'  \item \code{env}         : Environments.
#'  \item \code{geno}        : Genotypes (\code{accession_name} and \code{instn} are also valid).
#'  \item \code{type}        : Entry type (\code{clon}, \code{check}, \code{progeny}, \code{parent}).
#'  \item \code{rep}         : Replications (\code{rep_number} is also valid).
#'  \item \code{block}       : Blocks (\code{block_number} is also valid).
#'  \item \code{treat}       : Treatments.
#'  }
#' ----------------------------------- Plants ---------------------------------
#'  \itemize{
#'  \item \code{ntp}         : Number of tubers planted per plot (CO_330:0000265).
#'  \item \code{npe}         : Number of plants emerged per plot (CO_330:0000268).
#'  \item \code{npe_15dap}   : Number of plants emerged per plot, 15 days after planting (COMP:0000002).
#'  \item \code{npe_30dap}   : Number of plants emerged per plot, 30 days after planting (COMP:0000003).
#'  \item \code{nph}         : Number of plants harvested per plot (CO_330:0000287).
#'  \item \code{ppe}         : Proportion of plants emerged in percentage (CO_330:0000283).
#'  \item \code{pph}         : Proportion of plants harvested in percentage (CO_330:0000290).
#'  }
#' ----------------------------- Number of tubers -----------------------------
#'  \itemize{
#'  \item \code{tntp}        : Total number of tubers per plot (CO_330:0000304).
#'  \item \code{tntpl}       : Total number of tubers per plant (CO_330:0000305).
#'  \item \code{nmtp}        : Number of marketable tubers per plot (CO_330:0000293).
#'  \item \code{nmtpl}       : Number of marketable tubers per plant (CO_330:0000297).
#'  \item \code{nnomtp}      : Number of non-marketable tubers per plot (CO_330:0000300).
#'  \item \code{nmtci}       : Number of marketable tubers category I per plot.
#'  \item \code{nmtcii}      : Number of marketable tubers category II per plot.
#'  }
#' ----------------------------- Weight of tubers -----------------------------
#'  \itemize{
#'  \item \code{ttwp}        : Total tuber weight per plot in kg (CO_330:0000317).
#'  \item \code{ttwpl}       : Total tuber weight per plant in kg (CO_330:0000321).
#'  \item \code{mtwp}        : Marketable tuber weight per plot in kg (CO_330:0000308).
#'  \item \code{mtwpl}       : Marketable tuber weight per plant in kg (CO_330:0000311).
#'  \item \code{nomtwp}      : Non-marketable tuber weight per plot in kg (CO_330:0000314).
#'  \item \code{mtwci}       : Marketable tuber weight category I per plot in kg.
#'  \item \code{mtwcii}      : Marketable tuber weight category II per plot in kg.
#'  \item \code{ttya}        : Total tuber yield adjusted by number of plants harvested per plot in tons per ha (CO_330:0000323).
#'  \item \code{ttyna}       : Total tuber yield not adjusted in tons per ha (CO_330:0000324).
#'  \item \code{mtya}        : Marketable tuber yield adjusted by number of plants harvested per plot in tons per ha (CO_330:0000327).
#'  \item \code{mtyna}       : Marketable tuber yield not adjusted in tons per ha (CO_330:0000330).
#'  \item \code{atw}         : Average tuber weight in g (CO_330:0000333).
#'  \item \code{atmw}        : Average marketable tuber weight in g (CO_330:0000336).
#'  }
#' ------------------------ Plant and tuber appearance ------------------------
#'  \itemize{
#'  \item \code{flowering}         : Flowering degree estimating 0-7 (CO_330:0000224).
#'  \item \code{flowering_45dap}   : Flowering degree estimating 0-7, 45 days after planting (COMP:0000010).
#'  \item \code{flowering_60dap}   : Flowering degree estimating 0-7, 60 days after planting (COMP:0000011).
#'  \item \code{plant_unif}        : Plant uniformity estimating 1-9 (CO_330:0000272).
#'  \item \code{plant_unif_45dap}  : Plant uniformity estimating 1-9, 45 days after planting (COMP:0000008).
#'  \item \code{plant_unif_60dap}  : Plant uniformity estimating 1-9, 60 days after planting (COMP:0000009).
#'  \item \code{plant_vigor}       : Plant vigor estimating 1-9 (CO_330:0000276).
#'  \item \code{plant_vigor_30dap} : Plant vigor estimating 1-9, 30 days after planting (COMP:0000004).
#'  \item \code{plant_vigor_45dap} : Plant vigor estimating 1-9, 45 days after planting (COMP:0000005).
#'  \item \code{plant_vigor_60dap} : Plant vigor estimating 1-9, 60 days after planting (COMP:0000006).
#'  \item \code{se}                : Senescence estimating 1-9 (CO_330:0000280).
#'  \item \code{tuber_apper}       : Tuber appearance estimating 1-9 (CO_330:0000348).
#'  \item \code{tub_unif}          : Tuber uniformity estimating 1-9 (CO_330:0000352).
#'  \item \code{tub_size}          : Tuber size estimating 1-9 (CO_330:0000356).
#'  \item \code{tbskn1}            : Predominant tuber skin color (CO_330:0000249).
#'  \item \code{tbfsh1}            : Predominant tuber flesh color (CO_330:0000253).
#'  \item \code{tbshp1}            : Tuber shape (CO_330:0000256).
#'  \item \code{tbshp3}            : Tuber depth eyes (CO_330:0000258).
#'  \item \code{dorpd}             : Tuber dormancy period (CO_330:0000397).
#'  }
#' ------------------------------- Organoleptic -------------------------------
#'  \itemize{
#'  \item \code{chipping}    : Chips color 1-5 (CO_330:0000384).
#'  \item \code{ffr}         : French fries color 1-5 (CO_330:0000388).
#'  \item \code{aocp}        : Chips oil absorption rate (CO_330:0000395).
#'  \item \code{flavour}     : Tuber flavor after cooking (CO_330:0000379).
#'  \item \code{textac}      : Tuber texture after cooking (CO_330:0000380).
#'  \item \code{cookqu}      : Tuber cooking quality (CO_330:0000381).
#'  \item \code{cootime}     : Tuber cooking time (CO_330:0000389).
#'  } 
#' ---------------------------------- Biomass ---------------------------------
#'  \itemize{
#'  \item \code{tafw}        : Total aerial fresh weight per plot in g.
#'  \item \code{tbfwp}       : Total biomass fresh weight per plot in g. (CO_330:0000829)
#'  \item \code{tbfwpl}      : Total biomass fresh weight per plant in g (CO_330:0000799).
#'  \item \code{hi_fw}       : Harvest index fresh weight in percentage (CO_330:0000738).
#'  \item \code{tadw}        : Total aerial dry weight per plot in g.
#'  \item \code{tbdwp}       : Total biomass dry weight per plot in g. (CO_330:0000828)
#'  \item \code{tbdwpl}      : Total biomass dry weight per plant in g (CO_330:0000798).
#'  \item \code{hi_dw}       : Harvest index dry weight in percentage (CO_330:0000735).
#'  }
#' -------------------------------- Dry weight --------------------------------
#'  \itemize{
#'  \item \code{fwts}        : Fresh weight of tuber for dry matter determination.
#'  \item \code{dwts}        : Dry weight of tuber for dry matter determination.
#'  \item \code{fwts1}       : Fresh weight of tuber sample 1.
#'  \item \code{fwts2}       : Fresh weight of tuber sample 2.
#'  \item \code{dwts1}       : Dry weight of tuber sample 1.
#'  \item \code{dwts2}       : Dry weight of tuber sample 2.
#'  \item \code{dm1}         : Tuber dry matter content sample 1.
#'  \item \code{dm2}         : Tuber dry matter content sample 2.
#'  \item \code{dm}          : Tuber dry matter content in percentage (CO_330:0000390).
#'  }
#' --------------------------- Nutrients evaluations --------------------------
#'  \itemize{
#'  \item \code{fedw}        : Tuber iron concentration on dry weight basis in mg per 100 g (CO_330:0000402).
#'  \item \code{fefw}        : Tuber iron concentration on fresh weight basis in mg per 100 g (CO_330:0000407).
#'  \item \code{fedw_xrf}    : Tuber iron concentration on dry weight basis in mg per kg by XRF (CO_330:0000404).
#'  \item \code{zndw}        : Tuber zinc concentration on dry weight basis in mg per 100 g (CO_330:0000408).
#'  \item \code{znfw}        : Tuber zinc concentration on fresh weight basis in mg per 100 g (CO_330:0000412).
#'  \item \code{zndw_xrf}    : Tuber zinc concentration on dry weight basis in mg per kg by XRF (CO_330:0000409).
#'  \item \code{antho_dw}    : Tuber anthocyanin concentration on dry weight basis in mg per 100 g (CO_330:0000414).
#'  \item \code{antho_fw}    : Tuber anthocyanin concentration on fresh weight basis in mg per 100 g (CO_330:0000417).
#'  \item \code{aah_dw}      : Tuber hydrophilic antioxidant activity on dry weight basis in ugTroloxequiv/g (CO_330:0000418).
#'  \item \code{aah_fw}      : Tuber hydrophilic antioxidant activity on fresh weight basis in ugTroloxequiv/g (CO_330:0000421).
#'  \item \code{aal_dw}      : Tuber lipophilic antioxidant activity on dry weight basis in ugTroloxequiv/g (CO_330:0000422).
#'  \item \code{aal_fw}      : Tuber lipophilic antioxidant activity on fresh weight basis in ugTroloxequiv/g (CO_330:0000425).
#'  \item \code{asc_dw}      : Tuber ascorbic acid concentration on dry weight basis in mg per 100 g (CO_330:0000426).
#'  \item \code{asc_fw}      : Tuber ascorbic acid concentration on fresh weight basis in mg per 100 g (CO_330:0000429).
#'  \item \code{pro}         : Tuber protein content in percentage (CO_330:0000432).
#'  \item \code{star}        : Tuber starch content in percentage (CO_330:0000435).
#'  \item \code{fruc}        : Tuber fructose content in percentage (CO_330:0000438).
#'  \item \code{gluc}        : Tuber glucose content in percentage (CO_330:0000441).
#'  \item \code{sucr}        : Tuber sucrose content in percentage (CO_330:0000444).
#'  \item \code{malt}        : Tuber maltose content in percentage (CO_330:0000447).
#'  \item \code{fiber}       : Tuber fiber content in percentage (CO_330:0000450).
#'  \item \code{tgly_fw}     : Tuber glycoalkaloid concentration on fresh weight basis in mg per 100 g (CO_330:0000661).
#'  }
#' ---------------------------------- Deseases --------------------------------
#'  \itemize{
#'  \item \code{lb}          : Late Blight severity estimating percentage (CO_330:0000359).
#'  \item \code{lbi}         : Late Blight severity estimating percentage evaluation i (i = 1, 2,... 8).
#'  \item \code{rlb}         : Late Blight resistance estimating 1-6 (CO_330:0000372).
#'  \item \code{rlb_30dap}   : Late Blight resistance estimating 1-6, 30 days after planting (COMP:0000012).
#'  \item \code{rlb_45dap}   : Late Blight resistance estimating 1-6, 45 days after planting (COMP:0000013).
#'  \item \code{rlb_60dap}   : Late Blight resistance estimating 1-6, 60 days after planting (COMP:0000014).
#'  \item \code{rlb_75dap}   : Late Blight resistance estimating 1-6, 75 days after planting (COMP:0000015).
#'  \item \code{audpc}       : Late blight AUDPC computed (CO_330:0000363).
#'  \item \code{raudpc}      : Late blight relative AUDPC computed (CO_330:0000367).
#'  \item \code{saudpc}      : Late blight susceptibility computed (CO_330:0000371).
#'  \item \code{pvx}         : Potato virus X resistance (CO_330:0000373).
#'  \item \code{pvy}         : Potato virus Y resistance (CO_330:0000374).
#'  \item \code{prlv}        : Potatoleaf roll virus resistance (CO_330:0000375).
#'  \item \code{bw}          : Bacterial wilt resistance (CO_330:0000376).
#'  \item \code{rkn}         : Root knot nematode resistance (CO_330:0000377).
#'  \item \code{rlmf}        : Resistance to leaf miner fly 1-5 (CO_330:0000378).
#'  \item \code{rlmf_45dap}  : Resistance to leaf miner fly 1-5 (COMP_0000016).
#'  \item \code{rlmf_60dap}  : Resistance to leaf miner fly 1-5 (COMP_0000017).
#'  \item \code{rlmf_75dap}  : Resistance to leaf miner fly 1-5 (COMP_0000018).
#'  }
#' ------------------------------- Other traits -------------------------------
#'  \itemize{
#'  \item \code{rd}          : Root density estimating 1-3 (CO_330:0000834).
#'  \item \code{rl}          : Root length in cm (CO_330:0000764).
#'  \item \code{sg}          : Tuber specific gravity (CO_330:0000394).
#'  \item \code{dsi}         : Drought susceptibility index computed (CO_330:0000776).
#'  \item \code{dti}         : Drought tolerance index computed (CO_330:0000779).
#'  \item \code{leaflet_tw}  : Leaflet turgid weight in mg (CO_330:0000818).
#'  \item \code{leaflet_twi} : Leaflet turgid weight evaluation i (i = 1,..., 5).
#'  \item \code{pw_ev}       : Degree of canopy wilting estimating 1-9 (CO_330:0000773).
#'  \item \code{pw_evi}      : Degree of canopy wilting estimating 1-9 evaluation i (i = 1,..., 5).
#'  \item \code{snpp}        : Stem number per plant (CO_330:0000851).
#'  \item \code{snppi}       : Stem number per plant evaluation i (i = 1,..., 5).
#'  \item \code{nipp}        : Number of inflorescences per plant estimating 0-7 (CO_330:0000226).
#'  \item \code{nfwp}        : Number of flowers per main inflorescence estimating 0-5 (CO_330:0000227).
#'  \item \code{Num_Stolon}  : Number of stolons estimating 1-9 (CO_330:0000340).
#'  \item \code{Leng_Stolon} : Length of stolons estimating 1-9 (CO_330:0000344).
#'  \item \code{plahe_ev}    : Plant height in cm (CO_330:0000758).
#'  \item \code{plahe_evi}   : Plant height evaluation i (i = 1,..., 5).
#'  \item \code{sd_ev}       : Stem diameter in mm (CO_330:0000847).
#'  \item \code{sd_evi}      : Stem diameter evaluation i (i = 1,..., 5)
#'  \item \code{chlspad_ev}  : Chlorophyll content index using SPAD-meter in SPAD units (CO_330:0000805).
#'  \item \code{chlspad_evi} : Chlorophyll content index evaluation i (i = 1,..., 5).
#'  \item \code{cr_ev}       : Canopy reflectance using Green-Seeker-meter in NDVI units (CO_330:0000860).
#'  \item \code{cr_evi}      : Canopy reflectance evaluation i (i = 1,..., 5).
#'  \item \code{lfa_ev}      : Leaflet area using Compu Eye, Leaf Symptom Area software in cm^2 (CO_330:0000766).
#'  \item \code{lfa_evi}     : Leaflet area using Compu Eye, Leaf Symptom Area software in cm^2 evaluation i (i = 1,..., 5).
#'  \item \code{rwc_ev}      : Relative water content in percentage (CO_330:0000741).
#'  \item \code{rwc_evi}     : Relative water content evaluation i (i = 1,..., 5).
#'  \item \code{sla_ev}      : Leaflet area using EasyLeafArea software in cm^2 per g (CO_330:0000812).
#'  \item \code{sla_evi}     : Leaflet area using EasyLeafArea software in cm^2 per g evaluation i (i = 1,..., 5).
#'  }
#' @return It returns:
#' \itemize{
#' \item The fieldbook data frame with all column names in lowercase and
#' with some possible modifications in the names. Traits that are stored
#' with a non-numeric class are transformed to numeric.
#' \item A list of warnings for all the column names that have been changed.
#' \item A list of warnings for all the column names not recognized.
#' \item A list of warnings for all the column traits that have been changed to numeric.
#' }
#' @author Raul Eyzaguirre, Johan Ninanya. 
#' @examples
#' check.names.pt(potatoyield)
#' @export

check.names.pt <- function(dfr, add = NULL) {
  
  # Valid names for factors
  
  plot.id <- c("plot", "row", "col")
  
  factors <- c("loc", "year", "season", "env", "geno", 'type',
               "rep", "block", "treat")
  
  # Repeated measures traits (i = 1...5)
  
  rm.traits <- c("leaflet_tw", "pw_ev", "snpp", "plahe_ev", "sd_ev", "chlspad_ev",
                 "cr_ev", "lfa_ev", "rwc_ev", "sla_ev")
  
  rm.traitsi <- c(rm.traits, sapply(rm.traits, function(x) paste0(x, 1:5)))
  
  # Late blight evaluations
  
  lb.ev <- paste0("lb", 1:8)
  
  # Valid names for traits

  traits <- c(ptont$Label, "nmtci", "nmtcii", "mtwci", "mtwcii", "tbskn1", "tbfsh1",
              "tbshp1", "tbshp3", "dorpd", "chipping", "ffr", "aocp", "flavour",
              "textac", "cookqu", "cootime", "tafw", "tadw", "fwts", "dwts",
              "fwts1", "fwts2", "dwts1", "dwts2", "dm1", "dm2", "avdm", "tgly_fw",
              "pvx", "pvy", "prlv", "bw", "rkn", lb.ev, rm.traits, rm.traitsi, tolower(add))
  
  # Valid names for factors and traits
  
  colnames.valid <- c(plot.id, factors, traits)
  
  # Factors and traits in field book (original names)
  
  colnames.fb <- colnames(dfr)
  
  # Convert all fieldbook names to lower case (except CO numbers)
  
  cond <- substring(colnames(dfr), 1, 7) != 'CO_330:' & substring(colnames(dfr), 1, 5) != 'COMP:'
  colnames(dfr)[cond] <- tolower(colnames(dfr))[cond]
  
  if (sum(colnames.fb != colnames(dfr)) > 0)
    warning("Some labels converted to lower case", call. = FALSE)
  
  # Solve synonyms for factors
  
  # All options for genotypes
  
  old.geno <- c("accession_name", "cipno", "cip.number", 'genotype', "instn")
  new.geno <- rep('geno', length(old.geno))
  
  change.names.f <- NULL 
  
  old.names.f <- c('location', "replication", 'rep_number', "block_number", "row_number", "col_number", old.geno)
  new.names.f <- c('loc',      "rep",         'rep',        "block",        "row",        "col",        new.geno)
  
  for (i in 1:length(old.names.f)) {
    
    if (exists(old.names.f[i], dfr) & !exists(new.names.f[i], dfr)) {
      change.names.f <- c(change.names.f, old.names.f[i])
      colnames(dfr)[colnames(dfr) == old.names.f[i]] <- new.names.f[i]
    }
    
  }  
  
  if (!is.null(change.names.f)) {
    change.names.list <- old.names.f %in% change.names.f
    warning("Factors' names ", list(old.names.f[change.names.list]), " changed to ", list(new.names.f[change.names.list]), call. = FALSE)
  }
  
  # Solve synonyms for traits
  
  change.names.t <- NULL 
  
  old.names.t <- c("mwt", "mwmt", "stfw", "stdw", "pdm", 'avdm', "protein")
  new.names.t <- c("atw", "atmw", "sfw",  "sdw",  "dm",  'dm',   "pro")
  
  for (i in 1:length(old.names.t)) {
    
    if (exists(old.names.t[i], dfr) & !exists(new.names.t[i], dfr)) {
      change.names.t <- c(change.names.t, old.names.t[i])
      colnames(dfr)[colnames(dfr) == old.names.t[i]] <- new.names.t[i]
    }
    
  }  

  if (!is.null(change.names.t)) {
    change.names.list <- old.names.t %in% change.names.t
    warning("Traits' names ", list(old.names.t[change.names.list]), " changed to ", list(new.names.t[change.names.list]), call. = FALSE)
  }
  
  # Names not valid
  
  names.not.valid <- !(colnames(dfr) %in% colnames.valid)
  
  if (max(names.not.valid) == 1)
    warning("Some columns with invalid names: ", list(colnames(dfr)[names.not.valid]), call. = FALSE)
  
  # Check traits are numeric
  
  nonumeric.list <- NULL
  
  column.class <- unlist(lapply(dfr, class))
  for(i in colnames(dfr)) {
    if(i %in% traits & column.class[i] != "numeric") {
      dfr[, i] <- suppressWarnings(as.numeric(as.character(dfr[, i])))
      nonumeric.list <- c(nonumeric.list, i)
    }
  }

  if (!is.null(nonumeric.list))
    warning("Some traits converted to numeric: ", list(nonumeric.list), call. = FALSE)
  
  # Return
  
  dfr
  
}
