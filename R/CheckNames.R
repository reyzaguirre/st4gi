#' Check fieldbook traits names
#'
#' Check that fieldbook traits names correspond with the names defined in the document
#' "PROCEDURES  FOR THE EVALUATION AND ANALYSIS OF SWEETPOTATO TRIALS".
#' @param fb The name of the fieldbook data frame.
#' @details The data frame must use the following labels (lower or upper case):
#' \itemize{
#'  \item L       : Locations (LOC is also valid)
#'  \item Y       : Years
#'  \item S       : Seasons
#'  \item G       : Genotypes (GENO is also valid)
#'  \item NAME    : Names for genotypes
#'  \item E       : Environments (ENV is also valid)
#'  \item R       : Replications (REP is also valid)
#'  \item NOPS    : Number of plants sowed
#'  \item NOPE    : Number of plants established
#'  \item VIR1    : Virus symptoms (1-9), first evaluation
#'  \item VIR2    : Virus symptoms (1-9), second evaluation
#'  \item VIR3    : Virus symptoms (1-9), third evaluation
#'  \item ALT1    : Alternaria symptoms (1-9), first evaluation
#'  \item ALT2    : Alternaria symptoms (1-9), second evaluation
#'  \item VV1     : Vine vigor (1-9), first evaluation
#'  \item VV2     : Vine vigor2 (1-9), second evaluation
#'  \item VW      : Vine weight
#'  \item NOPH    : Number of plants harvested
#'  \item NOPR    : Number of plants with roots
#'  \item NOCR    : Number of commercial roots
#'  \item NONC    : Number of non commercial roots
#'  \item CRW     : Commercial root weight
#'  \item NCRW    : Non commercial root weight
#'  \item RFCP    : Root primary flesh color using CIP color charts
#'  \item RFCS    : Root secondary flesh color using CIP color charts
#'  \item SCOL    : Storage root skin color (1-9)
#'  \item FCOL    : Storage root flesh color (1-9)
#'  \item RFCP    : Storage root primary flesh color (1-9)
#'  \item RFCS    : Storage root secondary flesh color (1-9)
#'  \item RS      : Root size (1-9)
#'  \item RF      : Root form (1-9)
#'  \item DAMR    : Root defects (1-9)
#'  \item RSPR    : Root sprouting (1-9)
#'  \item WED1    : Weevil damage (1-9), first evaluation
#'  \item WED2    : Weevil damage2 (1-9), second evaluation
#'  \item DMF     : Fresh weight of roots for dry matter assessment
#'  \item DMD     : Dry weight of DMF samples
#'  \item DM      : Storage root dry matter content (\%)
#'  \item DMRY    : Dry matter root yield
#'  \item DMVF    : Fresh weight vines for dry matter assessment
#'  \item DMVD    : Dry weight of DMVF samples
#'  \item DMV     : Vines dry matter content (\%)
#'  \item DMFY    : Dry matter foliage yield
#'  \item FRAW1   : Root fiber (1-9), first determination
#'  \item SURAW1  : Root sugar (1-9), first determination
#'  \item STRAW1  : Root starch (1-9), first determination
#'  \item COOF1   : Cooked fiber (1-9), first evaluation
#'  \item COOSU1  : Cooked sugars (1-9), first evaluation
#'  \item COOST1  : Cooked starch (1-9), first evaluation
#'  \item COOT1   : Cooked taste (1-9), first evaluation
#'  \item COOAP1  : Cooked appearance (1-9), first evaluation
#'  \item FRAW2   : Root fiber (1-9), second determination
#'  \item SURAW2  : Root sugar (1-9), second determination
#'  \item STRAW2  : Root starch (1-9), second determination
#'  \item COOF2   : Cooked fiber (1-9), second evaluation
#'  \item COOSU2  : Cooked sugars (1-9), second evaluation
#'  \item COOST2  : Cooked starch (1-9), second evaluation
#'  \item COOT2   : Cooked taste (1-9), second evaluation
#'  \item COOAP2  : Cooked appearance (1-9), second evaluation
#'  \item PROT    : Protein (\%)
#'  \item FE      : Iron (mg/100 g dry weight)
#'  \item ZN      : Zinc (mg/100 g dry weight)
#'  \item CA      : Calcium (mg/100 g dry weight)
#'  \item MG      : Magnesium (mg/100 g dry weight)
#'  \item BC      : Beta-carotene (mg/100 g dry weight)
#'  \item BC.CC   : Beta-carotene with color charts
#'  \item TC      : Total carotenoids (mg/100 g dry weight)
#'  \item STAR    : Starch (\%)
#'  \item FRUC    : Fructose (\%)
#'  \item GLUC    : Glucose (\%)
#'  \item SUCR    : Sucrose (\%)
#'  \item MALT    : Maltose (\%)
#'  \item TRW     : Total root weight
#'  \item CYTHA   : Commercial root yield t/ha
#'  \item RYTHA   : Total root yield t/ha
#'  \item ACRW    : Average commercial root weight = CRW / NOCR
#'  \item NRPP    : Number of roots per plant
#'  \item YPP     : Yield per plant Kg
#'  \item CI      : Percent marketable roots (commercial index)
#'  \item HI      : Harvest index
#'  \item SHI     : Harvest sowing index (survival)
#'  \item BIOM    : Biomass yield
#'  \item FYTHA   : Foliage total yield t/ha
#'  \item RFR     : Root foliage ratio
#'  }
#' @return It returns a data frame with all traits names in upper case, and a list of the traits
#' with names not included in the list shown above.
#' @author Raul Eyzaguirre.
#' @examples
#'  # The data
#'  head(pjpz09)
#'  str(pjpz09)
#'
#'  # Check the trait names
#'  checknames(pjpz09)
#' @export

checknames <- function(fb) {


  colnames.valid <- c("L", "LOC", "Y", "S", "G", "GENO", "NAME", "E", "ENV", "R", "REP", "NOPS",
                      "NOPE", "VIR1", "VIR2", "VIR3", "ALT1", "ALT2", "VV1", "VV2", "VW", "NOPH",
                      "NOPR", "NOCR", "NONC", "CRW", "NCRW", "RFCP", "RFCS", "SCOL", "FCOL", "RFCP",
                      "RFCS", "RS", "RF", "DAMR", "RSPR", "WED1", "WED2", "DMF", "DMD", "DM", "DMRY",
                      "DMVF", "DMVD", "DMV", "DMFY", "FRAW1", "SURAW1", "STRAW1", "COOF1", "COOSU1",
                      "COOST1", "COOT1", "COOAP1", "FRAW2", "SURAW2", "STRAW2", "COOF2", "COOSU2",
                      "COOST2", "COOT2", "COOAP2", "PROT", "FE", "ZN", "CA", "MG", "BC", "BC.CC",
                      "TC", "STAR", "FRUC", "GLUC", "SUCR", "MALT", "TRW", "CYTHA", "RYTHA", "ACRW",
                      "NRPP", "YPP", "CI", "HI", "SHI", "BIOM", "FYTHA", "RFR")
    
  colnames.list <- colnames(fb)
  
  check.list.1 <- !(toupper(colnames.list) %in% colnames.valid) # which are not valid
  temp <- colnames.list[!check.list.1]                          # list of valid names
  check.list.2 <- !(temp %in% colnames.valid)                   # which are valid but lower case
  
  colnames(fb) <- toupper(colnames(fb))
    
  # Warnings
  
  if (max(check.list.1) == 1)
    warning("Invalid labels not included for checking: ", list(colnames.list[check.list.1]), call. = FALSE)
  
  if (max(check.list.2) == 1)
    warning("Some labels converted to upper case: ", list(temp[check.list.2]), call. = FALSE)
  
  fb
}
