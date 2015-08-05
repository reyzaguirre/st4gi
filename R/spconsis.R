#' Check consistency for sweetpotato experimental data
#'
#' Set of rules to check for consistency of sweetpotato experimental data.
#' Data labels must be defined as specified in the PROCEDURES FOR THE EVALUATION
#' AND ANALYSIS OF SWEETPOTATO TRIALS document.
#' @param data The name of the data frame.
#' @param plot.size Plot size in square meters.
#' @param width Number of columns for the output file.
#' @details The data frame must use the following labels:
#' \itemize{
#'  \item L       : Locations
#'  \item Y       : Year
#'  \item S       : Season
#'  \item G       : Genotypes
#'  \item NAME    : Names for genotypes
#'  \item E       : Environments
#'  \item R       : Replications or blocks
#'  \item NOPS    : Number of plants planted
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
#'  \item DMVF    : Fresh weight vines for dry matter assessment
#'  \item DMVD    : Dry weight of DMVF samples
#'  \item DM      : Storage root dry matter content (%)
#'  \item DMFY    : Dry matter foliage yield
#'  \item DMRY    : Dry matter root yield
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
#'  \item PROT    : Protein
#'  \item FE      : Iron in dry weight
#'  \item ZN      : Zinc in dry weight
#'  \item CA      : Calcium in dry weight
#'  \item MG      : Magnesium in dry weight
#'  \item BC      : Beta-carotene in dry weight (NIRS)
#'  \item BC.CC   : Beta-carotene with color charts
#'  \item TC      : Total carotenoids  in dry weight (NIRS)
#'  \item STAR    : Starch
#'  \item FRUC    : Fructose
#'  \item GLUC    : Glucose
#'  \item SUCR    : Sucrose
#'  \item MALT    : Maltose
#'  \item TRW     : Total root weight
#'  \item CYTHA   : Commercial root yield t/ ha
#'  \item RYTHA   : Total root yield t/ha
#'  \item ACRW    : Average commercial root weight = CRW/NOCR
#'  \item NRPP    : Number of roots per plant
#'  \item YPP     : Yield per plant Kg
#'  \item CI      : Percent marketable roots (commercial index)
#'  \item HI      : Harvest index
#'  \item SHI     : Harvest sowing index  (survival)
#'  \item BIOM    : Biomass yield
#'  \item FYTHA   : Foliage total yield t/ha
#'  \item RFR     : Root foliage ratio
#'  }
#' @return It returns a file with name checks.txt with a list of
#' all rows with some kind of inconsistency and all rows with outliers.
#' @author Raul Eyzaguirre.
#' @examples
#'  # The data
#'  head(pjpz09)
#'  str(pjpz09)
#'
#'  # Check the data
#'  spconsis(pjpz09, 4.5)
#' @export

spconsis <- function(data, plot.size, width = 240){

  options(width = width)
  
  colnames.valid <- c("L", "Y", "S", "G", "NAME", "E", "R", "NOPS", "NOPE", "VIR1", "VIR2", "VIR3",
                      "ALT1", "ALT2", "VV1", "VV2", "VW", "NOPH", "NOPR", "NOCR", "NONC", "CRW",
                      "NCRW", "RFCP", "RFCS", "SCOL", "FCOL", "RFCP", "RFCS", "RS", "RF", "DAMR",
                      "RSPR", "WED1", "WED2", "DMF", "DMD", "DMVF", "DMVD", "DM", "DMFY", "DMRY",
                      "FRAW1", "SURAW1", "STRAW1", "COOF1", "COOSU1", "COOST1", "COOT1", "COOAP1",
                      "FRAW2", "SURAW2", "STRAW2", "COOF2", "COOSU2", "COOST2", "COOT2", "COOAP2",
                      "PROT", "FE", "ZN", "CA", "MG", "BC", "BC.CC", "TC", "STAR", "FRUC", "GLUC",
                      "SUCR", "MALT", "TRW", "CYTHA", "RYTHA", "ACRW", "NRPP", "YPP", "CI", "HI",
                      "SHI", "BIOM", "FYTHA", "RFR")
    
  colnames.list <- colnames(data)
  
  check.list.1 <- !(toupper(colnames.list) %in% colnames.valid)
  temp <- colnames.list[!check.list.1]
  check.list.2 <- !(temp %in% colnames.valid)
  
  colnames(data) <- toupper(colnames(data))
    
  # Warnings
  
  if (max(check.list.1) == 1)
    warning("Invalid labels not included for checking :", list(colnames.list[check.list.1]),
            call. = FALSE)
  
  if (max(check.list.2) == 1)
    warning("Some labels converted to upper case: ", list(temp[check.list.2]), call. = FALSE)
  
  sink("checks.txt")

  spconsis01(data) # NOPS > NOPE > NOPH > NOPR
  spconsis02(data) # NOPE and dependencies
  spconsis03(data) # NOPH and VW
  spconsis04(data) # VW and dependencies
  spconsis05(data) # NOPR and number of roots
  spconsis06(data) # Number of roots and root weight
  spconsis07(data) # TRW, CRW+NCRW, NOCR+NONC, NOPR
  spconsis08(data) # Roots and dependencies
  spconsis09(data, plot.size) # Calculated variables
  spconsis10(data) # Outliers detection and values out of range for field data
  spconsis11(data) # Outliers detection and values out of range for DM data
  spconsis12(data) # Outliers detection and values out of range for cooked traits
  spconsis13(data) # Outliers detection and values out of range for lab data
  spconsis14(data) # Outliers detection and values out of range for derived variables

  sink()
}

# Check consistency for sweetpotato experimental data, part 1.
# Inconsistencies for NOPS > NOPE > NOPH > NOPR.

spconsis01 <- function(data){

  if (exists("NOPE", where = data) & exists("NOPS", where=data))
    if (dim(subset(data, NOPE>NOPS))[1]>0){
      cat("\n", "-  Number of plants established (NOPE) greater than number of plants sowed (NOPS):", "\n")
      print(subset(data, NOPE>NOPS))
    }

  if (exists("NOPH", where=data) & exists("NOPE", where=data))
    if (dim(subset(data, NOPH>NOPE))[1]>0){
      cat("\n", "-  Number of plants harvested (NOPH) greater than number of plants established (NOPE):", "\n")
      print(subset(data, NOPH>NOPE))
    }

  if (exists("NOPH", where=data) & exists("NOPS", where=data)){
    if (exists("NOPE", where=data)){
      if (dim(subset(data, is.na(NOPE) & NOPH>NOPS))[1]>0){
        cat("\n", "- Number of plants harvested (NOPH) greater than number of plants sowed (NOPS):", "\n")
        print(subset(data, is.na(NOPE) & NOPH>NOPS))
      }
    } else {
      if (dim(subset(data, NOPH>NOPS))[1]>0){
        cat("\n", "- Number of plants harvested (NOPH) greater than number of plants sowed (NOPS):", "\n")
        print(subset(data, NOPH>NOPS))
      }
    }
  }

  if (exists("NOPR", where=data) & exists("NOPH", where=data))
    if (dim(subset(data, NOPR>NOPH))[1]>0){
      cat("\n", "- Number of plants with roots (NOPR) greater than number of plants harvested (NOPH):", "\n")
      print(subset(data, NOPR>NOPH))
    }

  if (exists("NOPR", where=data) & exists("NOPE", where=data)){
    if (exists("NOPH", where=data)){
      if (dim(subset(data, is.na(NOPH) & NOPR>NOPE))[1]>0){
        cat("\n", "- Number of plants with roots (NOPR) greater than number of plants established (NOPE):", "\n")
        print(subset(data, is.na(NOPH) & NOPR>NOPE))
      }
    } else {
      if (dim(subset(data, NOPR>NOPE))[1]>0){
        cat("\n", "- Number of plants with roots (NOPR) greater than number of plants established (NOPE):", "\n")
        print(subset(data, NOPR>NOPE))
      }
    }
  }

  if (exists("NOPR", where=data) & exists("NOPS", where=data)){
    if (exists("NOPH", where=data) & exists("NOPE", where=data)){
      if (dim(subset(data, is.na(NOPH) & is.na(NOPE) & NOPR>NOPS))[1]>0){
        cat("\n", "- Number of plants with roots (NOPR) greater than number of plants sowed (NOPS):", "\n")
        print(subset(data, is.na(NOPH) & is.na(NOPE) & NOPR>NOPS))
      }
    } else {
      if (dim(subset(data, NOPR>NOPS))[1]>0){
        cat("\n", "- Number of plants with roots (NOPR) greater than number of plants sowed (NOPS):", "\n")
        print(subset(data, NOPR>NOPS))
      }
    }
  }
}

# Check consistency for sweetpotato experimental data, part 2.
# Inconsistencies for NOPE and dependencies.

spconsis02 <- function(data){

  if (exists("NOPE", where=data)==1 & exists("VIR1", where=data)==1)
    if (dim(subset(data, (NOPE==0 | is.na(NOPE)==1) & is.na(VIR1)==0))[1]>0){
      cat("\n", "- Number of plants established (NOPE) is zero or NA but there is data for virus symptoms first evaluation (VIR1):", "\n")
      print(subset(data, (NOPE==0 | is.na(NOPE)==1) & is.na(VIR1)==0))
    }

  if (exists("NOPE", where=data)==1 & exists("VIR2", where=data)==1)
    if (dim(subset(data, (NOPE==0 | is.na(NOPE)==1) & is.na(VIR2)==0))[1]>0){
      cat("\n", "- Number of plants established (NOPE) is zero or NA but there is data for virus symptoms second evaluation (VIR2):", "\n")
      print(subset(data, (NOPE==0 | is.na(NOPE)==1) & is.na(VIR2)==0))
    }

  if (exists("NOPE", where=data)==1 & exists("ALT1", where=data)==1)
    if (dim(subset(data, (NOPE==0 | is.na(NOPE)==1) & is.na(ALT1)==0))[1]>0){
      cat("\n", "- Number of plants established (NOPE) is zero or NA but there is data for alternaria symptoms first evaluation (ALT1):", "\n")
      print(subset(data, (NOPE==0 | is.na(NOPE)==1) & is.na(ALT1)==0))
    }

  if (exists("NOPE", where=data)==1 & exists("ALT2", where=data)==1)
    if (dim(subset(data, (NOPE==0 | is.na(NOPE)==1) & is.na(ALT2)==0))[1]>0){
      cat("\n", "- Number of plants established (NOPE) is zero or NA but there is data for alternaria symptoms second evaluation (ALT2):", "\n")
      print(subset(data, (NOPE==0 | is.na(NOPE)==1) & is.na(ALT2)==0))
    }

  if (exists("NOPE", where=data)==1 & exists("VV1", where=data)==1)
    if (dim(subset(data, (NOPE==0 | is.na(NOPE)==1) & is.na(VV1)==0))[1]>0){
      cat("\n", "- Number of plants established (NOPE) is zero or NA but there is data for vine vigor first evaluation (VV1):", "\n")
      print(subset(data, (NOPE==0 | is.na(NOPE)==1) & is.na(VV1)==0))
    }
}

# Check consistency for sweetpotato experimental data, part 3.
# Inconsistencies for NOPH and VW.

spconsis03 <- function(data){

  if (exists("NOPH", where=data)==1 & exists("VW", where=data)==1)
    if (dim(subset(data, (NOPH==0 | is.na(NOPH)==1) & VW>0))[1]>0){
      cat("\n", "- Number of plants harvested (NOPH) is zero or NA but vine weight (VW) is greater than zero:", "\n")
      print(subset(data, (NOPH==0 | is.na(NOPH)==1) & VW>0))
    }

  if (exists("NOPH", where=data)==1 & exists("VW", where=data)==1)
    if (dim(subset(data, NOPH>0 & (VW==0 | is.na(VW)==1)))[1]>0){
      cat("\n", "- Vine weight (VW) is zero or NA but the number of plants harvested (NOPH) is greater than zero:", "\n")
      print(subset(data, NOPH>0 & (VW==0 | is.na(VW)==1)))
    }
}

# Check consistency for sweetpotato experimental data, part 4.
# Inconsistencies for VW and dependencies.

spconsis04 <- function(data){

  if (exists("VW", where=data)==1 & exists("DMVF", where=data)==1)
    if (dim(subset(data, (VW==0 | is.na(VW)==1) & DMVF>0))[1]>0){
      cat("\n", "- Vine weight (VW) is zero or NA but there is fresh weight vines for dry matter assessment (DMVF):", "\n")
      print(subset(data, (VW==0 | is.na(VW)==1) & DMVF>0))
    }

  if (exists("VW", where=data)==1 & exists("DMVD", where=data)==1)
    if (dim(subset(data, (VW==0 | is.na(VW)==1) & DMVD>0))[1]>0){
      cat("\n", "- Vine weight (VW) is zero or NA but there is dry weight vines for dry matter assessment (DMVD):", "\n")
      print(subset(data, (VW==0 | is.na(VW)==1) & DMVD>0))
    }

  if (exists("DMVF", where=data)==1 & exists("DMVD", where=data)==1)
    if (dim(subset(data, DMVD>DMVF))[1]>0){
      cat("\n", "- Dry weight vines for dry matter assessment (DMVD) is greater than fresh weight vines for dry matter assessment (DBVF):", "\n")
      print(subset(data, DMVD>DMVF))
    }

  if (exists("VW", where=data)==1 & exists("VV2", where=data)==1)
    if (dim(subset(data, (VW==0 | is.na(VW)==1) & is.na(VV2)==0))[1]>0){
      cat("\n", "- Vine weight (VW) is zero or NA but there is data for vine vigor second evaluation (VV2):", "\n")
      print(subset(data, (VW==0 | is.na(VW)==1) & is.na(VV2)==0))
    }

  if (exists("VW", where=data)==1 & exists("VIR3", where=data)==1)
    if (dim(subset(data, (VW==0 | is.na(VW)==1) & is.na(VIR3)==0))[1]>0){
      cat("\n", "- Vine weight (VW) is zero or NA but there is data for virus symptoms third evaluation (VIR3):", "\n")
      print(subset(data, (VW==0 | is.na(VW)==1) & is.na(VIR3)==0))
    }
}

# Check consistency for sweetpotato experimental data, part 5.
# Inconsistencies for NOPR and number of roots.

spconsis05 <- function(data){

  if (exists("NOPR", where=data)==1 & exists("NOCR", where=data)==1 & exists("NONC", where=data)==1)
    if (dim(subset(data, (NOPR==0 | is.na(NOPR)==1) & (NOCR>0 | NONC>0)))[1]>0){
      cat("\n", "- Number of plants with roots (NOPR) is zero or NA but number of roots (NOCR+NONC) is greater than zero:", "\n")
      print(subset(data, (NOPR==0 | is.na(NOPR)==1) & (NOCR>0 | NONC>0)))
    }

  if (exists("NOPR", where=data)==1 & exists("NOCR", where=data)==1 & exists("NONC", where=data)==1)
    if (dim(subset(data, NOPR>0 & ((NOCR+NONC)==0 | (NOCR==0 & is.na(NONC)==1) | (is.na(NOCR)==1 & NONC==0) |
                                     (is.na(NOCR)==1 & is.na(NONC)==1))))[1]>0){
      cat("\n", "- Number of roots (NOCR+NONC) is zero or NA but number of plants with roots (NOPR) is greater than zero:", "\n")
      print(subset(data, NOPR>0 & ((NOCR+NONC)==0 | (NOCR==0 & is.na(NONC)==1) | (is.na(NOCR)==1 & NONC==0) |
                                     (is.na(NOCR)==1 & is.na(NONC)==1))))
    }
}

# Check consistency for sweetpotato experimental data, part 6.
# Inconsistencies for number of roots and root weight.

spconsis06 <- function(data){

  if (exists("NOCR", where=data)==1 & exists("CRW", where=data)==1)
    if (dim(subset(data, (NOCR==0 | is.na(NOCR)==1) & CRW>0))[1]>0){
      cat("\n", "- Number of commercial roots (NOCR) is zero or NA but the commercial root weight (CRW) is greater than zero:", "\n")
      print(subset(data, (NOCR==0 | is.na(NOCR)==1) & CRW>0))
    }

  if (exists("NOCR", where=data)==1 & exists("CRW", where=data)==1)
    if (dim(subset(data, NOCR>0 & (CRW==0 | is.na(CRW)==1)))[1]>0){
      cat("\n", "- Commercial root weight (CRW) is zero or NA but the number of commercial roots (NOCR) is greater than zero:", "\n")
      print(subset(data, NOCR>0 & (CRW==0 | is.na(CRW)==1)))
    }

  if (exists("NONC", where=data)==1 & exists("NCRW", where=data)==1)
    if (dim(subset(data, (NONC==0 | is.na(NONC)==1) & NCRW>0))[1]>0){
      cat("\n", "- Number of non commercial roots (NONC) is zero or NA but the non commercial root weight (NCRW) is greater than zero:", "\n")
      print(subset(data, (NONC==0 | is.na(NONC)==1) & NCRW>0))
    }

  if (exists("NONC", where=data)==1 & exists("NCRW", where=data)==1)
    if (dim(subset(data, NONC>0 & (NCRW==0 | is.na(NCRW)==1)))[1]>0){
      cat("\n", "- Non commercial root weight (NCRW) is zero or NA but the number of non commercial roots (NONC) is greater than zero:", "\n")
      print(subset(data, NONC>0 & (NCRW==0 | is.na(NCRW)==1)))
    }
}

# Check consistency for sweetpotato experimental data, part 7.
# Inconsistencies for TRW, CRW+NCRW, NOCR+NONC, NOPR.

spconsis07 <- function(data){

  if (exists("TRW", where=data)==1 & exists("NOPR", where=data)==1)
    if (dim(subset(data, (TRW==0 | is.na(TRW)==1) & NOPR>0))[1]>0){
      cat("\n", "- Total root weight (TRW) is zero or NA but number of plants with roots (NOPR) is greater than zero:", "\n")
      print(subset(data, (TRW==0 | is.na(TRW)==1) & NOPR>0))
    }

  if (exists("TRW", where=data)==1 & exists("CRW", where=data)==1 & exists("NCRW", where=data)==1)
    if (dim(subset(data, (TRW==0 | is.na(TRW)==1) & (CRW>0 | NCRW>0)))[1]>0){
      cat("\n", "- Total root weight (TRW) is zero or NA but root weight (CRW+NCRW) is greater than zero:", "\n")
      print(subset(data, (TRW==0 | is.na(TRW)==1) & (CRW>0 | NCRW>0)))
    }

  if (exists("TRW", where=data)==1 & exists("NOCR", where=data)==1 & exists("NONC", where=data)==1)
    if (dim(subset(data, (TRW==0 | is.na(TRW)==1) & (NOCR>0 | NONC>0)))[1]>0){
      cat("\n", "- Total root weight (TRW) is zero or NA but number of roots (NOCR+NONC) is greater than zero:", "\n")
      print(subset(data, (TRW==0 | is.na(TRW)==1) & (NOCR>0 | NONC>0)))
    }

  if (exists("CRW", where=data)==1 & exists("NCRW", where=data)==1 & exists("NOPR", where=data)==1)
    if (dim(subset(data, NOPR>0 & ((CRW+NCRW)==0 | (CRW==0 & is.na(NCRW)==1) | (is.na(CRW)==1 & NCRW==0) |
                                     (is.na(CRW)==1 & is.na(NCRW)==1))))[1]>0){
      cat("\n", "- Root weight (CRW+NCRW) is zero or NA but number of plants with roots (NOPR) is greater than zero:", "\n")
      print(subset(data, NOPR>0 & ((CRW+NCRW)==0 | (CRW==0 & is.na(NCRW)==1) | (is.na(CRW)==1 & NCRW==0) |
                                     (is.na(CRW)==1 & is.na(NCRW)==1))))
    }
}

# Check consistency for sweetpotato experimental data, part 8.
# Inconsistencies for roots and dependencies.

spconsis08 <- function(data){

  if (exists("NOPR", where=data)==1)
    data$RAUX <- data$NOPR else
      if (exists("NOCR", where=data)==1 & exists("NONC", where=data)==1)
        data$RAUX <- apply(cbind(data$NOCR, data$NONC), 1, sum, na.rm=T) else
          if (exists("CRW", where=data)==1 & exists("NCRW", where=data)==1)
            data$RAUX <- apply(cbind(data$CRW, data$NCRW), 1, sum, na.rm=T) else
              if (exists("TRW", where=data)==1)
                data$RAUX <- data$TRW else
                  if (exists("RYTHA", where=data)==1)
                    data$RAUX <- data$RYTHA else
                      if (exists("CRW", where=data)==1)
                        data$RAUX <- data$CRW else
                          if (exists("CYTHA", where=data)==1)
                            data$RAUX <- data$CYTHA

  if (exists("RAUX", where=data)==1 & exists("RFCP", where=data)==1)
    if (dim(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(RFCP)==0))[1]>0){
      cat("\n", "- There are no roots but there is data for root primary flesh color (RFCP):", "\n")
      print(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(RFCP)==0, select = -RAUX))
    }

  if (exists("RAUX", where=data)==1 & exists("RFCS", where=data)==1)
    if (dim(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(RFCS)==0))[1]>0){
      cat("\n", "- There are no roots but there is data for root secondary flesh color (RFCS):", "\n")
      print(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(RFCS)==0, select = -RAUX))
    }

  if (exists("RAUX", where=data)==1 & exists("SCOL", where=data)==1)
    if (dim(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(SCOL)==0))[1]>0){
      cat("\n", "- There are no roots but there is data for storage root skin color (SCOL):", "\n")
      print(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(SCOL)==0, select = -RAUX))
    }

  if (exists("RAUX", where=data)==1 & exists("FCOL", where=data)==1)
    if (dim(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(FCOL)==0))[1]>0){
      cat("\n", "- There are no roots but there is data for storage root flesh color (FCOL):", "\n")
      print(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(FCOL)==0, select = -RAUX))
    }

  if (exists("RAUX", where=data)==1 & exists("RS", where=data)==1)
    if (dim(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(RS)==0))[1]>0){
      cat("\n", "- There are no roots but there is data for root size (RS):", "\n")
      print(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(RS)==0, select = -RAUX))
    }

  if (exists("RAUX", where=data)==1 & exists("RF", where=data)==1)
    if (dim(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(RF)==0))[1]>0){
      cat("\n", "- There are no roots but there is data for root form (RF):", "\n")
      print(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(RF)==0, select = -RAUX))
    }

  if (exists("RAUX", where=data)==1 & exists("DAMR", where=data)==1)
    if (dim(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(DAMR)==0))[1]>0){
      cat("\n", "- There are no roots but there is data for root defects (DAMR):", "\n")
      print(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(DAMR)==0, select = -RAUX))
    }

  if (exists("RAUX", where=data)==1 & exists("RSPR", where=data)==1)
    if (dim(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(RSPR)==0))[1]>0){
      cat("\n", "- There are no roots but there is data for root sprouting (RSPR):", "\n")
      print(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(RSPR)==0, select = -RAUX))
    }

  if (exists("RAUX", where=data)==1 & exists("WED1", where=data)==1)
    if (dim(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(WED1)==0))[1]>0){
      cat("\n", "- There are no roots but there is data for weevil damage first evaluation (WED1):", "\n")
      print(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(WED1)==0, select = -RAUX))
    }

  if (exists("RAUX", where=data)==1 & exists("WED2", where=data)==1)
    if (dim(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(WED2)==0))[1]>0){
      cat("\n", "- There are no roots but there is data for weevil damage second evaluation (WED2):", "\n")
      print(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(WED2)==0, select = -RAUX))
    }

  if (exists("RAUX", where=data)==1 & exists("DMF", where=data)==1)
    if (dim(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(DMF)==0))[1]>0){
      cat("\n", "- There are no roots but there is data for fresh weight of roots for dry matter assessment (DMF):", "\n")
      print(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(DMF)==0, select = -RAUX))
    }

  if (exists("RAUX", where=data)==1 & exists("DMD", where=data)==1)
    if (dim(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(DMD)==0))[1]>0){
      cat("\n", "- There are no roots but there is data for dry weight of roots for dry matter assessment (DMD):", "\n")
      print(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(DMD)==0, select = -RAUX))
    }

  if (exists("DMF", where=data)==1 & exists("DMD", where=data)==1)
    if (dim(subset(data, DMF<DMD))[1]>0){
      cat("\n", "- Dry weight of roots for dry matter assessment (DMD) is greater than fresh weight of roots for dry matter assessment (DMF):", "\n")
      print(subset(data, DMF<DMD, select = -RAUX))
    }

  if (exists("RAUX", where=data)==1 & exists("FRAW1", where=data)==1)
    if (dim(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(FRAW1)==0))[1]>0){
      cat("\n", "- There are no roots but there is data for root fiber first determination (FRAW1):", "\n")
      print(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(FRAW1)==0, select = -RAUX))
    }

  if (exists("RAUX", where=data)==1 & exists("SURAW1", where=data)==1)
    if (dim(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(SURAW1)==0))[1]>0){
      cat("\n", "- There are no roots but there is data for root sugar first determination (SURAW1):", "\n")
      print(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(SURAW1)==0, select = -RAUX))
    }

  if (exists("RAUX", where=data)==1 & exists("STRAW1", where=data)==1)
    if (dim(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(STRAW1)==0))[1]>0){
      cat("\n", "- There are no roots but there is data for root starch first determination (STRAW1):", "\n")
      print(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(STRAW1)==0, select = -RAUX))
    }

  if (exists("RAUX", where=data)==1 & exists("COOF1", where=data)==1)
    if (dim(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(COOF1)==0))[1]>0){
      cat("\n", "- There are no roots but there is data for cooked fiber first evaluation (COOF1):", "\n")
      print(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(COOF1)==0, select = -RAUX))
    }

  if (exists("RAUX", where=data)==1 & exists("COOSU1", where=data)==1)
    if (dim(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(COOSU1)==0))[1]>0){
      cat("\n", "- There are no roots but there is data for cooked sugars first evaluation (COOSU1):", "\n")
      print(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(COOSU1)==0, select = -RAUX))
    }

  if (exists("RAUX", where=data)==1 & exists("COOST1", where=data)==1)
    if (dim(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(COOST1)==0))[1]>0){
      cat("\n", "- There are no roots but there is data for cooked starch first evaluation (COOST1):", "\n")
      print(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(COOST1)==0, select = -RAUX))
    }

  if (exists("RAUX", where=data)==1 & exists("COOT1", where=data)==1)
    if (dim(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(COOT1)==0))[1]>0){
      cat("\n", "- There are no roots but there is data for cooked taste first evaluation (COOT1):", "\n")
      print(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(COOT1)==0, select = -RAUX))
    }

  if (exists("RAUX", where=data)==1 & exists("COOAP1", where=data)==1)
    if (dim(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(COOAP1)==0))[1]>0){
      cat("\n", "- There are no roots but there is data for cooked appearance first evaluation (COOAP1):", "\n")
      print(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(COOAP1)==0, select = -RAUX))
    }

  if (exists("RAUX", where=data)==1 & exists("FRAW2", where=data)==1)
    if (dim(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(FRAW2)==0))[1]>0){
      cat("\n", "- There are no roots but there is data for root fiber second determination (FRAW2):", "\n")
      print(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(FRAW2)==0, select = -RAUX))
    }

  if (exists("RAUX", where=data)==1 & exists("SURAW2", where=data)==1)
    if (dim(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(SURAW2)==0))[1]>0){
      cat("\n", "- There are no roots but there is data for root sugar second determination (SURAW2):", "\n")
      print(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(SURAW2)==0, select = -RAUX))
    }

  if (exists("RAUX", where=data)==1 & exists("STRAW2", where=data)==1)
    if (dim(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(STRAW2)==0))[1]>0){
      cat("\n", "- There are no roots but there is data for root starch second determination (STRAW2):", "\n")
      print(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(STRAW2)==0, select = -RAUX))
    }

  if (exists("RAUX", where=data)==1 & exists("COOF2", where=data)==1)
    if (dim(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(COOF2)==0))[1]>0){
      cat("\n", "- There are no roots but there is data for cooked fiber second evaluation (COOF2):", "\n")
      print(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(COOF2)==0, select = -RAUX))
    }

  if (exists("RAUX", where=data)==1 & exists("COOSU2", where=data)==1)
    if (dim(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(COOSU2)==0))[1]>0){
      cat("\n", "- There are no roots but there is data for cooked sugars second evaluation (COOSU2):", "\n")
      print(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(COOSU2)==0, select = -RAUX))
    }

  if (exists("RAUX", where=data)==1 & exists("COOST2", where=data)==1)
    if (dim(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(COOST2)==0))[1]>0){
      cat("\n", "- There are no roots but there is data for cooked starch second evaluation (COOST2):", "\n")
      print(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(COOST2)==0, select = -RAUX))
    }

  if (exists("RAUX", where=data)==1 & exists("COOT2", where=data)==1)
    if (dim(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(COOT2)==0))[1]>0){
      cat("\n", "- There are no roots but there is data for cooked taste second evaluation (COOT2):", "\n")
      print(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(COOT2)==0, select = -RAUX))
    }

  if (exists("RAUX", where=data)==1 & exists("COOAP2", where=data)==1)
    if (dim(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(COOAP2)==0))[1]>0){
      cat("\n", "- There are no roots but there is data for cooked appearance second evaluation (COOAP2):", "\n")
      print(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(COOAP2)==0, select = -RAUX))
    }

  if (exists("RAUX", where=data)==1 & exists("PROT", where=data)==1)
    if (dim(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(PROT)==0))[1]>0){
      cat("\n", "- There are no roots but there is data for protein (PROT):", "\n")
      print(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(PROT)==0, select = -RAUX))
    }

  if (exists("RAUX", where=data)==1 & exists("FE", where=data)==1)
    if (dim(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(FE)==0))[1]>0){
      cat("\n", "- There are no roots but there is data for iron in dry weight (FE):", "\n")
      print(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(FE)==0, select = -RAUX))
    }

  if (exists("RAUX", where=data)==1 & exists("ZN", where=data)==1)
    if (dim(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(ZN)==0))[1]>0){
      cat("\n", "- There are no roots but there is data for zinc in dry weight (ZN):", "\n")
      print(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(ZN)==0, select = -RAUX))
    }

  if (exists("RAUX", where=data)==1 & exists("CA", where=data)==1)
    if (dim(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(CA)==0))[1]>0){
      cat("\n", "- There are no roots but there is data for calcium in dry weight (CA):", "\n")
      print(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(CA)==0, select = -RAUX))
    }

  if (exists("RAUX", where=data)==1 & exists("MG", where=data)==1)
    if (dim(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(MG)==0))[1]>0){
      cat("\n", "- There are no roots but there is data for magnesium in dry weight (MG):", "\n")
      print(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(MG)==0, select = -RAUX))
    }

  if (exists("RAUX", where=data)==1 & exists("BC", where=data)==1)
    if (dim(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(BC)==0))[1]>0){
      cat("\n", "- There are no roots but there is data for beta-carotene in dry weight (BC):", "\n")
      print(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(BC)==0, select = -RAUX))
    }

  if (exists("RAUX", where=data)==1 & exists("TC", where=data)==1)
    if (dim(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(TC)==0))[1]>0){
      cat("\n", "- There are no roots but there is data for total carotenoids in dry weight (TC):", "\n")
      print(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(TC)==0, select = -RAUX))
    }

  if (exists("RAUX", where=data)==1 & exists("STAR", where=data)==1)
    if (dim(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(STAR)==0))[1]>0){
      cat("\n", "- There are no roots but there is data for starch (STAR):", "\n")
      print(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(STAR)==0, select = -RAUX))
    }

  if (exists("RAUX", where=data)==1 & exists("FRUC", where=data)==1)
    if (dim(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(FRUC)==0))[1]>0){
      cat("\n", "- There are no roots but there is data for fructose (FRUC):", "\n")
      print(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(FRUC)==0, select = -RAUX))
    }

  if (exists("RAUX", where=data)==1 & exists("GLUC", where=data)==1)
    if (dim(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(GLUC)==0))[1]>0){
      cat("\n", "- There are no roots but there is data for glucose (GLUC):", "\n")
      print(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(GLUC)==0, select = -RAUX))
    }

  if (exists("RAUX", where=data)==1 & exists("SUCR", where=data)==1)
    if (dim(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(SUCR)==0))[1]>0){
      cat("\n", "- There are no roots but there is data for sucrose (SUCR):", "\n")
      print(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(SUCR)==0, select = -RAUX))
    }

  if (exists("RAUX", where=data)==1 & exists("MALT", where=data)==1)
    if (dim(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(MALT)==0))[1]>0){
      cat("\n", "- There are no roots but there is data for maltose (MALT):", "\n")
      print(subset(data, (RAUX==0 | is.na(RAUX)==1) & is.na(MALT)==0, select = -RAUX))
    }

  if (exists("RAUX", where=data)==1 & exists("TRW", where=data)==1)
    if (dim(subset(data, (RAUX==0 | is.na(RAUX)==1) & TRW>0))[1]>0){
      cat("\n", "- There are no roots but total root weight (TRW) is greater than zero:", "\n")
      print(subset(data, (RAUX==0 | is.na(RAUX)==1) & TRW>0, select = -RAUX))
    }

  if (exists("RAUX", where=data)==1 & exists("RYTHA", where=data)==1)
    if (dim(subset(data, (RAUX==0 | is.na(RAUX)==1) & RYTHA>0))[1]>0){
      cat("\n", "- There are no roots but total root yield (RYTHA) is greater than zero:", "\n")
      print(subset(data, (RAUX==0 | is.na(RAUX)==1) & RYTHA>0, select = -RAUX))
    }

  data$RAUX <- NULL
}

# Check consistency for sweetpotato experimental data, part 9.
# Inconsistencies for calculated variables.

spconsis09 <- function(data, plot.size){

  if (exists("TRW", where=data)==1 & exists("CRW", where=data)==1 & exists("NCRW", where=data)==1)
    if (dim(subset(data, abs(TRW-apply(cbind(data$CRW, data$NCRW), 1, sum, na.rm=T))>1e-10))[1]>0){
      cat("\n", "- Total root weight (TRW) different from CRW+NCRW:", "\n")
      print(subset(data, abs(TRW-apply(cbind(data$CRW, data$NCRW), 1, sum, na.rm=T))>1e-10))
    }

  if (exists("CYTHA", where=data)==1 & exists("CRW", where=data)==1)
    if (dim(subset(data, abs(CYTHA-CRW*10/plot.size)>1e-10))[1]>0){
      cat("\n", "- Commercial root yield in tons per hectare (CYTHA) is different from CRW*10/plot.size:", "\n")
      print(subset(data, abs(CYTHA-CRW*10/plot.size)>1e-10))
    }

  if (exists("RYTHA", where=data)==1 & exists("CRW", where=data)==1 & exists("NCRW", where=data)==1)
    if (dim(subset(data, abs(RYTHA-apply(cbind(CRW,NCRW), 1, sum, na.rm=T)*10/plot.size)>1e-10))[1]>0){
      cat("\n", "- Total root yield in tons per hectare (RYTHA) is different from (CRW+NCRW)*10/plot.size:", "\n")
      print(subset(data, abs(RYTHA-apply(cbind(CRW,NCRW), 1, sum, na.rm=T)*10/plot.size)>1e-10))
    }

  if (exists("ACRW", where=data)==1 & exists("CRW", where=data)==1 & exists("NOCR", where=data)==1)
    if (dim(subset(data, abs(ACRW-CRW/NOCR)>1e-10))[1]>0){
      cat("\n", "- Average commercial root weight (ACRW) is different from CRW/NOCR:", "\n")
      print(subset(data, abs(ACRW-CRW/NOCR)>1e-10))
    }

  if (exists("NRPP", where=data)==1 & exists("NOCR", where=data)==1 & exists("NONC", where=data)==1 & exists("NOPH", where=data)==1)
    if (dim(subset(data, abs(NRPP-apply(cbind(NOCR,NONC), 1, sum, na.rm=T)/NOPH)>1e-10))[1]>0){
      cat("\n", "- Number of roots per plant (NRPP) is different from (NOCR+NONC)/NOPH:", "\n")
      print(subset(data, abs(NRPP-apply(cbind(NOCR,NONC), 1, sum, na.rm=T)/NOPH)>1e-10))
    }

  if (exists("YPP", where=data)==1 & exists("CRW", where=data)==1 & exists("NCRW", where=data)==1 & exists("NOPH", where=data)==1)
    if (dim(subset(data, abs(YPP-apply(cbind(CRW, NCRW), 1, sum, na.rm=T)/NOPH)>1e-10))[1]>0){
      cat("\n", "- Yield per plant (YPP) is different from (CRW+NCRW)/NOPH:", "\n")
      print(subset(data, abs(YPP-apply(cbind(CRW, NCRW), 1, sum, na.rm=T)/NOPH)>1e-10))
    }

  if (exists("CI", where=data)==1 & exists("NOCR", where=data)==1 & exists("NONC", where=data)==1)
    if (dim(subset(data, abs(CI-NOCR/apply(cbind(NOCR,NONC), 1, sum, na.rm=T)*100)>1e-10))[1]>0){
      cat("\n", "- Commercial index (CI) is different from NOCR/(NOCR+NONC)*100:", "\n")
      print(subset(data, abs(CI-NOCR/apply(cbind(NOCR,NONC), 1, sum, na.rm=T)*100)>1e-10))
    }

  if (exists("HI", where=data)==1 & exists("CRW", where=data)==1 & exists("NCRW", where=data)==1 & exists("VW", where=data)==1)
    if (dim(subset(data, abs(HI-apply(cbind(CRW, NCRW), 1, sum, na.rm=T)/apply(cbind(VW, CRW, NCRW), 1, sum, na.rm=T)*100)>1e-10))[1]>0){
      cat("\n", "- Harvest index (HI) is different from (CRW+NCRW)/(VW+CRW+NCRW)*100:", "\n")
      print(subset(data, abs(HI-apply(cbind(CRW, NCRW), 1, sum, na.rm=T)/apply(cbind(VW, CRW, NCRW), 1, sum, na.rm=T)*100)>1e-10))
    }

  if (exists("SHI", where=data)==1 & exists("NOPH", where=data)==1 & exists("NOPS", where=data)==1)
    if (dim(subset(data, abs(SHI-NOPH/NOPS*100)>1e-10))[1]>0){
      cat("\n", "- Harvest sowing index (SHI) is different from NOPH/NOPS*100:", "\n")
      print(subset(data, abs(SHI-NOPH/NOPS*100)>1e-10))
    }

  if (exists("BIOM", where=data)==1 & exists("CRW", where=data)==1 & exists("NCRW", where=data)==1 & exists("VW", where=data)==1)
    if (dim(subset(data, abs(BIOM-apply(cbind(VW, CRW, NCRW), 1, sum, na.rm=T)*10/plot.size)>1e-10))[1]>0){
      cat("\n", "- Biomass yield (BIOM) is different from (CRW+NCRW+VW)*10/plot.size:", "\n")
      print(subset(data, abs(BIOM-apply(cbind(VW, CRW, NCRW), 1, sum, na.rm=T)*10/plot.size)>1e-10))
    }

  if (exists("FYTHA", where=data)==1 & exists("VW", where=data)==1)
    if (dim(subset(data, abs(FYTHA-VW*10/plot.size)>1e-10))[1]>0){
      cat("\n", "- Foliage total yield in tons per hectare (FYTHA) is different from VW*10/plot.size:", "\n")
      print(subset(data, abs(FYTHA-VW*10/plot.size)>1e-10))
    }

  if (exists("DM", where=data)==1 & exists("DMD", where=data)==1 & exists("DMF", where=data)==1)
    if (dim(subset(data, abs(DM-DMD/DMF*100)>1e-10))[1]>0){
      cat("\n", "- Storage root dry matter content (DM) is different from DMD/DMF*100:", "\n")
      print(subset(data, abs(DM-DMD/DMF*100)>1e-10))
    }

  if (exists("DMFY", where=data)==1 & exists("VW", where=data)==1 & exists("DMVD", where=data)==1 & exists("DMVF", where=data)==1)
    if (dim(subset(data, abs(DMFY-VW*10/plot.size*DMVD/DMVF)>1e-10))[1]>0){
      cat("\n", "- Dry matter foliage yield (DMFY) is different from VW*10/plot.size*DMVD/DMVF:", "\n")
      print(subset(data, abs(DMFY-VW*10/plot.size*DMVD/DMVF)>1e-10))
    }

  if (exists("DMRY", where=data)==1 & exists("CRW", where=data)==1 & exists("NCRW", where=data)==1
      & exists("DMD", where=data)==1 & exists("DMF", where=data)==1)
    if (dim(subset(data, DMRY!=apply(cbind(CRW, NCRW), 1, sum, na.rm=T)*10/plot.size*DMD/DMF))[1]>0){
      cat("\n", "- Dry matter root yield (DMRY) is different from (CRW+NCRW)*10/plot.size*DMD/DMF:", "\n")
      print(subset(data, DMRY!=apply(cbind(CRW, NCRW), 1, sum, na.rm=T)*10/plot.size*DMD/DMF))
    }

  if (exists("RFR", where=data)==1 & exists("CRW", where=data)==1 & exists("NCRW", where=data)==1 & exists("DMD", where=data)==1
      & exists("DMF", where=data)==1 & exists("VW", where=data)==1 & exists("DMVD", where=data)==1 & exists("DMVF", where=data)==1)
    if (dim(subset(data, abs(RFR-apply(cbind(CRW, NCRW), 1, sum, na.rm=T)*(DMD/DMF)/(VW*DMVD/DMVF))>1e-10))[1]>0){
      cat("\n", "- Root foliage ratio (RFR) is different from (CRW+NCRW)*(DMD/DMF)/(VW*DMVD/DMVF)*100:", "\n")
      print(subset(data, abs(RFR-apply(cbind(CRW, NCRW), 1, sum, na.rm=T)*(DMD/DMF)/(VW*DMVD/DMVF))>1e-10))
    }
}

# Check consistency for sweetpotato experimental data, part 10.
# Outliers detection based on interquartile range and values out of range for field data.

spconsis10 <- function(data){
  
  if (exists("NOPE", where=data)==1)
    if (dim(subset(data, !(NOPE %in% c(0:100, NA))))[1]>0){
      cat("\n", "- Out of range values for number of plants established (NOPE):", "\n")
      print(subset(data, !(NOPE %in% c(0:100, NA))))
    }

  if (exists("VIR1", where=data)==1)
    if (dim(subset(data, !(VIR1 %in% c(1:9, NA))))[1]>0){
      cat("\n", "- Out of range values for virus symptoms first evaluation (VIR1):", "\n")
      print(subset(data, !(VIR1 %in% c(1:9, NA))))
    }

  if (exists("VIR2", where=data)==1)
    if (dim(subset(data, !(VIR2 %in% c(1:9, NA))))[1]>0){
      cat("\n", "- Out of range values for virus symptoms second evaluation (VIR2):", "\n")
      print(subset(data, !(VIR2 %in% c(1:9, NA))))
    }

  if (exists("VIR3", where=data)==1)
    if (dim(subset(data, !(VIR3 %in% c(1:9, NA))))[1]>0){
      cat("\n", "- Out of range values for virus symptoms third evaluation (VIR3):", "\n")
      print(subset(data, !(VIR3 %in% c(1:9, NA))))
    }

  if (exists("ALT1", where=data)==1)
    if (dim(subset(data, !(ALT1 %in% c(1:9, NA))))[1]>0){
      cat("\n", "- Out of range values for alternaria symptoms first evaluation (ALT1):", "\n")
      print(subset(data, !(ALT1 %in% c(1:9, NA))))
    }

  if (exists("ALT2", where=data)==1)
    if (dim(subset(data, !(ALT2 %in% c(1:9, NA))))[1]>0){
      cat("\n", "- Out of range values for alternaria symptoms second evaluation (ALT2):", "\n")
      print(subset(data, !(ALT2 %in% c(1:9, NA))))
    }

  if (exists("VV1", where=data)==1)
    if (dim(subset(data, !(VV1 %in% c(1:9, NA))))[1]>0){
      cat("\n", "- Out of range values for vine vigor first evaluation (VV1):", "\n")
      print(subset(data, !(VV1 %in% c(1:9, NA))))
    }

  if (exists("VV2", where=data)==1)
    if (dim(subset(data, !(VV2 %in% c(1:9, NA))))[1]>0){
      cat("\n", "- Out of range values for vine vigor second evaluation (VV2):", "\n")
      print(subset(data, !(VV2 %in% c(1:9, NA))))
    }

  if (exists("VW", where=data)==1)
    if (dim(subset(data, VW < 0))[1]>0){
      cat("\n", "- Out of range values for vine weight (VW):", "\n")
      print(subset(data, VW < 0))
    }

  if (exists("VW", where=data)==1)
    if (dim(subset(data, VW < quantile(VW, 0.25, na.rm=T)-3*IQR(VW, na.rm=T)))[1]>0){
      cat("\n", "- Extreme low values for vine weight (VW):", "\n")
      print(subset(data, VW < quantile(VW, 0.25, na.rm=T)-3*IQR(VW, na.rm=T)))
    }

  if (exists("VW", where=data)==1)
    if (dim(subset(data, VW > quantile(VW, 0.75, na.rm=T)+3*IQR(VW, na.rm=T)))[1]>0){
      cat("\n", "- Extreme high values for vine weight (VW):", "\n")
      print(subset(data, VW > quantile(VW, 0.75, na.rm=T)+3*IQR(VW, na.rm=T)))
    }

  if (exists("NOPH", where=data)==1)
    if (dim(subset(data, !(NOPH %in% c(0:100, NA))))[1]>0){
      cat("\n", "- Out of range values for number of plants harvested (NOPH):", "\n")
      print(subset(data, !(NOPE %in% c(0:100, NA))))
    }

  if (exists("NOPR", where=data)==1)
    if (dim(subset(data, !(NOPR %in% c(0:100, NA))))[1]>0){
      cat("\n", "- Out of range values for number of plants with roots (NOPR):", "\n")
      print(subset(data, !(NOPR %in% c(0:100, NA))))
    }

  if (exists("NOCR", where=data)==1)
    if (dim(subset(data, !(NOCR %in% c(0:1000, NA))))[1]>0){
      cat("\n", "- Out of range values for number of commercial roots (NOCR):", "\n")
      print(subset(data, !(NOCR %in% c(0:1000, NA))))
    }

  if (exists("NOCR", where=data)==1)
    if (dim(subset(data, NOCR < quantile(NOCR, 0.25, na.rm=T)-3*IQR(NOCR, na.rm=T)))[1]>0){
      cat("\n", "- Extreme low values for number of commercial roots (NOCR):", "\n")
      print(subset(data, NOCR < quantile(NOCR, 0.25, na.rm=T)-3*IQR(NOCR, na.rm=T)))
    }

  if (exists("NOCR", where=data)==1)
    if (dim(subset(data, NOCR > quantile(NOCR, 0.75, na.rm=T)+3*IQR(NOCR, na.rm=T)))[1]>0){
      cat("\n", "- Extreme high values for number of commercial roots (NOCR):", "\n")
      print(subset(data, NOCR > quantile(NOCR, 0.75, na.rm=T)+3*IQR(NOCR, na.rm=T)))
    }

  if (exists("NONC", where=data)==1)
    if (dim(subset(data, !(NONC %in% c(0:1000, NA))))[1]>0){
      cat("\n", "- Out of range values for number of non commercial roots (NONC):", "\n")
      print(subset(data, !(NONC %in% c(0:1000, NA))))
    }

  if (exists("NONC", where=data)==1)
    if (dim(subset(data, NONC < quantile(NONC, 0.25, na.rm=T)-3*IQR(NONC, na.rm=T)))[1]>0){
      cat("\n", "- Extreme low values for number of non commercial roots (NONC):", "\n")
      print(subset(data, NONC < quantile(NONC, 0.25, na.rm=T)-3*IQR(NONC, na.rm=T)))
    }

  if (exists("NONC", where=data)==1)
    if (dim(subset(data, NONC > quantile(NONC, 0.75, na.rm=T)+3*IQR(NONC, na.rm=T)))[1]>0){
      cat("\n", "- Extreme high values for number of non commercial roots (NONC):", "\n")
      print(subset(data, NONC > quantile(NONC, 0.75, na.rm=T)+3*IQR(NONC, na.rm=T)))
    }

  if (exists("CRW", where=data)==1)
    if (dim(subset(data, CRW < 0))[1]>0){
      cat("\n", "- Out of range values for commercial root weight (CRW):", "\n")
      print(subset(data, CRW < 0))
    }

  if (exists("CRW", where=data)==1)
    if (dim(subset(data, CRW < quantile(CRW, 0.25, na.rm=T)-3*IQR(CRW, na.rm=T)))[1]>0){
      cat("\n", "- Extreme low values for commercial root weight (CRW):", "\n")
      print(subset(data, CRW < quantile(CRW, 0.25, na.rm=T)-3*IQR(CRW, na.rm=T)))
    }

  if (exists("CRW", where=data)==1)
    if (dim(subset(data, CRW > quantile(CRW, 0.75, na.rm=T)+3*IQR(CRW, na.rm=T)))[1]>0){
      cat("\n", "- Extreme high values for commercial root weight (CRW):", "\n")
      print(subset(data, CRW > quantile(CRW, 0.75, na.rm=T)+3*IQR(CRW, na.rm=T)))
    }

  if (exists("NCRW", where=data)==1)
    if (dim(subset(data, NCRW < 0))[1]>0){
      cat("\n", "- Out of range values for non commercial root weight (NCRW):", "\n")
      print(subset(data, NCRW < 0))
    }

  if (exists("NCRW", where=data)==1)
    if (dim(subset(data, NCRW < quantile(NCRW, 0.25, na.rm=T)-3*IQR(NCRW, na.rm=T)))[1]>0){
      cat("\n", "- Extreme low values for non commercial root weight (NCRW):", "\n")
      print(subset(data, NCRW < quantile(NCRW, 0.25, na.rm=T)-3*IQR(NCRW, na.rm=T)))
    }

  if (exists("NCRW", where=data)==1)
    if (dim(subset(data, NCRW > quantile(NCRW, 0.75, na.rm=T)+3*IQR(NCRW, na.rm=T)))[1]>0){
      cat("\n", "- Extreme high values for non commercial root weight (NCRW):", "\n")
      print(subset(data, NCRW > quantile(NCRW, 0.75, na.rm=T)+3*IQR(NCRW, na.rm=T)))
    }

  if (exists("SCOL", where=data)==1)
    if (dim(subset(data, !(SCOL %in% c(1:9, NA))))[1]>0){
      cat("\n", "- Out of range values for storage root skin color (SCOL):", "\n")
      print(subset(data, !(SCOL %in% c(1:9, NA))))
    }

  if (exists("FCOL", where=data)==1)
    if (dim(subset(data, !(FCOL %in% c(1:9, NA))))[1]>0){
      cat("\n", "- Out of range values for storage root flesh color (FCOL):", "\n")
      print(subset(data, !(FCOL %in% c(1:9, NA))))
    }

  if (exists("RS", where=data)==1)
    if (dim(subset(data, !(RS %in% c(1:9, NA))))[1]>0){
      cat("\n", "- Out of range values for root size (RS):", "\n")
      print(subset(data, !(RS %in% c(1:9, NA))))
    }

  if (exists("RF", where=data)==1)
    if (dim(subset(data, !(RF %in% c(1:9, NA))))[1]>0){
      cat("\n", "- Out of range values for root form (RF):", "\n")
      print(subset(data, !(RF %in% c(1:9, NA))))
    }

  if (exists("DAMR", where=data)==1)
    if (dim(subset(data, !(DAMR %in% c(1:9, NA))))[1]>0){
      cat("\n", "- Out of range values for root defects (DAMR):", "\n")
      print(subset(data, !(DAMR %in% c(1:9, NA))))
    }

  if (exists("RSPR", where=data)==1)
    if (dim(subset(data, !(RSPR %in% c(1:9, NA))))[1]>0){
      cat("\n", "- Out of range values for root sprouting (RSPR):", "\n")
      print(subset(data, !(RSPR %in% c(1:9, NA))))
    }

  if (exists("WED1", where=data)==1)
    if (dim(subset(data, !(WED1 %in% c(1:9, NA))))[1]>0){
      cat("\n", "- Out of range values for weevil damage first evaluation (WED1):", "\n")
      print(subset(data, !(WED1 %in% c(1:9, NA))))
    }

  if (exists("WED2", where=data)==1)
    if (dim(subset(data, !(WED2 %in% c(1:9, NA))))[1]>0){
      cat("\n", "- Out of range values for weevil damage second evaluation (WED2):", "\n")
      print(subset(data, !(WED2 %in% c(1:9, NA))))
    }
}

# Check consistency for sweetpotato experimental data, part 11.
# Outliers detection based on interquartile range and values out of range for DM data.

spconsis11 <- function(data){

  if (exists("DMF", where=data)==1)
    if (dim(subset(data, DMF < 0))[1]>0){
      cat("\n", "- Out of range values for fresh weight of roots for dry matter assessment (DMF):", "\n")
      print(subset(data, DMF < 0))
    }

  if (exists("DMF", where=data)==1)
    if (dim(subset(data, DMF < quantile(DMF, 0.25, na.rm=T)-3*IQR(DMF, na.rm=T)))[1]>0){
      cat("\n", "- Extreme low values for fresh weight of roots for dry matter assessment (DMF):", "\n")
      print(subset(data, DMF < quantile(DMF, 0.25, na.rm=T)-3*IQR(DMF, na.rm=T)))
    }

  if (exists("DMF", where=data)==1)
    if (dim(subset(data, DMF > quantile(DMF, 0.75, na.rm=T)+3*IQR(DMF, na.rm=T)))[1]>0){
      cat("\n", "- Extreme high values for fresh weight of roots for dry matter assessment (DMF):", "\n")
      print(subset(data, DMF > quantile(DMF, 0.75, na.rm=T)+3*IQR(DMF, na.rm=T)))
    }

  if (exists("DMD", where=data)==1)
    if (dim(subset(data, DMD < 0))[1]>0){
      cat("\n", "- Out of range values for dry weight of roots for dry matter assessment (DMD):", "\n")
      print(subset(data, DMD < 0))
    }

  if (exists("DMD", where=data)==1)
    if (dim(subset(data, DMD < quantile(DMD, 0.25, na.rm=T)-3*IQR(DMD, na.rm=T)))[1]>0){
      cat("\n", "- Extreme low values for dry weight of roots for dry matter assessment (DMD):", "\n")
      print(subset(data, DMD < quantile(DMD, 0.25, na.rm=T)-3*IQR(DMD, na.rm=T)))
    }

  if (exists("DMD", where=data)==1)
    if (dim(subset(data, DMD > quantile(DMD, 0.75, na.rm=T)+3*IQR(DMD, na.rm=T)))[1]>0){
      cat("\n", "- Extreme high values for dry weight of roots for dry matter assessment (DMD):", "\n")
      print(subset(data, DMD > quantile(DMD, 0.75, na.rm=T)+3*IQR(DMD, na.rm=T)))
    }

  if (exists("DMVF", where=data)==1)
    if (dim(subset(data, DMVF < 0))[1]>0){
      cat("\n", "- Out of range values for fresh weight vines for dry matter assessment (DMVF):", "\n")
      print(subset(data, DMVF < 0))
    }

  if (exists("DMVF", where=data)==1)
    if (dim(subset(data, DMVF < quantile(DMVF, 0.25, na.rm=T)-3*IQR(DMVF, na.rm=T)))[1]>0){
      cat("\n", "- Extreme low values for fresh weight of vines for dry matter assessment (DMVF):", "\n")
      print(subset(data, DMVF < quantile(DMVF, 0.25, na.rm=T)-3*IQR(DMVF, na.rm=T)))
    }

  if (exists("DMVF", where=data)==1)
    if (dim(subset(data, DMVF > quantile(DMVF, 0.75, na.rm=T)+3*IQR(DMVF, na.rm=T)))[1]>0){
      cat("\n", "- Extreme high values for fresh weight of vines for dry matter assessment (DMVF):", "\n")
      print(subset(data, DMVF > quantile(DMVF, 0.75, na.rm=T)+3*IQR(DMVF, na.rm=T)))
    }

  if (exists("DMVD", where=data)==1)
    if (dim(subset(data, DMVD < 0))[1]>0){
      cat("\n", "- Out of range values for dry weight of vines for dry matter assessment (DMVD):", "\n")
      print(subset(data, DMVD < 0))
    }

  if (exists("DMVD", where=data)==1)
    if (dim(subset(data, DMVD < quantile(DMVD, 0.25, na.rm=T)-3*IQR(DMVD, na.rm=T)))[1]>0){
      cat("\n", "- Extreme low values for dry weight of vines for dry matter assessment (DMVD):", "\n")
      print(subset(data, DMVD < quantile(DMVD, 0.25, na.rm=T)-3*IQR(DMVD, na.rm=T)))
    }

  if (exists("DMVD", where=data)==1)
    if (dim(subset(data, DMVD > quantile(DMVD, 0.75, na.rm=T)+3*IQR(DMVD, na.rm=T)))[1]>0){
      cat("\n", "- Extreme high values for dry weight of vines for dry matter assessment (DMVD):", "\n")
      print(subset(data, DMVD > quantile(DMVD, 0.75, na.rm=T)+3*IQR(DMVD, na.rm=T)))
    }

  if (exists("DM", where=data)==1)
    if (dim(subset(data, DM < 0 | DM > 100))[1]>0){
      cat("\n", "- Out of range values for storage root dry matter content (DM):", "\n")
      print(subset(data, DM < 0 | DM > 100))
    }

  if (exists("DM", where=data)==1)
    if (dim(subset(data, DM < quantile(DM, 0.25, na.rm=T)-3*IQR(DM, na.rm=T)))[1]>0){
      cat("\n", "- Extreme low values for storage root dry matter content (DM):", "\n")
      print(subset(data, DM < quantile(DM, 0.25, na.rm=T)-3*IQR(DM, na.rm=T)))
    }

  if (exists("DM", where=data)==1)
    if (dim(subset(data, DM > quantile(DM, 0.75, na.rm=T)+3*IQR(DM, na.rm=T)))[1]>0){
      cat("\n", "- Extreme high values for storage root dry matter content (DM):", "\n")
      print(subset(data, DM > quantile(DM, 0.75, na.rm=T)+3*IQR(DM, na.rm=T)))
    }

  if (exists("DMFY", where=data)==1)
    if (dim(subset(data, DMFY < 0))[1]>0){
      cat("\n", "- Out of range values for dry matter foliage yield (DMFY):", "\n")
      print(subset(data, DMFY < 0))
    }

  if (exists("DMFY", where=data)==1)
    if (dim(subset(data, DMFY < quantile(DMFY, 0.25, na.rm=T)-3*IQR(DMFY, na.rm=T)))[1]>0){
      cat("\n", "- Extreme low values for dry matter foliage yield (DMFY):", "\n")
      print(subset(data, DMFY < quantile(DMFY, 0.25, na.rm=T)-3*IQR(DMFY, na.rm=T)))
    }

  if (exists("DMFY", where=data)==1)
    if (dim(subset(data, DMFY > quantile(DMFY, 0.75, na.rm=T)+3*IQR(DMFY, na.rm=T)))[1]>0){
      cat("\n", "- Extreme high values for dry matter foliage yield (DMFY):", "\n")
      print(subset(data, DMFY > quantile(DMFY, 0.75, na.rm=T)+3*IQR(DMFY, na.rm=T)))
    }

  if (exists("DMRY", where=data)==1)
    if (dim(subset(data, DMRY < 0))[1]>0){
      cat("\n", "- Out of range values for dry matter root yield (DMRY):", "\n")
      print(subset(data, DMRY < 0))
    }

  if (exists("DMRY", where=data)==1)
    if (dim(subset(data, DMRY < quantile(DMRY, 0.25, na.rm=T)-3*IQR(DMRY, na.rm=T)))[1]>0){
      cat("\n", "- Extreme low values for dry matter root yield (DMRY):", "\n")
      print(subset(data, DMRY < quantile(DMRY, 0.25, na.rm=T)-3*IQR(DMRY, na.rm=T)))
    }

  if (exists("DMRY", where=data)==1)
    if (dim(subset(data, DMRY > quantile(DMRY, 0.75, na.rm=T)+3*IQR(DMRY, na.rm=T)))[1]>0){
      cat("\n", "- Extreme high values for dry matter root yield (DMRY):", "\n")
      print(subset(data, DMRY > quantile(DMRY, 0.75, na.rm=T)+3*IQR(DMRY, na.rm=T)))
    }
}

# Check consistency for sweetpotato experimental data, part 12.
# Outliers detection based on interquartile range and values out of range for cooked traits.

spconsis12 <- function(data){

  if (exists("FRAW1", where=data)==1)
    if (dim(subset(data, !(FRAW1 %in% c(1:9, NA))))[1]>0){
      cat("\n", "- Out of range values for root fiber first determination (FRAW1):", "\n")
      print(subset(data, !(FRAW1 %in% c(1:9, NA))))
    }

  if (exists("SURAW1", where=data)==1)
    if (dim(subset(data, !(SURAW1 %in% c(1:9, NA))))[1]>0){
      cat("\n", "- Out of range values for root sugar first determination (SURAW1):", "\n")
      print(subset(data, !(SURAW1 %in% c(1:9, NA))))
    }

  if (exists("STRAW1", where=data)==1)
    if (dim(subset(data, !(STRAW1 %in% c(1:9, NA))))[1]>0){
      cat("\n", "- Out of range values for root starch first determination (STRAW1):", "\n")
      print(subset(data, !(STRAW1 %in% c(1:9, NA))))
    }

  if (exists("COOF1", where=data)==1)
    if (dim(subset(data, !(COOF1 %in% c(1:9, NA))))[1]>0){
      cat("\n", "- Out of range values for cooked fiber first evaluation (COOF1):", "\n")
      print(subset(data, !(COOF1 %in% c(1:9, NA))))
    }

  if (exists("COOSU1", where=data)==1)
    if (dim(subset(data, !(COOSU1 %in% c(1:9, NA))))[1]>0){
      cat("\n", "- Out of range values for cooked sugars first evaluation (COOSU1):", "\n")
      print(subset(data, !(COOSU1 %in% c(1:9, NA))))
    }

  if (exists("COOST1", where=data)==1)
    if (dim(subset(data, !(COOST1 %in% c(1:9, NA))))[1]>0){
      cat("\n", "- Out of range values for cooked starch first evaluation (COOST1):", "\n")
      print(subset(data, !(COOST1 %in% c(1:9, NA))))
    }

  if (exists("COOT1", where=data)==1)
    if (dim(subset(data, !(COOT1 %in% c(1:9, NA))))[1]>0){
      cat("\n", "- Out of range values for cooked taste first evaluation (COOT1):", "\n")
      print(subset(data, !(COOT1 %in% c(1:9, NA))))
    }

  if (exists("COOAP1", where=data)==1)
    if (dim(subset(data, !(COOAP1 %in% c(1:9, NA))))[1]>0){
      cat("\n", "- Out of range values for cooked appearance first evaluation (COOAP1):", "\n")
      print(subset(data, !(COOAP1 %in% c(1:9, NA))))
    }

  if (exists("FRAW2", where=data)==1)
    if (dim(subset(data, !(FRAW2 %in% c(1:9, NA))))[1]>0){
      cat("\n", "- Out of range values for root fiber second determination (FRAW2):", "\n")
      print(subset(data, !(FRAW2 %in% c(1:9, NA))))
    }

  if (exists("SURAW2", where=data)==1)
    if (dim(subset(data, !(SURAW2 %in% c(1:9, NA))))[1]>0){
      cat("\n", "- Out of range values for root sugar second determination (SURAW2):", "\n")
      print(subset(data, !(SURAW2 %in% c(1:9, NA))))
    }

  if (exists("STRAW2", where=data)==1)
    if (dim(subset(data, !(STRAW2 %in% c(1:9, NA))))[1]>0){
      cat("\n", "- Out of range values for root starch second determination (STRAW2):", "\n")
      print(subset(data, !(STRAW2 %in% c(1:9, NA))))
    }

  if (exists("COOF2", where=data)==1)
    if (dim(subset(data, !(COOF2 %in% c(1:9, NA))))[1]>0){
      cat("\n", "- Out of range values for cooked fiber second evaluation (COOF2):", "\n")
      print(subset(data, !(COOF2 %in% c(1:9, NA))))
    }

  if (exists("COOSU2", where=data)==1)
    if (dim(subset(data, !(COOSU2 %in% c(1:9, NA))))[1]>0){
      cat("\n", "- Out of range values for cooked sugars second evaluation (COOSU2):", "\n")
      print(subset(data, !(COOSU2 %in% c(1:9, NA))))
    }

  if (exists("COOST2", where=data)==1)
    if (dim(subset(data, !(COOST2 %in% c(1:9, NA))))[1]>0){
      cat("\n", "- Out of range values for cooked starch second evaluation (COOST2):", "\n")
      print(subset(data, !(COOST2 %in% c(1:9, NA))))
    }

  if (exists("COOT2", where=data)==1)
    if (dim(subset(data, !(COOT2 %in% c(1:9, NA))))[1]>0){
      cat("\n", "- Out of range values for cooked taste second evaluation (COOT2):", "\n")
      print(subset(data, !(COOT2 %in% c(1:9, NA))))
    }

  if (exists("COOAP2", where=data)==1)
    if (dim(subset(data, !(COOAP2 %in% c(1:9, NA))))[1]>0){
      cat("\n", "- Out of range values for cooked appearance second evaluation (COOAP2):", "\n")
      print(subset(data, !(COOAP2 %in% c(1:9, NA))))
    }
}

# Check consistency for sweetpotato experimental data, part 13.
# Outliers detection based on interquartile range and values out of range for lab data.

spconsis13 <- function(data){

  bc.cc.values <- c(0.03, 0, 0.12, 0.02, 0.15, 1.38, 1.65, 1.5, 1.74, 1.76, 0.69, 1.17, 1.32,
                    1.04, 4.41, 4.92, 6.12, 5.46, 3.96, 5.49, 3.03, 3.76, 4.61, 7.23, 7.76,
                    10.5, 11.03, 12.39, 14.37)

  if (exists("PROT", where=data)==1)
    if (dim(subset(data, PROT < 0))[1]>0){
      cat("\n", "- Out of range values for protein (PROT):", "\n")
      print(subset(data, PROT < 0))
    }

  if (exists("PROT", where=data)==1)
    if (dim(subset(data, PROT < quantile(PROT, 0.25, na.rm=T)-3*IQR(PROT, na.rm=T)))[1]>0){
      cat("\n", "- Extreme low values for protein (PROT):", "\n")
      print(subset(data, PROT < quantile(PROT, 0.25, na.rm=T)-3*IQR(PROT, na.rm=T)))
    }

  if (exists("PROT", where=data)==1)
    if (dim(subset(data, PROT > quantile(PROT, 0.75, na.rm=T)+3*IQR(PROT, na.rm=T)))[1]>0){
      cat("\n", "- Extreme high values for protein (PROT):", "\n")
      print(subset(data, PROT > quantile(PROT, 0.75, na.rm=T)+3*IQR(PROT, na.rm=T)))
    }

  if (exists("FE", where=data)==1)
    if (dim(subset(data, FE < 0))[1]>0){
      cat("\n", "- Out of range values for iron in dry weight (FE):", "\n")
      print(subset(data, FE < 0))
    }

  if (exists("FE", where=data)==1)
    if (dim(subset(data, FE < quantile(FE, 0.25, na.rm=T)-3*IQR(FE, na.rm=T)))[1]>0){
      cat("\n", "- Extreme low values for iron in dry weight (FE):", "\n")
      print(subset(data, FE < quantile(FE, 0.25, na.rm=T)-3*IQR(FE, na.rm=T)))
    }

  if (exists("FE", where=data)==1)
    if (dim(subset(data, FE > quantile(FE, 0.75, na.rm=T)+3*IQR(FE, na.rm=T)))[1]>0){
      cat("\n", "- Extreme high values for iron in dry weight (FE):", "\n")
      print(subset(data, FE > quantile(FE, 0.75, na.rm=T)+3*IQR(FE, na.rm=T)))
    }

  if (exists("ZN", where=data)==1)
    if (dim(subset(data, ZN < 0))[1]>0){
      cat("\n", "- Out of range values for zinc in dry weight (ZN):", "\n")
      print(subset(data, ZN < 0))
    }

  if (exists("ZN", where=data)==1)
    if (dim(subset(data, ZN < quantile(ZN, 0.25, na.rm=T)-3*IQR(ZN, na.rm=T)))[1]>0){
      cat("\n", "- Extreme low values for zinc in dry weight (ZN):", "\n")
      print(subset(data, ZN < quantile(ZN, 0.25, na.rm=T)-3*IQR(ZN, na.rm=T)))
    }

  if (exists("ZN", where=data)==1)
    if (dim(subset(data, ZN > quantile(ZN, 0.75, na.rm=T)+3*IQR(ZN, na.rm=T)))[1]>0){
      cat("\n", "- Extreme high values for zinc in dry weight (ZN):", "\n")
      print(subset(data, ZN > quantile(ZN, 0.75, na.rm=T)+3*IQR(ZN, na.rm=T)))
    }

  if (exists("CA", where=data)==1)
    if (dim(subset(data, CA < 0))[1]>0){
      cat("\n", "- Out of range values for calcium in dry weight (CA):", "\n")
      print(subset(data, CA < 0))
    }

  if (exists("CA", where=data)==1)
    if (dim(subset(data, CA < quantile(CA, 0.25, na.rm=T)-3*IQR(CA, na.rm=T)))[1]>0){
      cat("\n", "- Extreme low values for calcium in dry weight (CA):", "\n")
      print(subset(data, CA < quantile(CA, 0.25, na.rm=T)-3*IQR(CA, na.rm=T)))
    }

  if (exists("CA", where=data)==1)
    if (dim(subset(data, CA > quantile(CA, 0.75, na.rm=T)+3*IQR(CA, na.rm=T)))[1]>0){
      cat("\n", "- Extreme high values for calcium in dry weight (CA):", "\n")
      print(subset(data, CA > quantile(CA, 0.75, na.rm=T)+3*IQR(CA, na.rm=T)))
    }

  if (exists("MG", where=data)==1)
    if (dim(subset(data, MG < 0))[1]>0){
      cat("\n", "- Out of range values for magnesium in dry weight (MG):", "\n")
      print(subset(data, MG < 0))
    }

  if (exists("MG", where=data)==1)
    if (dim(subset(data, MG < quantile(MG, 0.25, na.rm=T)-3*IQR(MG, na.rm=T)))[1]>0){
      cat("\n", "- Extreme low values for magnesium in dry weight (MG):", "\n")
      print(subset(data, MG < quantile(MG, 0.25, na.rm=T)-3*IQR(MG, na.rm=T)))
    }

  if (exists("MG", where=data)==1)
    if (dim(subset(data, MG > quantile(MG, 0.75, na.rm=T)+3*IQR(MG, na.rm=T)))[1]>0){
      cat("\n", "- Extreme high values for magnesium in dry weight (MG):", "\n")
      print(subset(data, MG > quantile(MG, 0.75, na.rm=T)+3*IQR(MG, na.rm=T)))
    }

  if (exists("BC", where=data)==1)
    if (dim(subset(data, BC < 0))[1]>0){
      cat("\n", "- Out of range values for beta-carotene in dry weight (BC):", "\n")
      print(subset(data, BC < 0))
    }

  if (exists("BC", where=data)==1)
    if (dim(subset(data, BC < quantile(BC, 0.25, na.rm=T)-3*IQR(BC, na.rm=T)))[1]>0){
      cat("\n", "- Extreme low values for beta-carotene in dry weight (BC):", "\n")
      print(subset(data, BC < quantile(BC, 0.25, na.rm=T)-3*IQR(BC, na.rm=T)))
    }

  if (exists("BC", where=data)==1)
    if (dim(subset(data, BC > quantile(BC, 0.75, na.rm=T)+3*IQR(BC, na.rm=T)))[1]>0){
      cat("\n", "- Extreme high values for beta-carotene in dry weight (BC):", "\n")
      print(subset(data, BC > quantile(BC, 0.75, na.rm=T)+3*IQR(BC, na.rm=T)))
    }

  if (exists("BC.CC", where=data)==1)
    if (dim(subset(data, !(BC.CC %in% c(bc.cc.values, NA))))[1]>0){
      cat("\n", "- Out of range values for beta-carotene with color chart (BC.CC):", "\n")
      print(subset(data, !(BC.CC %in% c(bc.cc.values, NA))))
    }

  if (exists("TC", where=data)==1)
    if (dim(subset(data, TC < 0))[1]>0){
      cat("\n", "- Out of range values for total carotenoids in dry weight (TC):", "\n")
      print(subset(data, TC < 0))
    }

  if (exists("TC", where=data)==1)
    if (dim(subset(data, TC < quantile(TC, 0.25, na.rm=T)-3*IQR(TC, na.rm=T)))[1]>0){
      cat("\n", "- Extreme low values for total carotenoids in dry weight (TC):", "\n")
      print(subset(data, TC < quantile(TC, 0.25, na.rm=T)-3*IQR(TC, na.rm=T)))
    }

  if (exists("TC", where=data)==1)
    if (dim(subset(data, TC > quantile(TC, 0.75, na.rm=T)+3*IQR(TC, na.rm=T)))[1]>0){
      cat("\n", "- Extreme high values for total carotenoids in dry weight (TC):", "\n")
      print(subset(data, TC > quantile(TC, 0.75, na.rm=T)+3*IQR(TC, na.rm=T)))
    }

  if (exists("STAR", where=data)==1)
    if (dim(subset(data, STAR < 0))[1]>0){
      cat("\n", "- Out of range values for starch (STAR):", "\n")
      print(subset(data, STAR < 0))
    }

  if (exists("STAR", where=data)==1)
    if (dim(subset(data, STAR < quantile(STAR, 0.25, na.rm=T)-3*IQR(STAR, na.rm=T)))[1]>0){
      cat("\n", "- Extreme low values for starch (STAR):", "\n")
      print(subset(data, STAR < quantile(STAR, 0.25, na.rm=T)-3*IQR(STAR, na.rm=T)))
    }

  if (exists("STAR", where=data)==1)
    if (dim(subset(data, STAR > quantile(STAR, 0.75, na.rm=T)+3*IQR(STAR, na.rm=T)))[1]>0){
      cat("\n", "- Extreme high values for starch (STAR):", "\n")
      print(subset(data, STAR > quantile(STAR, 0.75, na.rm=T)+3*IQR(STAR, na.rm=T)))
    }

  if (exists("FRUC", where=data)==1)
    if (dim(subset(data, FRUC < 0))[1]>0){
      cat("\n", "- Out of range values for fructose (FRUC):", "\n")
      print(subset(data, FRUC < 0))
    }

  if (exists("FRUC", where=data)==1)
    if (dim(subset(data, FRUC < quantile(FRUC, 0.25, na.rm=T)-3*IQR(FRUC, na.rm=T)))[1]>0){
      cat("\n", "- Extreme low values for fructose (FRUC):", "\n")
      print(subset(data, FRUC < quantile(FRUC, 0.25, na.rm=T)-3*IQR(FRUC, na.rm=T)))
    }

  if (exists("FRUC", where=data)==1)
    if (dim(subset(data, FRUC > quantile(FRUC, 0.75, na.rm=T)+3*IQR(FRUC, na.rm=T)))[1]>0){
      cat("\n", "- Extreme high values for fructose (FRUC):", "\n")
      print(subset(data, FRUC > quantile(FRUC, 0.75, na.rm=T)+3*IQR(FRUC, na.rm=T)))
    }

  if (exists("GLUC", where=data)==1)
    if (dim(subset(data, GLUC < 0))[1]>0){
      cat("\n", "- Out of range values for glucose (GLUC):", "\n")
      print(subset(data, GLUC < 0))
    }

  if (exists("GLUC", where=data)==1)
    if (dim(subset(data, GLUC < quantile(GLUC, 0.25, na.rm=T)-3*IQR(GLUC, na.rm=T)))[1]>0){
      cat("\n", "- Extreme low values for glucose (GLUC):", "\n")
      print(subset(data, GLUC < quantile(GLUC, 0.25, na.rm=T)-3*IQR(GLUC, na.rm=T)))
    }

  if (exists("GLUC", where=data)==1)
    if (dim(subset(data, GLUC > quantile(GLUC, 0.75, na.rm=T)+3*IQR(GLUC, na.rm=T)))[1]>0){
      cat("\n", "- Extreme high values for glucose (GLUC):", "\n")
      print(subset(data, GLUC > quantile(GLUC, 0.75, na.rm=T)+3*IQR(GLUC, na.rm=T)))
    }

  if (exists("SUCR", where=data)==1)
    if (dim(subset(data, SUCR < 0))[1]>0){
      cat("\n", "- Out of range values for sucrose (SUCR):", "\n")
      print(subset(data, SUCR < 0))
    }

  if (exists("SUCR", where=data)==1)
    if (dim(subset(data, SUCR < quantile(SUCR, 0.25, na.rm=T)-3*IQR(SUCR, na.rm=T)))[1]>0){
      cat("\n", "- Extreme low values for sucrose (SUCR):", "\n")
      print(subset(data, SUCR < quantile(SUCR, 0.25, na.rm=T)-3*IQR(SUCR, na.rm=T)))
    }

  if (exists("SUCR", where=data)==1)
    if (dim(subset(data, SUCR > quantile(SUCR, 0.75, na.rm=T)+3*IQR(SUCR, na.rm=T)))[1]>0){
      cat("\n", "- Extreme high values for sucrose (SUCR):", "\n")
      print(subset(data, SUCR > quantile(SUCR, 0.75, na.rm=T)+3*IQR(SUCR, na.rm=T)))
    }

  if (exists("MALT", where=data)==1)
    if (dim(subset(data, MALT < 0))[1]>0){
      cat("\n", "- Out of range values for maltose (MALT):", "\n")
      print(subset(data, MALT < 0))
    }

  if (exists("MALT", where=data)==1)
    if (dim(subset(data, MALT < quantile(MALT, 0.25, na.rm=T)-3*IQR(MALT, na.rm=T)))[1]>0){
      cat("\n", "- Extreme low values for maltose (MALT):", "\n")
      print(subset(data, MALT < quantile(MALT, 0.25, na.rm=T)-3*IQR(MALT, na.rm=T)))
    }

  if (exists("MALT", where=data)==1)
    if (dim(subset(data, MALT > quantile(MALT, 0.75, na.rm=T)+3*IQR(MALT, na.rm=T)))[1]>0){
      cat("\n", "- Extreme high values for maltose (MALT):", "\n")
      print(subset(data, MALT > quantile(MALT, 0.75, na.rm=T)+3*IQR(MALT, na.rm=T)))
    }
}

# Check consistency for sweetpotato experimental data, part 14.
# Outliers detection based on interquartile range and values out of range for
# derived variables.

spconsis14 <- function(data){

  if (exists("TRW", where=data)==1)
    if (dim(subset(data, TRW < 0))[1]>0){
      cat("\n", "- Out of range values for total root weight (TRW):", "\n")
      print(subset(data, TRW < 0))
    }

  if (exists("TRW", where=data)==1)
    if (dim(subset(data, TRW < quantile(TRW, 0.25, na.rm=T)-3*IQR(TRW, na.rm=T)))[1]>0){
      cat("\n", "- Extreme low values for total root weight (TRW):", "\n")
      print(subset(data, TRW < quantile(TRW, 0.25, na.rm=T)-3*IQR(TRW, na.rm=T)))
    }

  if (exists("TRW", where=data)==1)
    if (dim(subset(data, TRW > quantile(TRW, 0.75, na.rm=T)+3*IQR(TRW, na.rm=T)))[1]>0){
      cat("\n", "- Extreme high values for total root weight (TRW):", "\n")
      print(subset(data, TRW > quantile(TRW, 0.75, na.rm=T)+3*IQR(TRW, na.rm=T)))
    }

  if (exists("CYTHA", where=data)==1)
    if (dim(subset(data, CYTHA < 0))[1]>0){
      cat("\n", "- Out of range values for commercial root yield in tons per hectare (CYTHA):", "\n")
      print(subset(data, CYTHA < 0))
    }

  if (exists("CYTHA", where=data)==1)
    if (dim(subset(data, CYTHA < quantile(CYTHA, 0.25, na.rm=T)-3*IQR(CYTHA, na.rm=T)))[1]>0){
      cat("\n", "- Extreme low values for commercial root yield in tons per hectare (CYTHA):", "\n")
      print(subset(data, CYTHA < quantile(CYTHA, 0.25, na.rm=T)-3*IQR(CYTHA, na.rm=T)))
    }

  if (exists("CYTHA", where=data)==1)
    if (dim(subset(data, CYTHA > quantile(CYTHA, 0.75, na.rm=T)+3*IQR(CYTHA, na.rm=T)))[1]>0){
      cat("\n", "- Extreme high values for commercial root yield in tons per hectare (CYTHA):", "\n")
      print(subset(data, CYTHA > quantile(CYTHA, 0.75, na.rm=T)+3*IQR(CYTHA, na.rm=T)))
    }

  if (exists("RYTHA", where=data)==1)
    if (dim(subset(data, RYTHA < 0))[1]>0){
      cat("\n", "- Out of range values for total root yield in tons per hectare (RYTHA):", "\n")
      print(subset(data, RYTHA < 0))
    }

  if (exists("RYTHA", where=data)==1)
    if (dim(subset(data, RYTHA < quantile(RYTHA, 0.25, na.rm=T)-3*IQR(RYTHA, na.rm=T)))[1]>0){
      cat("\n", "- Extreme low values for total root yield in tons per hectare (RYTHA):", "\n")
      print(subset(data, RYTHA < quantile(RYTHA, 0.25, na.rm=T)-3*IQR(RYTHA, na.rm=T)))
    }

  if (exists("RYTHA", where=data)==1)
    if (dim(subset(data, RYTHA > quantile(RYTHA, 0.75, na.rm=T)+3*IQR(RYTHA, na.rm=T)))[1]>0){
      cat("\n", "- Extreme high values for total root yield in tons per hectare (RYTHA):", "\n")
      print(subset(data, RYTHA > quantile(RYTHA, 0.75, na.rm=T)+3*IQR(RYTHA, na.rm=T)))
    }

  if (exists("ACRW", where=data)==1)
    if (dim(subset(data, ACRW < 0))[1]>0){
      cat("\n", "- Out of range values for average commercial root weight (ACRW):", "\n")
      print(subset(data, ACRW < 0))
    }

  if (exists("ACRW", where=data)==1)
    if (dim(subset(data, ACRW < quantile(ACRW, 0.25, na.rm=T)-3*IQR(ACRW, na.rm=T)))[1]>0){
      cat("\n", "- Extreme low values for average commercial root weight (ACRW):", "\n")
      print(subset(data, ACRW < quantile(ACRW, 0.25, na.rm=T)-3*IQR(ACRW, na.rm=T)))
    }

  if (exists("ACRW", where=data)==1)
    if (dim(subset(data, ACRW > quantile(ACRW, 0.75, na.rm=T)+3*IQR(ACRW, na.rm=T)))[1]>0){
      cat("\n", "- Extreme high values for average commercial root weight (ACRW):", "\n")
      print(subset(data, ACRW > quantile(ACRW, 0.75, na.rm=T)+3*IQR(ACRW, na.rm=T)))
    }

  if (exists("NRPP", where=data)==1)
    if (dim(subset(data, NRPP < 0))[1]>0){
      cat("\n", "- Out of range values for number of roots per plant (NRPP):", "\n")
      print(subset(data, NRPP < 0))
    }

  if (exists("NRPP", where=data)==1)
    if (dim(subset(data, NRPP < quantile(NRPP, 0.25, na.rm=T)-3*IQR(NRPP, na.rm=T)))[1]>0){
      cat("\n", "- Extreme low values for number of roots per plant (NRPP):", "\n")
      print(subset(data, NRPP < quantile(NRPP, 0.25, na.rm=T)-3*IQR(NRPP, na.rm=T)))
    }

  if (exists("NRPP", where=data)==1)
    if (dim(subset(data, NRPP > quantile(NRPP, 0.75, na.rm=T)+3*IQR(NRPP, na.rm=T)))[1]>0){
      cat("\n", "- Extreme high values for number of roots per plant (NRPP):", "\n")
      print(subset(data, NRPP > quantile(NRPP, 0.75, na.rm=T)+3*IQR(NRPP, na.rm=T)))
    }

  if (exists("YPP", where=data)==1)
    if (dim(subset(data, YPP < 0))[1]>0){
      cat("\n", "- Out of range values for yield per plant (YPP):", "\n")
      print(subset(data, YPP < 0))
    }

  if (exists("YPP", where=data)==1)
    if (dim(subset(data, YPP < quantile(YPP, 0.25, na.rm=T)-3*IQR(YPP, na.rm=T)))[1]>0){
      cat("\n", "- Extreme low values for yield per plant (YPP):", "\n")
      print(subset(data, YPP < quantile(YPP, 0.25, na.rm=T)-3*IQR(YPP, na.rm=T)))
    }

  if (exists("YPP", where=data)==1)
    if (dim(subset(data, YPP > quantile(YPP, 0.75, na.rm=T)+3*IQR(YPP, na.rm=T)))[1]>0){
      cat("\n", "- Extreme high values for yield per plant (YPP):", "\n")
      print(subset(data, YPP > quantile(YPP, 0.75, na.rm=T)+3*IQR(YPP, na.rm=T)))
    }

  if (exists("CI", where=data)==1)
    if (dim(subset(data, CI < 0 | CI > 100))[1]>0){
      cat("\n", "- Out of range values for commercial index (CI):", "\n")
      print(subset(data, CI < 0 | CI > 100))
    }

  if (exists("CI", where=data)==1)
    if (dim(subset(data, CI < quantile(CI, 0.25, na.rm=T)-3*IQR(CI, na.rm=T)))[1]>0){
      cat("\n", "- Extreme low values for commercial index (CI):", "\n")
      print(subset(data, CI < quantile(CI, 0.25, na.rm=T)-3*IQR(CI, na.rm=T)))
    }

  if (exists("CI", where=data)==1)
    if (dim(subset(data, CI > quantile(CI, 0.75, na.rm=T)+3*IQR(CI, na.rm=T)))[1]>0){
      cat("\n", "- Extreme high values for commercial index (CI):", "\n")
      print(subset(data, CI > quantile(CI, 0.75, na.rm=T)+3*IQR(CI, na.rm=T)))
    }

  if (exists("HI", where=data)==1)
    if (dim(subset(data, HI < 0 | HI > 100))[1]>0){
      cat("\n", "- Out of range values for harvest index (HI):", "\n")
      print(subset(data, HI < 0 | HI > 100))
    }

  if (exists("HI", where=data)==1)
    if (dim(subset(data, HI < quantile(HI, 0.25, na.rm=T)-3*IQR(HI, na.rm=T)))[1]>0){
      cat("\n", "- Extreme low values for harvest index (HI):", "\n")
      print(subset(data, HI < quantile(HI, 0.25, na.rm=T)-3*IQR(HI, na.rm=T)))
    }

  if (exists("HI", where=data)==1)
    if (dim(subset(data, HI > quantile(HI, 0.75, na.rm=T)+3*IQR(HI, na.rm=T)))[1]>0){
      cat("\n", "- Extreme high values for harvest index (HI):", "\n")
      print(subset(data, HI > quantile(HI, 0.75, na.rm=T)+3*IQR(HI, na.rm=T)))
    }

  if (exists("SHI", where=data)==1)
    if (dim(subset(data, SHI < 0 | SHI > 100))[1]>0){
      cat("\n", "- Out of range values for harvest sowing index (SHI):", "\n")
      print(subset(data, SHI < 0 | SHI > 100))
    }

  if (exists("SHI", where=data)==1)
    if (dim(subset(data, SHI < quantile(SHI, 0.25, na.rm=T)-3*IQR(SHI, na.rm=T)))[1]>0){
      cat("\n", "- Extreme low values for harvest sowing index (SHI):", "\n")
      print(subset(data, SHI < quantile(SHI, 0.25, na.rm=T)-3*IQR(SHI, na.rm=T)))
    }

  if (exists("SHI", where=data)==1)
    if (dim(subset(data, SHI > quantile(SHI, 0.75, na.rm=T)+3*IQR(SHI, na.rm=T)))[1]>0){
      cat("\n", "- Extreme high values for harvest sowing index (SHI):", "\n")
      print(subset(data, SHI > quantile(SHI, 0.75, na.rm=T)+3*IQR(SHI, na.rm=T)))
    }

  if (exists("BIOM", where=data)==1)
    if (dim(subset(data, BIOM < 0))[1]>0){
      cat("\n", "- Out of range values for biomass yield (BIOM):", "\n")
      print(subset(data, BIOM < 0))
    }

  if (exists("BIOM", where=data)==1)
    if (dim(subset(data, BIOM < quantile(BIOM, 0.25, na.rm=T)-3*IQR(BIOM, na.rm=T)))[1]>0){
      cat("\n", "- Extreme low values for biomass yield (BIOM):", "\n")
      print(subset(data, BIOM < quantile(BIOM, 0.25, na.rm=T)-3*IQR(BIOM, na.rm=T)))
    }

  if (exists("BIOM", where=data)==1)
    if (dim(subset(data, BIOM > quantile(BIOM, 0.75, na.rm=T)+3*IQR(BIOM, na.rm=T)))[1]>0){
      cat("\n", "- Extreme high values for biomass yield (BIOM):", "\n")
      print(subset(data, BIOM > quantile(BIOM, 0.75, na.rm=T)+3*IQR(BIOM, na.rm=T)))
    }

  if (exists("FYTHA", where=data)==1)
    if (dim(subset(data, FYTHA < 0))[1]>0){
      cat("\n", "- Out of range values for foliage total yield in tons per hectare (FYTHA):", "\n")
      print(subset(data, FYTHA < 0))
    }

  if (exists("FYTHA", where=data)==1)
    if (dim(subset(data, FYTHA < quantile(FYTHA, 0.25, na.rm=T)-3*IQR(FYTHA, na.rm=T)))[1]>0){
      cat("\n", "- Extreme low values for foliage total yield in tons per hectare (FYTHA):", "\n")
      print(subset(data, FYTHA < quantile(FYTHA, 0.25, na.rm=T)-3*IQR(FYTHA, na.rm=T)))
    }

  if (exists("FYTHA", where=data)==1)
    if (dim(subset(data, FYTHA > quantile(FYTHA, 0.75, na.rm=T)+3*IQR(FYTHA, na.rm=T)))[1]>0){
      cat("\n", "- Extreme high values for foliage total yield in tons per hectare (FYTHA):", "\n")
      print(subset(data, FYTHA > quantile(FYTHA, 0.75, na.rm=T)+3*IQR(FYTHA, na.rm=T)))
    }

  if (exists("RFR", where=data)==1)
    if (dim(subset(data, RFR < 0 | RFR > 100))[1]>0){
      cat("\n", "- Out of range values for root foliage ratio (RFR):", "\n")
      print(subset(data, RFR < 0 | RFR > 100))
    }

  if (exists("RFR", where=data)==1)
    if (dim(subset(data, RFR < quantile(RFR, 0.25, na.rm=T)-3*IQR(RFR, na.rm=T)))[1]>0){
      cat("\n", "- Extreme low values for root foliage ratio (RFR):", "\n")
      print(subset(data, RFR < quantile(RFR, 0.25, na.rm=T)-3*IQR(RFR, na.rm=T)))
    }

  if (exists("RFR", where=data)==1)
    if (dim(subset(data, RFR > quantile(RFR, 0.75, na.rm=T)+3*IQR(RFR, na.rm=T)))[1]>0){
      cat("\n", "- Extreme high values for root foliage ratio (RFR):", "\n")
      print(subset(data, RFR > quantile(RFR, 0.75, na.rm=T)+3*IQR(RFR, na.rm=T)))
    }
}
