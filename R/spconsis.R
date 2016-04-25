#' Check consistency for sweetpotato experimental data
#'
#' Set of rules to check for consistency of sweetpotato experimental data.
#' Data labels must be defined as specified in the PROCEDURES FOR THE EVALUATION
#' AND ANALYSIS OF SWEETPOTATO TRIALS document.
#' @param fb The name of the fieldbook data frame.
#' @param plot.size Plot size in square meters.
#' @param f Factor for extreme values detection. See details.
#' @param width Number of columns for the output.
#' @param file Logigal, if TRUE the output goes to a file.
#' @details The data frame must use the labels (lower or upper case) listed in function \code{checknames}.
#' See \code{?checknames} for details.
#' Extreme values are detected using the interquartile range.
#' The rule is to detect any value out of the interval 
#' \eqn{[Q_1 - f \times IQR; Q_3 + f \times IQR]}. By default \code{f = 3}.
#' @return If \code{file = TRUE} it returns a file with name checks.txt with a list of
#' all rows with some kind of inconsistency and all rows with outliers. If \code{file = FALSE}
#' the output is shown in the R console.
#' @author Raul Eyzaguirre.
#' @examples
#'  # The data
#'  head(pjpz09)
#'  str(pjpz09)
#'
#'  # Check the data
#'  spconsis(pjpz09, 4.5)
#' @export

spconsis <- function(fb, plot.size, f = 3, width = 240, file = TRUE) {

  options(width = width)
  
  fb <- checknames(fb)

  if (file == TRUE) sink("checks.txt")

  spconsis01(fb) # NOPS > NOPE > NOPH > NOPR
  spconsis02(fb) # NOPE and dependencies
  spconsis03(fb) # NOPH and VW
  spconsis04(fb) # VW and dependencies
  spconsis05(fb) # NOPR and number of roots
  spconsis06(fb) # Number of roots and root weight
  spconsis07(fb) # TRW, CRW + NCRW, NOCR + NONC, NOPR
  spconsis08(fb) # Roots and dependencies
  spconsis09(fb, plot.size) # Calculated variables
  spconsis10(fb, f) # Outliers detection and values out of range for field data
  spconsis11(fb, f) # Outliers detection and values out of range for DM data
  spconsis12(fb, f) # Outliers detection and values out of range for cooked traits
  spconsis13(fb, f) # Outliers detection and values out of range for lab data
  spconsis14(fb, f) # Outliers detection and values out of range for derived variables

  if (file == TRUE) sink()
}

# Check consistency for sweetpotato experimental data, part 1.
# Inconsistencies for NOPS > NOPE > NOPH > NOPR.

spconsis01 <- function(fb) {

  if (exists("NOPE", fb) & exists("NOPS", fb))
    if (dim(subset(fb, NOPE > NOPS))[1] > 0) {
      cat("\n", "- Number of plants established (NOPE) greater than number of plants sowed (NOPS):", "\n")
      print(subset(fb, NOPE > NOPS))
    }

  if (exists("NOPH", fb) & exists("NOPE", fb))
    if (dim(subset(fb, NOPH > NOPE))[1] > 0) {
      cat("\n", "- Number of plants harvested (NOPH) greater than number of plants established (NOPE):", "\n")
      print(subset(fb, NOPH > NOPE))
    }

  if (exists("NOPH", fb) & exists("NOPS", fb)) {
    if (exists("NOPE", fb)) {
      if (dim(subset(fb, is.na(NOPE) & NOPH > NOPS))[1] > 0) {
        cat("\n", "- Number of plants harvested (NOPH) greater than number of plants sowed (NOPS):", "\n")
        print(subset(fb, is.na(NOPE) & NOPH > NOPS))
      }
    } else {
      if (dim(subset(fb, NOPH > NOPS))[1] > 0) {
        cat("\n", "- Number of plants harvested (NOPH) greater than number of plants sowed (NOPS):", "\n")
        print(subset(fb, NOPH > NOPS))
      }
    }
  }

  if (exists("NOPR", fb) & exists("NOPH", fb))
    if (dim(subset(fb, NOPR > NOPH))[1] > 0) {
      cat("\n", "- Number of plants with roots (NOPR) greater than number of plants harvested (NOPH):", "\n")
      print(subset(fb, NOPR > NOPH))
    }

  if (exists("NOPR", fb) & exists("NOPE", fb)) {
    if (exists("NOPH", fb)) {
      if (dim(subset(fb, is.na(NOPH) & NOPR > NOPE))[1] > 0) {
        cat("\n", "- Number of plants with roots (NOPR) greater than number of plants established (NOPE):", "\n")
        print(subset(fb, is.na(NOPH) & NOPR > NOPE))
      }
    } else {
      if (dim(subset(fb, NOPR > NOPE))[1] > 0) {
        cat("\n", "- Number of plants with roots (NOPR) greater than number of plants established (NOPE):", "\n")
        print(subset(fb, NOPR > NOPE))
      }
    }
  }

  if (exists("NOPR", fb) & exists("NOPS", fb)) {
    if (exists("NOPH", fb) & exists("NOPE", fb)) {
      if (dim(subset(fb, is.na(NOPH) & is.na(NOPE) & NOPR > NOPS))[1] > 0) {
        cat("\n", "- Number of plants with roots (NOPR) greater than number of plants sowed (NOPS):", "\n")
        print(subset(fb, is.na(NOPH) & is.na(NOPE) & NOPR > NOPS))
      }
    } else {
      if (dim(subset(fb, NOPR > NOPS))[1] > 0) {
        cat("\n", "- Number of plants with roots (NOPR) greater than number of plants sowed (NOPS):", "\n")
        print(subset(fb, NOPR > NOPS))
      }
    }
  }
}

# Check consistency for sweetpotato experimental data, part 2.
# Inconsistencies for NOPE and dependencies.

spconsis02 <- function(fb) {

  if (exists("NOPE", fb) & exists("VIR1", fb))
    if (dim(subset(fb, (NOPE == 0 | is.na(NOPE)) & !is.na(VIR1)))[1] > 0) {
      cat("\n", "- Number of plants established (NOPE) is zero or NA but there is data for virus symptoms first evaluation (VIR1):", "\n")
      print(subset(fb, (NOPE == 0 | is.na(NOPE)) & !is.na(VIR1)))
    }

  if (exists("NOPE", fb) & exists("VIR2", fb))
    if (dim(subset(fb, (NOPE == 0 | is.na(NOPE)) & !is.na(VIR2)))[1] > 0) {
      cat("\n", "- Number of plants established (NOPE) is zero or NA but there is data for virus symptoms second evaluation (VIR2):", "\n")
      print(subset(fb, (NOPE == 0 | is.na(NOPE)) & !is.na(VIR2)))
    }

  if (exists("NOPE", fb) & exists("ALT1", fb))
    if (dim(subset(fb, (NOPE == 0 | is.na(NOPE)) & !is.na(ALT1)))[1] > 0) {
      cat("\n", "- Number of plants established (NOPE) is zero or NA but there is data for alternaria symptoms first evaluation (ALT1):", "\n")
      print(subset(fb, (NOPE == 0 | is.na(NOPE)) & !is.na(ALT1)))
    }

  if (exists("NOPE", fb) & exists("ALT2", fb))
    if (dim(subset(fb, (NOPE == 0 | is.na(NOPE)) & !is.na(ALT2)))[1] > 0) {
      cat("\n", "- Number of plants established (NOPE) is zero or NA but there is data for alternaria symptoms second evaluation (ALT2):", "\n")
      print(subset(fb, (NOPE == 0 | is.na(NOPE)) & !is.na(ALT2)))
    }

  if (exists("NOPE", fb) & exists("VV1", fb))
    if (dim(subset(fb, (NOPE == 0 | is.na(NOPE)) & !is.na(VV1)))[1] > 0) {
      cat("\n", "- Number of plants established (NOPE) is zero or NA but there is data for vine vigor first evaluation (VV1):", "\n")
      print(subset(fb, (NOPE == 0 | is.na(NOPE)) & !is.na(VV1)))
    }
}

# Check consistency for sweetpotato experimental data, part 3.
# Inconsistencies for NOPH and VW.

spconsis03 <- function(fb) {

  if (exists("NOPH", fb) & exists("VW", fb))
    if (dim(subset(fb, (NOPH == 0 | is.na(NOPH)) & VW > 0))[1] > 0) {
      cat("\n", "- Number of plants harvested (NOPH) is zero or NA but vine weight (VW) is greater than zero:", "\n")
      print(subset(fb, (NOPH == 0 | is.na(NOPH)) & VW > 0))
    }

  if (exists("NOPH", fb) & exists("VW", fb))
    if (dim(subset(fb, NOPH > 0 & (VW == 0 | is.na(VW))))[1] > 0) {
      cat("\n", "- Vine weight (VW) is zero or NA but the number of plants harvested (NOPH) is greater than zero:", "\n")
      print(subset(fb, NOPH > 0 & (VW == 0 | is.na(VW))))
    }
}

# Check consistency for sweetpotato experimental data, part 4.
# Inconsistencies for VW and dependencies.

spconsis04 <- function(fb) {

  if (exists("VW", fb) & exists("DMVF", fb))
    if (dim(subset(fb, (VW == 0 | is.na(VW)) & DMVF > 0))[1] > 0) {
      cat("\n", "- Vine weight (VW) is zero or NA but there is fresh weight vines for dry matter assessment (DMVF):", "\n")
      print(subset(fb, (VW == 0 | is.na(VW)) & DMVF > 0))
    }

  if (exists("VW", fb) & exists("DMVD", fb))
    if (dim(subset(fb, (VW == 0 | is.na(VW)) & DMVD > 0))[1] > 0) {
      cat("\n", "- Vine weight (VW) is zero or NA but there is dry weight vines for dry matter assessment (DMVD):", "\n")
      print(subset(fb, (VW == 0 | is.na(VW)) & DMVD > 0))
    }

  if (exists("DMVF", fb) & exists("DMVD", fb))
    if (dim(subset(fb, DMVD > DMVF))[1] > 0) {
      cat("\n", "- Dry weight vines for dry matter assessment (DMVD) is greater than fresh weight vines for dry matter assessment (DBVF):", "\n")
      print(subset(fb, DMVD > DMVF))
    }

  if (exists("VW", fb) & exists("VV2", fb))
    if (dim(subset(fb, (VW == 0 | is.na(VW)) & !is.na(VV2)))[1] > 0) {
      cat("\n", "- Vine weight (VW) is zero or NA but there is data for vine vigor second evaluation (VV2):", "\n")
      print(subset(fb, (VW == 0 | is.na(VW)) & !is.na(VV2)))
    }

  if (exists("VW", fb) & exists("VIR3", fb))
    if (dim(subset(fb, (VW == 0 | is.na(VW)) & !is.na(VIR3)))[1] > 0) {
      cat("\n", "- Vine weight (VW) is zero or NA but there is data for virus symptoms third evaluation (VIR3):", "\n")
      print(subset(fb, (VW == 0 | is.na(VW)) & !is.na(VIR3)))
    }
}

# Check consistency for sweetpotato experimental data, part 5.
# Inconsistencies for NOPR and number of roots.

spconsis05 <- function(fb) {

  if (exists("NOPR", fb) & exists("NOCR", fb) & exists("NONC", fb))
    if (dim(subset(fb, (NOPR == 0 | is.na(NOPR)) & (NOCR > 0 | NONC > 0)))[1] > 0) {
      cat("\n", "- Number of plants with roots (NOPR) is zero or NA but number of roots (NOCR + NONC) is greater than zero:", "\n")
      print(subset(fb, (NOPR == 0 | is.na(NOPR)) & (NOCR > 0 | NONC > 0)))
    }

  if (exists("NOPR", fb) & exists("NOCR", fb) & exists("NONC", fb))
    if (dim(subset(fb, NOPR > 0 & (suma(NOCR, NONC) == 0 | is.na(suma(NOCR, NONC)))))[1] > 0) {
      cat("\n", "- Number of roots (NOCR + NONC) is zero or NA but number of plants with roots (NOPR) is greater than zero:", "\n")
      print(subset(fb, NOPR > 0 & (suma(NOCR, NONC) == 0 | is.na(suma(NOCR, NONC)))))
    }
}

# Check consistency for sweetpotato experimental data, part 6.
# Inconsistencies for number of roots and root weight.

spconsis06 <- function(fb) {

  if (exists("NOCR", fb) & exists("CRW", fb))
    if (dim(subset(fb, (NOCR == 0 | is.na(NOCR)) & CRW > 0))[1] > 0) {
      cat("\n", "- Number of commercial roots (NOCR) is zero or NA but the commercial root weight (CRW) is greater than zero:", "\n")
      print(subset(fb, (NOCR == 0 | is.na(NOCR)) & CRW > 0))
    }

  if (exists("NOCR", fb) & exists("CRW", fb))
    if (dim(subset(fb, NOCR > 0 & (CRW == 0 | is.na(CRW))))[1] > 0) {
      cat("\n", "- Commercial root weight (CRW) is zero or NA but the number of commercial roots (NOCR) is greater than zero:", "\n")
      print(subset(fb, NOCR > 0 & (CRW == 0 | is.na(CRW))))
    }

  if (exists("NONC", fb) & exists("NCRW", fb))
    if (dim(subset(fb, (NONC == 0 | is.na(NONC)) & NCRW > 0))[1] > 0) {
      cat("\n", "- Number of non commercial roots (NONC) is zero or NA but the non commercial root weight (NCRW) is greater than zero:", "\n")
      print(subset(fb, (NONC == 0 | is.na(NONC)) & NCRW > 0))
    }

  if (exists("NONC", fb) & exists("NCRW", fb))
    if (dim(subset(fb, NONC > 0 & (NCRW == 0 | is.na(NCRW))))[1] > 0) {
      cat("\n", "- Non commercial root weight (NCRW) is zero or NA but the number of non commercial roots (NONC) is greater than zero:", "\n")
      print(subset(fb, NONC > 0 & (NCRW == 0 | is.na(NCRW))))
    }
}

# Check consistency for sweetpotato experimental data, part 7.
# Inconsistencies for TRW, CRW + NCRW, NOCR + NONC, NOPR.

spconsis07 <- function(fb) {

  if (exists("TRW", fb) & exists("NOPR", fb))
    if (dim(subset(fb, (TRW == 0 | is.na(TRW)) & NOPR > 0))[1] > 0) {
      cat("\n", "- Total root weight (TRW) is zero or NA but number of plants with roots (NOPR) is greater than zero:", "\n")
      print(subset(fb, (TRW == 0 | is.na(TRW)) & NOPR > 0))
    }

  if (exists("TRW", fb) & exists("CRW", fb) & exists("NCRW", fb))
    if (dim(subset(fb, (TRW == 0 | is.na(TRW)) & suma(CRW, NCRW) > 0))[1] > 0) {
      cat("\n", "- Total root weight (TRW) is zero or NA but root weight (CRW + NCRW) is greater than zero:", "\n")
      print(subset(fb, (TRW == 0 | is.na(TRW)) & suma(CRW, NCRW) > 0))
    }

  if (exists("TRW", fb) & exists("NOCR", fb) & exists("NONC", fb))
    if (dim(subset(fb, (TRW == 0 | is.na(TRW)) & suma(NOCR, NONC) > 0))[1] > 0) {
      cat("\n", "- Total root weight (TRW) is zero or NA but number of roots (NOCR + NONC) is greater than zero:", "\n")
      print(subset(fb, (TRW == 0 | is.na(TRW)) & suma(NOCR, NONC) > 0))
    }

  if (exists("CRW", fb) & exists("NCRW", fb) & exists("NOPR", fb))
    if (dim(subset(fb, NOPR > 0 & (suma(CRW, NCRW) == 0 | is.na(suma(CRW, NCRW)))))[1] > 0) {
      cat("\n", "- Root weight (CRW + NCRW) is zero or NA but number of plants with roots (NOPR) is greater than zero:", "\n")
      print(subset(fb, NOPR > 0 & (suma(CRW, NCRW) == 0 | is.na(suma(CRW, NCRW)))))
    }
}

# Check consistency for sweetpotato experimental data, part 8.
# Inconsistencies for roots and dependencies.

spconsis08 <- function(fb) {

  if (exists("NOPR", fb))
    fb$RAUX <- fb$NOPR else
      if (exists("NOCR", fb) & exists("NONC", fb))
        fb$RAUX <- suma(fb$NOCR, fb$NONC) else
          if (exists("CRW", fb) & exists("NCRW", fb))
            fb$RAUX <- suma(fb$CRW, fb$NCRW) else
              if (exists("TRW", fb))
                fb$RAUX <- fb$TRW else
                  if (exists("RYTHA", fb))
                    fb$RAUX <- fb$RYTHA else
                      if (exists("CRW", fb))
                        fb$RAUX <- fb$CRW else
                          if (exists("CYTHA", fb))
                            fb$RAUX <- fb$CYTHA
                          
  if (exists("RAUX", fb) & exists("RFCP", fb))
    if (dim(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(RFCP)))[1] > 0) {
      cat("\n", "- There are no roots but there is data for root primary flesh color (RFCP):", "\n")
      print(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(RFCP), select = -RAUX))
    }

  if (exists("RAUX", fb) & exists("RFCS", fb))
    if (dim(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(RFCS)))[1] > 0) {
      cat("\n", "- There are no roots but there is data for root secondary flesh color (RFCS):", "\n")
      print(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(RFCS), select = -RAUX))
    }

  if (exists("RAUX", fb) & exists("SCOL", fb))
    if (dim(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(SCOL)))[1] > 0) {
      cat("\n", "- There are no roots but there is data for storage root skin color (SCOL):", "\n")
      print(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(SCOL), select = -RAUX))
    }

  if (exists("RAUX", fb) & exists("FCOL", fb))
    if (dim(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(FCOL)))[1] > 0) {
      cat("\n", "- There are no roots but there is data for storage root flesh color (FCOL):", "\n")
      print(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(FCOL), select = -RAUX))
    }

  if (exists("RAUX", fb) & exists("RS", fb))
    if (dim(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(RS)))[1] > 0) {
      cat("\n", "- There are no roots but there is data for root size (RS):", "\n")
      print(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(RS), select = -RAUX))
    }

  if (exists("RAUX", fb) & exists("RF", fb))
    if (dim(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(RF)))[1] > 0) {
      cat("\n", "- There are no roots but there is data for root form (RF):", "\n")
      print(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(RF), select = -RAUX))
    }

  if (exists("RAUX", fb) & exists("DAMR", fb))
    if (dim(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(DAMR)))[1] > 0) {
      cat("\n", "- There are no roots but there is data for root defects (DAMR):", "\n")
      print(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(DAMR), select = -RAUX))
    }

  if (exists("RAUX", fb) & exists("RSPR", fb))
    if (dim(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(RSPR)))[1] > 0) {
      cat("\n", "- There are no roots but there is data for root sprouting (RSPR):", "\n")
      print(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(RSPR), select = -RAUX))
    }

  if (exists("RAUX", fb) & exists("WED1", fb))
    if (dim(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(WED1)))[1] > 0) {
      cat("\n", "- There are no roots but there is data for weevil damage first evaluation (WED1):", "\n")
      print(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(WED1), select = -RAUX))
    }

  if (exists("RAUX", fb) & exists("WED2", fb))
    if (dim(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(WED2)))[1] > 0) {
      cat("\n", "- There are no roots but there is data for weevil damage second evaluation (WED2):", "\n")
      print(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(WED2), select = -RAUX))
    }

  if (exists("RAUX", fb) & exists("DMF", fb))
    if (dim(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(DMF)))[1] > 0) {
      cat("\n", "- There are no roots but there is data for fresh weight of roots for dry matter assessment (DMF):", "\n")
      print(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(DMF), select = -RAUX))
    }

  if (exists("RAUX", fb) & exists("DMD", fb))
    if (dim(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(DMD)))[1] > 0) {
      cat("\n", "- There are no roots but there is data for dry weight of roots for dry matter assessment (DMD):", "\n")
      print(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(DMD), select = -RAUX))
    }

  if (exists("DMF", fb) & exists("DMD", fb))
    if (dim(subset(fb, DMF < DMD))[1] > 0) {
      cat("\n", "- Dry weight of roots for dry matter assessment (DMD) is greater than fresh weight of roots for dry matter assessment (DMF):", "\n")
      print(subset(fb, DMF < DMD, select = -RAUX))
    }

  if (exists("RAUX", fb) & exists("FRAW1", fb))
    if (dim(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(FRAW1)))[1] > 0) {
      cat("\n", "- There are no roots but there is data for root fiber first determination (FRAW1):", "\n")
      print(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(FRAW1), select = -RAUX))
    }

  if (exists("RAUX", fb) & exists("SURAW1", fb))
    if (dim(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(SURAW1)))[1] > 0) {
      cat("\n", "- There are no roots but there is data for root sugar first determination (SURAW1):", "\n")
      print(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(SURAW1), select = -RAUX))
    }

  if (exists("RAUX", fb) & exists("STRAW1", fb))
    if (dim(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(STRAW1)))[1] > 0) {
      cat("\n", "- There are no roots but there is data for root starch first determination (STRAW1):", "\n")
      print(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(STRAW1), select = -RAUX))
    }

  if (exists("RAUX", fb) & exists("COOF1", fb))
    if (dim(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(COOF1)))[1] > 0) {
      cat("\n", "- There are no roots but there is data for cooked fiber first evaluation (COOF1):", "\n")
      print(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(COOF1), select = -RAUX))
    }

  if (exists("RAUX", fb) & exists("COOSU1", fb))
    if (dim(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(COOSU1)))[1] > 0) {
      cat("\n", "- There are no roots but there is data for cooked sugars first evaluation (COOSU1):", "\n")
      print(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(COOSU1), select = -RAUX))
    }

  if (exists("RAUX", fb) & exists("COOST1", fb))
    if (dim(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(COOST1)))[1] > 0) {
      cat("\n", "- There are no roots but there is data for cooked starch first evaluation (COOST1):", "\n")
      print(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(COOST1), select = -RAUX))
    }

  if (exists("RAUX", fb) & exists("COOT1", fb))
    if (dim(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(COOT1)))[1] > 0) {
      cat("\n", "- There are no roots but there is data for cooked taste first evaluation (COOT1):", "\n")
      print(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(COOT1), select = -RAUX))
    }

  if (exists("RAUX", fb) & exists("COOAP1", fb))
    if (dim(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(COOAP1)))[1] > 0) {
      cat("\n", "- There are no roots but there is data for cooked appearance first evaluation (COOAP1):", "\n")
      print(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(COOAP1), select = -RAUX))
    }

  if (exists("RAUX", fb) & exists("FRAW2", fb))
    if (dim(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(FRAW2)))[1] > 0) {
      cat("\n", "- There are no roots but there is data for root fiber second determination (FRAW2):", "\n")
      print(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(FRAW2), select = -RAUX))
    }

  if (exists("RAUX", fb) & exists("SURAW2", fb))
    if (dim(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(SURAW2)))[1] > 0) {
      cat("\n", "- There are no roots but there is data for root sugar second determination (SURAW2):", "\n")
      print(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(SURAW2), select = -RAUX))
    }

  if (exists("RAUX", fb) & exists("STRAW2", fb))
    if (dim(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(STRAW2)))[1] > 0) {
      cat("\n", "- There are no roots but there is data for root starch second determination (STRAW2):", "\n")
      print(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(STRAW2), select = -RAUX))
    }

  if (exists("RAUX", fb) & exists("COOF2", fb))
    if (dim(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(COOF2)))[1] > 0) {
      cat("\n", "- There are no roots but there is data for cooked fiber second evaluation (COOF2):", "\n")
      print(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(COOF2), select = -RAUX))
    }

  if (exists("RAUX", fb) & exists("COOSU2", fb))
    if (dim(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(COOSU2)))[1] > 0) {
      cat("\n", "- There are no roots but there is data for cooked sugars second evaluation (COOSU2):", "\n")
      print(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(COOSU2), select = -RAUX))
    }

  if (exists("RAUX", fb) & exists("COOST2", fb))
    if (dim(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(COOST2)))[1] > 0) {
      cat("\n", "- There are no roots but there is data for cooked starch second evaluation (COOST2):", "\n")
      print(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(COOST2), select = -RAUX))
    }

  if (exists("RAUX", fb) & exists("COOT2", fb))
    if (dim(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(COOT2)))[1] > 0) {
      cat("\n", "- There are no roots but there is data for cooked taste second evaluation (COOT2):", "\n")
      print(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(COOT2), select = -RAUX))
    }

  if (exists("RAUX", fb) & exists("COOAP2", fb))
    if (dim(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(COOAP2)))[1] > 0) {
      cat("\n", "- There are no roots but there is data for cooked appearance second evaluation (COOAP2):", "\n")
      print(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(COOAP2), select = -RAUX))
    }

  if (exists("RAUX", fb) & exists("PROT", fb))
    if (dim(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(PROT)))[1] > 0) {
      cat("\n", "- There are no roots but there is data for protein (PROT):", "\n")
      print(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(PROT), select = -RAUX))
    }

  if (exists("RAUX", fb) & exists("FE", fb))
    if (dim(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(FE)))[1] > 0) {
      cat("\n", "- There are no roots but there is data for iron in dry weight (FE):", "\n")
      print(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(FE), select = -RAUX))
    }

  if (exists("RAUX", fb) & exists("ZN", fb))
    if (dim(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(ZN)))[1] > 0) {
      cat("\n", "- There are no roots but there is data for zinc in dry weight (ZN):", "\n")
      print(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(ZN), select = -RAUX))
    }

  if (exists("RAUX", fb) & exists("CA", fb))
    if (dim(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(CA)))[1] > 0) {
      cat("\n", "- There are no roots but there is data for calcium in dry weight (CA):", "\n")
      print(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(CA), select = -RAUX))
    }

  if (exists("RAUX", fb) & exists("MG", fb))
    if (dim(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(MG)))[1] > 0) {
      cat("\n", "- There are no roots but there is data for magnesium in dry weight (MG):", "\n")
      print(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(MG), select = -RAUX))
    }

  if (exists("RAUX", fb) & exists("BC", fb))
    if (dim(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(BC)))[1] > 0) {
      cat("\n", "- There are no roots but there is data for beta-carotene in dry weight (BC):", "\n")
      print(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(BC), select = -RAUX))
    }

  if (exists("RAUX", fb) & exists("TC", fb))
    if (dim(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(TC)))[1] > 0) {
      cat("\n", "- There are no roots but there is data for total carotenoids in dry weight (TC):", "\n")
      print(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(TC), select = -RAUX))
    }

  if (exists("RAUX", fb) & exists("STAR", fb))
    if (dim(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(STAR)))[1] > 0) {
      cat("\n", "- There are no roots but there is data for starch (STAR):", "\n")
      print(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(STAR), select = -RAUX))
    }

  if (exists("RAUX", fb) & exists("FRUC", fb))
    if (dim(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(FRUC)))[1] > 0) {
      cat("\n", "- There are no roots but there is data for fructose (FRUC):", "\n")
      print(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(FRUC), select = -RAUX))
    }

  if (exists("RAUX", fb) & exists("GLUC", fb))
    if (dim(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(GLUC)))[1] > 0) {
      cat("\n", "- There are no roots but there is data for glucose (GLUC):", "\n")
      print(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(GLUC), select = -RAUX))
    }

  if (exists("RAUX", fb) & exists("SUCR", fb))
    if (dim(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(SUCR)))[1] > 0) {
      cat("\n", "- There are no roots but there is data for sucrose (SUCR):", "\n")
      print(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(SUCR), select = -RAUX))
    }

  if (exists("RAUX", fb) & exists("MALT", fb))
    if (dim(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(MALT)))[1] > 0) {
      cat("\n", "- There are no roots but there is data for maltose (MALT):", "\n")
      print(subset(fb, (RAUX == 0 | is.na(RAUX)) & !is.na(MALT), select = -RAUX))
    }

  if (exists("RAUX", fb) & exists("TRW", fb))
    if (dim(subset(fb, (RAUX == 0 | is.na(RAUX)) & TRW > 0))[1] > 0) {
      cat("\n", "- There are no roots but total root weight (TRW) is greater than zero:", "\n")
      print(subset(fb, (RAUX == 0 | is.na(RAUX)) & TRW > 0, select = -RAUX))
    }

  if (exists("RAUX", fb) & exists("RYTHA", fb))
    if (dim(subset(fb, (RAUX == 0 | is.na(RAUX)) & RYTHA > 0))[1] > 0) {
      cat("\n", "- There are no roots but total root yield (RYTHA) is greater than zero:", "\n")
      print(subset(fb, (RAUX == 0 | is.na(RAUX)) & RYTHA > 0, select = -RAUX))
    }

  fb$RAUX <- NULL
}

# Check consistency for sweetpotato experimental data, part 9.
# Inconsistencies for calculated variables.

spconsis09 <- function(fb, plot.size) {

  if (exists("TRW", fb) & exists("CRW", fb) & exists("NCRW", fb))
    if (dim(subset(fb, abs(TRW - suma(CRW, NCRW)) > 1e-10))[1] > 0) {
      cat("\n", "- Total root weight (TRW) different from CRW + NCRW:", "\n")
      print(subset(fb, abs(TRW - suma(CRW, NCRW)) > 1e-10))
    }

  if (exists("CYTHA", fb) & exists("CRW", fb))
    if (dim(subset(fb, abs(CYTHA - CRW * 10 / plot.size) > 1e-10))[1] > 0) {
      cat("\n", "- Commercial root yield in tons per hectare (CYTHA) is different from CRW * 10 / plot.size:", "\n")
      print(subset(fb, abs(CYTHA - CRW * 10 / plot.size) > 1e-10))
    }

  if (exists("RYTHA", fb) & exists("CRW", fb) & exists("NCRW", fb))
    if (dim(subset(fb, abs(RYTHA - suma(CRW, NCRW) * 10 / plot.size) > 1e-10))[1] > 0) {
      cat("\n", "- Total root yield in tons per hectare (RYTHA) is different from (CRW + NCRW) * 10 / plot.size:", "\n")
      print(subset(fb, abs(RYTHA - suma(CRW, NCRW) * 10 / plot.size) > 1e-10))
    }

  if (exists("ACRW", fb) & exists("CRW", fb) & exists("NOCR", fb))
    if (dim(subset(fb, abs(ACRW - CRW / NOCR) > 1e-10))[1] > 0) {
      cat("\n", "- Average commercial root weight (ACRW) is different from CRW / NOCR:", "\n")
      print(subset(fb, abs(ACRW - CRW / NOCR) > 1e-10))
    }

  if (exists("NRPP", fb) & exists("NOCR", fb) & exists("NONC", fb) & exists("NOPH", fb))
    if (dim(subset(fb, abs(NRPP - suma(NOCR, NONC) / NOPH) > 1e-10))[1] > 0) {
      cat("\n", "- Number of roots per plant (NRPP) is different from (NOCR + NONC) / NOPH:", "\n")
      print(subset(fb, abs(NRPP - suma(NOCR, NONC) / NOPH) > 1e-10))
    }

  if (exists("YPP", fb) & exists("CRW", fb) & exists("NCRW", fb) & exists("NOPH", fb))
    if (dim(subset(fb, abs(YPP - suma(CRW, NCRW) / NOPH) > 1e-10))[1] > 0) {
      cat("\n", "- Yield per plant (YPP) is different from (CRW + NCRW) / NOPH:", "\n")
      print(subset(fb, abs(YPP - suma(CRW, NCRW) / NOPH) > 1e-10))
    }

  if (exists("CI", fb) & exists("NOCR", fb) & exists("NONC", fb))
    if (dim(subset(fb, abs(CI - NOCR / suma(NOCR, NONC) * 100) > 1e-10))[1] > 0) {
      cat("\n", "- Commercial index (CI) is different from NOCR / (NOCR + NONC) * 100:", "\n")
      print(subset(fb, abs(CI - NOCR / suma(NOCR, NONC) * 100) > 1e-10))
    }

  if (exists("HI", fb) & exists("CRW", fb) & exists("NCRW", fb) & exists("VW", fb))
    if (dim(subset(fb, abs(HI - suma(CRW, NCRW) / suma(suma(VW, CRW), NCRW) * 100) > 1e-10))[1] > 0) {
      cat("\n", "- Harvest index (HI) is different from (CRW + NCRW) / (VW + CRW + NCRW) * 100:", "\n")
      print(subset(fb, abs(HI - suma(CRW, NCRW) / suma(suma(VW, CRW), NCRW) * 100) > 1e-10))
    }

  if (exists("SHI", fb) & exists("NOPH", fb) & exists("NOPS", fb))
    if (dim(subset(fb, abs(SHI - NOPH / NOPS * 100) > 1e-10))[1] > 0) {
      cat("\n", "- Harvest sowing index (SHI) is different from NOPH / NOPS * 100:", "\n")
      print(subset(fb, abs(SHI - NOPH / NOPS * 100) > 1e-10))
    }

  if (exists("BIOM", fb) & exists("CRW", fb) & exists("NCRW", fb) & exists("VW", fb))
    if (dim(subset(fb, abs(BIOM - suma(suma(VW, CRW), NCRW) * 10 / plot.size) > 1e-10))[1] > 0) {
      cat("\n", "- Biomass yield (BIOM) is different from (CRW + NCRW + VW) * 10 / plot.size:", "\n")
      print(subset(fb, abs(BIOM - suma(suma(VW, CRW), NCRW) * 10 / plot.size) > 1e-10))
    }

  if (exists("FYTHA", fb) & exists("VW", fb))
    if (dim(subset(fb, abs(FYTHA - VW * 10 / plot.size) > 1e-10))[1] > 0) {
      cat("\n", "- Foliage total yield in tons per hectare (FYTHA) is different from VW * 10 / plot.size:", "\n")
      print(subset(fb, abs(FYTHA - VW * 10 / plot.size) > 1e-10))
    }

  if (exists("DM", fb) & exists("DMD", fb) & exists("DMF", fb))
    if (dim(subset(fb, abs(DM - DMD / DMF * 100) > 1e-10))[1] > 0) {
      cat("\n", "- Storage root dry matter content (DM) is different from DMD / DMF * 100:", "\n")
      print(subset(fb, abs(DM - DMD / DMF * 100) > 1e-10))
    }

  if (exists("DMV", fb) & exists("DMVD", fb) & exists("DMVF", fb))
    if (dim(subset(fb, abs(DMV - DMVD / DMVF * 100) > 1e-10))[1] > 0) {
      cat("\n", "- Vine dry matter content (DMV) is different from DMVD / DMVF * 100:", "\n")
      print(subset(fb, abs(DMV - DMVD / DMVF * 100) > 1e-10))
    }

    if (exists("DMFY", fb) & exists("VW", fb) & exists("DMVD", fb) & exists("DMVF", fb))
    if (dim(subset(fb, abs(DMFY - VW * 10 / plot.size * DMVD / DMVF) > 1e-10))[1] > 0) {
      cat("\n", "- Dry matter foliage yield (DMFY) is different from VW * 10 / plot.size * DMVD / DMVF:", "\n")
      print(subset(fb, abs(DMFY - VW * 10 / plot.size * DMVD / DMVF) > 1e-10))
    }

  if (exists("DMRY", fb) & exists("CRW", fb) & exists("NCRW", fb) & exists("DMD", fb) & exists("DMF", fb))
    if (dim(subset(fb, DMRY != suma(CRW, NCRW) * 10 / plot.size * DMD / DMF))[1] > 0) {
      cat("\n", "- Dry matter root yield (DMRY) is different from (CRW + NCRW) * 10 / plot.size * DMD / DMF:", "\n")
      print(subset(fb, DMRY != suma(CRW, NCRW) * 10 / plot.size * DMD / DMF))
    }

  if (exists("RFR", fb) & exists("CRW", fb) & exists("NCRW", fb) & exists("DMD", fb)
      & exists("DMF", fb) & exists("VW", fb) & exists("DMVD", fb) & exists("DMVF", fb))
    if (dim(subset(fb, abs(RFR - suma(CRW, NCRW) * (DMD / DMF) / (VW * DMVD / DMVF)) > 1e-10))[1] > 0) {
      cat("\n", "- Root foliage ratio (RFR) is different from (CRW + NCRW) * (DMD / DMF) / (VW * DMVD / DMVF) * 100:", "\n")
      print(subset(fb, abs(RFR - suma(CRW, NCRW) * (DMD / DMF) / (VW * DMVD / DMVF)) > 1e-10))
    }
}

# Check consistency for sweetpotato experimental data, part 10.
# Outliers detection based on interquartile range and values out of range for field data.

spconsis10 <- function(fb, f) {
  
  if (exists("NOPE", fb))
    if (dim(subset(fb, !(NOPE %in% c(0:100, NA))))[1] > 0) {
      cat("\n", "- Out of range values for number of plants established (NOPE):", "\n")
      print(subset(fb, !(NOPE %in% c(0:100, NA))))
    }

  if (exists("VIR1", fb))
    if (dim(subset(fb, !(VIR1 %in% c(1:9, NA))))[1] > 0) {
      cat("\n", "- Out of range values for virus symptoms first evaluation (VIR1):", "\n")
      print(subset(fb, !(VIR1 %in% c(1:9, NA))))
    }

  if (exists("VIR2", fb))
    if (dim(subset(fb, !(VIR2 %in% c(1:9, NA))))[1] > 0) {
      cat("\n", "- Out of range values for virus symptoms second evaluation (VIR2):", "\n")
      print(subset(fb, !(VIR2 %in% c(1:9, NA))))
    }

  if (exists("VIR3", fb))
    if (dim(subset(fb, !(VIR3 %in% c(1:9, NA))))[1] > 0) {
      cat("\n", "- Out of range values for virus symptoms third evaluation (VIR3):", "\n")
      print(subset(fb, !(VIR3 %in% c(1:9, NA))))
    }

  if (exists("ALT1", fb))
    if (dim(subset(fb, !(ALT1 %in% c(1:9, NA))))[1] > 0) {
      cat("\n", "- Out of range values for alternaria symptoms first evaluation (ALT1):", "\n")
      print(subset(fb, !(ALT1 %in% c(1:9, NA))))
    }

  if (exists("ALT2", fb))
    if (dim(subset(fb, !(ALT2 %in% c(1:9, NA))))[1] > 0) {
      cat("\n", "- Out of range values for alternaria symptoms second evaluation (ALT2):", "\n")
      print(subset(fb, !(ALT2 %in% c(1:9, NA))))
    }

  if (exists("VV1", fb))
    if (dim(subset(fb, !(VV1 %in% c(1:9, NA))))[1] > 0) {
      cat("\n", "- Out of range values for vine vigor first evaluation (VV1):", "\n")
      print(subset(fb, !(VV1 %in% c(1:9, NA))))
    }

  if (exists("VV2", fb))
    if (dim(subset(fb, !(VV2 %in% c(1:9, NA))))[1] > 0) {
      cat("\n", "- Out of range values for vine vigor second evaluation (VV2):", "\n")
      print(subset(fb, !(VV2 %in% c(1:9, NA))))
    }

  if (exists("VW", fb))
    if (dim(subset(fb, VW < 0))[1] > 0) {
      cat("\n", "- Out of range values for vine weight (VW):", "\n")
      print(subset(fb, VW < 0))
    }

  if (exists("VW", fb))
    if (dim(subset(fb, VW < quantile(VW, 0.25, na.rm = TRUE) - f * IQR(VW, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme low values for vine weight (VW):", "\n")
      print(subset(fb, VW < quantile(VW, 0.25, na.rm = TRUE) - f * IQR(VW, na.rm = TRUE)))
    }

  if (exists("VW", fb))
    if (dim(subset(fb, VW > quantile(VW, 0.75, na.rm = TRUE) + f * IQR(VW, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme high values for vine weight (VW):", "\n")
      print(subset(fb, VW > quantile(VW, 0.75, na.rm = TRUE) + f * IQR(VW, na.rm = TRUE)))
    }

  if (exists("NOPH", fb))
    if (dim(subset(fb, !(NOPH %in% c(0:100, NA))))[1] > 0) {
      cat("\n", "- Out of range values for number of plants harvested (NOPH):", "\n")
      print(subset(fb, !(NOPE %in% c(0:100, NA))))
    }

  if (exists("NOPR", fb))
    if (dim(subset(fb, !(NOPR %in% c(0:100, NA))))[1] > 0) {
      cat("\n", "- Out of range values for number of plants with roots (NOPR):", "\n")
      print(subset(fb, !(NOPR %in% c(0:100, NA))))
    }

  if (exists("NOCR", fb))
    if (dim(subset(fb, !(NOCR %in% c(0:1000, NA))))[1] > 0) {
      cat("\n", "- Out of range values for number of commercial roots (NOCR):", "\n")
      print(subset(fb, !(NOCR %in% c(0:1000, NA))))
    }

  if (exists("NOCR", fb))
    if (dim(subset(fb, NOCR < quantile(NOCR, 0.25, na.rm = TRUE) - f * IQR(NOCR, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme low values for number of commercial roots (NOCR):", "\n")
      print(subset(fb, NOCR < quantile(NOCR, 0.25, na.rm = TRUE) - f * IQR(NOCR, na.rm = TRUE)))
    }

  if (exists("NOCR", fb))
    if (dim(subset(fb, NOCR > quantile(NOCR, 0.75, na.rm = TRUE) + f * IQR(NOCR, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme high values for number of commercial roots (NOCR):", "\n")
      print(subset(fb, NOCR > quantile(NOCR, 0.75, na.rm = TRUE) + f * IQR(NOCR, na.rm = TRUE)))
    }

  if (exists("NONC", fb))
    if (dim(subset(fb, !(NONC %in% c(0:1000, NA))))[1] > 0) {
      cat("\n", "- Out of range values for number of non commercial roots (NONC):", "\n")
      print(subset(fb, !(NONC %in% c(0:1000, NA))))
    }

  if (exists("NONC", fb))
    if (dim(subset(fb, NONC < quantile(NONC, 0.25, na.rm = TRUE) - f * IQR(NONC, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme low values for number of non commercial roots (NONC):", "\n")
      print(subset(fb, NONC < quantile(NONC, 0.25, na.rm = TRUE) - f * IQR(NONC, na.rm = TRUE)))
    }

  if (exists("NONC", fb))
    if (dim(subset(fb, NONC > quantile(NONC, 0.75, na.rm = TRUE) + f * IQR(NONC, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme high values for number of non commercial roots (NONC):", "\n")
      print(subset(fb, NONC > quantile(NONC, 0.75, na.rm = TRUE) + f * IQR(NONC, na.rm = TRUE)))
    }

  if (exists("CRW", fb))
    if (dim(subset(fb, CRW < 0))[1] > 0) {
      cat("\n", "- Out of range values for commercial root weight (CRW):", "\n")
      print(subset(fb, CRW < 0))
    }

  if (exists("CRW", fb))
    if (dim(subset(fb, CRW < quantile(CRW, 0.25, na.rm = TRUE) - f * IQR(CRW, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme low values for commercial root weight (CRW):", "\n")
      print(subset(fb, CRW < quantile(CRW, 0.25, na.rm = TRUE) - f * IQR(CRW, na.rm = TRUE)))
    }

  if (exists("CRW", fb))
    if (dim(subset(fb, CRW > quantile(CRW, 0.75, na.rm = TRUE) + f * IQR(CRW, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme high values for commercial root weight (CRW):", "\n")
      print(subset(fb, CRW > quantile(CRW, 0.75, na.rm = TRUE) + f * IQR(CRW, na.rm = TRUE)))
    }

  if (exists("NCRW", fb))
    if (dim(subset(fb, NCRW < 0))[1] > 0) {
      cat("\n", "- Out of range values for non commercial root weight (NCRW):", "\n")
      print(subset(fb, NCRW < 0))
    }

  if (exists("NCRW", fb))
    if (dim(subset(fb, NCRW < quantile(NCRW, 0.25, na.rm = TRUE) - f * IQR(NCRW, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme low values for non commercial root weight (NCRW):", "\n")
      print(subset(fb, NCRW < quantile(NCRW, 0.25, na.rm = TRUE) - f * IQR(NCRW, na.rm = TRUE)))
    }

  if (exists("NCRW", fb))
    if (dim(subset(fb, NCRW > quantile(NCRW, 0.75, na.rm = TRUE) + f * IQR(NCRW, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme high values for non commercial root weight (NCRW):", "\n")
      print(subset(fb, NCRW > quantile(NCRW, 0.75, na.rm = TRUE) + f * IQR(NCRW, na.rm = TRUE)))
    }

  if (exists("SCOL", fb))
    if (dim(subset(fb, !(SCOL %in% c(1:9, NA))))[1] > 0) {
      cat("\n", "- Out of range values for storage root skin color (SCOL):", "\n")
      print(subset(fb, !(SCOL %in% c(1:9, NA))))
    }

  if (exists("FCOL", fb))
    if (dim(subset(fb, !(FCOL %in% c(1:9, NA))))[1] > 0) {
      cat("\n", "- Out of range values for storage root flesh color (FCOL):", "\n")
      print(subset(fb, !(FCOL %in% c(1:9, NA))))
    }

  if (exists("RS", fb))
    if (dim(subset(fb, !(RS %in% c(1:9, NA))))[1] > 0) {
      cat("\n", "- Out of range values for root size (RS):", "\n")
      print(subset(fb, !(RS %in% c(1:9, NA))))
    }

  if (exists("RF", fb))
    if (dim(subset(fb, !(RF %in% c(1:9, NA))))[1] > 0) {
      cat("\n", "- Out of range values for root form (RF):", "\n")
      print(subset(fb, !(RF %in% c(1:9, NA))))
    }

  if (exists("DAMR", fb))
    if (dim(subset(fb, !(DAMR %in% c(1:9, NA))))[1] > 0) {
      cat("\n", "- Out of range values for root defects (DAMR):", "\n")
      print(subset(fb, !(DAMR %in% c(1:9, NA))))
    }

  if (exists("RSPR", fb))
    if (dim(subset(fb, !(RSPR %in% c(1:9, NA))))[1] > 0) {
      cat("\n", "- Out of range values for root sprouting (RSPR):", "\n")
      print(subset(fb, !(RSPR %in% c(1:9, NA))))
    }

  if (exists("WED1", fb))
    if (dim(subset(fb, !(WED1 %in% c(1:9, NA))))[1] > 0) {
      cat("\n", "- Out of range values for weevil damage first evaluation (WED1):", "\n")
      print(subset(fb, !(WED1 %in% c(1:9, NA))))
    }

  if (exists("WED2", fb))
    if (dim(subset(fb, !(WED2 %in% c(1:9, NA))))[1] > 0) {
      cat("\n", "- Out of range values for weevil damage second evaluation (WED2):", "\n")
      print(subset(fb, !(WED2 %in% c(1:9, NA))))
    }
}

# Check consistency for sweetpotato experimental data, part 11.
# Outliers detection based on interquartile range and values out of range for DM data.

spconsis11 <- function(fb, f) {

  if (exists("DMF", fb))
    if (dim(subset(fb, DMF < 0))[1] > 0) {
      cat("\n", "- Out of range values for fresh weight of roots for dry matter assessment (DMF):", "\n")
      print(subset(fb, DMF < 0))
    }

  if (exists("DMF", fb))
    if (dim(subset(fb, DMF < quantile(DMF, 0.25, na.rm = TRUE) - f * IQR(DMF, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme low values for fresh weight of roots for dry matter assessment (DMF):", "\n")
      print(subset(fb, DMF < quantile(DMF, 0.25, na.rm = TRUE) - f * IQR(DMF, na.rm = TRUE)))
    }

  if (exists("DMF", fb))
    if (dim(subset(fb, DMF > quantile(DMF, 0.75, na.rm = TRUE) + f * IQR(DMF, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme high values for fresh weight of roots for dry matter assessment (DMF):", "\n")
      print(subset(fb, DMF > quantile(DMF, 0.75, na.rm = TRUE) + f * IQR(DMF, na.rm = TRUE)))
    }

  if (exists("DMD", fb))
    if (dim(subset(fb, DMD < 0))[1] > 0) {
      cat("\n", "- Out of range values for dry weight of roots for dry matter assessment (DMD):", "\n")
      print(subset(fb, DMD < 0))
    }

  if (exists("DMD", fb))
    if (dim(subset(fb, DMD < quantile(DMD, 0.25, na.rm = TRUE) - f * IQR(DMD, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme low values for dry weight of roots for dry matter assessment (DMD):", "\n")
      print(subset(fb, DMD < quantile(DMD, 0.25, na.rm = TRUE) - f * IQR(DMD, na.rm = TRUE)))
    }

  if (exists("DMD", fb))
    if (dim(subset(fb, DMD > quantile(DMD, 0.75, na.rm = TRUE) + f * IQR(DMD, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme high values for dry weight of roots for dry matter assessment (DMD):", "\n")
      print(subset(fb, DMD > quantile(DMD, 0.75, na.rm = TRUE) + f * IQR(DMD, na.rm = TRUE)))
    }

  if (exists("DMVF", fb))
    if (dim(subset(fb, DMVF < 0))[1] > 0) {
      cat("\n", "- Out of range values for fresh weight vines for dry matter assessment (DMVF):", "\n")
      print(subset(fb, DMVF < 0))
    }

  if (exists("DMVF", fb))
    if (dim(subset(fb, DMVF < quantile(DMVF, 0.25, na.rm = TRUE) - f * IQR(DMVF, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme low values for fresh weight of vines for dry matter assessment (DMVF):", "\n")
      print(subset(fb, DMVF < quantile(DMVF, 0.25, na.rm = TRUE) - f * IQR(DMVF, na.rm = TRUE)))
    }

  if (exists("DMVF", fb))
    if (dim(subset(fb, DMVF > quantile(DMVF, 0.75, na.rm = TRUE) + f * IQR(DMVF, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme high values for fresh weight of vines for dry matter assessment (DMVF):", "\n")
      print(subset(fb, DMVF > quantile(DMVF, 0.75, na.rm = TRUE) + f * IQR(DMVF, na.rm = TRUE)))
    }

  if (exists("DMVD", fb))
    if (dim(subset(fb, DMVD < 0))[1] > 0) {
      cat("\n", "- Out of range values for dry weight of vines for dry matter assessment (DMVD):", "\n")
      print(subset(fb, DMVD < 0))
    }

  if (exists("DMVD", fb))
    if (dim(subset(fb, DMVD < quantile(DMVD, 0.25, na.rm = TRUE) - f * IQR(DMVD, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme low values for dry weight of vines for dry matter assessment (DMVD):", "\n")
      print(subset(fb, DMVD < quantile(DMVD, 0.25, na.rm = TRUE) - f * IQR(DMVD, na.rm = TRUE)))
    }

  if (exists("DMVD", fb))
    if (dim(subset(fb, DMVD > quantile(DMVD, 0.75, na.rm = TRUE) + f * IQR(DMVD, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme high values for dry weight of vines for dry matter assessment (DMVD):", "\n")
      print(subset(fb, DMVD > quantile(DMVD, 0.75, na.rm = TRUE) + f * IQR(DMVD, na.rm = TRUE)))
    }

  if (exists("DM", fb))
    if (dim(subset(fb, DM < 0 | DM > 100))[1] > 0) {
      cat("\n", "- Out of range values for storage root dry matter content (DM):", "\n")
      print(subset(fb, DM < 0 | DM > 100))
    }

  if (exists("DM", fb))
    if (dim(subset(fb, DM < quantile(DM, 0.25, na.rm = TRUE) - f * IQR(DM, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme low values for storage root dry matter content (DM):", "\n")
      print(subset(fb, DM < quantile(DM, 0.25, na.rm = TRUE) - f * IQR(DM, na.rm = TRUE)))
    }

  if (exists("DM", fb))
    if (dim(subset(fb, DM > quantile(DM, 0.75, na.rm = TRUE) + f * IQR(DM, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme high values for storage root dry matter content (DM):", "\n")
      print(subset(fb, DM > quantile(DM, 0.75, na.rm = TRUE) + f * IQR(DM, na.rm = TRUE)))
    }

  if (exists("DMV", fb))
    if (dim(subset(fb, DMV < 0 | DMV > 100))[1] > 0) {
      cat("\n", "- Out of range values for vine dry matter content (DMV):", "\n")
      print(subset(fb, DMV < 0 | DMV > 100))
    }
  
  if (exists("DMV", fb))
    if (dim(subset(fb, DMV < quantile(DMV, 0.25, na.rm = TRUE) - f * IQR(DMV, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme low values for vine dry matter content (DMV):", "\n")
      print(subset(fb, DMV < quantile(DMV, 0.25, na.rm = TRUE) - f * IQR(DMV, na.rm = TRUE)))
    }
  
  if (exists("DMV", fb))
    if (dim(subset(fb, DMV > quantile(DMV, 0.75, na.rm = TRUE) + f * IQR(DMV, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme high values for vine dry matter content (DMV):", "\n")
      print(subset(fb, DMV > quantile(DMV, 0.75, na.rm = TRUE) + f * IQR(DMV, na.rm = TRUE)))
    }

    if (exists("DMFY", fb))
    if (dim(subset(fb, DMFY < 0))[1] > 0) {
      cat("\n", "- Out of range values for dry matter foliage yield (DMFY):", "\n")
      print(subset(fb, DMFY < 0))
    }

  if (exists("DMFY", fb))
    if (dim(subset(fb, DMFY < quantile(DMFY, 0.25, na.rm = TRUE) - f * IQR(DMFY, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme low values for dry matter foliage yield (DMFY):", "\n")
      print(subset(fb, DMFY < quantile(DMFY, 0.25, na.rm = TRUE) - f * IQR(DMFY, na.rm = TRUE)))
    }

  if (exists("DMFY", fb))
    if (dim(subset(fb, DMFY > quantile(DMFY, 0.75, na.rm = TRUE) + f * IQR(DMFY, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme high values for dry matter foliage yield (DMFY):", "\n")
      print(subset(fb, DMFY > quantile(DMFY, 0.75, na.rm = TRUE) + f * IQR(DMFY, na.rm = TRUE)))
    }

  if (exists("DMRY", fb))
    if (dim(subset(fb, DMRY < 0))[1] > 0) {
      cat("\n", "- Out of range values for dry matter root yield (DMRY):", "\n")
      print(subset(fb, DMRY < 0))
    }

  if (exists("DMRY", fb))
    if (dim(subset(fb, DMRY < quantile(DMRY, 0.25, na.rm = TRUE) - f * IQR(DMRY, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme low values for dry matter root yield (DMRY):", "\n")
      print(subset(fb, DMRY < quantile(DMRY, 0.25, na.rm = TRUE) - f * IQR(DMRY, na.rm = TRUE)))
    }

  if (exists("DMRY", fb))
    if (dim(subset(fb, DMRY > quantile(DMRY, 0.75, na.rm = TRUE) + f * IQR(DMRY, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme high values for dry matter root yield (DMRY):", "\n")
      print(subset(fb, DMRY > quantile(DMRY, 0.75, na.rm = TRUE) + f * IQR(DMRY, na.rm = TRUE)))
    }
}

# Check consistency for sweetpotato experimental data, part 12.
# Outliers detection based on interquartile range and values out of range for cooked traits.

spconsis12 <- function(fb, f) {

  if (exists("FRAW1", fb))
    if (dim(subset(fb, !(FRAW1 %in% c(1:9, NA))))[1] > 0) {
      cat("\n", "- Out of range values for root fiber first determination (FRAW1):", "\n")
      print(subset(fb, !(FRAW1 %in% c(1:9, NA))))
    }

  if (exists("SURAW1", fb))
    if (dim(subset(fb, !(SURAW1 %in% c(1:9, NA))))[1] > 0) {
      cat("\n", "- Out of range values for root sugar first determination (SURAW1):", "\n")
      print(subset(fb, !(SURAW1 %in% c(1:9, NA))))
    }

  if (exists("STRAW1", fb))
    if (dim(subset(fb, !(STRAW1 %in% c(1:9, NA))))[1] > 0) {
      cat("\n", "- Out of range values for root starch first determination (STRAW1):", "\n")
      print(subset(fb, !(STRAW1 %in% c(1:9, NA))))
    }

  if (exists("COOF1", fb))
    if (dim(subset(fb, !(COOF1 %in% c(1:9, NA))))[1] > 0) {
      cat("\n", "- Out of range values for cooked fiber first evaluation (COOF1):", "\n")
      print(subset(fb, !(COOF1 %in% c(1:9, NA))))
    }

  if (exists("COOSU1", fb))
    if (dim(subset(fb, !(COOSU1 %in% c(1:9, NA))))[1] > 0) {
      cat("\n", "- Out of range values for cooked sugars first evaluation (COOSU1):", "\n")
      print(subset(fb, !(COOSU1 %in% c(1:9, NA))))
    }

  if (exists("COOST1", fb))
    if (dim(subset(fb, !(COOST1 %in% c(1:9, NA))))[1] > 0) {
      cat("\n", "- Out of range values for cooked starch first evaluation (COOST1):", "\n")
      print(subset(fb, !(COOST1 %in% c(1:9, NA))))
    }

  if (exists("COOT1", fb))
    if (dim(subset(fb, !(COOT1 %in% c(1:9, NA))))[1] > 0) {
      cat("\n", "- Out of range values for cooked taste first evaluation (COOT1):", "\n")
      print(subset(fb, !(COOT1 %in% c(1:9, NA))))
    }

  if (exists("COOAP1", fb))
    if (dim(subset(fb, !(COOAP1 %in% c(1:9, NA))))[1] > 0) {
      cat("\n", "- Out of range values for cooked appearance first evaluation (COOAP1):", "\n")
      print(subset(fb, !(COOAP1 %in% c(1:9, NA))))
    }

  if (exists("FRAW2", fb))
    if (dim(subset(fb, !(FRAW2 %in% c(1:9, NA))))[1] > 0) {
      cat("\n", "- Out of range values for root fiber second determination (FRAW2):", "\n")
      print(subset(fb, !(FRAW2 %in% c(1:9, NA))))
    }

  if (exists("SURAW2", fb))
    if (dim(subset(fb, !(SURAW2 %in% c(1:9, NA))))[1] > 0) {
      cat("\n", "- Out of range values for root sugar second determination (SURAW2):", "\n")
      print(subset(fb, !(SURAW2 %in% c(1:9, NA))))
    }

  if (exists("STRAW2", fb))
    if (dim(subset(fb, !(STRAW2 %in% c(1:9, NA))))[1] > 0) {
      cat("\n", "- Out of range values for root starch second determination (STRAW2):", "\n")
      print(subset(fb, !(STRAW2 %in% c(1:9, NA))))
    }

  if (exists("COOF2", fb))
    if (dim(subset(fb, !(COOF2 %in% c(1:9, NA))))[1] > 0) {
      cat("\n", "- Out of range values for cooked fiber second evaluation (COOF2):", "\n")
      print(subset(fb, !(COOF2 %in% c(1:9, NA))))
    }

  if (exists("COOSU2", fb))
    if (dim(subset(fb, !(COOSU2 %in% c(1:9, NA))))[1] > 0) {
      cat("\n", "- Out of range values for cooked sugars second evaluation (COOSU2):", "\n")
      print(subset(fb, !(COOSU2 %in% c(1:9, NA))))
    }

  if (exists("COOST2", fb))
    if (dim(subset(fb, !(COOST2 %in% c(1:9, NA))))[1] > 0) {
      cat("\n", "- Out of range values for cooked starch second evaluation (COOST2):", "\n")
      print(subset(fb, !(COOST2 %in% c(1:9, NA))))
    }

  if (exists("COOT2", fb))
    if (dim(subset(fb, !(COOT2 %in% c(1:9, NA))))[1] > 0) {
      cat("\n", "- Out of range values for cooked taste second evaluation (COOT2):", "\n")
      print(subset(fb, !(COOT2 %in% c(1:9, NA))))
    }

  if (exists("COOAP2", fb))
    if (dim(subset(fb, !(COOAP2 %in% c(1:9, NA))))[1] > 0) {
      cat("\n", "- Out of range values for cooked appearance second evaluation (COOAP2):", "\n")
      print(subset(fb, !(COOAP2 %in% c(1:9, NA))))
    }
}

# Check consistency for sweetpotato experimental data, part 13.
# Outliers detection based on interquartile range and values out of range for lab data.

spconsis13 <- function(fb, f) {

  bc.cc.values <- c(0.03, 0, 0.12, 0.02, 0.15, 1.38, 1.65, 1.5, 1.74, 1.76, 0.69, 1.17, 1.32,
                    1.04, 4.41, 4.92, 6.12, 5.46, 3.96, 5.49, 3.03, 3.76, 4.61, 7.23, 7.76,
                    10.5, 11.03, 12.39, 14.37)

  if (exists("PROT", fb))
    if (dim(subset(fb, PROT < 0))[1] > 0) {
      cat("\n", "- Out of range values for protein (PROT):", "\n")
      print(subset(fb, PROT < 0))
    }

  if (exists("PROT", fb))
    if (dim(subset(fb, PROT < quantile(PROT, 0.25, na.rm = TRUE) - f * IQR(PROT, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme low values for protein (PROT):", "\n")
      print(subset(fb, PROT < quantile(PROT, 0.25, na.rm = TRUE) - f * IQR(PROT, na.rm = TRUE)))
    }

  if (exists("PROT", fb))
    if (dim(subset(fb, PROT > quantile(PROT, 0.75, na.rm = TRUE) + f * IQR(PROT, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme high values for protein (PROT):", "\n")
      print(subset(fb, PROT > quantile(PROT, 0.75, na.rm = TRUE) + f * IQR(PROT, na.rm = TRUE)))
    }

  if (exists("FE", fb))
    if (dim(subset(fb, FE < 0))[1] > 0) {
      cat("\n", "- Out of range values for iron in dry weight (FE):", "\n")
      print(subset(fb, FE < 0))
    }

  if (exists("FE", fb))
    if (dim(subset(fb, FE < quantile(FE, 0.25, na.rm = TRUE) - f * IQR(FE, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme low values for iron in dry weight (FE):", "\n")
      print(subset(fb, FE < quantile(FE, 0.25, na.rm = TRUE) - f * IQR(FE, na.rm = TRUE)))
    }

  if (exists("FE", fb))
    if (dim(subset(fb, FE > quantile(FE, 0.75, na.rm = TRUE) + f * IQR(FE, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme high values for iron in dry weight (FE):", "\n")
      print(subset(fb, FE > quantile(FE, 0.75, na.rm = TRUE) + f * IQR(FE, na.rm = TRUE)))
    }

  if (exists("ZN", fb))
    if (dim(subset(fb, ZN < 0))[1] > 0) {
      cat("\n", "- Out of range values for zinc in dry weight (ZN):", "\n")
      print(subset(fb, ZN < 0))
    }

  if (exists("ZN", fb))
    if (dim(subset(fb, ZN < quantile(ZN, 0.25, na.rm = TRUE) - f * IQR(ZN, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme low values for zinc in dry weight (ZN):", "\n")
      print(subset(fb, ZN < quantile(ZN, 0.25, na.rm = TRUE) - f * IQR(ZN, na.rm = TRUE)))
    }

  if (exists("ZN", fb))
    if (dim(subset(fb, ZN > quantile(ZN, 0.75, na.rm = TRUE) + f * IQR(ZN, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme high values for zinc in dry weight (ZN):", "\n")
      print(subset(fb, ZN > quantile(ZN, 0.75, na.rm = TRUE) + f * IQR(ZN, na.rm = TRUE)))
    }

  if (exists("CA", fb))
    if (dim(subset(fb, CA < 0))[1] > 0) {
      cat("\n", "- Out of range values for calcium in dry weight (CA):", "\n")
      print(subset(fb, CA < 0))
    }

  if (exists("CA", fb))
    if (dim(subset(fb, CA < quantile(CA, 0.25, na.rm = TRUE) - f * IQR(CA, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme low values for calcium in dry weight (CA):", "\n")
      print(subset(fb, CA < quantile(CA, 0.25, na.rm = TRUE) - f * IQR(CA, na.rm = TRUE)))
    }

  if (exists("CA", fb))
    if (dim(subset(fb, CA > quantile(CA, 0.75, na.rm = TRUE) + f * IQR(CA, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme high values for calcium in dry weight (CA):", "\n")
      print(subset(fb, CA > quantile(CA, 0.75, na.rm = TRUE) + f * IQR(CA, na.rm = TRUE)))
    }

  if (exists("MG", fb))
    if (dim(subset(fb, MG < 0))[1] > 0) {
      cat("\n", "- Out of range values for magnesium in dry weight (MG):", "\n")
      print(subset(fb, MG < 0))
    }

  if (exists("MG", fb))
    if (dim(subset(fb, MG < quantile(MG, 0.25, na.rm = TRUE) - f * IQR(MG, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme low values for magnesium in dry weight (MG):", "\n")
      print(subset(fb, MG < quantile(MG, 0.25, na.rm = TRUE) - f * IQR(MG, na.rm = TRUE)))
    }

  if (exists("MG", fb))
    if (dim(subset(fb, MG > quantile(MG, 0.75, na.rm = TRUE) + f * IQR(MG, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme high values for magnesium in dry weight (MG):", "\n")
      print(subset(fb, MG > quantile(MG, 0.75, na.rm = TRUE) + f * IQR(MG, na.rm = TRUE)))
    }

  if (exists("BC", fb))
    if (dim(subset(fb, BC < 0))[1] > 0) {
      cat("\n", "- Out of range values for beta-carotene in dry weight (BC):", "\n")
      print(subset(fb, BC < 0))
    }

  if (exists("BC", fb))
    if (dim(subset(fb, BC < quantile(BC, 0.25, na.rm = TRUE) - f * IQR(BC, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme low values for beta-carotene in dry weight (BC):", "\n")
      print(subset(fb, BC < quantile(BC, 0.25, na.rm = TRUE) - f * IQR(BC, na.rm = TRUE)))
    }

  if (exists("BC", fb))
    if (dim(subset(fb, BC > quantile(BC, 0.75, na.rm = TRUE) + f * IQR(BC, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme high values for beta-carotene in dry weight (BC):", "\n")
      print(subset(fb, BC > quantile(BC, 0.75, na.rm = TRUE) + f * IQR(BC, na.rm = TRUE)))
    }

  if (exists("BC.CC", fb))
    if (dim(subset(fb, !(BC.CC %in% c(bc.cc.values, NA))))[1] > 0) {
      cat("\n", "- Out of range values for beta-carotene with color chart (BC.CC):", "\n")
      print(subset(fb, !(BC.CC %in% c(bc.cc.values, NA))))
    }

  if (exists("TC", fb))
    if (dim(subset(fb, TC < 0))[1] > 0) {
      cat("\n", "- Out of range values for total carotenoids in dry weight (TC):", "\n")
      print(subset(fb, TC < 0))
    }

  if (exists("TC", fb))
    if (dim(subset(fb, TC < quantile(TC, 0.25, na.rm = TRUE) - f * IQR(TC, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme low values for total carotenoids in dry weight (TC):", "\n")
      print(subset(fb, TC < quantile(TC, 0.25, na.rm = TRUE) - f * IQR(TC, na.rm = TRUE)))
    }

  if (exists("TC", fb))
    if (dim(subset(fb, TC > quantile(TC, 0.75, na.rm = TRUE) + f * IQR(TC, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme high values for total carotenoids in dry weight (TC):", "\n")
      print(subset(fb, TC > quantile(TC, 0.75, na.rm = TRUE) + f * IQR(TC, na.rm = TRUE)))
    }

  if (exists("STAR", fb))
    if (dim(subset(fb, STAR < 0))[1] > 0) {
      cat("\n", "- Out of range values for starch (STAR):", "\n")
      print(subset(fb, STAR < 0))
    }

  if (exists("STAR", fb))
    if (dim(subset(fb, STAR < quantile(STAR, 0.25, na.rm = TRUE) - f * IQR(STAR, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme low values for starch (STAR):", "\n")
      print(subset(fb, STAR < quantile(STAR, 0.25, na.rm = TRUE) - f * IQR(STAR, na.rm = TRUE)))
    }

  if (exists("STAR", fb))
    if (dim(subset(fb, STAR > quantile(STAR, 0.75, na.rm = TRUE) + f * IQR(STAR, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme high values for starch (STAR):", "\n")
      print(subset(fb, STAR > quantile(STAR, 0.75, na.rm = TRUE) + f * IQR(STAR, na.rm = TRUE)))
    }

  if (exists("FRUC", fb))
    if (dim(subset(fb, FRUC < 0))[1] > 0) {
      cat("\n", "- Out of range values for fructose (FRUC):", "\n")
      print(subset(fb, FRUC < 0))
    }

  if (exists("FRUC", fb))
    if (dim(subset(fb, FRUC < quantile(FRUC, 0.25, na.rm = TRUE) - f * IQR(FRUC, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme low values for fructose (FRUC):", "\n")
      print(subset(fb, FRUC < quantile(FRUC, 0.25, na.rm = TRUE) - f * IQR(FRUC, na.rm = TRUE)))
    }

  if (exists("FRUC", fb))
    if (dim(subset(fb, FRUC > quantile(FRUC, 0.75, na.rm = TRUE) + f * IQR(FRUC, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme high values for fructose (FRUC):", "\n")
      print(subset(fb, FRUC > quantile(FRUC, 0.75, na.rm = TRUE) + f * IQR(FRUC, na.rm = TRUE)))
    }

  if (exists("GLUC", fb))
    if (dim(subset(fb, GLUC < 0))[1] > 0) {
      cat("\n", "- Out of range values for glucose (GLUC):", "\n")
      print(subset(fb, GLUC < 0))
    }

  if (exists("GLUC", fb))
    if (dim(subset(fb, GLUC < quantile(GLUC, 0.25, na.rm = TRUE) - f * IQR(GLUC, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme low values for glucose (GLUC):", "\n")
      print(subset(fb, GLUC < quantile(GLUC, 0.25, na.rm = TRUE) - f * IQR(GLUC, na.rm = TRUE)))
    }

  if (exists("GLUC", fb))
    if (dim(subset(fb, GLUC > quantile(GLUC, 0.75, na.rm = TRUE) + f * IQR(GLUC, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme high values for glucose (GLUC):", "\n")
      print(subset(fb, GLUC > quantile(GLUC, 0.75, na.rm = TRUE) + f * IQR(GLUC, na.rm = TRUE)))
    }

  if (exists("SUCR", fb))
    if (dim(subset(fb, SUCR < 0))[1] > 0) {
      cat("\n", "- Out of range values for sucrose (SUCR):", "\n")
      print(subset(fb, SUCR < 0))
    }

  if (exists("SUCR", fb))
    if (dim(subset(fb, SUCR < quantile(SUCR, 0.25, na.rm = TRUE) - f * IQR(SUCR, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme low values for sucrose (SUCR):", "\n")
      print(subset(fb, SUCR < quantile(SUCR, 0.25, na.rm = TRUE) - f * IQR(SUCR, na.rm = TRUE)))
    }

  if (exists("SUCR", fb))
    if (dim(subset(fb, SUCR > quantile(SUCR, 0.75, na.rm = TRUE) + f * IQR(SUCR, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme high values for sucrose (SUCR):", "\n")
      print(subset(fb, SUCR > quantile(SUCR, 0.75, na.rm = TRUE) + f * IQR(SUCR, na.rm = TRUE)))
    }

  if (exists("MALT", fb))
    if (dim(subset(fb, MALT < 0))[1] > 0) {
      cat("\n", "- Out of range values for maltose (MALT):", "\n")
      print(subset(fb, MALT < 0))
    }

  if (exists("MALT", fb))
    if (dim(subset(fb, MALT < quantile(MALT, 0.25, na.rm = TRUE) - f * IQR(MALT, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme low values for maltose (MALT):", "\n")
      print(subset(fb, MALT < quantile(MALT, 0.25, na.rm = TRUE) - f * IQR(MALT, na.rm = TRUE)))
    }

  if (exists("MALT", fb))
    if (dim(subset(fb, MALT > quantile(MALT, 0.75, na.rm = TRUE) + f * IQR(MALT, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme high values for maltose (MALT):", "\n")
      print(subset(fb, MALT > quantile(MALT, 0.75, na.rm = TRUE) + f * IQR(MALT, na.rm = TRUE)))
    }
}

# Check consistency for sweetpotato experimental data, part 14.
# Outliers detection based on interquartile range and values out of range for
# derived variables.

spconsis14 <- function(fb, f) {

  if (exists("TRW", fb))
    if (dim(subset(fb, TRW < 0))[1] > 0) {
      cat("\n", "- Out of range values for total root weight (TRW):", "\n")
      print(subset(fb, TRW < 0))
    }

  if (exists("TRW", fb))
    if (dim(subset(fb, TRW < quantile(TRW, 0.25, na.rm = TRUE) - f * IQR(TRW, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme low values for total root weight (TRW):", "\n")
      print(subset(fb, TRW < quantile(TRW, 0.25, na.rm = TRUE) - f * IQR(TRW, na.rm = TRUE)))
    }

  if (exists("TRW", fb))
    if (dim(subset(fb, TRW > quantile(TRW, 0.75, na.rm = TRUE) + f * IQR(TRW, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme high values for total root weight (TRW):", "\n")
      print(subset(fb, TRW > quantile(TRW, 0.75, na.rm = TRUE) + f * IQR(TRW, na.rm = TRUE)))
    }

  if (exists("CYTHA", fb))
    if (dim(subset(fb, CYTHA < 0))[1] > 0) {
      cat("\n", "- Out of range values for commercial root yield in tons per hectare (CYTHA):", "\n")
      print(subset(fb, CYTHA < 0))
    }

  if (exists("CYTHA", fb))
    if (dim(subset(fb, CYTHA < quantile(CYTHA, 0.25, na.rm = TRUE) - f * IQR(CYTHA, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme low values for commercial root yield in tons per hectare (CYTHA):", "\n")
      print(subset(fb, CYTHA < quantile(CYTHA, 0.25, na.rm = TRUE) - f * IQR(CYTHA, na.rm = TRUE)))
    }

  if (exists("CYTHA", fb))
    if (dim(subset(fb, CYTHA > quantile(CYTHA, 0.75, na.rm = TRUE) + f * IQR(CYTHA, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme high values for commercial root yield in tons per hectare (CYTHA):", "\n")
      print(subset(fb, CYTHA > quantile(CYTHA, 0.75, na.rm = TRUE) + f * IQR(CYTHA, na.rm = TRUE)))
    }

  if (exists("RYTHA", fb))
    if (dim(subset(fb, RYTHA < 0))[1] > 0) {
      cat("\n", "- Out of range values for total root yield in tons per hectare (RYTHA):", "\n")
      print(subset(fb, RYTHA < 0))
    }

  if (exists("RYTHA", fb))
    if (dim(subset(fb, RYTHA < quantile(RYTHA, 0.25, na.rm = TRUE) - f * IQR(RYTHA, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme low values for total root yield in tons per hectare (RYTHA):", "\n")
      print(subset(fb, RYTHA < quantile(RYTHA, 0.25, na.rm = TRUE) - f * IQR(RYTHA, na.rm = TRUE)))
    }

  if (exists("RYTHA", fb))
    if (dim(subset(fb, RYTHA > quantile(RYTHA, 0.75, na.rm = TRUE) + f * IQR(RYTHA, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme high values for total root yield in tons per hectare (RYTHA):", "\n")
      print(subset(fb, RYTHA > quantile(RYTHA, 0.75, na.rm = TRUE) + f * IQR(RYTHA, na.rm = TRUE)))
    }

  if (exists("ACRW", fb))
    if (dim(subset(fb, ACRW < 0))[1] > 0) {
      cat("\n", "- Out of range values for average commercial root weight (ACRW):", "\n")
      print(subset(fb, ACRW < 0))
    }

  if (exists("ACRW", fb))
    if (dim(subset(fb, ACRW < quantile(ACRW, 0.25, na.rm = TRUE) - f * IQR(ACRW, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme low values for average commercial root weight (ACRW):", "\n")
      print(subset(fb, ACRW < quantile(ACRW, 0.25, na.rm = TRUE) - f * IQR(ACRW, na.rm = TRUE)))
    }

  if (exists("ACRW", fb))
    if (dim(subset(fb, ACRW > quantile(ACRW, 0.75, na.rm = TRUE) + f * IQR(ACRW, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme high values for average commercial root weight (ACRW):", "\n")
      print(subset(fb, ACRW > quantile(ACRW, 0.75, na.rm = TRUE) + f * IQR(ACRW, na.rm = TRUE)))
    }

  if (exists("NRPP", fb))
    if (dim(subset(fb, NRPP < 0))[1] > 0) {
      cat("\n", "- Out of range values for number of roots per plant (NRPP):", "\n")
      print(subset(fb, NRPP < 0))
    }

  if (exists("NRPP", fb))
    if (dim(subset(fb, NRPP < quantile(NRPP, 0.25, na.rm = TRUE) - f * IQR(NRPP, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme low values for number of roots per plant (NRPP):", "\n")
      print(subset(fb, NRPP < quantile(NRPP, 0.25, na.rm = TRUE) - f * IQR(NRPP, na.rm = TRUE)))
    }

  if (exists("NRPP", fb))
    if (dim(subset(fb, NRPP > quantile(NRPP, 0.75, na.rm = TRUE) + f * IQR(NRPP, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme high values for number of roots per plant (NRPP):", "\n")
      print(subset(fb, NRPP > quantile(NRPP, 0.75, na.rm = TRUE) + f * IQR(NRPP, na.rm = TRUE)))
    }

  if (exists("YPP", fb))
    if (dim(subset(fb, YPP < 0))[1] > 0) {
      cat("\n", "- Out of range values for yield per plant (YPP):", "\n")
      print(subset(fb, YPP < 0))
    }

  if (exists("YPP", fb))
    if (dim(subset(fb, YPP < quantile(YPP, 0.25, na.rm = TRUE) - f * IQR(YPP, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme low values for yield per plant (YPP):", "\n")
      print(subset(fb, YPP < quantile(YPP, 0.25, na.rm = TRUE) - f * IQR(YPP, na.rm = TRUE)))
    }

  if (exists("YPP", fb))
    if (dim(subset(fb, YPP > quantile(YPP, 0.75, na.rm = TRUE) + f * IQR(YPP, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme high values for yield per plant (YPP):", "\n")
      print(subset(fb, YPP > quantile(YPP, 0.75, na.rm = TRUE) + f * IQR(YPP, na.rm = TRUE)))
    }

  if (exists("CI", fb))
    if (dim(subset(fb, CI < 0 | CI > 100))[1] > 0) {
      cat("\n", "- Out of range values for commercial index (CI):", "\n")
      print(subset(fb, CI < 0 | CI > 100))
    }

  if (exists("CI", fb))
    if (dim(subset(fb, CI < quantile(CI, 0.25, na.rm = TRUE) - f * IQR(CI, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme low values for commercial index (CI):", "\n")
      print(subset(fb, CI < quantile(CI, 0.25, na.rm = TRUE) - f * IQR(CI, na.rm = TRUE)))
    }

  if (exists("CI", fb))
    if (dim(subset(fb, CI > quantile(CI, 0.75, na.rm = TRUE) + f * IQR(CI, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme high values for commercial index (CI):", "\n")
      print(subset(fb, CI > quantile(CI, 0.75, na.rm = TRUE) + f * IQR(CI, na.rm = TRUE)))
    }

  if (exists("HI", fb))
    if (dim(subset(fb, HI < 0 | HI > 100))[1] > 0) {
      cat("\n", "- Out of range values for harvest index (HI):", "\n")
      print(subset(fb, HI < 0 | HI > 100))
    }

  if (exists("HI", fb))
    if (dim(subset(fb, HI < quantile(HI, 0.25, na.rm = TRUE) - f * IQR(HI, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme low values for harvest index (HI):", "\n")
      print(subset(fb, HI < quantile(HI, 0.25, na.rm = TRUE) - f * IQR(HI, na.rm = TRUE)))
    }

  if (exists("HI", fb))
    if (dim(subset(fb, HI > quantile(HI, 0.75, na.rm = TRUE) + f * IQR(HI, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme high values for harvest index (HI):", "\n")
      print(subset(fb, HI > quantile(HI, 0.75, na.rm = TRUE) + f * IQR(HI, na.rm = TRUE)))
    }

  if (exists("SHI", fb))
    if (dim(subset(fb, SHI < 0 | SHI > 100))[1] > 0) {
      cat("\n", "- Out of range values for harvest sowing index (SHI):", "\n")
      print(subset(fb, SHI < 0 | SHI > 100))
    }

  if (exists("SHI", fb))
    if (dim(subset(fb, SHI < quantile(SHI, 0.25, na.rm = TRUE) - f * IQR(SHI, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme low values for harvest sowing index (SHI):", "\n")
      print(subset(fb, SHI < quantile(SHI, 0.25, na.rm = TRUE) - f * IQR(SHI, na.rm = TRUE)))
    }

  if (exists("SHI", fb))
    if (dim(subset(fb, SHI > quantile(SHI, 0.75, na.rm = TRUE) + f * IQR(SHI, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme high values for harvest sowing index (SHI):", "\n")
      print(subset(fb, SHI > quantile(SHI, 0.75, na.rm = TRUE) + f * IQR(SHI, na.rm = TRUE)))
    }

  if (exists("BIOM", fb))
    if (dim(subset(fb, BIOM < 0))[1] > 0) {
      cat("\n", "- Out of range values for biomass yield (BIOM):", "\n")
      print(subset(fb, BIOM < 0))
    }

  if (exists("BIOM", fb))
    if (dim(subset(fb, BIOM < quantile(BIOM, 0.25, na.rm = TRUE) - f * IQR(BIOM, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme low values for biomass yield (BIOM):", "\n")
      print(subset(fb, BIOM < quantile(BIOM, 0.25, na.rm = TRUE) - f * IQR(BIOM, na.rm = TRUE)))
    }

  if (exists("BIOM", fb))
    if (dim(subset(fb, BIOM > quantile(BIOM, 0.75, na.rm = TRUE) + f * IQR(BIOM, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme high values for biomass yield (BIOM):", "\n")
      print(subset(fb, BIOM > quantile(BIOM, 0.75, na.rm = TRUE) + f * IQR(BIOM, na.rm = TRUE)))
    }

  if (exists("FYTHA", fb))
    if (dim(subset(fb, FYTHA < 0))[1] > 0) {
      cat("\n", "- Out of range values for foliage total yield in tons per hectare (FYTHA):", "\n")
      print(subset(fb, FYTHA < 0))
    }

  if (exists("FYTHA", fb))
    if (dim(subset(fb, FYTHA < quantile(FYTHA, 0.25, na.rm = TRUE) - f * IQR(FYTHA, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme low values for foliage total yield in tons per hectare (FYTHA):", "\n")
      print(subset(fb, FYTHA < quantile(FYTHA, 0.25, na.rm = TRUE) - f * IQR(FYTHA, na.rm = TRUE)))
    }

  if (exists("FYTHA", fb))
    if (dim(subset(fb, FYTHA > quantile(FYTHA, 0.75, na.rm = TRUE) + f * IQR(FYTHA, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme high values for foliage total yield in tons per hectare (FYTHA):", "\n")
      print(subset(fb, FYTHA > quantile(FYTHA, 0.75, na.rm = TRUE) + f * IQR(FYTHA, na.rm = TRUE)))
    }

  if (exists("RFR", fb))
    if (dim(subset(fb, RFR < 0 | RFR > 100))[1] > 0) {
      cat("\n", "- Out of range values for root foliage ratio (RFR):", "\n")
      print(subset(fb, RFR < 0 | RFR > 100))
    }

  if (exists("RFR", fb))
    if (dim(subset(fb, RFR < quantile(RFR, 0.25, na.rm = TRUE) - f * IQR(RFR, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme low values for root foliage ratio (RFR):", "\n")
      print(subset(fb, RFR < quantile(RFR, 0.25, na.rm = TRUE) - f * IQR(RFR, na.rm = TRUE)))
    }

  if (exists("RFR", fb))
    if (dim(subset(fb, RFR > quantile(RFR, 0.75, na.rm = TRUE) + f * IQR(RFR, na.rm = TRUE)))[1] > 0) {
      cat("\n", "- Extreme high values for root foliage ratio (RFR):", "\n")
      print(subset(fb, RFR > quantile(RFR, 0.75, na.rm = TRUE) + f * IQR(RFR, na.rm = TRUE)))
    }
}
