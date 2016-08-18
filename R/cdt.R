#' Compute derived traits
#'
#' Compute derived traits for a given fieldbook.
#' @param fb The name of the fieldbook data frame.
#' @param method Method to scale data from plot to hectare level. Options are plot size
#' \code{"ps"} and number of plants for a full hectare \code{"np"}. See details.
#' @param value Value for the method selected in square meters if \code{method = "ps"}
#' and in number of plants per hectare if \code{method = "np"}.
#' @param nops Number of plants sowed per plot.
#' @details The data frame must use the labels (lower or upper case) listed in function
#' \code{checknames}. See \code{?checknames} for details.
#' 
#' Conversion from kilograms per plot to tons per hectare can be done using \code{ps}, the
#' plot size, or \code{np}, the total number of plants that are expected to be allocated in a
#' full hectare. In both cases computations can be adjusted by the number of harvested plants,
#' if available in the fieldbook. For \code{method = "np"}, \code{nops} must be specified to
#' compute non adjusted values.
#' @return It returns a data frame with the original and derived traits.
#' @author Raul Eyzaguirre.
#' @examples
#' cdt(pjpz09)
#' @export

cdt <- function(fb, method = c("none", "ps", "np"), value = NULL, nops = NULL) {
  
  # Check names
  
  fb <- checknames(fb)

  # Warnings
  
  method = match.arg(method)

  if (method == "ps" & is.null(value))
    warning("Plot size value is missing.", call. = FALSE)
  
  if (method == "np" & is.null(value))
    warning("Total number of plants per hectare value is missing.", call. = FALSE)

  if (!is.null(nops)) {
    if (exists("NOPS", fb)) {
      fb$NOPS <- nops
      warning("NOPS has been replace in the fieldbook by nops = ", nops, call. = FALSE)
    } else {
      fb$NOPS <- nops
    }
  }
  
  if (method == "np" & !exists("NOPS", fb))
    warning("Number of plants sowed, nops, is missing.", call. = FALSE)

  # List to overwrite

  ow <- NULL 
  
  # General computations a priori
  
  if (exists("CRW", fb) & exists("NCRW", fb)) {
    if (exists("TRW", fb))
      ow <- c(ow, "TRW")
    fb$TRW <- suma(fb$CRW, fb$NCRW)
  }

  if (exists("CRW", fb) & exists("NOCR", fb)) {
    if (exists("ACRW", fb))
      ow <- c(ow, "ACRW")
    fb$ACRW <- fb$CRW / fb$NOCR
    fb$ACRW[fb$NOCR == 0] <- NA
  }
  
  if (exists("NOCR", fb) & exists("NONC", fb) & exists("NOPH", fb)) {
    if (exists("NRPP", fb))
      ow <- c(ow, "NRPP")
    fb$NRPP <- suma(fb$NOCR, fb$NONC) / fb$NOPH
    fb$NRPP[fb$NOPH == 0] <- NA
  }

  if (exists("TRW", fb) & exists("NOPH", fb)) {
    if (exists("YPP", fb))
      ow <- c(ow, "YPP")
    fb$YPP <- fb$TRW / fb$NOPH
    fb$YPP[fb$NOPH == 0] <- NA
  }

  if (exists("NOCR", fb) & exists("NONC", fb)) {
    if (exists("CI", fb))
      ow <- c(ow, "CI")
    temp <- suma(fb$NOCR, fb$NONC)
    fb$CI <- fb$NOCR / temp * 100
    fb$CI[temp == 0] <- NA
  }

  if (exists("TRW", fb) & exists("VW", fb)) {
    if (exists("HI", fb))
      ow <- c(ow, "HI")
    temp <- suma(fb$VW, fb$TRW)
    fb$HI <- fb$TRW / temp * 100
    fb$HI[temp == 0] <- NA
  }
  
  if (exists("NOPH", fb) & exists("NOPS", fb)) {
    if (exists("SHI", fb))
      ow <- c(ow, "SHI")
    fb$SHI <- fb$NOPH / fb$NOPS * 100
    fb$SHI[fb$NOPS == 0] <- NA
  }
  
  if (exists("DMD", fb) & exists("DMF", fb)) {
    if (exists("DM", fb))
      ow <- c(ow, "DM")
    fb$DM <- fb$DMD / fb$DMF * 100
    fb$DM[fb$DMF == 0] <- NA
  }
  
  if (exists("DMVD", fb) & exists("DMVF", fb)) {
    if (exists("DMV", fb))
      ow <- c(ow, "DMV")
    fb$DMV <- fb$DMVD / fb$DMVF * 100
    fb$DMV[fb$DMVF == 0] <- NA
  }
  
  # Computations based on plot size
  
  if (method == "ps" & !is.null(value)) {

    if (exists("CRW", fb)) {
      if (exists("CYTHA", fb))
        ow <- c(ow, "CYTHA")
      fb$CYTHA <- fb$CRW * 10 / value
      if (exists("NOPH", fb) & exists("NOPS", fb)) {
        if (exists("CYTHA.AJ", fb))
          ow <- c(ow, "CYTHA.AJ")
        fb$CYTHA.AJ <- fb$CRW / fb$NOPH * fb$NOPS * 10 / value
        fb$CYTHA.AJ[fb$NOPH == 0] <- NA
      }
    }

    if (exists("TRW", fb)) {
      if (exists("RYTHA", fb))
        ow <- c(ow, "RYTHA")
      fb$RYTHA <- fb$TRW * 10 / value
      if (exists("NOPH", fb) & exists("NOPS", fb)) {
        if (exists("RYTHA.AJ", fb))
          ow <- c(ow, "RYTHA.AJ")
        fb$RYTHA.AJ <- fb$TRW / fb$NOPH * fb$NOPS * 10 / value
        fb$RYTHA.AJ[fb$NOPH == 0] <- NA
      }
    }
    
    if (exists("VW", fb)) {
      if (exists("FYTHA", fb))
        ow <- c(ow, "FYTHA")
      fb$FYTHA <- fb$VW * 10 / value
      if (exists("NOPH", fb) & exists("NOPS", fb)) {
        if (exists("FYTHA.AJ", fb))
          ow <- c(ow, "FYTHA.AJ")
        fb$FYTHA.AJ <- fb$VW / fb$NOPH * fb$NOPS * 10 / value
        fb$FYTHA.AJ[fb$NOPH == 0] <- NA
      }
    }

    if (exists("TRW", fb) & exists("DM", fb)) {
      temp1 <- fb$TRW * fb$DM / 100
      if (exists("DMRY", fb))
        ow <- c(ow, "DMRY")
      fb$DMRY <- temp1 * 10 / value
      if (exists("NOPH", fb) & exists("NOPS", fb)) {
        if (exists("DMRY.AJ", fb))
          ow <- c(ow, "DMRY.AJ")
        fb$DMRY.AJ <- temp1 / fb$NOPH * fb$NOPS * 10 / value
        fb$DMRY.AJ[fb$NOPH == 0] <- NA
      }
    }
   
    if (exists("VW", fb) & exists("DMV", fb)) {
      temp2 <- fb$VW * fb$DMV / 100
      if (exists("DMFY", fb))
        ow <- c(ow, "DMFY")
      fb$DMFY <- temp2 * 10 / value
      if (exists("NOPH", fb) & exists("NOPS", fb)) {
        if (exists("DMFY.AJ", fb))
          ow <- c(ow, "DMFY.AJ")
        fb$DMFY.AJ <- temp2 / fb$NOPH * fb$NOPS * 10 / value
        fb$DMFY.AJ[fb$NOPH == 0] <- NA
      }
    }

  }
  
  # Computations based on number of plants
  
  if (method == "np" & !is.null(value)) {

    if (exists("CRW", fb)) {
      if (exists("NOPS", fb)) {
        if (exists("CYTHA", fb))
          ow <- c(ow, "CYTHA")
        fb$CYTHA <- fb$CRW / fb$NOPS * value / 1000
        fb$CYTHA[fb$NOPS == 0] <- NA
      }
      if (exists("NOPH", fb)) {
        if (exists("CYTHA.AJ", fb))
          ow <- c(ow, "CYTHA.AJ")
        fb$CYTHA.AJ <- fb$CRW / fb$NOPH * value / 1000
        fb$CYTHA.AJ[fb$NOPH == 0] <- NA
      }
    }

    if (exists("TRW", fb)) {
      if (exists("NOPS", fb)) {
        if (exists("RYTHA", fb))
          ow <- c(ow, "RYTHA")
        fb$RYTHA <- fb$TRW / fb$NOPS * value / 1000
        fb$RYTHA[fb$NOPS == 0] <- NA
      }
      if (exists("NOPH", fb)) {
        if (exists("RYTHA.AJ", fb))
          ow <- c(ow, "RYTHA.AJ")
        fb$RYTHA.AJ <- fb$TRW / fb$NOPH * value / 1000
        fb$RYTHA.AJ[fb$NOPH == 0] <- NA
      }
    }

    if (exists("VW", fb)) {
      if (exists("NOPS", fb)) {
        if (exists("FYTHA", fb))
          ow <- c(ow, "FYTHA")
        fb$FYTHA <- fb$VW / fb$NOPS * value / 1000
        fb$FYTHA[fb$NOPS == 0] <- NA
      }
      if (exists("NOPH", fb)) {
        if (exists("FYTHA.AJ", fb))
          ow <- c(ow, "FYTHA.AJ")
        fb$FYTHA.AJ <- fb$VW / fb$NOPH * value / 1000
        fb$FYTHA.AJ[fb$NOPH == 0] <- NA
      }
    }

    if (exists("TRW", fb) & exists("DM", fb)) {
      temp1 <- fb$TRW * fb$DM / 100
      if (exists("NOPS", fb)) {
        if (exists("DMRY", fb))
          ow <- c(ow, "DMRY")
        fb$DMRY <- temp1 / fb$NOPS * value / 1000
        fb$DMRY[fb$NOPS == 0] <- NA
      }
      if (exists("NOPH", fb)) {
        if (exists("DMRY.AJ", fb))
          ow <- c(ow, "DMRY.AJ")
        fb$DMRY.AJ <- temp1 / fb$NOPH * value / 1000
        fb$DMRY.AJ[fb$NOPH == 0] <- NA
      }
    }
    
    if (exists("VW", fb) & exists("DMV", fb)) {
      temp2 <- fb$VW * fb$DMV / 100
      if (exists("NOPS", fb)) {
        if (exists("DMFY", fb))
          ow <- c(ow, "DMFY")
        fb$DMFY <- temp2 / fb$NOPS * value / 1000
        fb$DMFY[fb$NOPS == 0] <- NA
      }
      if (exists("NOPH", fb)) {
        if (exists("DMFY.AJ", fb))
          ow <- c(ow, "DMFY.AJ")
        fb$DMFY.AJ <- temp2 / fb$NOPH * value / 1000
        fb$DMFY.AJ[fb$NOPH == 0] <- NA
      }
    }
    
  }
  
  # General computations a posteriori
  
  if (exists("RYTHA", fb) & exists("FYTHA", fb)) {
    if (exists("BIOM", fb))
      ow <- c(ow, "BIOM")
    fb$BIOM <- suma(fb$RYTHA, fb$FYTHA)
  }
  
  if (exists("RYTHA.AJ", fb) & exists("FYTHA.AJ", fb)) {
    if (exists("BIOM.AJ", fb))
      ow <- c(ow, "BIOM.AJ")
    fb$BIOM.AJ <- suma(fb$RYTHA.AJ, fb$FYTHA.AJ)
  }
  
  if (exists("temp1") & exists("temp2")) {
    if (exists("RFR", fb))
      ow <- c(ow, "RFR")
      fb$RFR <- temp1 / temp2
      fb$RFR[temp2 == 0] <- NA
    }
  

  # Warning: Overwritten traits
  
  if (length(ow) > 0)
    warning("Some traits have been overwritten: ", list(ow), call. = FALSE)

  # Return
  
  fb

}
