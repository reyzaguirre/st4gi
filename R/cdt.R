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
    if (exists("nops", fb)) {
      fb$nops <- nops
      warning("nops has been replace in the fieldbook by ", nops, call. = FALSE)
    } else {
      fb$nops <- nops
    }
  }
  
  if (method == "np" & !exists("nops", fb))
    warning("Number of plants sowed, nops, is missing.", call. = FALSE)

  # List to overwrite

  ow <- NULL 
  
  # General computations a priori
  
  if (exists("crw", fb) & exists("ncrw", fb)) {
    if (exists("trw", fb))
      ow <- c(ow, "trw")
    fb$trw <- suma(fb$crw, fb$ncrw)
  }

  if (exists("crw", fb) & exists("nocr", fb)) {
    if (exists("acrw", fb))
      ow <- c(ow, "acrw")
    fb$acrw <- fb$crw / fb$nocr
    fb$acrw[fb$nocr == 0] <- NA
  }
  
  if (exists("nocr", fb) & exists("nonc", fb)) {
    if (exists("tnr", fb))
      ow <- c(ow, "tnr")
    fb$tnr <- suma(fb$nocr, fb$nonc)
  }

  if (exists("tnr", fb) & exists("noph", fb)) {
    if (exists("nrpp", fb))
      ow <- c(ow, "nrpp")
    fb$nrpp <- fb$tnr / fb$noph
    fb$nrpp[fb$noph == 0] <- NA
  }

  if (exists("tnr", fb) & exists("nops", fb)) {
    if (exists("nrpsp", fb))
      ow <- c(ow, "nrpsp")
    fb$nrpsp <- fb$tnr / fb$nops
  }

  if (exists("nocr", fb) & exists("noph", fb)) {
    if (exists("ncrpp", fb))
      ow <- c(ow, "ncrpp")
    fb$ncrpp <- fb$nocr / fb$noph
    fb$ncrpp[fb$noph == 0] <- NA
  }

  if (exists("nocr", fb) & exists("nops", fb)) {
    if (exists("ncrpsp", fb))
      ow <- c(ow, "ncrpsp")
    fb$ncrpsp <- fb$nocr / fb$nops
  }
  
  if (exists("trw", fb) & exists("noph", fb)) {
    if (exists("ypp", fb))
      ow <- c(ow, "ypp")
    fb$ypp <- fb$trw / fb$noph
    fb$ypp[fb$noph == 0] <- NA
  }

  if (exists("trw", fb) & exists("nops", fb)) {
    if (exists("ypsp", fb))
      ow <- c(ow, "ypsp")
    fb$ypsp <- fb$trw / fb$nops
  }

  if (exists("vw", fb) & exists("noph", fb)) {
    if (exists("vpp", fb))
      ow <- c(ow, "vpp")
    fb$vpp <- fb$vw / fb$noph
    fb$vpp[fb$noph == 0] <- NA
  }

  if (exists("vw", fb) & exists("nops", fb)) {
    if (exists("vpsp", fb))
      ow <- c(ow, "vpsp")
    fb$vpsp <- fb$vw / fb$nops
  }
  
  if (exists("nocr", fb) & exists("nonc", fb)) {
    if (exists("ci", fb))
      ow <- c(ow, "ci")
    fb$ci <- fb$nocr / fb$tnr * 100
    fb$ci[fb$tnr == 0] <- NA
  }

  if (exists("trw", fb) & exists("vw", fb)) {
    if (exists("hi", fb))
      ow <- c(ow, "hi")
    temp <- suma(fb$vw, fb$trw)
    fb$hi <- fb$trw / temp * 100
    fb$hi[temp == 0] <- NA
  }
  
  if (exists("noph", fb) & exists("nops", fb)) {
    if (exists("shi", fb))
      ow <- c(ow, "shi")
    fb$shi <- fb$noph / fb$nops * 100
    fb$shi[fb$nops == 0] <- NA
  }
  
  if (exists("dmd", fb) & exists("dmf", fb)) {
    if (exists("dm", fb))
      ow <- c(ow, "dm")
    fb$dm <- fb$dmd / fb$dmf * 100
    fb$dm[fb$dmf == 0] <- NA
  }
  
  if (exists("dmvd", fb) & exists("dmvf", fb)) {
    if (exists("dmv", fb))
      ow <- c(ow, "dmv")
    fb$dmv <- fb$dmvd / fb$dmvf * 100
    fb$dmv[fb$dmvf == 0] <- NA
  }
  
  if (exists("trw", fb) & exists("dm", fb)) {
    if (exists("trw.d", fb))
      ow <- c(ow, "trw.d")
    fb$trw.d <- fb$trw * fb$dm / 100
  }

  if (exists("vw", fb) & exists("dmv", fb)) {
    if (exists("vw.d", fb))
      ow <- c(ow, "vw.d")
    fb$vw.d <- fb$vw * fb$dmv / 100
  }

  # Computations based on plot size
  
  if (method == "ps" & !is.null(value)) {

    if (exists("crw", fb)) {
      if (exists("cytha", fb))
        ow <- c(ow, "cytha")
      fb$cytha <- fb$crw * 10 / value
      if (exists("noph", fb) & exists("nops", fb)) {
        if (exists("cytha.aj", fb))
          ow <- c(ow, "cytha.aj")
        fb$cytha.aj <- fb$crw / fb$noph * fb$nops * 10 / value
        fb$cytha.aj[fb$noph == 0] <- NA
      }
    }

    if (exists("trw", fb)) {
      if (exists("rytha", fb))
        ow <- c(ow, "rytha")
      fb$rytha <- fb$trw * 10 / value
      if (exists("noph", fb) & exists("nops", fb)) {
        if (exists("rytha.aj", fb))
          ow <- c(ow, "rytha.aj")
        fb$rytha.aj <- fb$trw / fb$noph * fb$nops * 10 / value
        fb$rytha.aj[fb$noph == 0] <- NA
      }
    }
    
    if (exists("vw", fb)) {
      if (exists("fytha", fb))
        ow <- c(ow, "fytha")
      fb$fytha <- fb$vw * 10 / value
      if (exists("noph", fb) & exists("nops", fb)) {
        if (exists("fytha.aj", fb))
          ow <- c(ow, "fytha.aj")
        fb$fytha.aj <- fb$vw / fb$noph * fb$nops * 10 / value
        fb$fytha.aj[fb$noph == 0] <- NA
      }
    }

    if (exists("trw.d", fb)) {
      if (exists("dmry", fb))
        ow <- c(ow, "dmry")
      fb$dmry <- fb$trw.d * 10 / value
      if (exists("noph", fb) & exists("nops", fb)) {
        if (exists("dmry.aj", fb))
          ow <- c(ow, "dmry.aj")
        fb$dmry.aj <- fb$trw.d / fb$noph * fb$nops * 10 / value
        fb$dmry.aj[fb$noph == 0] <- NA
      }
    }
   
    if (exists("vw.d", fb)) {
      if (exists("dmvy", fb))
        ow <- c(ow, "dmvy")
      fb$dmvy <- fb$vw.d * 10 / value
      if (exists("noph", fb) & exists("nops", fb)) {
        if (exists("dmvy.aj", fb))
          ow <- c(ow, "dmvy.aj")
        fb$dmvy.aj <- fb$vw.d / fb$noph * fb$nops * 10 / value
        fb$dmvy.aj[fb$noph == 0] <- NA
      }
    }
  
  }
  
  # Computations based on number of plants
  
  if (method == "np" & !is.null(value)) {

    if (exists("crw", fb)) {
      if (exists("nops", fb)) {
        if (exists("cytha", fb))
          ow <- c(ow, "cytha")
        fb$cytha <- fb$crw / fb$nops * value / 1000
        fb$cytha[fb$nops == 0] <- NA
      }
      if (exists("noph", fb)) {
        if (exists("cytha.aj", fb))
          ow <- c(ow, "cytha.aj")
        fb$cytha.aj <- fb$crw / fb$noph * value / 1000
        fb$cytha.aj[fb$noph == 0] <- NA
      }
    }

    if (exists("trw", fb)) {
      if (exists("nops", fb)) {
        if (exists("rytha", fb))
          ow <- c(ow, "rytha")
        fb$rytha <- fb$trw / fb$nops * value / 1000
        fb$rytha[fb$nops == 0] <- NA
      }
      if (exists("noph", fb)) {
        if (exists("rytha.aj", fb))
          ow <- c(ow, "rytha.aj")
        fb$rytha.aj <- fb$trw / fb$noph * value / 1000
        fb$rytha.aj[fb$noph == 0] <- NA
      }
    }

    if (exists("vw", fb)) {
      if (exists("nops", fb)) {
        if (exists("fytha", fb))
          ow <- c(ow, "fytha")
        fb$fytha <- fb$vw / fb$nops * value / 1000
        fb$fytha[fb$nops == 0] <- NA
      }
      if (exists("noph", fb)) {
        if (exists("fytha.aj", fb))
          ow <- c(ow, "fytha.aj")
        fb$fytha.aj <- fb$vw / fb$noph * value / 1000
        fb$fytha.aj[fb$noph == 0] <- NA
      }
    }

    if (exists("trw.d", fb)) {
      if (exists("nops", fb)) {
        if (exists("dmry", fb))
          ow <- c(ow, "dmry")
        fb$dmry <- fb$trw.d / fb$nops * value / 1000
        fb$dmry[fb$nops == 0] <- NA
      }
      if (exists("noph", fb)) {
        if (exists("dmry.aj", fb))
          ow <- c(ow, "dmry.aj")
        fb$dmry.aj <- fb$trw.d / fb$noph * value / 1000
        fb$dmry.aj[fb$noph == 0] <- NA
      }
    }
    
    if (exists("vw.d", fb)) {
      if (exists("nops", fb)) {
        if (exists("dmvy", fb))
          ow <- c(ow, "dmvy")
        fb$dmvy <- fb$vw.d / fb$nops * value / 1000
        fb$dmvy[fb$nops == 0] <- NA
      }
      if (exists("noph", fb)) {
        if (exists("dmvy.aj", fb))
          ow <- c(ow, "dmvy.aj")
        fb$dmvy.aj <- fb$vw.d / fb$noph * value / 1000
        fb$dmvy.aj[fb$noph == 0] <- NA
      }
    }
  }
  
  # General computations a posteriori
  
  if (exists("rytha", fb) & exists("fytha", fb)) {
    if (exists("biom", fb))
      ow <- c(ow, "biom")
    fb$biom <- suma(fb$rytha, fb$fytha)
  }
  
  if (exists("rytha.aj", fb) & exists("fytha.aj", fb)) {
    if (exists("biom.aj", fb))
      ow <- c(ow, "biom.aj")
    fb$biom.aj <- suma(fb$rytha.aj, fb$fytha.aj)
  }
  
  if (exists("rytha", fb) & exists("fytha", fb) & exists("dm", fb) & exists("dmv", fb)) {
    if (exists("dmbiom", fb))
      ow <- c(ow, "dmbiom")
    fb$dmbiom <- suma(fb$rytha * fb$dm / 100, fb$fytha * fb$dmv / 100)
  }
  
  if (exists("rytha.aj", fb) & exists("fytha.aj", fb) & exists("dm", fb) & exists("dmv", fb)) {
    if (exists("dmbiom.aj", fb))
      ow <- c(ow, "dmbiom.aj")
    fb$dmbiom.aj <- suma(fb$rytha.aj * fb$dm / 100, fb$fytha.aj * fb$dmv / 100)
  }

  if (exists("trw.d", fb) & exists("vw.d", fb)) {
    if (exists("rfr", fb))
      ow <- c(ow, "rfr")
    fb$rfr <- fb$trw.d / fb$vw.d
    fb$rfr[fb$vw.d == 0] <- NA
  }

  # Betacarotene from color chart
  
  if (exists("rfc.cc", fb)) {
    if (exists("bc.cc", fb))
      ow <- c(ow, "bc.cc")
    fb$bc_cc[fb$rfc.cc == "1"] <- 0.03
    fb$bc_cc[fb$rfc.cc == "2"] <- 0
    fb$bc_cc[fb$rfc.cc == "3"] <- 0.12
    fb$bc_cc[fb$rfc.cc == "4"] <- 0.02
    fb$bc_cc[fb$rfc.cc == "5"] <- 0
    fb$bc_cc[fb$rfc.cc == "6"] <- 0.15
    fb$bc_cc[fb$rfc.cc == "7"] <- 1.38
    fb$bc_cc[fb$rfc.cc == "8"] <- 1.65
    fb$bc_cc[fb$rfc.cc == "9"] <- 1.5
    fb$bc_cc[fb$rfc.cc == "10"] <- 1.74
    fb$bc_cc[fb$rfc.cc == "11"] <- 1.76
    fb$bc_cc[fb$rfc.cc == "12"] <- 0.69
    fb$bc_cc[fb$rfc.cc == "13"] <- 1.17
    fb$bc_cc[fb$rfc.cc == "14"] <- 1.32
    fb$bc_cc[fb$rfc.cc == "15"] <- 1.04
    fb$bc_cc[fb$rfc.cc == "16"] <- 4.41
    fb$bc_cc[fb$rfc.cc == "17"] <- 4.92
    fb$bc_cc[fb$rfc.cc == "18"] <- 6.12
    fb$bc_cc[fb$rfc.cc == "19"] <- 5.46
    fb$bc_cc[fb$rfc.cc == "20"] <- 3.96
    fb$bc_cc[fb$rfc.cc == "21"] <- 5.49
    fb$bc_cc[fb$rfc.cc == "22"] <- 3.03
    fb$bc_cc[fb$rfc.cc == "23"] <- 3.76
    fb$bc_cc[fb$rfc.cc == "24"] <- 4.61
    fb$bc_cc[fb$rfc.cc == "25"] <- 7.23
    fb$bc_cc[fb$rfc.cc == "26"] <- 7.76
    fb$bc_cc[fb$rfc.cc == "27"] <- 10.5
    fb$bc_cc[fb$rfc.cc == "28"] <- 11.03
    fb$bc_cc[fb$rfc.cc == "29"] <- 12.39
    fb$bc_cc[fb$rfc.cc == "30"] <- 14.37
  }
  
  # Warning: Overwritten traits
  
  if (length(ow) > 0)
    warning("Some traits have been overwritten: ", list(ow), call. = FALSE)

  # Return
  
  fb

}
