#' Compute derived traits for sweetpotato
#'
#' Compute derived traits for a given fieldbook data frame.
#' @param dfr The name of the data frame.
#' @param method Method to scale data from plot to hectare level. Options are plot size
#' \code{"ps"} and number of plants for a full hectare \code{"np"}. See details.
#' @param value Value for the method selected in square meters if \code{method = "ps"}
#' and in number of plants per hectare if \code{method = "np"}.
#' @param nops Number of plants sowed per plot.
#' @details The data frame must use the labels (lower or upper case) listed in function
#' \code{check.names.sp}. See \code{?check.names.sp} for details.
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

cdt <- function(dfr, method = c("none", "ps", "np"), value = NULL, nops = NULL) {
  
  # Check names
  
  dfr <- check.names.sp(dfr)
  
  # Original trait names
  
  on <- names(dfr)
  
  # Warnings
  
  method = match.arg(method)

  if (method == "ps" & is.null(value))
    warning("Plot size value is missing.", call. = FALSE)
  
  if (method == "np" & is.null(value))
    warning("Total number of plants per hectare value is missing.", call. = FALSE)

  if (!is.null(nops)) {
    if (exists("nops", dfr)) {
      dfr$nops <- nops
      warning("nops has been replace in the fieldbook by ", nops, call. = FALSE)
    } else {
      dfr$nops <- nops
    }
  }
  
  if (method == "np" & !exists("nops", dfr))
    warning("Number of plants sowed, nops, is missing.", call. = FALSE)

  # List of traits to overwrite

  ow <- NULL 
  
  # General computations a priori
  
  if (exists("crw", dfr) & exists("ncrw", dfr)) {
    if ("trw" %in% on)
      ow <- c(ow, "trw")
    dfr$trw <- suma(dfr$crw, dfr$ncrw)
  }

  if (exists("trw", dfr) & exists("vw", dfr)) {
    if ("biom" %in% on)
      ow <- c(ow, "biom")
    dfr$biom <- suma(dfr$trw, dfr$vw)
  }

  if (exists("crw", dfr) & exists("nocr", dfr)) {
    if ("acrw" %in% on)
      ow <- c(ow, "acrw")
    dfr$acrw <- dfr$crw / dfr$nocr
    dfr$acrw[dfr$nocr == 0] <- NA
  }
  
  if (exists("ncrw", dfr) & exists("nonc", dfr)) {
    if ("ancrw" %in% on)
      ow <- c(ow, "ancrw")
    dfr$ancrw <- dfr$ncrw / dfr$nonc
    dfr$ancrw[dfr$nonc == 0] <- NA
  }

  if (exists("trw", dfr) & exists("tnr", dfr)) {
    if ("atrw" %in% on)
      ow <- c(ow, "atrw")
    dfr$atrw <- dfr$trw / dfr$tnr
    dfr$atrw[dfr$tnr == 0] <- NA
  }

  if (exists("nocr", dfr) & exists("nonc", dfr)) {
    if ("tnr" %in% on)
      ow <- c(ow, "tnr")
    dfr$tnr <- suma(dfr$nocr, dfr$nonc)
  }

  if (exists("tnr", dfr) & exists("noph", dfr)) {
    if ("nrpp" %in% on)
      ow <- c(ow, "nrpp")
    dfr$nrpp <- dfr$tnr / dfr$noph
    dfr$nrpp[dfr$noph == 0] <- NA
  }

  if (exists("tnr", dfr) & exists("nops", dfr)) {
    if ("nrpsp" %in% on)
      ow <- c(ow, "nrpsp")
    dfr$nrpsp <- dfr$tnr / dfr$nops
  }

  if (exists("nocr", dfr) & exists("noph", dfr)) {
    if ("ncrpp" %in% on)
      ow <- c(ow, "ncrpp")
    dfr$ncrpp <- dfr$nocr / dfr$noph
    dfr$ncrpp[dfr$noph == 0] <- NA
  }

  if (exists("nocr", dfr) & exists("nops", dfr)) {
    if ("ncrpsp" %in% on)
      ow <- c(ow, "ncrpsp")
    dfr$ncrpsp <- dfr$nocr / dfr$nops
  }
  
  if (exists("trw", dfr) & exists("noph", dfr)) {
    if ("ypp" %in% on)
      ow <- c(ow, "ypp")
    dfr$ypp <- dfr$trw / dfr$noph
    dfr$ypp[dfr$noph == 0] <- NA
  }

  if (exists("trw", dfr) & exists("nops", dfr)) {
    if ("ypsp" %in% on)
      ow <- c(ow, "ypsp")
    dfr$ypsp <- dfr$trw / dfr$nops
  }

  if (exists("vw", dfr) & exists("noph", dfr)) {
    if ("vpp" %in% on)
      ow <- c(ow, "vpp")
    dfr$vpp <- dfr$vw / dfr$noph
    dfr$vpp[dfr$noph == 0] <- NA
  }

  if (exists("vw", dfr) & exists("nops", dfr)) {
    if ("vpsp" %in% on)
      ow <- c(ow, "vpsp")
    dfr$vpsp <- dfr$vw / dfr$nops
  }
  
  if (exists("nocr", dfr) & exists("nonc", dfr)) {
    if ("ci" %in% on)
      ow <- c(ow, "ci")
    dfr$ci <- dfr$nocr / dfr$tnr * 100
    dfr$ci[dfr$tnr == 0] <- NA
  }

  if (exists("trw", dfr) & exists("vw", dfr)) {
    if ("hi" %in% on)
      ow <- c(ow, "hi")
    dfr$hi <- dfr$trw / dfr$biom * 100
    dfr$hi[dfr$biom == 0] <- NA
  }
  
  if (exists("noph", dfr) & exists("nops", dfr)) {
    if ("shi" %in% on)
      ow <- c(ow, "shi")
    dfr$shi <- dfr$noph / dfr$nops * 100
    dfr$shi[dfr$nops == 0] <- NA
  }
  
  if (exists("dmd", dfr) & exists("dmf", dfr)) {
    if ("dm" %in% on)
      ow <- c(ow, "dm")
    dfr$dm <- dfr$dmd / dfr$dmf * 100
    dfr$dm[dfr$dmf == 0] <- NA
  }
  
  if (exists("dmvd", dfr) & exists("dmvf", dfr)) {
    if ("dmv" %in% on)
      ow <- c(ow, "dmv")
    dfr$dmv <- dfr$dmvd / dfr$dmvf * 100
    dfr$dmv[dfr$dmvf == 0] <- NA
  }
  
  if (exists("trw", dfr) & exists("dm", dfr)) {
    if ("trw.d" %in% on)
      ow <- c(ow, "trw.d")
    dfr$trw.d <- dfr$trw * dfr$dm / 100
  }

  if (exists("vw", dfr) & exists("dmv", dfr)) {
    if ("vw.d" %in% on)
      ow <- c(ow, "vw.d")
    dfr$vw.d <- dfr$vw * dfr$dmv / 100
  }

  if (exists("trw.d", dfr) & exists("vw.d", dfr)) {
    if ("biom.d" %in% on)
      ow <- c(ow, "biom.d")
    dfr$biom.d <- suma(dfr$trw.d, dfr$vw.d)
  }
  
  # Computations based on plot size
  
  if (method == "ps" & !is.null(value)) {

    if (exists("crw", dfr)) {
      if ("cytha" %in% on)
        ow <- c(ow, "cytha")
      dfr$cytha <- dfr$crw * 10 / value
      if (exists("noph", dfr) & exists("nops", dfr)) {
        if ("cytha.aj" %in% on)
          ow <- c(ow, "cytha.aj")
        dfr$cytha.aj <- dfr$crw / dfr$noph * dfr$nops * 10 / value
        dfr$cytha.aj[dfr$noph == 0] <- NA
      }
    }

    if (exists("trw", dfr)) {
      if ("rytha" %in% on)
        ow <- c(ow, "rytha")
      dfr$rytha <- dfr$trw * 10 / value
      if (exists("noph", dfr) & exists("nops", dfr)) {
        if ("rytha.aj" %in% on)
          ow <- c(ow, "rytha.aj")
        dfr$rytha.aj <- dfr$trw / dfr$noph * dfr$nops * 10 / value
        dfr$rytha.aj[dfr$noph == 0] <- NA
      }
    }
    
    if (exists("vw", dfr)) {
      if ("fytha" %in% on)
        ow <- c(ow, "fytha")
      dfr$fytha <- dfr$vw * 10 / value
      if (exists("noph", dfr) & exists("nops", dfr)) {
        if ("fytha.aj" %in% on)
          ow <- c(ow, "fytha.aj")
        dfr$fytha.aj <- dfr$vw / dfr$noph * dfr$nops * 10 / value
        dfr$fytha.aj[dfr$noph == 0] <- NA
      }
    }

    if (exists("trw.d", dfr)) {
      if ("dmry" %in% on)
        ow <- c(ow, "dmry")
      dfr$dmry <- dfr$trw.d * 10 / value
      if (exists("noph", dfr) & exists("nops", dfr)) {
        if ("dmry.aj" %in% on)
          ow <- c(ow, "dmry.aj")
        dfr$dmry.aj <- dfr$trw.d / dfr$noph * dfr$nops * 10 / value
        dfr$dmry.aj[dfr$noph == 0] <- NA
      }
    }
   
    if (exists("vw.d", dfr)) {
      if ("dmvy" %in% on)
        ow <- c(ow, "dmvy")
      dfr$dmvy <- dfr$vw.d * 10 / value
      if (exists("noph", dfr) & exists("nops", dfr)) {
        if ("dmvy.aj" %in% on)
          ow <- c(ow, "dmvy.aj")
        dfr$dmvy.aj <- dfr$vw.d / dfr$noph * dfr$nops * 10 / value
        dfr$dmvy.aj[dfr$noph == 0] <- NA
      }
    }
  
  }
  
  # Computations based on number of plants
  
  if (method == "np" & !is.null(value)) {

    if (exists("crw", dfr)) {
      if (exists("nops", dfr)) {
        if ("cytha" %in% on)
          ow <- c(ow, "cytha")
        dfr$cytha <- dfr$crw / dfr$nops * value / 1000
        dfr$cytha[dfr$nops == 0] <- NA
      }
      if (exists("noph", dfr)) {
        if ("cytha.aj" %in% on)
          ow <- c(ow, "cytha.aj")
        dfr$cytha.aj <- dfr$crw / dfr$noph * value / 1000
        dfr$cytha.aj[dfr$noph == 0] <- NA
      }
    }

    if (exists("trw", dfr)) {
      if (exists("nops", dfr)) {
        if ("rytha" %in% on)
          ow <- c(ow, "rytha")
        dfr$rytha <- dfr$trw / dfr$nops * value / 1000
        dfr$rytha[dfr$nops == 0] <- NA
      }
      if (exists("noph", dfr)) {
        if ("rytha.aj" %in% on)
          ow <- c(ow, "rytha.aj")
        dfr$rytha.aj <- dfr$trw / dfr$noph * value / 1000
        dfr$rytha.aj[dfr$noph == 0] <- NA
      }
    }

    if (exists("vw", dfr)) {
      if (exists("nops", dfr)) {
        if ("fytha" %in% on)
          ow <- c(ow, "fytha")
        dfr$fytha <- dfr$vw / dfr$nops * value / 1000
        dfr$fytha[dfr$nops == 0] <- NA
      }
      if (exists("noph", dfr)) {
        if ("fytha.aj" %in% on)
          ow <- c(ow, "fytha.aj")
        dfr$fytha.aj <- dfr$vw / dfr$noph * value / 1000
        dfr$fytha.aj[dfr$noph == 0] <- NA
      }
    }

    if (exists("trw.d", dfr)) {
      if (exists("nops", dfr)) {
        if ("dmry" %in% on)
          ow <- c(ow, "dmry")
        dfr$dmry <- dfr$trw.d / dfr$nops * value / 1000
        dfr$dmry[dfr$nops == 0] <- NA
      }
      if (exists("noph", dfr)) {
        if ("dmry.aj" %in% on)
          ow <- c(ow, "dmry.aj")
        dfr$dmry.aj <- dfr$trw.d / dfr$noph * value / 1000
        dfr$dmry.aj[dfr$noph == 0] <- NA
      }
    }
    
    if (exists("vw.d", dfr)) {
      if (exists("nops", dfr)) {
        if ("dmvy" %in% on)
          ow <- c(ow, "dmvy")
        dfr$dmvy <- dfr$vw.d / dfr$nops * value / 1000
        dfr$dmvy[dfr$nops == 0] <- NA
      }
      if (exists("noph", dfr)) {
        if ("dmvy.aj" %in% on)
          ow <- c(ow, "dmvy.aj")
        dfr$dmvy.aj <- dfr$vw.d / dfr$noph * value / 1000
        dfr$dmvy.aj[dfr$noph == 0] <- NA
      }
    }
  }
  
  # General computations a posteriori
  
  if (exists("rytha", dfr) & exists("fytha", dfr)) {
    if ("bytha" %in% on)
      ow <- c(ow, "bytha")
    dfr$bytha <- suma(dfr$rytha, dfr$fytha)
  }
  
  if (exists("rytha.aj", dfr) & exists("fytha.aj", dfr)) {
    if ("bytha.aj" %in% on)
      ow <- c(ow, "bytha.aj")
    dfr$bytha.aj <- suma(dfr$rytha.aj, dfr$fytha.aj)
  }
  
  if (exists("rytha", dfr) & exists("fytha", dfr) & exists("dm", dfr) & exists("dmv", dfr)) {
    if ("dmby" %in% on)
      ow <- c(ow, "dmby")
    dfr$dmby <- suma(dfr$rytha * dfr$dm / 100, dfr$fytha * dfr$dmv / 100)
  }
  
  if (exists("rytha.aj", dfr) & exists("fytha.aj", dfr) & exists("dm", dfr) & exists("dmv", dfr)) {
    if ("dmby.aj" %in% on)
      ow <- c(ow, "dmby.aj")
    dfr$dmby.aj <- suma(dfr$rytha.aj * dfr$dm / 100, dfr$fytha.aj * dfr$dmv / 100)
  }

  if (exists("trw.d", dfr) & exists("vw.d", dfr)) {
    if ("rfr" %in% on)
      ow <- c(ow, "rfr")
    dfr$rfr <- dfr$trw.d / dfr$vw.d * 100
    dfr$rfr[dfr$vw.d == 0] <- NA
  }

  # Betacarotene from color chart
  
  if (exists("fcol.cc", dfr)) {
    if ("bc.cc" %in% on)
      ow <- c(ow, "bc.cc")
    dfr$bc.cc[dfr$fcol.cc == "1"] <- 0.03
    dfr$bc.cc[dfr$fcol.cc == "2"] <- 0
    dfr$bc.cc[dfr$fcol.cc == "3"] <- 0.12
    dfr$bc.cc[dfr$fcol.cc == "4"] <- 0.02
    dfr$bc.cc[dfr$fcol.cc == "5"] <- 0
    dfr$bc.cc[dfr$fcol.cc == "6"] <- 0.15
    dfr$bc.cc[dfr$fcol.cc == "7"] <- 1.38
    dfr$bc.cc[dfr$fcol.cc == "8"] <- 1.65
    dfr$bc.cc[dfr$fcol.cc == "9"] <- 1.5
    dfr$bc.cc[dfr$fcol.cc == "10"] <- 1.74
    dfr$bc.cc[dfr$fcol.cc == "11"] <- 1.76
    dfr$bc.cc[dfr$fcol.cc == "12"] <- 0.69
    dfr$bc.cc[dfr$fcol.cc == "13"] <- 1.17
    dfr$bc.cc[dfr$fcol.cc == "14"] <- 1.32
    dfr$bc.cc[dfr$fcol.cc == "15"] <- 1.04
    dfr$bc.cc[dfr$fcol.cc == "16"] <- 4.41
    dfr$bc.cc[dfr$fcol.cc == "17"] <- 4.92
    dfr$bc.cc[dfr$fcol.cc == "18"] <- 6.12
    dfr$bc.cc[dfr$fcol.cc == "19"] <- 5.46
    dfr$bc.cc[dfr$fcol.cc == "20"] <- 3.96
    dfr$bc.cc[dfr$fcol.cc == "21"] <- 5.49
    dfr$bc.cc[dfr$fcol.cc == "22"] <- 3.03
    dfr$bc.cc[dfr$fcol.cc == "23"] <- 3.76
    dfr$bc.cc[dfr$fcol.cc == "24"] <- 4.61
    dfr$bc.cc[dfr$fcol.cc == "25"] <- 7.23
    dfr$bc.cc[dfr$fcol.cc == "26"] <- 7.76
    dfr$bc.cc[dfr$fcol.cc == "27"] <- 10.5
    dfr$bc.cc[dfr$fcol.cc == "28"] <- 11.03
    dfr$bc.cc[dfr$fcol.cc == "29"] <- 12.39
    dfr$bc.cc[dfr$fcol.cc == "30"] <- 14.37
  }
  
  # Warning: Overwritten traits
  
  if (length(ow) > 0)
    warning("Some traits have been overwritten: ", list(ow), call. = FALSE)

  # Return
  
  dfr

}
