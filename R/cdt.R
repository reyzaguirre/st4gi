#' Compute derived variables for potato and sweetpotato
#'
#' Compute derived variables for a given fieldbook data frame.
#' @param dfr The name of the data frame.
#' @param method Method to scale data from plot to hectare level. Options are
#' plot size \code{"ps"} and number of plants for a full hectare \code{"np"}.
#' See details.
#' @param value Value for the method selected in square meters if \code{method = "ps"}
#' and in number of plants per hectare if \code{method = "np"}.
#' @param crop \code{"auto"} for autodetection or \code{"pt"} for potato and \code{"sp"} for sweetpotato.
#' @details The data frame must use the labels (lower or upper case) listed in
#' functions \code{pt.ont()} and \code{sp.ont()}. 
#' Conversion from kilograms per plot to tons per hectare can be done using
#' \code{ps}, the plot size, or \code{np}, the total number of plants that can be
#' allocated in a full hectare. In both cases computations can be adjusted by the
#' number of harvested plants if available in the fieldbook. For \code{method = "np"},
#' number of plants sowed must be specified to compute non adjusted values.
#' @return It returns a data frame with the original and derived variables.
#' @author Raul Eyzaguirre.
#' @examples
#' cdt(potatoyield)
#' cdt(pjpz09)
#' @export

cdt <- function(dfr, method = c("none", "ps", "np"), value = NULL,
                   crop = c('auto', 'pt', 'sp')) {
  
  # Match arguments
  
  method = match.arg(method)
  crop = match.arg(crop)
  
  # Check names
  
  dfr <- check.names(dfr, crop = crop)
  if (crop == 'auto')
    crop <- detect.crop(dfr)
  
  # Original variable names
  
  on <- names(dfr)
  
  # Warnings
  
  if (method == "ps" & is.null(value))
    warning("Plot size value is missing.", call. = FALSE)
  
  if (method == "np" & is.null(value))
    warning("Total number of plants per hectare value is missing.", call. = FALSE)

  if (method == "np" & crop == 'pt' & !exists("ntp", dfr))
    warning("Number of tubers planted, ntp, is missing.", call. = FALSE)
  
  if (method == "np" & crop == 'sp' & !exists("nops", dfr))
    warning("Number of plants sowed, nops, is missing.", call. = FALSE)

  # List of variables to overwrite
  
  ow <- NULL
  
  # -------------------------------
  # Compute variables for potato
  # -------------------------------
  
  if (crop == 'pt') {
    
    # General computations for marketable and nonmarketable tubers
    
    if (exists("mtwci", dfr) & exists("mtwcii", dfr)) {
      if ("mtwp" %in% on)
        ow <- c(ow, "mtwp")
      dfr$mtwp <- dfr$mtwci + dfr$mtwcii
    }
    
    if (exists("mtwp", dfr) & exists("nomtwp", dfr)) {
      if ("ttwp" %in% on)
        ow <- c(ow, "ttwp")
      dfr$ttwp <- dfr$mtwp + dfr$nomtwp
    }
    
    if (exists("nmtci", dfr) & exists("nmtcii", dfr)) {
      if ("nmtp" %in% on)
        ow <- c(ow, "nmtp")
      dfr$nmtp <- dfr$nmtci + dfr$nmtcii
    }
    
    if (exists("nmtp", dfr) & exists("nnomtp", dfr)) {
      if ("tntp" %in% on)
        ow <- c(ow, "tntp")
      dfr$tntp <- dfr$nmtp + dfr$nnomtp
    }
    
    if (exists("mtwp", dfr) & exists("nmtp", dfr)) {
      if ("atmw" %in% on)
        ow <- c(ow, "atmw")
      dfr$atmw <- dfr$mtwp / dfr$nmtp * 1000
      dfr$atmw[dfr$nmtp == 0] <- NA
    }
    
    if (exists("ttwp", dfr) & exists("tntp", dfr)) {
      if ("atw" %in% on)
        ow <- c(ow, "atw")
      dfr$atw <- dfr$ttwp / dfr$tntp * 1000
      dfr$atw[dfr$tntp == 0] <- NA
    }
    
    if (exists("tntp", dfr) & exists("nph", dfr)) {
      if ("tntpl" %in% on)
        ow <- c(ow, "tntpl")
      dfr$tntpl <- dfr$tntp / dfr$nph
      dfr$tntpl[dfr$nph == 0] <- NA
    }
    
    if (exists("nmtp", dfr) & exists("nph", dfr)) {
      if ("nmtpl" %in% on)
        ow <- c(ow, "nmtpl")
      dfr$nmtpl <- dfr$nmtp / dfr$nph
      dfr$nmtpl[dfr$nph == 0] <- NA
    }
    
    if (exists("ttwp", dfr) & exists("nph", dfr)) {
      if ("ttwpl" %in% on)
        ow <- c(ow, "ttwpl")
      dfr$ttwpl <- dfr$ttwp / dfr$nph
      dfr$ttwpl[dfr$nph == 0] <- NA
    }
    
    if (exists("mtwp", dfr) & exists("nph", dfr)) {
      if ("mtwpl" %in% on)
        ow <- c(ow, "mtwpl")
      dfr$mtwpl <- dfr$mtwp / dfr$nph
      dfr$mtwpl[dfr$nph == 0] <- NA
    }
    
    # General computations for dry weight of tubers
    
    if (exists("dwts1", dfr) & exists("fwts1", dfr)) {
      if ("dm1" %in% on)
        ow <- c(ow, "dm1")
      dfr$dm1 <- dfr$dwts1 / dfr$fwts1 * 100
      dfr$dm1[dfr$fwts1 == 0] <- NA
    }
    
    if (exists("dwts2", dfr) & exists("fwts2", dfr)) {
      if ("dm2" %in% on)
        ow <- c(ow, "dm2")
      dfr$dm2 <- dfr$dwts2 / dfr$fwts2 * 100
      dfr$dm2[dfr$fwts2 == 0] <- NA
    }
    
    if (exists("dm1", dfr) & exists("dm2", dfr)) {
      if ("dm" %in% on)
        ow <- c(ow, "dm")
      dfr$dm <- apply(dfr[, c("dm1", "dm2")], 1, mean, na.rm = TRUE)
    }
    
    if (exists("dwts", dfr) & exists("fwts", dfr)) {
      if ("dm" %in% on)
        ow <- c(ow, "dm")
      dfr$dm <- dfr$dwts / dfr$fwts * 100
      dfr$dm[dfr$fwts == 0] <- NA
    }
    
    # General computation for harvest index
    
    if (exists("ttwp", dfr) & exists("tbfwp", dfr)) {
      if ("hi_fw" %in% on)
        ow <- c(ow, "hi_fw")
      dfr$hi_fw <- dfr$ttwp * 1000 / dfr$tbfwp
    }
    
    # Percentages for plants emerged and harvested
    
    if (exists("npe", dfr) & exists("ntp", dfr)) {
      if ("ppe" %in% on)
        ow <- c(ow, "ppe")
      dfr$ppe <- dfr$npe / dfr$ntp * 100
    }
    
    if (exists("nph", dfr) & exists("ntp", dfr)) {
      if ("pph" %in% on)
        ow <- c(ow, "pph")
      dfr$pph <- dfr$nph / dfr$ntp * 100
    }
    
    # Computations based on plot size
    
    if (method == "ps" & !is.null(value)) {
      
      if (exists("mtwp", dfr)) {
        if ("mtyna" %in% on)
          ow <- c(ow, "mtyna")
        dfr$mtyna <- dfr$mtwp * 10 / value
        if (exists("nph", dfr) & exists("ntp", dfr)) {
          if ("mtya" %in% on)
            ow <- c(ow, "mtya")
          dfr$mtya <- dfr$mtwp / dfr$nph * dfr$ntp * 10 / value
          dfr$mtya[dfr$nph == 0] <- NA
        }
      }
      
      if (exists("ttwp", dfr)) {
        if ("ttyna" %in% on)
          ow <- c(ow, "ttyna")
        dfr$ttyna <- dfr$ttwp * 10 / value
        if (exists("nph", dfr) & exists("ntp", dfr)) {
          if ("ttya" %in% on)
            ow <- c(ow, "ttya")
          dfr$ttya <- dfr$ttwp / dfr$nph * dfr$ntp * 10 / value
          dfr$ttya[dfr$nph == 0] <- NA
        }
      }
      
    }
    
    # Computations based on number of plants
    
    if (method == "np" & !is.null(value)) {
      
      if (exists("mtwp", dfr)) {
        if (exists("ntp", dfr)) {
          if ("mtyna" %in% on)
            ow <- c(ow, "mtyna")
          dfr$mtyna <- dfr$mtwp / dfr$ntp * value / 1000
          dfr$mtyna[dfr$ntp == 0] <- NA
        }
        if (exists("nph", dfr)) {
          if ("mtya" %in% on)
            ow <- c(ow, "mtya")
          dfr$mtya <- dfr$mtwp / dfr$nph * value / 1000
          dfr$mtya[dfr$nph == 0] <- NA
        }
      }
      
      if (exists("ttwp", dfr)) {
        if (exists("ntp", dfr)) {
          if ("ttyna" %in% on)
            ow <- c(ow, "ttyna")
          dfr$ttyna <- dfr$ttwp / dfr$ntp * value / 1000
          dfr$ttyna[dfr$ntp == 0] <- NA
        }
        if (exists("nph", dfr)) {
          if ("ttya" %in% on)
            ow <- c(ow, "ttya")
          dfr$ttya <- dfr$ttwp / dfr$nph * value / 1000
          dfr$ttya[dfr$nph == 0] <- NA
        }
      }
      
    }
    
  }
  
  # ----------------------------------
  # Compute variables for sweetpotato
  # ----------------------------------
  
  if (crop == 'sp') {

    # General computations a priori
    
    if (exists("crw", dfr) & exists("ncrw", dfr)) {
      if ("trw" %in% on)
        ow <- c(ow, "trw")
      dfr$trw <- dfr$crw + dfr$ncrw
    }
    
    if (exists("trw", dfr) & exists("vw", dfr)) {
      if ("biom" %in% on)
        ow <- c(ow, "biom")
      dfr$biom <- dfr$trw + dfr$vw
    }
    
    if (exists("nocr", dfr) & exists("nonc", dfr)) {
      if ("tnr" %in% on)
        ow <- c(ow, "tnr")
      dfr$tnr <- dfr$nocr + dfr$nonc
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
      dfr$biom.d <- dfr$trw.d + dfr$vw.d
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
      dfr$bytha <- dfr$rytha + dfr$fytha
    }
    
    if (exists("rytha.aj", dfr) & exists("fytha.aj", dfr)) {
      if ("bytha.aj" %in% on)
        ow <- c(ow, "bytha.aj")
      dfr$bytha.aj <- dfr$rytha.aj + dfr$fytha.aj
    }
    
    if (exists("rytha", dfr) & exists("fytha", dfr) & exists("dm", dfr) & exists("dmv", dfr)) {
      if ("dmby" %in% on)
        ow <- c(ow, "dmby")
      dfr$dmby <- dfr$rytha * dfr$dm / 100 + dfr$fytha * dfr$dmv / 100
    }
    
    if (exists("rytha.aj", dfr) & exists("fytha.aj", dfr) & exists("dm", dfr) & exists("dmv", dfr)) {
      if ("dmby.aj" %in% on)
        ow <- c(ow, "dmby.aj")
      dfr$dmby.aj <- dfr$rytha.aj * dfr$dm / 100 + dfr$fytha.aj * dfr$dmv / 100
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
      dfr$bc.cc[dfr$fcol.cc == 1] <- 0.03
      dfr$bc.cc[dfr$fcol.cc == 2] <- 0
      dfr$bc.cc[dfr$fcol.cc == 3] <- 0.12
      dfr$bc.cc[dfr$fcol.cc == 4] <- 0.02
      dfr$bc.cc[dfr$fcol.cc == 5] <- 0
      dfr$bc.cc[dfr$fcol.cc == 6] <- 0.15
      dfr$bc.cc[dfr$fcol.cc == 7] <- 1.38
      dfr$bc.cc[dfr$fcol.cc == 8] <- 1.65
      dfr$bc.cc[dfr$fcol.cc == 9] <- 1.5
      dfr$bc.cc[dfr$fcol.cc == 10] <- 1.74
      dfr$bc.cc[dfr$fcol.cc == 11] <- 1.76
      dfr$bc.cc[dfr$fcol.cc == 12] <- 0.69
      dfr$bc.cc[dfr$fcol.cc == 13] <- 1.17
      dfr$bc.cc[dfr$fcol.cc == 14] <- 1.32
      dfr$bc.cc[dfr$fcol.cc == 15] <- 1.04
      dfr$bc.cc[dfr$fcol.cc == 16] <- 4.41
      dfr$bc.cc[dfr$fcol.cc == 17] <- 4.92
      dfr$bc.cc[dfr$fcol.cc == 18] <- 6.12
      dfr$bc.cc[dfr$fcol.cc == 19] <- 5.46
      dfr$bc.cc[dfr$fcol.cc == 20] <- 3.96
      dfr$bc.cc[dfr$fcol.cc == 21] <- 5.49
      dfr$bc.cc[dfr$fcol.cc == 22] <- 3.03
      dfr$bc.cc[dfr$fcol.cc == 23] <- 3.76
      dfr$bc.cc[dfr$fcol.cc == 24] <- 4.61
      dfr$bc.cc[dfr$fcol.cc == 25] <- 7.23
      dfr$bc.cc[dfr$fcol.cc == 26] <- 7.76
      dfr$bc.cc[dfr$fcol.cc == 27] <- 10.5
      dfr$bc.cc[dfr$fcol.cc == 28] <- 11.03
      dfr$bc.cc[dfr$fcol.cc == 29] <- 12.39
      dfr$bc.cc[dfr$fcol.cc == 30] <- 14.37
    }
    
  }
  
  # Warning: Overwritten variables
  
  if (length(ow) > 0)
    warning("Some variables have been overwritten: ", list(ow), call. = FALSE)
  
  # Return
  
  dfr
  
}
