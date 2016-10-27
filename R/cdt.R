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
  
  if (exists("nocr", fb) & exists("nonc", fb) & exists("noph", fb)) {
    if (exists("nrpp", fb))
      ow <- c(ow, "nrpp")
    fb$nrpp <- suma(fb$nocr, fb$nonc) / fb$noph
    fb$nrpp[fb$noph == 0] <- NA
  }

  if (exists("trw", fb) & exists("noph", fb)) {
    if (exists("ypp", fb))
      ow <- c(ow, "ypp")
    fb$ypp <- fb$trw / fb$noph
    fb$ypp[fb$noph == 0] <- NA
  }

  if (exists("nocr", fb) & exists("nonc", fb)) {
    if (exists("ci", fb))
      ow <- c(ow, "ci")
    temp <- suma(fb$nocr, fb$nonc)
    fb$ci <- fb$nocr / temp * 100
    fb$ci[temp == 0] <- NA
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

    if (exists("trw", fb) & exists("dm", fb)) {
      temp1 <- fb$trw * fb$dm / 100
      if (exists("dmry", fb))
        ow <- c(ow, "dmry")
      fb$dmry <- temp1 * 10 / value
      if (exists("noph", fb) & exists("nops", fb)) {
        if (exists("dmry.aj", fb))
          ow <- c(ow, "dmry.aj")
        fb$dmry.aj <- temp1 / fb$noph * fb$nops * 10 / value
        fb$dmry.aj[fb$noph == 0] <- NA
      }
    }
   
    if (exists("vw", fb) & exists("dmv", fb)) {
      temp2 <- fb$vw * fb$dmv / 100
      if (exists("dmvy", fb))
        ow <- c(ow, "dmvy")
      fb$dmvy <- temp2 * 10 / value
      if (exists("noph", fb) & exists("nops", fb)) {
        if (exists("dmvy.aj", fb))
          ow <- c(ow, "dmvy.aj")
        fb$dmvy.aj <- temp2 / fb$noph * fb$nops * 10 / value
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

    if (exists("trw", fb) & exists("dm", fb)) {
      temp1 <- fb$trw * fb$dm / 100
      if (exists("nops", fb)) {
        if (exists("dmry", fb))
          ow <- c(ow, "dmry")
        fb$dmry <- temp1 / fb$nops * value / 1000
        fb$dmry[fb$nops == 0] <- NA
      }
      if (exists("noph", fb)) {
        if (exists("dmry.aj", fb))
          ow <- c(ow, "dmry.aj")
        fb$dmry.aj <- temp1 / fb$noph * value / 1000
        fb$dmry.aj[fb$noph == 0] <- NA
      }
    }
    
    if (exists("vw", fb) & exists("dmv", fb)) {
      temp2 <- fb$vw * fb$dmv / 100
      if (exists("nops", fb)) {
        if (exists("dmvy", fb))
          ow <- c(ow, "dmvy")
        fb$dmvy <- temp2 / fb$nops * value / 1000
        fb$dmvy[fb$nops == 0] <- NA
      }
      if (exists("noph", fb)) {
        if (exists("dmvy.aj", fb))
          ow <- c(ow, "dmvy.aj")
        fb$dmvy.aj <- temp2 / fb$noph * value / 1000
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
  
  if (exists("temp1") & exists("temp2")) {
    if (exists("rfr", fb))
      ow <- c(ow, "rfr")
      fb$rfr <- temp1 / temp2
      fb$rfr[temp2 == 0] <- NA
    }

  # Warning: Overwritten traits
  
  if (length(ow) > 0)
    warning("Some traits have been overwritten: ", list(ow), call. = FALSE)

  # Return
  
  fb

}
