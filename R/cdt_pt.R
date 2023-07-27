#' Compute derived traits for potato
#'
#' Compute derived traits for a given fieldbook data frame.
#' @param dfr The name of the data frame.
#' @param method Method to scale data from plot to hectare level. Options are
#' plot size \code{"ps"} and number of plants for a full hectare \code{"np"}.
#' See details.
#' @param value Value for the method selected in square meters if \code{method = "ps"}
#' and in number of plants per hectare if \code{method = "np"}.
#' @param ntp Number of tubers planted per plot.
#' @details The data frame must use the labels (lower or upper case) listed in
#' function \code{check.names.pt}.
#' 
#' Conversion from kilograms per plot to tons per hectare can be done using
#' \code{ps}, the plot size, or \code{np}, the total number of plants that can be
#' allocated in a full hectare. In both cases computations can be adjusted by the
#' number of harvested plants if available in the fieldbook. For \code{method = "np"},
#' \code{ntp} must be specified to compute non adjusted values.
#' @return It returns a data frame with the original and derived traits.
#' @author Raul Eyzaguirre.
#' @examples
#' cdt.pt(potatoyield)
#' @export

cdt.pt <- function(dfr, method = c("none", "ps", "np"),
                   value = NULL, ntp = NULL) {
  
  # Match arguments
  
  method = match.arg(method)

  # Check names
  
  dfr <- check.names.pt(dfr)
  
  # Original trait names
  
  on <- names(dfr)
  
  # Warnings
  
  if (method == "ps" & is.null(value))
    warning("Plot size value is missing.", call. = FALSE)
  
  if (method == "np" & is.null(value))
    warning("Total number of plants per hectare value is missing.", call. = FALSE)

  if (!is.null(ntp)) {
    if (exists("ntp", dfr))
      warning("ntp has been replace in the fieldbook by ", ntp, call. = FALSE)
    dfr$ntp <- ntp
  }
  
  if (method == "np" & !exists("ntp", dfr))
    warning("Number of tubers planted, ntp, is missing.", call. = FALSE)

  # List of traits to overwrite

  ow <- NULL 
  
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

  # General computations for harvest index
  
  if (exists("ttwp", dfr) & exists("tbfw", dfr)) {
    if ("hi_fw" %in% on)
      ow <- c(ow, "hi_fw")
    dfr$hi_fw <- dfr$ttwp * 1000 / dfr$tbfw
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
  
  # Warning: Overwritten traits
  
  if (length(ow) > 0)
    warning("Some traits have been overwritten: ", list(ow), call. = FALSE)

  # Return
  
  dfr

}
