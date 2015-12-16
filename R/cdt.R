#' Compute derived traits
#'
#' Compute derived traits for a given fieldbook.
#' @param data The name of the data frame.
#' @param plot.size Plot size in square meters.
#' @details The data frame must use the labels (lower or upper case) listed in function \code{checknames}.
#' See \code{?checknames} for details.
#' @return It returns a data frame with the original and derived traits.
#' @author Raul Eyzaguirre.
#' @examples
#'  # The data
#'  head(pjpz09)
#'  str(pjpz09)
#'
#'  # Compute derived traits
#'  cdt(pjpz09, 4.5)
#' @export

cdt <- function(fb, plot.size = NULL) {
  
  # Check names
  
  fb <- checknames(fb)
  
  # Check derived traits already computed
  
  derived.traits <- c("TRW", "CYTHA", "RYTHA", "ACRW", "NRPP", "YPP", "CI", "HI",
                      "SHI", "BIOM", "FYTHA", "DM", "DMV", "DMFY", "DMRY", "RFR")
  colnames.list <- colnames(fb)
  to.overwrite <- colnames.list %in% derived.traits
  
  # Warning
  
  if (max(to.overwrite) == 1)
    warning("Some traits have been overwritten: ", list(colnames.list[to.overwrite]), call. = FALSE)
  
  # Compute derived traits
  
  if (exists("CRW", fb) & exists("NCRW", fb)) {
    fb$TRW <- apply(cbind(fb$CRW, fb$NCRW), 1, sum, na.rm = TRUE)
    fb$TRW[is.na(fb$CRW) & is.na(fb$NCRW)] <- NA
  }
    
  if (exists("CRW", fb) & !is.null(plot.size))
    fb$CYTHA <- fb$CRW * 10 / plot.size
  
  if (exists("CRW", fb) & exists("NCRW", fb) & !is.null(plot.size))
    fb$RYTHA <- fb$TRW * 10 / plot.size
  
  if (exists("CRW", fb) & exists("NOCR", fb)) {
    fb$ACRW <- fb$CRW / fb$NOCR
    fb$ACRW[fb$NOCR == 0] <- NA
  }
  
  if (exists("NOCR", fb) & exists("NONC", fb) & exists("NOPH", fb)) {
    fb$NRPP <- apply(cbind(fb$NOCR, fb$NONC), 1, sum, na.rm = TRUE) / fb$NOPH
    fb$NRPP[is.na(fb$NOCR) & is.na(fb$NONC)] <- NA
    fb$NRPP[fb$NOPH == 0] <- NA
  }
  
  if (exists("CRW", fb) & exists("NCRW", fb) & exists("NOPH", fb)) {
    fb$YPP <- fb$TRW / fb$NOPH
    fb$YPP[fb$NOPH == 0] <- NA
  }
  
  if (exists("NOCR", fb) & exists("NONC", fb)) {
    temp <- apply(cbind(fb$NOCR, fb$NONC), 1, sum, na.rm = TRUE)
    fb$CI <- fb$NOCR /  temp * 100
    fb$CI[temp == 0] <- NA
  }

  if (exists("CRW", fb) & exists("NCRW", fb) & exists("VW", fb)) {
    temp <- apply(cbind(fb$VW, fb$TRW), 1, sum, na.rm = TRUE)
    fb$HI <- fb$TRW / temp * 100
    fb$HI[temp == 0] <- NA
  }
  
  if (exists("NOPH", fb) & exists("NOPS", fb)) {
    fb$SHI <- fb$NOPH / fb$NOPS * 100
    fb$SHI[fb$NOPS == 0] <- NA
  }
  
  if (exists("CRW", fb) & exists("NCRW", fb) & exists("VW", fb) & !is.null(plot.size)) {
    fb$BIOM <- apply(cbind(fb$VW, fb$TRW), 1, sum, na.rm = TRUE) * 10 / plot.size
    fb$BIOM[is.na(fb$VW) & is.na(fb$TRW)] <- NA
  }
  
  if (exists("VW", fb) & !is.null(plot.size))
    fb$FYTHA <- fb$VW * 10 / plot.size
  
  if (exists("DMD", fb) & exists("DMF", fb)) {
    fb$DM <- fb$DMD / fb$DMF * 100
    fb$DM[fb$DMF == 0] <- NA
  }
  
  if (exists("DMVD", fb) & exists("DMVF", fb)) {
    fb$DMV <- fb$DMVD / fb$DMVF * 100
    fb$DMV[fb$DMVF == 0] <- NA
  }
  
  if (exists("VW", fb) & exists("DMVD", fb) & exists("DMVF", fb)) {
    temp1 <- fb$VW * fb$DMVD / fb$DMVF
    temp1[fb$DMVF == 0] <- NA
    if (!is.null(plot.size))
      fb$DMFY <- temp1 * 10 / plot.size
  }
  
  if (exists("CRW", fb) & exists("NCRW", fb) & exists("DMD", fb) & exists("DMF", fb)) {
    temp2 <- fb$TRW * fb$DMD / fb$DMF
    temp2[fb$DMF == 0] <- NA
    if (!is.null(plot.size))
      fb$DMRY <- temp2 * 10 / plot.size
  }

  if (exists("CRW", fb) & exists("NCRW", fb) & exists("DMD", fb) & exists("DMF", fb)
      & exists("VW", fb) & exists("DMVD", fb) & exists("DMVF", fb)) {
    fb$RFR <- temp2 / temp1
    fb$RFR[temp1 == 0] <- NA
  }
  
  # Return
  
  fb
}
