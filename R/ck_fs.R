#' Check factors structure
#' 
#' This function cheks the number of treatments and replications for different designs.
#' 
#' @param factors The factors
#' @param rep The replications
#' @param design The statistical design
#' @param dfr The name of the data frame
#' @return The number of factors (\code{nf}), the number of levels of the factors
#' (\code{nl}), the number of treatments (\code{nt}), the number of checks
#' (\code{nt.chk}), the list of treatments (\code{lt}), the list of checks
#' (\code{lt.chk}), and the number of replications (\code{nr}).
#' @author Raul Eyzaguirre.
#' @examples 
#' ## Example 1
#' # Create a design
#' dfr <- cr.rcbd(1:20, 3, 10)
#' dfr <- dfr$book
#' # Check the design
#' ck.fs("geno", "block", "rcbd", dfr)
#' 
#' ## Example 2
#' # Create a design
#' A <- paste0("a", 1:5)
#' B <- paste0("b", 1:3)
#' dfr <- cr.f(c("A", "B"), list(A, B), "rcbd", 3, 10)
#' dfr <- dfr$book
#' # Check the design
#' ck.fs(c("A", "B"), "block", "rcbd", dfr)
#' @export

ck.fs <- function(factors, rep = NULL, design = c('crd', 'rcbd', 'abd'), dfr) {

  # match arguments
  
  design <- match.arg(design)
  
  # Check and remove rows with missing values for factors
  
  dfr <- rm.fna(c(factors, rep), dfr)$dfr

  # Number of factors
  
  nf <- length(factors)
  
  # Number of levels for factors
  
  nl <- apply(data.frame(dfr[, factors]), 2, function(x) length(unique(x)))
  names(nl) <- factors
  
  # Define treatments
  
  treat <- factors[1]
  if (nf > 1)
    for (i in 2:nf)
      treat <- paste0(treat, factors[i])

  dfr[, treat] <- dfr[, factors[1]]
  
  if (nf > 1)
    for (i in 2:nf)
      dfr[, treat] <- paste(dfr[, treat], dfr[, factors[i]], sep = "_")

  # Number of treatments and checks
  
  if (design == 'abd') {
    tfreq <- data.frame(table(dfr[, treat]))
    lt.chk <- as.character(tfreq[tfreq$Freq > 1, 1])
    lt <- as.character(tfreq[tfreq$Freq == 1, 1])
    nt.chk <- length(lt.chk)
    nt <- length(lt)
  } else {
    lt <- as.character(unique(dfr[, treat]))
    nt <- length(lt)
    lt.chk <- NULL
    nt.chk <- NULL
  }
  
  # Number of replications
  
  if (design %in% c('abd', 'rcbd')) {
    nr <- length(unique(dfr[, rep]))
  } else {
    tfreq <- table(dfr[, treat])
    nr <- max(tfreq)
  }
  
  # Return
  
  list(nf = nf, nl = nl, nt = nt, nt.chk = nt.chk,
       lt = lt, lt.chk = lt.chk, nr = nr)
  
}
