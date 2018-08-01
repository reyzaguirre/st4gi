#' Check data for an ABD
#'
#' This function checks the frequencies of genotypes in an ABD.
#' @param trait The trait to analyze.
#' @param geno The genotypes including checks.
#' @param rep The replications.
#' @param dfr The name of the data frame.
#' @return The number of checks \code{ng.check}, the number of no checks \code{ng},
#' the number of missing values for checks \code{nmis.check}), the number of
#' missing values for no checks \code{nmis}, the number \code{ncheck.0}
#' and list \code{check.0} of checks without data, the number \code{ncheck.1}
#' and list \code{check.1} of checks with only one datum, the number of checks
#' with at least two data \code{ncheck.2}, the number of replications \code{nr},
#' and the number of rows in the data frame with missing values for factors
#' (\code{nmis.fac}).
#' @author Raul Eyzaguirre.
#' @examples
#' # Create design
#' dfr <- cr.abd(1:50, c('a', 'b', 'd'), 5, 10)
#' dfr <- dfr$book
#' # Create some random data
#' dfr$y <- rnorm(65)
#' # Delete some values
#' dfr[c(1, 5, 7, 56), 'y'] <- NA
#' # Delete some values for classification factors
#' dfr[64, 'geno'] <- NA
#' # Check the design
#' ck.abd('y', 'geno', 'block', dfr)
#' @export

ck.abd <- function(trait, geno, rep, dfr) {
  
  # Check and remove rows with missing values for factors
  
  out <- rm.fna(c(geno, rep), dfr)
  dfr <- out$dfr
  nmis.fac <- out$nmis.fac

  # Number and list of checks and nochecks, and number of replications
  
  out <- ck.fs(geno, rep, 'abd', dfr)
  ng <- out$nt
  lg <- out$lt
  ng.check <- out$nt.check
  lg.check <- out$lt.check
  nr <- out$nr
  
  # Number of missing values for no checks
  
  temp <- dfr[dfr[, geno] %in% lg, ]
  nmis <- sum(is.na(temp[, trait]))
  
  # Number of missing values for checks

  temp <- dfr[dfr[, geno] %in% lg.check, ]
  out <- ck.fq(trait, geno, rep, temp)
  tfr <- out$tfr
  nmis.check <- sum(tfr == 0)

  # Number of checks without data, with 1, and more data

  tf <- data.frame(out$tf)
  
  ncheck.0 <- sum(tf$Freq == 0)
  ncheck.1 <- sum(tf$Freq == 1)
  ncheck.2 <- sum(tf$Freq > 1)
  
  # List of checks witout data or only one datum
  
  check.0 <- NULL
  if (ncheck.0 > 0)
    check.0 <- tfreq[tfreq$Freq == 0, 1]
  
  check.1 <- NULL
  if (ncheck.1 > 0)
    check.1 <- tfreq[tfreq$Freq == 1, 1]

  # Return
  
  list(ng.check = ng.check, ng = ng, nmis.check = nmis.check, nmis = nmis,
       ncheck.0 = ncheck.0, check.0 = check.0, ncheck.1 = ncheck.1,
       check.1 = check.1, ncheck.2 = ncheck.2, nmis.fac = nmis.fac, nr = nr)
  
}
