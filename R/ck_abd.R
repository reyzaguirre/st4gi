#' Check data for an ABD
#'
#' This function checks the frequencies of genotypes in an ABD.
#' @param trait The trait to analyze.
#' @param geno The genotypes including checks.
#' @param rep The replications.
#' @param dfr The name of the data frame.
#' @return The number of checks \code{ng.ck}, the number of no checks \code{ng},
#' the number of missing values for checks \code{nmis.ck}), the number of missing
#' values for no checks \code{nmis}, the number \code{nck.0} and list \code{ck.0}
#' of checks without data, the number \code{nck.1} and list \code{ck.1} of checks
#' with only one datum, the number of checks with at least two data \code{nck.2},
#' the number of checks that appear more than once in a given replication
#' (\code{nck.mult}), the number of replications \code{nrep}, and the number of rows
#' in the data frame with missing values for factors (\code{nmis.fac}).
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
  
  out <- ck.fs(geno, rep, dfr)
  dfr <- out$dfr
  nmis.fac <- out$nmis.fac

  # Identify checks and no checks
  
  temp <- data.frame(table(dfr[, geno]))
  lg.ck <- temp[temp$Freq > 1, 1]
  lg <- temp[temp$Freq == 1, 1]
  
  # Number of checks and no checks
  
  ng.ck <- length(lg.ck)
  ng <- length(lg)

  # Number of replications
  
  nrep <- length(unique(dfr[, rep]))

  # Number of missing values for no checks
  
  temp <- dfr[dfr[, geno] %in% lg, ]
  nmis <- sum(is.na(temp[, trait]))
  
  # Evaluate checks
  
  temp <- dfr[dfr[, geno] %in% lg.ck, ]
  out <- ck.fq(trait, geno, rep, temp)
  
  # Number of missing values for checks
  
  nmis.ck <- out$nmis
  
  # Frequencies for checks

  tf <- data.frame(out$tf)
  tfr <- out$tfr
  
  # Number of checks without data, with 1, and more data

  nck.0 <- sum(tf$Freq == 0)
  nck.1 <- sum(tf$Freq == 1)
  nck.2 <- sum(tf$Freq > 1)
  
  # List of checks witout data or only one datum
  
  ck.0 <- NULL
  if (nck.0 > 0)
    ck.0 <- tf[tf$Freq == 0, 1]
  
  ck.1 <- NULL
  if (nck.1 > 0)
    ck.1 <- tf[tf$Freq == 1, 1]

  # Number of checks that appear more than once in a given replication
  
  nck.mult <- sum(tfr > 1)
  
  # Return
  
  list(ng.ck = ng.ck, ng = ng, nmis.ck = nmis.ck, nmis = nmis,
       nck.0 = nck.0, ck.0 = ck.0, nck.1 = nck.1, ck.1 = ck.1,
       nck.2 = nck.2, nck.mult = nck.mult, nrep = nrep,
       nmis.fac = nmis.fac)
  
}
