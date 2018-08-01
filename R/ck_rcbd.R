#' Check data for a RCBD
#'
#' This function checks the frequencies of genotypes in a RCBD.
#' @param trait The trait to analyze.
#' @param geno The genotypes.
#' @param rep The replications.
#' @param dfr The name of the data frame.
#' @return The number of genotypes without data (\code{ng.0}), the number of
#' genotypes with more than one plot in a given block (\code{ng.2}), the number
#' of missing values \code{nmis}, the proportion of missing values (\code{pmis}),
#' the number of genotypes (\code{ng}), the number of replications (\code{nr}),
#' a table with frequencies of valid cases for each genotype, and the number of
#' rows in the data frame with missing values for factors (\code{nmis.fac}).
#' @author Raul Eyzaguirre.
#' @examples 
#' # Create design
#' dfr <- cr.rcbd(1:20, 3, 10)
#' dfr <- dfr$book
#' # Create some random data
#' dfr$y <- rnorm(60)
#' # Delete some values
#' dfr[c(1, 5, 16, 17), 'y'] <- NA
#' # Check the design
#' ck.rcbd('y', 'geno', 'block', dfr)
#' @export

ck.rcbd <- function(trait, geno, rep, data) {
  
  # Check and remove rows with missing values for factors
  
  out <- rm.fna(c(geno, rep), dfr)
  dfr <- out$dfr
  nmis.fac <- out$nmis.fac
  
  # Number of genotypes and replications
  
  out <- ck.gr(geno, rep, 'rcbd', dfr)
  ng <- out$ng
  nr <- out$nr

  # Genotypes and replications as factors to preserve levels in the table of frequencies
  
  temp <- dfr
  temp[, geno] <- factor(temp[, geno])
  temp[, rep] <- factor(temp[, rep])

  # Frequencies by geno
  
  temp <- temp[!is.na(temp[, trait]), ]
  tfreq <- table(temp[, geno], temp[, rep])
  
  # Number of missing values
  
  nmis <- sum(tfreq == 0)
  pmis <- mean(tfreq == 0)
  
  # Number of genotypes with more than one plot in a given block
  
  ng.2 <- sum(tfreq > 1)
  
  # Number of genotypes without data
  
  tfreq <- table(temp[, geno])
  ng.0 <- sum(tfreq == 0)

  # Return
  
  list(ng.0 = ng.0, ng.2 = ng.2, nmis = nmis, pmis = pmis,
       ng = ng, nr = nr, tfreq = tfreq, nmis.fac = nmis.fac)
  
}
