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
#' the number of genotypes (\code{ng}), the number of replications (\code{nrep}),
#' a table with valid cases for each genotype (\code{tf}), a table with valid 
#' cases for each genotype in each replication (\code{tfr}), and the number of
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

ck.rcbd <- function(trait, geno, rep, dfr) {
  
  # Check factor structure
  
  out <- ck.fs(geno, rep, dfr)
  dfr <- out$dfr
  ng <- out$nt
  nrep <- out$nrep
  nmis.fac <- out$nmis.fac
  
  # Frequencies for genotypes and replications
  
  out <- ck.fq(trait, geno, rep, dfr)
  tf <- out$tf
  tfr <- out$tfr
  nmis <- out$nmis
  pmis <- out$pmis
  
  # Number of genotypes without data
  
  ng.0 <- sum(tf == 0)

  # Number of genotypes that appear more than once in a given replication

  ng.mult <- sum(tfr > 1)
  
  # Return
  
  list(ng.0 = ng.0, ng.mult = ng.mult, nmis = nmis, pmis = pmis, ng = ng,
       nrep = nrep, tf = tf, tfr = tfr, nmis.fac = nmis.fac)
  
}
