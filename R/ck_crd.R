#' Check data for a CRD
#'
#' This function checks the frequencies of genotypes in a CRD.
#' @param trait The trait to analyze.
#' @param geno The genotypes.
#' @param dfr The name of the data frame.
#' @return The number of genotypes (\code{ng}), the number of genotypes without
#' data (\code{ng.0}), the number of replications (\code{nrep}), and the number
#' of rows in the data frame with missing values for factors (\code{nmis.fac}).
#' @author Raul Eyzaguirre.
#' @examples
#' # Create design
#' dfr <- cr.crd(1:50, 3, 10)
#' dfr <- dfr$book
#' # Create some random data
#' dfr$y <- rnorm(150)
#' # Delete some values
#' dfr[c(1, 5, 56, 77, 111), 'y'] <- NA
#' # Delete some values for classification factors
#' dfr[c(27, 48), 'geno'] <- NA
#' # Check the design
#' ck.crd('y', 'geno', dfr)
#' @export

ck.crd <- function(trait, geno, dfr) {
  
  # Check factor structure
  
  out <- ck.fs(geno, NULL, dfr)
  dfr <- out$dfr
  ng <- out$nt
  nrep <- out$nrep
  nmis.fac <- out$nmis.fac
  
  # Frequencies for genotypes
  
  out <- ck.fq(trait, geno, NULL, dfr)
  tf <- out$tf

  # Number of genotypes without data
  
  ng.0 <- sum(tf == 0)
  
  # Return
  
  list(ng.0 = ng.0, ng = ng, nrep = nrep, tf = tf, nmis.fac = nmis.fac)
  
}
