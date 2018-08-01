#' Check data for a CRD
#'
#' This function checks the frequencies of genotypes in a CRD.
#' @param trait The trait to analyze.
#' @param geno The genotypes.
#' @param dfr The name of the data frame.
#' @return The number of genotypes (\code{ng}), the number of genotypes without
#' data (\code{ng.0}), the number of replications (\code{nr}), a table with
#' frequencies of valid cases for each genotype (\code{tfreq}), and the number
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
  
  # Check and remove rows with missing values for factors
  
  out <- rm.fna(geno, dfr)
  dfr <- out$dfr
  nmis.fac <- out$nmis.fac

  # Number of genotypes and replications
  
  out <- ck.gr(geno, NULL, 'crd', dfr)
  ng <- out$ng
  nr <- out$nr

  # Genotypes as factor to preserve levels in the table of frequencies
  
  temp <- dfr
  temp[, geno] <- factor(temp[, geno])
  
  # Frequencies for genotypes
  
  temp <- temp[!is.na(temp[, trait]), ]
  tfreq <- table(temp[, geno])

  # Number of genotypes without data
  
  ng.0 <- sum(tfreq == 0)
  
  # Return
  
  list(ng.0 = ng.0, ng = ng, nr = nr, tfreq = tfreq, nmis.fac = nmis.fac)
  
}
