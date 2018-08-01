#' Check number of genotypes and replications
#' 
#' This function cheks the number of genotypes and replications for different designs.
#' 
#' @param geno The genotypes
#' @param design The statistical design
#' @param dfr The name of the data frame
#' @return The number of genotypes (\code{ng}), the number of checks (\code{ng.check}),
#' the list of genotypes (\code{lg}), the list of checks (\code{lg.check}),
#' and the number of replications (\code{nr}).
#' @author Raul Eyzaguirre.
#' @examples 
#' # Create a design
#' dfr <- cr.rcbd(1:20, 3, 10)
#' dfr <- dfr$book
#' # Check the design
#' ck.fs('geno', 'block', 'rcbd', dfr = dfr)
#' @export

ck.fs <- function(geno, rep = NULL, design = c('crd', 'rcbd', 'abd'), dfr) {

  # match arguments
  
  design <- match.arg(design)
  
  # Check and remove rows with missing values for factors
  
  dfr <- rm.fna(c(geno, rep), dfr)$dfr

  # Number of genotypes and checks
  
  if (design == 'abd') {
    tfreq <- data.frame(table(dfr[, geno]))
    lg.check <- as.character(tfreq[tfreq$Freq > 1, 1])
    lg <- as.character(tfreq[tfreq$Freq == 1, 1])
    ng.check <- length(lg.check)
    ng <- length(lg)
  } else {
    lg <- as.character(unique(dfr[, geno]))
    ng <- length(lg)
    lg.check <- NULL
    ng.check <- NULL
  }
  
  # Number of replications
  
  if (design %in% c('abd', 'rcbd')) {
    nr <- length(unique(dfr[, rep]))
  } else {
    tfreq <- table(dfr[, geno])
    nr <- max(tfreq)
  }
  
  # Return
  
  list(ng = ng, ng.check = ng.check, lg = lg, lg.check = lg.check, nr = nr)
  
}
