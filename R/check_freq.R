#' Check frequencies for complete replications
#'
#' Check frequencies for designs with complete replications and one or several
#' environments. This is a wrapper for \code{ck.rcbd} and \code{ck.f} functions.
#' @param trait The trait to analyze.
#' @param geno Genotypes.
#' @param env Environments, \code{NULL} if there are no environments.
#' @param rep The replications.
#' @param dfr The name of the data frame.
#' @return Information about the balance, missing values, and replications of the design.
#' @author Raul Eyzaguirre.
#' @examples
#' check.freq("trw", "geno", NULL, "rep", pjpz09)
#' check.freq("rytha", "geno", "env", "rep", megaclones)
#' @export

check.freq <- function(trait, geno, env = NULL, rep, dfr) {
  
  # Levels for replications
  
  lrep <- sort(unique(dfr[, rep]))
  
  # Main text

  cat('----------------------------------------\n')
  cat('Check frequencies for trait', trait, '\n')
  cat('----------------------------------------\n')
  cat('\n')

  if (is.null(env)) {
    
    # Run check for rcbd
    
    lc <- ck.rcbd(trait, geno, rep, dfr)
    
    # Write warnings
    
    if (lc$nmis.fac > 0) {
      cat('There are missing values for classification factors. \n')
      cat('\n')
    }
  
    if (lc$ng.0 > 0) {
      lista <- names(lc$tf[lc$tf == 0])
      cat('There are genotypes without data: \n', lista, '\n')
      cat('\n')
    }

    if (lc$nrep == 1) {
      cat('There is only one replication. \n')
      cat('\n')
    }
    
    if (lc$ng.mult > 0) {
      temp <- lc$tfr > 1
      cat('There are genotypes that appear more than once in a given replication: \n')
      for (i in 1:lc$nrep) {
        if (sum(temp[, i]) > 0) {
          lista <- rownames(temp)[temp[, i]]
          cat(paste0('- Replication ', lrep[i], ':'), lista, '\n')
        }
      }
      cat('\n')
    }
    
    if (lc$nmis > 0) {
      cat("There are missing values:", format(lc$pmis * 100, digits = 3), '% \n')
      cat('\n')
    }
    
    # OK message
    
    if (lc$ng.0 == 0 & lc$nrep > 1 & lc$ng.mult == 0 & lc$nmis == 0 & lc$nmis.fac == 0) {
      cat('OK \n')
      cat('\n')
    }
    
  } else {
    
    # Levels for environments
    
    le <- sort(unique(dfr[, env]))
    
    # Run check for factorial

    lc <- ck.f(trait, c(geno, env), rep, dfr)
    
    # Write warnings
    
    if (lc$nmis.fac > 0) {
      cat('There are missing values for classification factors. \n')
      cat('\n')
    }

    if (lc$nt.0 > 0) {
      temp <- lc$tf == 0
      cat('There are genotypes without data in a given environment: \n')
      for (i in 1:lc$nl[2]) {
        if (sum(temp[, i]) > 0) {
          lista <- rownames(temp)[temp[, i]]
          cat(paste0('- Environment ', le[i], ':'), lista, '\n')
        }
      }
      cat('\n')
    }
    
    if (lc$nrep == 0) {
      cat('There is only one replication. \n')
      cat('\n')
    }
    
    if (lc$nt.mult > 0) {
      temp <- lc$tfr > 1
      cat('There are genotypes that appear more than once in a given replication: \n')
      for (i in 1:lc$nl[2]) {
        for (j in 1:lc$nrep) {
          if (sum(temp[, i, j]) > 0) {
            lista <- rownames(temp)[temp[, i, j]]
            cat(paste0('- Environment ', le[i], ', replication ', lrep[j], ':'), lista, '\n')
          }
        }
      }
      cat('\n')
    }
    
    if (lc$nmis > 0) {
      cat("There are missing values:", format(lc$pmis * 100, digits = 3), '% \n')
      cat('\n')
    }
    
    # OK message
    
    if (lc$nt.0 == 0 & lc$nrep > 1 & lc$nt.mult == 0 & lc$nmis == 0 & lc$nmis.fac == 0) {
      cat('OK \n')
      cat('\n')
    }    
  }
  
}
