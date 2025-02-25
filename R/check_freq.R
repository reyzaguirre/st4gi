#' Check frequencies for complete replications
#'
#' Check frequencies for designs with complete replications and one or several
#' environments. This is a wrapper for \code{ck.rcbd} and \code{ck.f} functions.
#' @param dfr The name of the data frame.
#' @param y The name of the column for the variable to analyze.
#' @param geno The name of the column that identifies the genotypes.
#' @param rep The name of the column that identifies the replications.
#' @param env The name of the column that identifies the environments,
#' \code{NULL} if there are no environments.
#' @return Information about the balance, missing values, and replications of the design.
#' @author Raul Eyzaguirre.
#' @examples
#' check.freq(pjpz09, "trw", "geno", "rep")
#' check.freq(megaclones, "rytha", "geno", "rep", "env")
#' @export

check.freq <- function(dfr, y, geno = 'geno', rep = 'rep', env = NULL) {
  
  # Levels for replications
  
  lrep <- sort(unique(dfr[, rep]))
  
  # Main text

  cat('----------------------------------------\n')
  cat('Check frequencies for variable', y, '\n')
  cat('----------------------------------------\n')
  cat('\n')

  if (is.null(env)) {
    
    # Run check for rcbd
    
    lc <- ck.rcbd(dfr, y, geno, rep)
    
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
      tmp <- lc$tfr > 1
      cat('There are genotypes that appear more than once in a given replication: \n')
      for (i in 1:lc$nrep) {
        if (sum(tmp[, i]) > 0) {
          lista <- rownames(tmp)[tmp[, i]]
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

    lc <- ck.f(dfr, y, c(geno, env), rep)
    
    # Write warnings
    
    if (lc$nmis.fac > 0) {
      cat('There are missing values for classification factors. \n')
      cat('\n')
    }

    if (lc$nt.0 > 0) {
      tmp <- lc$tf == 0
      cat('There are genotypes without data in a given environment: \n')
      for (i in 1:lc$nl[2]) {
        if (sum(tmp[, i]) > 0) {
          lista <- rownames(tmp)[tmp[, i]]
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
      tmp <- lc$tfr > 1
      cat('There are genotypes that appear more than once in a given replication: \n')
      for (i in 1:lc$nl[2]) {
        for (j in 1:lc$nrep) {
          if (sum(tmp[, i, j]) > 0) {
            lista <- rownames(tmp)[tmp[, i, j]]
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
