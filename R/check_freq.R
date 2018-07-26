#' Check frequencies
#'
#' Check frequencies for designs with complete replications and one or two factors.
#' This is a wrapper for \code{ck.rcbd} and \code{ck.f} functions.
#' @param trait The trait to analyze.
#' @param geno Genotypes.
#' @param env Environments.
#' @param rep The replications.
#' @param data The name of the data frame.
#' @return Information about the balance, missing values, and replications of the design.
#' @author Raul Eyzaguirre.
#' @export

check.freq <- function(trait, geno, env = NULL, rep, data) {
  
  # Levels for replications
  
  lr <- sort(unique(data[, rep]))
  
  # Main text

  cat('------------------------------\n')
  cat('Check design for trait', trait, '\n')
  cat('------------------------------\n')
  cat('\n')

  if (is.null(env)) {
    
    # Run check for rcbd
    
    lc <- ck.rcbd(trait, geno, rep, data)
    
    # Write warnings
    
    if (lc$nmis.fact > 0) {
      cat('There are missing values for classification factors. \n')
      cat('\n')
    }
  
    if (lc$c1 == 0) {
      lista <- apply(lc$tfreq, 1, sum)
      lista <- names(lista[lista == 0])
      cat('There are genotypes without data: \n', lista, '\n')
      cat('\n')
    }

    if (lc$c2 == 0) {
      cat('There is only one replication. \n')
      cat('\n')
    }
    
    if (lc$c3 == 0) {
      tf <- lc$tfreq > 1
      cat('There are genotypes that appear more than once in a given replication: \n')
      for (i in 1:lc$nr) {
        if (sum(tf[, i]) > 0) {
          lista <- rownames(tf)[tf[, i]]
          cat(paste0('- Replication ', lr[i], ':'), lista, '\n')
        }
      }
      cat('\n')
    }
    
    if (lc$c4 == 0) {
      cat("There are missing values:", format(lc$pmis * 100, digits = 3), '% \n')
      cat('\n')
    }
    
    # OK message
    
    if (lc$c1 == 1 & lc$c2 == 1 & lc$c3 == 1 & lc$c4 == 1 & lc$nmis.fact == 0) {
      cat('OK \n')
      cat('\n')
    }
    
  } else {
    
    # Levels for environments
    
    le <- sort(unique(data[, env]))
    
    # Run check for factorial

    lc <- ck.f(trait, c(geno, env), rep, data)
    
    # Write warnings
    
    if (lc$nmis.fact > 0) {
      cat('There are missing values for classification factors. \n')
      cat('\n')
    }

    if (lc$c1 == 0) {
      tf <- lc$tfreq == 0
      cat('There are genotypes without data in a given environment: \n')
      for (i in 1:lc$nl[2]) {
        if (sum(tf[, i]) > 0) {
          lista <- rownames(tf)[tf[, i]]
          cat(paste0('- Environment ', le[i], ':'), lista, '\n')
        }
      }
      cat('\n')
    }
    
    if (lc$c2 == 0) {
      cat('There is only one replication. \n')
      cat('\n')
    }
    
    if (lc$c3 == 0) {
      tf <- lc$tfreqr > 1
      cat('There are genotypes that appear more than once in a given replication: \n')
      for (i in 1:lc$nl[2]) {
        for (j in 1:lc$nr) {
          if (sum(tf[, i, j]) > 0) {
            lista <- rownames(tf)[tf[, i, j]]
            cat(paste0('- Environment ', le[i], ', replication ', lr[j], ':'), lista, '\n')
          }
        }
      }
      cat('\n')
    }
    
    if (lc$c4 == 0) {
      cat("There are missing values:", format(lc$pmis * 100, digits = 3), '% \n')
      cat('\n')
    }
    
    # OK message
    
    if (lc$c1 == 1 & lc$c2 == 1 & lc$c3 == 1 & lc$c4 == 1 & lc$nmis.fact == 0) {
      cat('OK \n')
      cat('\n')
    }    
  }
  
}
