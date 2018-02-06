#' Check design
#'
#' Check frequencies for designs with complete replications and one or two factors.
#' @param trait The trait to analyze.
#' @param geno Genotypes.
#' @param env Environments.
#' @param rep The replications.
#' @param data The name of the data frame.
#' @return Information about the balance, missing values, and replications of the design.
#' @author Raul Eyzaguirre.
#' @export

check.design <- function(trait, geno, env = NULL, rep, data) {
  
  data[, rep] <- factor(data[, rep])
  lr <- levels(data[, rep])

  cat('------------------------------\n')
  cat('Check design for trait', trait, '\n')
  cat('------------------------------\n')
  cat('\n')

  if (is.null(env)) {
    
    lc <- check.rcbd(trait, geno, rep, data)
  
    if (lc$c1 == 0) {
      lista <- apply(lc$tfreq, 1, sum)
      lista <- names(lista[lista == 0])
      cat('There are genotypes without data: \n', lista, '\n')
      cat('\n')
    }

    if (lc$c2 == 0) {
      cat('There is only one replication \n')
      cat('\n')
    }
    
    if (lc$c3 == 0) {
      tf <- lc$tfreq > 1
      cat('There are genotypes that appear more than once in a given replication: \n')
      for (i in 1:lc$nr) {
        if (sum(tf[, i]) > 0) {
          lista <- rownames(tf)[tf[, i]]
          cat(paste('- Replication ', lr[i], ':', sep = ""), lista, '\n')
        }
      }
      cat('\n')
    }
    
    if (lc$c4 == 0) {
      cat("There are missing values:", format(lc$pmis * 100, digits = 3), '% \n')
      cat('\n')
    }
    
    if (lc$c1 == 1 & lc$c2 == 1 & lc$c3 == 1 & lc$c4 == 1) {
      cat('OK \n')
      cat('\n')
    }
    
  } else {
    
    data[, env] <- factor(data[, env])
    le <- levels(data[, env])

    lc <- check.2f(trait, geno, env, rep, data)

    if (lc$c1 == 0) {
      tf <- lc$tfreq == 0
      cat('There are genotypes without data in a given environment: \n')
      for (i in 1:lc$nb) {
        if (sum(tf[, i]) > 0) {
          lista <- rownames(tf)[tf[, i]]
          cat(paste('- Environment ', le[i], ':', sep = ""), lista, '\n')
        }
      }
      cat('\n')
    }
    
    if (lc$c2 == 0) {
      cat('There is only one replication \n')
      cat('\n')
    }
    
    if (lc$c3 == 0) {
      tf <- lc$tfreqr > 1
      cat('There are genotypes that appear more than once in a given replication: \n')
      for (i in 1:lc$nb) {
        for (j in 1:lc$nr) {
          if (sum(tf[, i, j]) > 0) {
            lista <- rownames(tf)[tf[, i, j]]
            cat(paste('- Environment ', le[i], ', replication ', lr[j], ':', sep = ""),
                lista, '\n')
          }
        }
      }
      cat('\n')
    }
    
    if (lc$c4 == 0) {
      cat("There are missing values:", format(lc$pmis * 100, digits = 3), '% \n')
      cat('\n')
    }
    
    if (lc$c1 == 1 & lc$c2 == 1 & lc$c3 == 1 & lc$c4 == 1) {
      cat('OK \n')
      cat('\n')
    }    
  }
}
