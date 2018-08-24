#' Optimize Westcott adjustment
#'
#' This function optimizes the parameter values of the Westcott adjustment.
#' @param trait The trait to adjust.
#' @param geno The genotypes.
#' @param replicated A \code{yes/no} column to identify replicated genotypes.
#' @param ck1 Name of check 1.
#' @param ck2 Name of check 2.
#' @param row Label for rows.
#' @param col Label for columns.
#' @param ncb Number of columns between two check columns.
#' @param nrs Set of values for the \code{nrs} param on \code{aj.w}.
#' @param method Set of values for the \code{method} param on \code{aj.w}.
#' @param ind Set of values for the \code{ind} param on \code{aj.w}.
#' @param p Set of values for the \code{p} param on \code{aj.w}.
#' @param opt Criteria for optimization, \code{F.value}, \code{MSE} or \code{CV}.
#' See details.
#' @param dfr The name of the data frame.
#' @details If \code{optim = F.value} then the F value for ANOVA is maximized.
#' If \code{optim = MSE} then the MSE is minimized.
#' If \code{optim = CV} then the CV is minimized.
#' @return It returns optimal set of values under the selected criteria.
#' @author Raul Eyzaguirre.
#' @export

optim.w <- function(trait, geno, replicated, ck1, ck2, row, col, ncb, nrs, method, ind,
                    p, opt = c("F.value", "MSE", "CV"), dfr) {
  
  # Match arguments
  
  opt <- match.arg(opt)

  # All settings to run
  
  nrs.l <- length(nrs)
  method.l <- length(method)
  ind.l <- length(ind)
  p.l <- length(p)
  
  total <- nrs.l * method.l * ind.l * p.l
  
  # Anova without adjustment
  
  temp <- dfr[dfr[, replicated] == 'yes', ]
  ff <- as.formula(paste(trait, '~', geno))
  at <- anova(aov(ff, temp))
  MSE <- at[2, 3]
  F.value <- at[1, 4]
  CV <- sqrt(at[2, 3]) / mean(temp[, trait], na.rm = T) * 100
  
  # Original values
  
  best.method <- ''
  best.ind <- ''
  best.nrs <- ''
  best.p <- ''
  best.MSE <- MSE
  best.F <- F.value 
  best.CV <- CV
  
  # Run adjustments
  
  trait.aj <- paste0(trait, '.aj')
  i <- 0
  
  for (method.v in method[1:method.l])
    for (ind.v in ind[1:ind.l])
      for (nrs.v in nrs[1:nrs.l])
        for (p.v in p[1:p.l]) {

          dfr <- aj.w(trait, geno, ck1, ck2, row, col, ncb, nrs.v, method.v, ind.v, p.v, dfr)
          temp <- dfr[dfr[, replicated] == 'yes', ]
          ff <- as.formula(paste(trait.aj, '~', geno))
          at.aj <- anova(aov(ff, temp))
          
          new.F <- at.aj[1, 4]
          new.MSE <- at.aj[2, 3]
          new.CV <- sqrt(at.aj[2, 3]) / mean(temp[, trait.aj], na.rm = T) * 100
          
          cond1 <- opt == "F.value" & new.F > best.F
          cond2 <- opt == "MSE" & new.MSE < best.MSE
          cond3 <- opt == "CV" & new.CV < best.CV
          
          if (cond1 | cond2 | cond3) {
            best.method <- method.v
            best.ind <- ind.v
            best.nrs <- nrs.v
            best.p <- p.v
            best.MSE <- new.MSE
            best.F <- new.F
            best.CV <- new.CV
            cat('\n-------------------------------\n')
            cat("Method = ", best.method, '\n')
            cat("Ind    = ", best.ind, '\n')
            cat("nrs    = ", best.nrs, '\n')
            cat("p      = ", best.p, '\n')
            cat("MSE    = ", MSE, '\n')
            cat("MSE.aj = ", best.MSE, '\n')
            cat("F      = ", F.value, '\n')
            cat("F.aj   = ", best.F, '\n')
            cat("CV     = ", CV, '\n')
            cat("CV.aj  = ", best.CV, '\n')
            cat('-------------------------------\n')
          }
          
          i <- i + 1
          cat(round(i / total * 100, 0), '% ', sep = "")
          
        }
  
}
