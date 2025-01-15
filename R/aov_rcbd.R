#' ANOVA for a RCBD
#'
#' Fit an analysis of variance model for a RCBD.
#' @param dfr The name of the data frame.
#' @param y The name of the column for the variable to analyze.
#' @param geno The name of the column that identifies the genotypes.
#' @param rep The name of the column that identifies the replications.
#' @param maxp Maximum allowed proportion of missing values to estimate, default is 10\%.
#' @details If data is unbalanced, missing values are estimated up to an specified maximum
#' proportion, 10\% by default.
#' @return It returns an ANOVA table.
#' @author Raul Eyzaguirre.
#' @examples
#' # Get a copy with some missing values for trw and run ANOVA
#' tmp <- pjpz09
#' tmp[c(10, 20, 30), "trw"] <- NA
#' aov.rcbd(tmp, "trw", "geno", "rep")
#' @export

aov.rcbd <- function(dfr, y, geno, rep, maxp = 0.1) {

  # Everything as character

  dfr[, geno] <- as.character(dfr[, geno])
  dfr[, rep] <- as.character(dfr[, rep])

  # Check data
  
  lc <- ck.rcbd(dfr, y, geno, rep)

  # Estimate missing values and report errors from mve.rcbd
  
  y.est <- paste0(y, ".est")
  
  if (lc$ng.0 > 0 | lc$nrep == 1 | lc$ng.mult > 0 | lc$nmis > 0 | lc$nmis.fac > 0) {
    dfr[, y] <- mve.rcbd(dfr, y, geno, rep, maxp, tol = 1e-06)[, y.est]
    warning(paste0("The data set is unbalanced, ",
                   format(lc$pmis * 100, digits = 3),
                   "% missing values estimated."))
  }

  # ANOVA

  model <- aov(dfr[, y] ~ dfr[, geno] + dfr[, rep])
  model$terms[[2]] <- y
  
  at <- anova(model)
  
  rownames(at)[1:2] <- c(geno, rep)
  
  # Correction for missing values
  
  if (lc$nmis > 0) {
    at[3, 1] <- at[3, 1] - lc$nmis
    at[3, 3] <- at[3, 2] / at[3, 1]
    at[c(1, 2), 4] <- at[c(1, 2), 3] / at[3, 3]
    at[c(1, 2), 5] <- pf(at[c(1, 2), 4], at[c(1, 2), 1], at[3, 1], lower.tail = FALSE)
  }
  
  # Return

  at
  
}
