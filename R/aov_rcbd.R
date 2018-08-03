#' ANOVA for a RCBD
#'
#' Fit an analysis of variance model for a RCBD.
#' @param trait The trait to analyze.
#' @param geno The genotypes.
#' @param rep The replications.
#' @param dfr The name of the data frame.
#' @param maxp Maximum allowed proportion of missing values to estimate, default is 10\%.
#' @details If data is unbalanced, missing values are estimated up to an specified maximum
#' proportion, 10\% by default.
#' @return It returns an ANOVA table.
#' @author Raul Eyzaguirre.
#' @examples
#' # Get a copy with some missing values for trw and run ANOVA
#' temp <- pjpz09
#' temp[c(10, 20, 30), "trw"] <- NA
#' aov.rcbd("trw", "geno", "rep", temp)
#' @export

aov.rcbd <- function(trait, geno, rep, dfr, maxp = 0.1) {

  # Everything as character

  dfr[, geno] <- as.character(dfr[, geno])
  dfr[, rep] <- as.character(dfr[, rep])

  # Check data
  
  lc <- ck.rcbd(trait, geno, rep, dfr)

  # Estimate missing values and report errors from mve.rcbd
  
  trait.est <- paste0(trait, ".est")
  
  if (lc$ng.0 > 0 | lc$nrep == 1 | lc$ng.mult > 0 | lc$nmis > 0 | lc$nmis.fac > 0) {
    dfr[, trait] <- mve.rcbd(trait, geno, rep, dfr, maxp, tol = 1e-06)[, trait.est]
    warning(paste0("The data set is unbalanced, ",
                   format(lc$pmis * 100, digits = 3),
                   "% missing values estimated."))
  }

  # ANOVA

  model <- aov(dfr[, trait] ~ dfr[, geno] + dfr[, rep])
  model$terms[[2]] <- trait
  
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
