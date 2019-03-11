#' ANOVA for a factorial experiment
#'
#' Fit an analysis of variance model for a factorial experiment with a CRD or RCBD.
#' @param trait The trait to analyze.
#' @param factors The factors.
#' @param rep The replications or blocks, \code{NULL} for a CRD.
#' @param dfr The name of the data frame.
#' @param maxp Maximum allowed proportion of missing values to estimate, default is 10\%.
#' @return It returns the ANOVA table.
#' @author Raul Eyzaguirre.
#' @examples
#' aov.f("asc.dw", c("geno", "treat"), NULL, asc)
#' @importFrom stats anova
#' @export

aov.f <- function(trait, factors, rep, dfr, maxp = 0.1) {

  # Check data
  
  lc <- ck.f(trait, factors, rep, dfr)

  # Everything as character
  
  for (i in 1:lc$nf)
    dfr[, factors[i]] <- as.character(dfr[, factors[i]])

  dfr[, rep] <- as.character(dfr[, rep])

  # Estimate missing values and report errors from mve.f
  
  trait.est <- paste0(trait, ".est")
  
  if (lc$nt.0 > 0 | lc$nrep == 1 | lc$nt.mult > 0 | lc$nmis > 0 |
      lc$nmis.fac > 0 | sum(lc$nl < 2) > 0) {
    dfr[, trait] <- mve.f(trait, factors, rep, dfr, maxp)[, trait.est]
    warning(paste0("The data set is unbalanced, ",
                   format(lc$pmis * 100, digits = 3),
                   "% missing values estimated."))
  }

  # ANOVA
  
  expr <- paste(trait, '~', factors[1])
  for (i in 2:lc$nf)
    expr <- paste(expr, '*', factors[i])

  if (is.null(rep)) {
    ff <- as.formula(expr)
    model <- aov(ff, dfr)
  }
  
  if (!is.null(rep)) {
    expr <- paste(expr, '+', rep)
    ff <- as.formula(expr)
    model <- aov(ff, dfr)
  }
  
  at <- anova(model)
  
  # Residuals row
  
  rr <- dim(at)[1]
  
  # Correction for missing values
  
  if (lc$nmis > 0) {
    
    # Correction for Residuals row
    
    at[rr, 1] <- at[rr, 1] - lc$nmis
    at[rr, 3] <- at[rr, 2] / at[rr, 1]
    
    # Correction for F and p values
    
    for (i in 1:(rr - 1)) {
      at[i, 4] <- at[i, 3] / at[rr, 3]
      at[i, 5] <- pf(at[i, 4], at[i, 1], at[rr, 1], lower.tail = FALSE)
    }
  }
  
  # Return
  
  at
  
}
