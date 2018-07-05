#' ANOVA for a factorial experiment
#'
#' Fit an analysis of variance model for a factorial experiment with a CRD or RCBD
#' @param trait The trait to analyze.
#' @param factors The factors.
#' @param rep The replications or blocks.
#' @param design The statistical design, \code{crd} or \code{rcbd}.
#' @param data The name of the data frame.
#' @param maxp Maximum allowed proportion of missing values to estimate, default is 10\%.
#' @author Raul Eyzaguirre.
#' @details If data is unbalanced, missing values are estimated up to an specified maximum
#' proportion, 10\% by default.
#' @return It returns the ANOVA table.
#' @examples
#' aov.f("asc.dw", c("geno", "treat"), "rep", "crd", asc)
#' @importFrom stats anova
#' @export

aov.f <- function(trait, factors, rep, design = c("crd", "rcbd"), data, maxp = 0.1) {

  # match arguments
  
  design <- match.arg(design)

  # Number of factors
  
  nf <- length(factors)

  # Everything as factor
  
  for (i in 1:nf)
    data[, factors[i]] <- factor(data[, factors[i]])

  data[, rep] <- factor(data[, rep])

  # Check data

  lc <- ck.f(trait, factors, rep, data)

  # Estimate missing values and report errors from mve.f
  
  trait.est <- paste(trait, ".est", sep = "")
  
  if (lc$c1 == 0 | lc$c2 == 0 | lc$c3 == 0 | lc$c4 == 0) {
    data[, trait] <- mve.f(trait, factors, rep, design, data, maxp)[, trait.est]
    warning(paste("The data set is unbalanced, ",
                  format(lc$pmis * 100, digits = 3),
                  "% missing values estimated.", sep = ""))
  }

  # ANOVA
  
  expr <- paste(trait, '~', factors[1])
  for (i in 2:nf)
    expr <- paste(expr, '*', factors[i])

  if (design == "crd") {
    ff <- as.formula(expr)
    model <- aov(ff, data = data)
  }
  
  if (design == "rcbd") {
    expr <- paste(expr, '+', rep)
    ff <- as.formula(expr)
    model <- aov(ff, data = data)
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
